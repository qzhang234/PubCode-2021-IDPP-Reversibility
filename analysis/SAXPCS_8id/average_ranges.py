"""Reduce the raw 8-ID-I result files into the small set of files in ``data/``.

This is the ONE place where the beamline's reprocessed NeXus results are read.
Everything downstream -- saxpcs.py, saxs_evolution.py, g2_grid_SI.py,
thermal_cycle.py, contrast_calibration.py -- reads only what this script writes
into ``data/``, so the repository carries enough reduced data for the whole
analysis to be rerun without access to the beamline storage.

It produces two kinds of file:

1. AVERAGES over explicit frame-number ranges, one file per range, for every
   group listed in ``FILE_RANGES``.  This is the bulk of the output and the
   original purpose of the script.
2. The per-acquisition TEMPERATURE TRACE of the thermal-cycling sequence
   (``TRACE_GROUPS`` -> ``TRACE_CSV``).  Figure S4a plots the temperature of
   every one of its 2325 acquisitions, which is not an average and cannot be
   recovered from the averaged files.

For each range this script:
  * globs the group's result files and selects those whose frame number is in
    the range,
  * removes outliers in two passes: a cross-correlation threshold on the SHAPE
    of g2, g2_err and saxs_1d, then a robust cut on acquisitions sitting far
    above the group median in the 0.004-0.008 A^-1 band (``spike_removal``),
    and finally any frames listed for that group in ``MANUAL_EXCLUDE``,
  * averages g2, g2_err and saxs_1d over the surviving files,
  * ALSO averages the two scalar ion-chamber monitors -- the incident
    (upstream) and transmitted (downstream) beam intensities -- over the same
    surviving files.  These are needed downstream in saxpcs.py to compute the
    absolute scattering cross-section coefficient of each averaged dataset, so
    each frame range carries its own mean flux/transmission.  The scalars are
    averaged SEPARATELY from g2/saxs_1d and are deliberately NOT put through
    outlier_removal (that routine cross-correlates every field it is handed and
    is meaningless for single scalars),
  * writes the averaged data into a new HDF file that keeps the same structure
    as the originals -- a full copy of the first file in the range is used as
    the template, so every averaged file is a complete, self-describing NeXus
    file rather than a stripped extract.

The averaged file name starts with 'Average_' and carries the frame range so it
is easy to tell apart from the raw files, e.g.::

    Average_B0147_S3_7_300C10p_att00_Rq0_00951_01050_results.hdf

In the averaged file:
  * ``/entry/start_time`` is set to the acquisition time of the first file in
    the range, read from that file's own ``/entry/start_time``,
  * ``/xpcs/average/file_list`` is added, listing every file included in the
    average (i.e. the files that survived outlier removal).

Run it once, before any figure script::

    python average_ranges.py            # everything
    python average_ranges.py B0083      # just one group's averages
"""

import sys
import os
import re
import glob
import shutil

from datetime import datetime

import numpy as np
import h5py

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from common.utils import outlier_removal, average_datasets

# --- PARAMETERS ---
# BEAMLINE-ONLY: the raw reprocessed NeXus results.  /gdata is mounted on the
# 10.x hosts (amber) and not on the 164.x analysis machines, so this script has
# to run there; every other 8-ID-I script reads only the data/ folder it writes.
#
# The ion chambers, the ring current and the acquisition timestamps in this
# tree were restored from the 2022 Data Exchange archive by
# restore_dx_metadata.py, which see: the 2025 NeXus rewrite had replaced them
# with constants.
prefix = '/gdata/s8id-dmdtn/2022-1/babnigg202203_nexus/reprocess_results'

# Groups to average: file-name header -> list of (start, end) frame ranges
# (inclusive), one averaged file per range.  Usually one range per group.

FILE_RANGES = {
    # --- Figures 3, S8, S9, S10: the isothermal series and its two references
    'B0147': [(1, 200),
        (200, 350),
        (351, 500),
        (501, 650),
        (651, 800),
        (801, 950),
        (951, 1050),
        (1051, 1150),
        (1151, 1250),
        (1251, 1313)],
    'B0146': [(1, 50)],
    'D0138': [(1, 50)],

    # --- Figure S4: seven thermal cycles of sample S1, 2022-03-04 -----------
    # Each cycle is a 6 C group of 50 acquisitions followed by a 270-acquisition
    # ramp to 34 C, of which two ten-acquisition windows near the top are used
    # (they land near 31.9 and 33.8 C).  A range simply selects whatever files
    # exist inside it, so the nominal (1, 50) is right even where a few
    # acquisitions are missing (B0078 has 49, D0077 has 39).
    'B0075': [(1, 50)], 'B0076': [(241, 250), (261, 270)],
    'B0078': [(1, 50)], 'B0079': [(241, 250), (261, 270)],
    'B0081': [(1, 50)], 'B0082': [(241, 250), (261, 270)],
    'B0083': [(1, 50)], 'B0084': [(241, 250), (261, 270)],
    'B0085': [(1, 50)], 'B0086': [(241, 250), (261, 270)],
    'B0087': [(1, 50)], 'B0088': [(241, 250), (261, 270)],
    'B0089': [(1, 50)], 'B0090': [(241, 250), (261, 270)],
    # the buffer for those cycles, held cold and measured twice
    'D0077': [(1, 50)],
    'D0080': [(1, 50)],

    # --- Figure S7: the nano-porous glass speckle-contrast standard ---------
    'F0145': [(1, 50)],
}

# Groups whose results carry a boost_corr output suffix, set with its -u flag.
# F0145 was re-correlated against a qmap with three hot pixels masked out, and
# the tag is in the file name:
#     F0145_10nm_Glass_006C_att00_Rq0_00001_BadpixRm_results.hdf
# Both generations sit in the same directory, so the suffix has to be part of
# the selection: a bare F0145*_results.hdf glob matches all 100 files and would
# silently average the two reductions together.
GROUP_SUFFIX = {'F0145': 'BadpixRm'}

# --- SECOND-STAGE CUT: acquisitions far above the group at low q ---------
# outlier_removal() cuts on the SHAPE of log10 I(q) across the whole q range, so
# an acquisition that is normal everywhere except the lowest bins stays nearly
# parallel to the group mean and survives it.  That is what a large object
# drifting through the 10 x 10 um beam during one 2 s exposure looks like.  This
# second cut catches those.
#
# The statistic is the per-acquisition mean of I(q) over SPIKE_BAND, the same
# window Figure S4b reports its intensities over, so the cut introduces no new
# choice of q range.  Note that it tests the intensity in that band, not the
# low-q rise as such: an acquisition that is bright across the whole detector
# fails it too, which is intended.
#
# The test is a one-sided iterated modified z-score (Iglewicz & Hoaglin).  An
# acquisition is dropped when its band mean lies more than SPIKE_Z robust
# standard deviations ABOVE the median, the scale being 1.4826 x MAD, with both
# statistics recomputed from the survivors until the set stops changing.  Median
# and MAD rather than mean and standard deviation, because a handful of 8x
# acquisitions inflate the ordinary statistics enough to hide themselves.
# One-sided, because the artefact only ever ADDS scattering.
#
# SPIKE_Z = 3 is conventional and is not tuned: in the worst group the same
# acquisitions are selected for any threshold between 2.5 and 3.0, and for
# SPIKE_BAND moved to 0.0032-0.006 or 0.003-0.005.
#
# Across the paper the cut removes 35 of the 1749 acquisitions that survive
# outlier_removal(), 2.0 %.  It is nearly inert on the isothermal series behind
# Figures 3 and S8-S10, taking 2 acquisitions of B0147, which move that group's
# absolute-scale coefficient by 0.03 % and leave every fitted g2 parameter
# unchanged.  It takes none at all from the contrast standard F0145.  Its work
# is on the Figure S4 thermal cycles, where transient scatterers are common:
# 4, 2, 1, 11, 4, 1 and 1 from the seven 6 C groups and 2 and 7 from the two
# buffers, against none from any of the fourteen 34 C windows.
SPIKE_BAND = (0.004, 0.008)          # A^-1
SPIKE_Z = 3.0


def static_q(fname):
    """q of every column of saxs_1d, read from the file's own static q map.

    scattering_1d is stored as one value per (q, phi) bin, flattened; the index
    map gives each bin's position in the (dim0 = q, dim1 = phi) grid.
    """
    with h5py.File(fname, 'r') as hf:
        idx = hf['/xpcs/qmap/static_index_mapping'][()]
        q_list = hf['/xpcs/qmap/static_v_list_dim0'][()]
        n_phi = hf['/xpcs/qmap/static_v_list_dim1'].shape[0]
    return q_list[idx // n_phi]


def spike_removal(saxs_1d, q, mask, band=SPIKE_BAND, z=SPIKE_Z):
    """Flag acquisitions whose low-q band mean sits far above the group median.

    ``saxs_1d`` is (n_acquisitions, n_bins) as stored in the file, ``q`` the
    matching q of every bin from static_q(), and ``mask`` the acquisitions
    still standing after outlier_removal().  The statistic is built from the
    WHOLE group, not from ``mask``'s survivors, so the two cuts are independent
    tests of the same acquisitions rather than a chain in which the second one
    inherits the first one's percentile; ``mask`` is applied only at the end.
    Returns a keep-mask.
    """
    sel = (q >= band[0]) & (q < band[1])
    stat = np.nanmean(saxs_1d[:, sel], axis=1)
    keep = np.ones(len(stat), dtype=bool)
    for _ in range(20):
        med = np.median(stat[keep])
        scale = 1.4826 * np.median(np.abs(stat[keep] - med))
        if not scale > 0:
            break
        new = keep & (stat - med <= z * scale)
        if np.array_equal(new, keep):
            break
        keep = new
    return mask & keep


# Acquisitions dropped BY HAND, on top of the two automatic cuts:
#   group -> (start, end) -> [frame numbers to drop]
# Empty, and kept only as the escape hatch for a defect the two automatic cuts
# cannot see.  The eleven acquisitions of B0083 that were once listed here by
# eye are now removed by spike_removal(), which takes eleven acquisitions of
# that group, ten of them the same: it also takes frame 9, and frame 47 of the
# hand list survives.  Those two sit at the boundary of the group's own scatter
# and no ranking separates them, which is why the by-eye list was replaced.
MANUAL_EXCLUDE = {}

# Groups whose per-acquisition temperature and time are written to TRACE_CSV,
# in the order they were measured.  Figure S4a needs every acquisition, not the
# group averages, so this is an extraction rather than a reduction.  The
# temperature is that of the stage zone holding the SAMPLE (qnw1) for every
# group, buffer acquisitions included: panel (a) is the thermal history of the
# sample, and the buffer sat in a different zone that was held at 6 C
# throughout.
TRACE_GROUPS = ['B0075', 'B0076', 'D0077', 'B0078', 'B0079', 'D0080',
                'B0081', 'B0082', 'B0083', 'B0084', 'B0085', 'B0086',
                'B0087', 'B0088', 'B0089', 'B0090']
TRACE_CSV = 'thermal_cycle_temperature.csv'
TEMP_PATH = '/entry/sample/qnw1_temperature'


# Averaged files are written into the local 'data' folder next to this script
# (created if it does not already exist), so downstream analysis reads locally
# instead of from the beamline storage on /gdata.
out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

# Percentile cutoff for the cross-correlation outlier removal.
percentile = 5

# --- HDF field locations (same layout in the raw and averaged files) ---
G2_PATH     = '/xpcs/multitau/normalized_g2'
G2_ERR_PATH = '/xpcs/multitau/normalized_g2_err'
SAXS_PATH   = '/xpcs/temporal_mean/scattering_1d'
START_TIME  = '/entry/start_time'
FILE_LIST   = '/xpcs/average/file_list'
# Scalar ion-chamber monitors, populated from the raw cluster results.  These
# are averaged per range (separately from the array fields) and written back so
# saxpcs.py can derive an absolute cross-section coefficient per averaged file.
INCIDENT_PATH    = '/entry/instrument/incident_beam/incident_beam_intensity'
TRANSMITTED_PATH = '/entry/instrument/incident_beam/transmitted_beam_intensity'

_frame_re = re.compile(r'_(\d+)(?:_[A-Za-z0-9]+)?_results\.hdf$')
_suffix_re = re.compile(r'_\d+(?:_([A-Za-z0-9]+))?_results\.hdf$')
# A start_time dated 2023 or later is a reprocessing date, not an acquisition
# time, so it is rejected rather than used.
STALE_TIME_RE = re.compile(r'^20(2[3-9]|[3-9]\d)-')


def extract_frame(fname):
    """Return the integer frame number embedded in a result file name."""
    m = _frame_re.search(os.path.basename(fname))
    if m is None:
        return -1
    return int(m.group(1))


def get_start_time(fname):
    """Acquisition time of a result file as 'YYYY-MM-DD HH:MM:SS', or None.

    This is the timestamp the instrument recorded in 2022.  The 2025 NeXus
    rewrite overwrote /entry/start_time with its own run date;
    restore_dx_metadata.py put the recorded value back from the Data Exchange
    archive, so the field can be trusted again and the acquisition time no
    longer has to be recovered from a saved directory listing.
    """
    with h5py.File(fname, 'r') as hf:
        raw = np.asarray(hf[START_TIME][()]).ravel()[0]
    # HDF text comes back as bytes from some writers and as str from others.
    if isinstance(raw, bytes):
        raw = raw.decode('utf-8')
    stamp = str(raw).strip()
    if STALE_TIME_RE.match(stamp):
        return None
    return stamp


def file_suffix(fname):
    """The boost_corr output tag in a result file name, or '' if it has none."""
    m = _suffix_re.search(os.path.basename(fname))
    if m is None:
        return ''
    tag = m.group(1)
    # The tag group is optional in the pattern, so it is None when the name
    # ends in '_00001_results.hdf' with nothing between the two.
    if tag is None:
        return ''
    return tag


def group_files(group):
    """Raw result files of a group, in frame order, one generation only.

    Reductions of the same group can differ by a boost_corr output tag, and
    they live side by side in ``prefix``.  Selecting on GROUP_SUFFIX keeps the
    two apart; matching on ``_suffix_re`` also drops earlier averages, whose
    names carry a frame RANGE where a raw file carries a single frame number.
    """
    want = GROUP_SUFFIX.get(group, '')
    keep = []
    for f in sorted(glob.glob(os.path.join(prefix, f'{group}*_results.hdf'))):
        name = os.path.basename(f)
        if 'Average' in name:
            continue                      # an averaged file from an earlier run
        if _frame_re.search(name) is None:
            continue                      # not a raw per-acquisition result
        if file_suffix(f) != want:
            continue                      # the other reduction generation
        keep.append(f)
    assert keep, f'no {want or "untagged"} dataset in {prefix} for group {group}'
    return keep


def read_range_data(section_files):
    """Read g2, g2_err, saxs_1d and the two beam monitors for every file.

    Returns
    -------
    data_dict : dict of stacked arrays (first axis = file index) for the fields
        that go through outlier removal (g2, g2_err, saxs_1d).
    read_files : list of the files read successfully, in stack order.
    incident : 1d array of the incident (upstream) beam intensity per file.
    transmitted : 1d array of the transmitted (downstream) beam intensity per file.

    ``incident``/``transmitted`` are kept OUT of ``data_dict`` on purpose: they
    are scalars per file and must not be fed to ``outlier_removal`` (which
    cross-correlates every key it receives).  They stay aligned with
    ``read_files`` so the g2-based outlier mask can be applied to them too.
    """
    data_dict = {'g2': [], 'g2_err': [], 'saxs_1d': []}
    incident, transmitted = [], []
    read_files = []
    for f in section_files:
        try:
            with h5py.File(f, 'r') as hf:
                g2 = hf[G2_PATH][()]
                g2_err = hf[G2_ERR_PATH][()]
                saxs = np.asarray(hf[SAXS_PATH][()])
                # collapse a (1, N)/segmented SAXS to a single 1d curve
                if saxs.ndim == 2:
                    saxs = np.nanmean(saxs, axis=0)
                inc = float(np.asarray(hf[INCIDENT_PATH][()]).ravel()[0])
                trn = float(np.asarray(hf[TRANSMITTED_PATH][()]).ravel()[0])
        except Exception as exc:
            print(f'  failed to read {os.path.basename(f)}: {exc}')
            continue
        data_dict['g2'].append(g2)
        data_dict['g2_err'].append(g2_err)
        data_dict['saxs_1d'].append(saxs)
        incident.append(inc)
        transmitted.append(trn)
        read_files.append(f)

    for k in data_dict:
        data_dict[k] = np.array(data_dict[k])
    return data_dict, read_files, np.array(incident), np.array(transmitted)


def _write_scalar(hf, path, value):
    """Overwrite an existing scalar dataset, keeping its stored shape/dtype."""
    ds = hf[path]
    ds[...] = np.asarray(value, dtype=ds.dtype).reshape(ds.shape)


def save_average(template, out_path, avg_dict, start_time, included_files,
                 incident_avg, transmitted_avg):
    """Copy ``template`` to ``out_path`` and overwrite it with the averaged data."""
    shutil.copyfile(template, out_path)

    with h5py.File(out_path, 'r+') as hf:
        hf[G2_PATH][...]     = avg_dict['g2']
        hf[G2_ERR_PATH][...] = avg_dict['g2_err']
        # write back with the file's own SAXS shape (e.g. (1, N))
        hf[SAXS_PATH][...] = avg_dict['saxs_1d'].reshape(hf[SAXS_PATH].shape)

        # range-averaged ion-chamber monitors (upstream / downstream)
        _write_scalar(hf, INCIDENT_PATH, incident_avg)
        _write_scalar(hf, TRANSMITTED_PATH, transmitted_avg)

        if start_time is not None:
            hf[START_TIME][()] = start_time

        if FILE_LIST in hf:
            del hf[FILE_LIST]
        basenames = [os.path.basename(f) for f in included_files]
        hf.create_dataset(
            FILE_LIST,
            data=np.array(basenames, dtype=object),
            dtype=h5py.string_dtype(encoding='utf-8'),
        )


def process_group(group, file_ranges):
    """Average every frame range for a single group and write the output files."""
    flist_all = group_files(group)
    frame_numbers = np.array([extract_frame(f) for f in flist_all])

    print(f'\n=== {group}: {len(flist_all)} raw files in {prefix} ===')
    for start, end in file_ranges:
        sel = (frame_numbers >= start) & (frame_numbers <= end)
        section_files = [f for f, m in zip(flist_all, sel) if m]
        assert section_files, f'no files found for {group} range ({start}, {end})'

        label = f'{group}_{start:05d}_{end:05d}'
        print(f'\nrange {start}-{end}: {len(section_files)} files')

        data_dict, read_files, incident, transmitted = read_range_data(section_files)
        assert read_files, f'no readable files for {group} range ({start}, {end})'

        mask = outlier_removal(data_dict, label=label, percentile=percentile)
        print(f'  kept {np.sum(mask)} of {len(mask)} files after outlier removal')

        # transient low-q scatterers, invisible to the shape-based cut above
        spike = spike_removal(data_dict['saxs_1d'], static_q(read_files[0]), mask)
        n_spike = int(np.sum(mask & ~spike))
        if n_spike:
            hit = sorted(extract_frame(f) for f, m in zip(read_files, mask & ~spike) if m)
            print(f'  dropped {n_spike} more as low-q spikes {hit} '
                  f'-> {np.sum(spike)} remain')
        mask = spike

        # hand-listed frames, dropped on top of the automatic cuts
        drop = MANUAL_EXCLUDE.get(group, {}).get((start, end), [])
        if drop:
            hit = np.array([extract_frame(f) in drop for f in read_files])
            n_hit = int(np.sum(mask & hit))          # only those still standing
            mask = mask & ~hit
            print(f'  dropped {n_hit} more from the hand list {sorted(drop)} '
                  f'-> {np.sum(mask)} remain')

        included_files = [f for f, keep in zip(read_files, mask) if keep]

        avg_dict = average_datasets(data_dict=data_dict, mask=mask)

        # Average the two beam monitors over the SAME surviving files (the mask
        # from the g2/saxs outlier removal), so the averaged file's flux and
        # transmission are consistent with its averaged g2/saxs curves.
        incident_avg = float(np.nanmean(incident[mask]))
        transmitted_avg = float(np.nanmean(transmitted[mask]))

        # start time comes from the first (lowest frame number) file in the range
        start_time = get_start_time(section_files[0])
        if start_time is None:
            print(f'  WARNING: no acquisition time in '
                  f'{os.path.basename(section_files[0])}')

        template = section_files[0]
        # replace the single frame number with the range, then prepend 'Average_'
        # 'B0147_..._00801_results.hdf' -> 'Average_B0147_..._00801_00950_results.hdf'.
        # The tail is whatever follows the frame number, so an output tag such
        # as '_BadpixRm' is carried through to the averaged file name.
        name = os.path.basename(template)
        tail = '_results.hdf'
        tag = file_suffix(template)
        if tag:
            tail = f'_{tag}_results.hdf'
        stem = name[:name.index(f'_{extract_frame(template):05d}{tail}')]
        out_name = f'Average_{stem}_{start:05d}_{end:05d}{tail}'
        out_path = os.path.join(out_dir, out_name)

        save_average(template, out_path, avg_dict, start_time, included_files,
                     incident_avg, transmitted_avg)
        print(f'  saved -> {out_name}  (start_time={start_time}, '
              f'incident={incident_avg:.6g}, transmitted={transmitted_avg:.6g})')


def write_temperature_trace(groups, out_path):
    """One row per acquisition: dataset, frame, elapsed seconds, temperature.

    Elapsed time is measured from the first acquisition of the sequence.
    """
    rows = []
    for group in groups:
        for f in group_files(group):
            key = os.path.basename(f).replace('_results.hdf', '')
            stamp = get_start_time(f)
            if stamp is None:
                print(f'  WARNING: no acquisition time in {key}')
                continue
            with h5py.File(f, 'r') as hf:
                T = float(np.asarray(hf[TEMP_PATH][()]).ravel()[0])
            rows.append([key, extract_frame(f),
                         datetime.strptime(stamp, '%Y-%m-%d %H:%M:%S'), T])
    rows.sort(key=lambda r: r[2])
    t0 = rows[0][2]
    with open(out_path, 'w') as fh:
        fh.write('dataset,frame,elapsed_s,qnw1_temperature_C\n')
        for key, frame, t, T in rows:
            fh.write(f'{key},{frame},{(t - t0).total_seconds():.0f},{T:.4f}\n')
    print(f'  {len(rows)} acquisitions, {(rows[-1][2] - t0).total_seconds() / 3600:.2f} h, '
          f'T = {min(r[3] for r in rows):.2f}-{max(r[3] for r in rows):.2f} C '
          f'-> {os.path.basename(out_path)}')


def main(only=None):
    """Re-reduce every group, or only the ones named on the command line.

    ``python average_ranges.py`` rebuilds all of data/.  ``python
    average_ranges.py B0083`` rebuilds just that group's averages and leaves
    every other file and the temperature trace untouched -- useful when only
    one group's frame selection has changed.
    """
    os.makedirs(out_dir, exist_ok=True)

    print(f'output dir: {out_dir}')
    for group, file_ranges in FILE_RANGES.items():
        if only and group not in only:
            continue
        process_group(group, file_ranges)
    if only:
        print('\nsubset run: temperature trace left as it was')
        return

    print(f'\n=== temperature trace of the cycling sequence ===')
    write_temperature_trace(TRACE_GROUPS, os.path.join(out_dir, TRACE_CSV))


if __name__ == '__main__':
    main(only=set(sys.argv[1:]) or None)
