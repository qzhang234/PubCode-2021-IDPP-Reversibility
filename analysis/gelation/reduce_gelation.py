"""Reduce the 8-ID-I temperature ramps behind the gelation figures into data/.

READS /gdata, SO IT RUNS ON AMBER (10.54.116.10; see CLAUDE.md).  Its input is
the per-acquisition correlation results, one HDF file per 2 s acquisition,
all correlated from the raw frames on /gdata, against one q-map per beamtime:
babnigg202107_2_nexus_Sq270_Dq27_Sphi8_Dphi1_q0p0032_Lin_2_BadPixRm.hdf for
2021-2 and babnigg202203_nexus_Sq270_Dq27_Sphi8_Dphi1_Q0p0032_Lin.hdf for 2022-1:

    /gdata/s8id-dmdtn/2021-2/babnigg202107_2_nexus/reprocess_results
    /gdata/s8id-dmdtn/2022-1/babnigg202203_nexus/reprocess_results
        the October 2026 correlation of every acquisition of the ramps and of
        the 2022-1 buffer, all 100,000 frames;
    /home/8-id-i/2024-3/2024_1228_qz_llps_analysis/2026_10_gdata_reanalysis/
        2021-2_D029, 2021-2_B039_part1, 2021-2_B039_part2
        what the October run left out -- the 2021-2 water blank, and B039
        correlated over frames 0-49,999 and 50,000-99,999 -- correlated on
        Amber by run_gelation_correlations.sh with the same settings.

Every figure script reads only what this writes into data/, so `make figures`
never runs it.

THE OVERLEAF FIGURES were made from an earlier correlation of the same frames,
the August 2025 reanalysis under /home/8-id-i/2024-3/2024_1228_qz_llps_analysis/
2025_08_analysis/, which for 2021-2 used the q-map without the two cleared
pixels.  The version of this script in commit 3e5ce76 read those files and
reproduced the original section averages exactly in g2 and g2_err for all
eleven datasets, which is how the port was checked.

THE DATA
--------
Three chain lengths, each heated from 20 C at 0.1 C/min until g2 stopped
decaying, then cooled at the same rate:

    (VPAVG)10   B039 up, B040 down    2021-2, blank D029 (water, 25 C)
    (VPAVG)20   T214 up, T215 down    2021-2, blank D029
    (VPAVG)30   E0110 up, E0111 down  2022-1, buffer D0138 (34 C)

THE BINNING is the one the Overleaf figures used, by acquisition count: a run
is cut into N consecutive sections of equal length (np.array_split rounding),
optionally after skipping its first few acquisitions.  On a constant-rate ramp
equal counts are equal temperature steps.

    key                 run    sections  skip   used by
    B039_n10            B039   10        0      Figure 3, Figure S2, S12
    B040_n10            B040   10        0      Figure S2
    T214_n20            T214   20        0      Figure 4, Figure S3
    T215_n20            T215   20        0      Figure S3
    E0110_n10           E0110  10        0      Figure 2, Figure S4
    E0111_n10           E0111  10        0      Figure S4
    B039_n20_skip70     B039   20        70     Figure 4
    E0110_n20           E0110  20        0      Figure 4
    B039_part1_n10      B039   10        0      Figure S12 (frames 1-50,000)
    B039_part2_n10      B039   10        0      Figure S12 (frames 50,001-100,000)
    D029, D0138         blank / buffer, one section each

OUTLIERS.  Within each section, the cross-correlation cut of
common.utils.outlier_removal: every acquisition's g2, g2_err and log I(Q) is
compared with every other's, and the least similar 5 % on any of the three is
dropped.  This is the only cut the Overleaf figures received.  The routine in
common.utils is the corrected one, which works on a copy; the original edited
the arrays in place and could write 1.0 into I(Q) wherever a value was missing
or non-positive.  On these runs the two give identical keep-masks and identical
averages (checked against the original pickles, Oct 2026).

PER ACQUISITION, as the original pyXpcsViewer reader did:
  * I(Q) is the azimuthal average of /xpcs/temporal_mean/scattering_1d over the
    eight phi sectors (photons per pixel per frame).
  * g2_err entries at or below 1e-6 are replaced by the mean of the rest of
    that Q column (pyXpcsViewer's correct_g2_err), which keeps a fit from
    treating a numerically zero error bar as infinitely precise.
  * temperature is the RTD of the stage zone holding the cell.

MONITORS AND TIMES come from each result file: the ion-chamber readings and
the start time that restore_dx_metadata(_2021-2).py recovered from the Data
Exchange archive into the NeXus metadata in October 2026, before these
correlations were run.  2021-2 never recorded the downstream chamber, so its
transmitted value is a placeholder and is stored as NaN.  Section monitors are
means over the acquisitions that survived the outlier cut, the same
acquisitions the curves are averaged over.

FIGURE S11 needs one acquisition split into the twelve 8,333-frame subsets the
correlator stores in /xpcs/temporal_mean/scattering_1d_segments.  The stored
subsets are each divided by their own mean detector intensity: across the
twelve, their level runs opposite to that intensity (correlation 0.8-0.9 over
the ramp).  Each is multiplied back by its mean intensity, read from
/xpcs/spatial_mean/intensity_vs_time, and one constant per acquisition puts the
twelve on the scale of the time-averaged I(Q).  The script prints how constant
that factor is across Q, which tests the description.

    python reduce_gelation.py              # everything
    python reduce_gelation.py D029 segments  # just those
"""

import glob
import os
import sys
from datetime import datetime
from multiprocessing import Pool

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
from common.utils import average_datasets, outlier_removal       # noqa: E402

GDATA = '/gdata/s8id-dmdtn'
REANALYSIS = '/home/8-id-i/2024-3/2024_1228_qz_llps_analysis/2026_10_gdata_reanalysis'
SOURCES = {
    'ramps21': f'{GDATA}/2021-2/babnigg202107_2_nexus/reprocess_results',
    'ramps22': f'{GDATA}/2022-1/babnigg202203_nexus/reprocess_results',
    'D029':    f'{REANALYSIS}/2021-2_D029',
    'part1':   f'{REANALYSIS}/2021-2_B039_part1',
    'part2':   f'{REANALYSIS}/2021-2_B039_part2',
}
OUT_DIR = os.path.join(HERE, 'data')

INCIDENT = '/entry/instrument/incident_beam/incident_beam_intensity'
TRANSMITTED = '/entry/instrument/incident_beam/transmitted_beam_intensity'
START = '/entry/start_time'

# key -> (beamtime, source, run, sections, skipped acquisitions, QNW zone)
DATASETS = {
    'B039_n10':        ('2021-2', 'ramps21', 'B039',  10, 0,  1),
    'B040_n10':        ('2021-2', 'ramps21', 'B040',  10, 0,  1),
    'T214_n20':        ('2021-2', 'ramps21', 'T214',  20, 0,  1),
    'T215_n20':        ('2021-2', 'ramps21', 'T215',  20, 0,  1),
    'D029':            ('2021-2', 'D029',    'D029',   1, 0,  2),
    'E0110_n10':       ('2022-1', 'ramps22', 'E0110', 10, 0,  2),
    'E0111_n10':       ('2022-1', 'ramps22', 'E0111', 10, 0,  2),
    'D0138':           ('2022-1', 'ramps22', 'D0138',  1, 0,  2),
    'B039_n20_skip70': ('2021-2', 'ramps21', 'B039',  20, 70, 1),
    'E0110_n20':       ('2022-1', 'ramps22', 'E0110', 20, 0,  2),
    'B039_part1_n10':  ('2021-2', 'part1',   'B039',  10, 0,  1),
    'B039_part2_n10':  ('2021-2', 'part2',   'B039',  10, 0,  1),
}

# Figure S11: the first and the last acquisition of the (VPAVG)10 heating ramp
SEGMENT_FILES = ['B039_VPAVG_H02_500C10p_att00_Rq0_00001',
                 'B039_VPAVG_H02_500C10p_att00_Rq0_02780']

PERCENTILE = 5
NUM_CORES = 16

G2, G2_ERR = '/xpcs/multitau/normalized_g2', '/xpcs/multitau/normalized_g2_err'
DELAY, FRAME_TIME = '/xpcs/multitau/delay_list', '/entry/instrument/detector_1/frame_time'
SAXS, SEGMENTS = '/xpcs/temporal_mean/scattering_1d', '/xpcs/temporal_mean/scattering_1d_segments'
INT_T = '/xpcs/spatial_mean/intensity_vs_time'
S_MAP, S_Q, S_PHI = ('/xpcs/qmap/static_index_mapping', '/xpcs/qmap/static_v_list_dim0',
                     '/xpcs/qmap/static_v_list_dim1')
D_Q = '/xpcs/qmap/dynamic_v_list_dim0'


def split(arr, n):
    """n consecutive pieces of near-equal length, the larger ones first."""
    k, m = divmod(len(arr), n)
    return [arr[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n)]


def phi_average(compressed, hf):
    """(rows, n_bins) -> (rows, n_q): mean over the phi sectors present at each Q."""
    idx = hf[S_MAP][()]
    nq, nphi = hf[S_Q].shape[0], hf[S_PHI].shape[0]
    compressed = np.asarray(compressed, dtype=np.float64).reshape(-1, idx.size)
    full = np.full((compressed.shape[0], nq * nphi), np.nan)
    full[:, idx] = compressed
    return np.nanmean(full.reshape(-1, nq, nphi), axis=2)


def correct_g2_err(g2_err, threshold=1e-6):
    """pyXpcsViewer's fix for numerically zero error bars, column by column."""
    out = np.copy(g2_err)
    for n in range(g2_err.shape[1]):
        ok = g2_err[:, n] > threshold
        out[~ok, n] = np.mean(g2_err[ok, n]) if ok.any() else threshold
    return out


def scalar(hf, path):
    raw = np.ravel(hf[path][()])[0]
    return raw.decode() if isinstance(raw, bytes) else raw


def seconds(stamp):
    return datetime.strptime(str(stamp).strip()[:19], '%Y-%m-%d %H:%M:%S').timestamp()


def read_one(args):
    """g2, g2_err, I(Q), temperature, incident and transmitted monitor and start
    time (s) of one acquisition, or None if the file is unreadable."""
    fname, zone, beam = args
    try:
        with h5py.File(fname, 'r') as hf:
            trn = float(scalar(hf, TRANSMITTED)) if beam == '2022-1' else np.nan
            return (hf[G2][()], correct_g2_err(hf[G2_ERR][()]),
                    phi_average(hf[SAXS][()], hf)[0],
                    float(scalar(hf, f'/entry/sample/qnw{zone}_temperature')),
                    float(scalar(hf, INCIDENT)), trn, seconds(scalar(hf, START)))
    except (OSError, KeyError) as err:
        print(f'  skipped unreadable {os.path.basename(fname)}: {err}')
        return None


def acquisition(fname):
    return os.path.basename(fname).replace('_results.hdf', '')


def reduce(key):
    beam, source, run, nsec, skip, zone = DATASETS[key]
    src = SOURCES[source]
    files = sorted(glob.glob(os.path.join(src, f'{run}*_results.hdf')))[skip:]
    assert files, f'no {run} files under {src}'
    with Pool(NUM_CORES) as pool:
        per_file = pool.map(read_one, [(f, zone, beam) for f in files])

    acq = [acquisition(f) for f in files]
    meta = np.array([[np.nan] * 4 if r is None else r[3:] for r in per_file])
    temp, inc, trn, t_abs = meta.T
    kept = np.zeros(len(files), dtype=bool)
    section = np.zeros(len(files), dtype=int)

    out = {k: [] for k in ('g2', 'g2_err', 'saxs_1d', 'T_mean', 'T_min', 'T_max',
                           'n_files', 'n_kept', 'incident', 'transmitted')}
    print(f'{key}: {run} from {src}, {len(files)} acquisitions after skipping {skip}')
    print('  sec  T-min   T-max   T-mean   kept')
    for s, idx in enumerate(split(np.arange(len(files)), nsec)):
        section[idx] = s
        good = [i for i in idx if per_file[i] is not None]
        data = {'g2': np.array([per_file[i][0] for i in good]),
                'g2_err': np.array([per_file[i][1] for i in good]),
                'saxs_1d': np.array([per_file[i][2] for i in good])}
        mask = outlier_removal(data, label=f'{key}_{s:02d}', percentile=PERCENTILE)
        avg = average_datasets(data_dict=data, mask=mask)
        keep = np.array(good)[mask]
        kept[keep] = True
        for k in ('g2', 'g2_err', 'saxs_1d'):
            out[k].append(avg[k])
        T = temp[idx]
        out['T_mean'].append(np.nanmean(T)); out['T_min'].append(np.nanmin(T))
        out['T_max'].append(np.nanmax(T))
        out['n_files'].append(len(idx)); out['n_kept'].append(len(keep))
        out['incident'].append(np.mean(inc[keep])); out['transmitted'].append(np.mean(trn[keep]))
        print(f'  {s:3d}  {np.nanmin(T):6.2f}  {np.nanmax(T):6.2f}  {np.nanmean(T):7.3f}  '
              f'{len(keep)}/{len(idx)}')

    with h5py.File(files[0], 'r') as hf:
        t0 = float(np.ravel(hf[FRAME_TIME][()])[0])
        tau = np.ravel(hf[DELAY][()]) * t0
        q_static, q_dynamic = hf[S_Q][()], hf[D_Q][()]

    path = os.path.join(OUT_DIR, f'{key}.h5')
    with h5py.File(path, 'w') as h:
        h.attrs.update(run=run, beamtime=beam, n_sections=nsec,
                       skipped=skip, qnw_zone=zone, source=src,
                       outlier_cut=f'cross-correlation, lowest {PERCENTILE}% of g2, g2_err, log I(Q)',
                       written_by='reduce_gelation.py')
        h['tau'] = tau
        h['q_static'] = q_static
        h['q_dynamic'] = q_dynamic
        h['g2'] = np.array(out['g2'], dtype=np.float32)
        h['g2_err'] = np.array(out['g2_err'], dtype=np.float32)
        h['saxs_1d'] = np.array(out['saxs_1d'])
        for k in ('T_mean', 'T_min', 'T_max', 'incident', 'transmitted'):
            h[f'section/{k}'] = np.array(out[k])
        for k in ('n_files', 'n_kept'):
            h[f'section/{k}'] = np.array(out[k], dtype=int)
        h['acquisition/name'] = np.array(acq, dtype='S')
        h['acquisition/temperature'] = temp
        h['acquisition/elapsed_s'] = t_abs - np.nanmin(t_abs)
        h['acquisition/section'] = section
        h['acquisition/kept'] = kept
        h['acquisition/incident'] = inc
        h['acquisition/transmitted'] = trn
        for k, v in (('tau', 's'), ('q_static', 'A^-1'), ('q_dynamic', 'A^-1'),
                     ('saxs_1d', 'photons/pixel/frame')):
            h[k].attrs['units'] = v
    print(f'  -> {os.path.relpath(path, HERE)}\n')


def segments():
    """Figure S11: the twelve frame subsets of two acquisitions, in photons/pixel/frame."""
    src = SOURCES['ramps21']
    path = os.path.join(OUT_DIR, 'B039_segments.h5')
    with h5py.File(path, 'w') as h:
        h.attrs.update(source=src, written_by='reduce_gelation.py')
        for name in SEGMENT_FILES:
            with h5py.File(os.path.join(src, name + '_results.hdf'), 'r') as hf:
                stored = phi_average(hf[SEGMENTS][()], hf)
                mean = phi_average(hf[SAXS][()], hf)[0]
                intt = np.ravel(hf[INT_T][()][1])
                q = hf[S_Q][()]
                T = float(np.ravel(hf['/entry/sample/qnw1_temperature'][()])[0])
            n = stored.shape[0]
            level = np.array([seg.mean() for seg in np.array_split(intt, n)])
            raw = stored * level[:, None]
            factor = np.nanmean(raw, axis=0) / mean          # one number if the model holds
            c = np.nanmedian(factor)
            print(f'{name}: T = {T:.2f} C, {n} subsets; scale factor across Q: '
                  f'median {c:.4g}, relative spread {np.nanstd(factor) / c:.2%}')
            g = h.create_group(name)
            g.attrs['temperature'] = T
            g['q'] = q
            g['saxs_1d_subsets'] = raw / c
            g['saxs_1d'] = mean
            g['subset_mean_intensity'] = level
    print(f'  -> {os.path.relpath(path, HERE)}\n')


def main(only=None):
    """Reduce every dataset and the Figure S11 subsets, or only those named."""
    os.makedirs(OUT_DIR, exist_ok=True)
    for key in DATASETS:
        if not only or key in only:
            reduce(key)
    if not only or 'segments' in only:
        segments()


if __name__ == '__main__':
    main(set(sys.argv[1:]) or None)
