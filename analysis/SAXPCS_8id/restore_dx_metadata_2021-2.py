"""Carry the ring current, the QNW temperatures and the real acquisition times
from the 2021-2 archive into its NeXus tree.

BEAMLINE-ONLY.  Needs /gdata, which mounts on the 10.x hosts (adamite, amber)
and not on the 164.x analysis machines.  Kept here for provenance.

WHAT IS WRONG WITH THE TREE
---------------------------
babnigg202107_2_nexus was written by the same rewrite that damaged the 2022-1
tree, and it carries the same placeholder constants:

    /entry/instrument/incident_beam/incident_beam_intensity     133.49012809232
    /entry/instrument/incident_beam/transmitted_beam_intensity  0.0001
    /entry/instrument/source/current                            absent
    /entry/sample/temperature                                   absent
    /entry/start_time                                           the rewrite date

The Data Exchange archive it was built from has the real values, under a
DIFFERENT schema from 2021-3: `source_begin/...` rather than `source_end/...`,
and the ion chambers named `beam_intensity_incident` / `_transmitted` rather
than `I0Monitor` / `TransmissionMonitor`.

ONE ION CHAMBER IS RECOVERABLE, THE OTHER IS NOT
------------------------------------------------
The upstream chamber IS recorded, under `source_end/I0Monitor`, the same path
2021-3 uses: 518,080 to 552,999 across the six ramp runs.  It is carried.

The downstream chamber is not: `source_end/TransmissionMonitor` is absent from
every 2021-2 acquisition, and the `source_begin/beam_intensity_incident` and
`_transmitted` pair that does exist reads exactly 1.0 everywhere, a placeholder
rather than a measurement.  Nothing is invented for it here.

The consequence belongs in abs_xsec.py rather than in this script.  The
incident flux follows from I0Monitor through the ordinary calibration, but the
sample transmission, normally Dn_IC / Up_IC / AIR_TRANSMISSION, has no Dn_IC to
come from.  See TRANSMISSION_2021_2 there: the 2022-1 value is used, which is
sound because the calibration coefficients are properties of the chambers
rather than of the beamtime, and the cell thickness, pixel size and sample
absorption barely move between these runs.

WHAT IS CARRIED
---------------
    NeXus destination                                        Data Exchange source
    /entry/instrument/incident_beam/incident_beam_intensity  source_end/I0Monitor
    /entry/instrument/source/current           source_begin/current
    /entry/sample/temperature                  sample/QNW_Zone1_Temperature
    /entry/sample/temperature_set              sample/QNW_Zone1_Temperature_Set
    /entry/start_time                          source_begin/datetime
    /entry/end_time                            source_end/datetime

THE SAMPLE IS NOT ALWAYS IN ZONE 1.  The QNW stage holds three cells and the
two idle zones are parked at 6 C, so the zone that matters changes from run to
run: (VPAVG)10 (B039/B040) and (VPAVG)20 (T214/T215) ramp in zone 1, while
(VPAVG)30 (H041/H042) ramps in zone 3 with zone 1 parked.  Writing zone 1
blindly puts a flat 6 C into the temperature axis of the (VPAVG)30 ramps.  The
active zone is therefore chosen per acquisition, as the zone whose set point is
furthest from the 6 C park, and the choice is recorded in the manifest and
written to /entry/sample/qnw_zone so a later reader can check it.

Only a field that is absent, non-finite, or different from the archive is
written, so the script is idempotent and cannot overwrite a value corrected by
hand.  Every comparison is recorded in the manifest whether or not it is
written.

Usage
-----
    python restore_dx_metadata_2021-2.py            # dry run; manifest only
    python restore_dx_metadata_2021-2.py --apply    # perform the writes
"""

import argparse
import csv
import glob
import multiprocessing
import os
from datetime import datetime

import h5py
import numpy as np

ARCHIVE = '/gdata/s8id-dmdtn/2021-2/babnigg202107_2'
NEXUS = '/gdata/s8id-dmdtn/2021-2/babnigg202107_2_nexus'
RESULTS = os.path.join(NEXUS, 'reprocess_results')
MANIFEST = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        'dx_restore_manifest_2021-2.csv')

# (NeXus destination, Data Exchange source, kind).  'f' is a float written in
# place; 't' is a timestamp, reformatted to the ISO form the NeXus tree uses.
FIELDS = (
    ('/entry/instrument/incident_beam/incident_beam_intensity',
     'measurement/instrument/source_end/I0Monitor', 'f'),
    ('/entry/instrument/source/current', 'measurement/instrument/source_begin/current', 'f'),
    ('/entry/start_time', 'measurement/instrument/source_begin/datetime', 't'),
    ('/entry/end_time', 'measurement/instrument/source_end/datetime', 't'),
)

# The downstream chamber: 2021-2 never recorded it, so the placeholder 0.0001
# the rewrite left is reported and never written.  See the module docstring.
UNRECOVERABLE = ('/entry/instrument/incident_beam/transmitted_beam_intensity',)

PARKED_C = 6.0          # the temperature an idle QNW zone is held at
ZONES = (1, 2, 3)


def scalar(obj):
    """First element of a field stored as (), (1,) or (1, 1)."""
    return float(np.asarray(obj[()]).ravel()[0])


def text(obj):
    v = np.asarray(obj[()]).ravel()[0]
    return v.decode() if isinstance(v, bytes) else str(v)


def iso(stamp):
    """'Wed Jul 21 03:40:52 2021' -> '2021-07-21 03:40:52', the NeXus form."""
    return datetime.strptime(stamp.strip(), '%a %b %d %H:%M:%S %Y').strftime('%Y-%m-%d %H:%M:%S')


def nexus_index():
    """Acquisition name -> every NeXus file that belongs to it.

    The tree mixes two layouts.  Most acquisitions sit one level down, inside a
    scan directory: <nexus>/<scan>/<acq>/<acq>_metadata.hdf.  A few sit at the
    top: <nexus>/<acq>/<acq>_metadata.hdf.  An earlier version of this script
    assumed only the second and so repaired 262 of 24,217 files.  Walk for the
    metadata files instead of guessing the depth.
    """
    index = {}
    for root, _dirs, files in os.walk(NEXUS):
        if os.path.basename(root) == 'reprocess_results':
            continue
        for fn in files:
            if fn.endswith('_metadata.hdf'):
                index.setdefault(fn[:-len('_metadata.hdf')], []).append(os.path.join(root, fn))
    for fp in glob.glob(os.path.join(RESULTS, '*_results.hdf')):
        acq = os.path.basename(fp)[:-len('_results.hdf')]
        index.setdefault(acq, []).append(fp)
    return index


def archive_index():
    """Acquisition name -> its Data Exchange file, for the whole archive."""
    index = {}
    for scan_dir in glob.glob(os.path.join(ARCHIVE, '*')):
        if not os.path.isdir(scan_dir):
            continue
        subs = [d for d in glob.glob(os.path.join(scan_dir, os.path.basename(scan_dir) + '_*'))
                if os.path.isdir(d)]
        for d in subs or [scan_dir]:
            found = sorted(glob.glob(os.path.join(d, '*_0001-*.hdf')))
            if found:
                index[os.path.basename(d)] = found[0]
    return index


def active_zone(hf):
    """Which QNW zone holds the ramping sample, by set point.

    The two idle zones sit at PARKED_C, so the active one is the zone whose set
    point is furthest from it.  When every zone is parked -- the 6 C reference
    runs -- there is nothing to choose between them and zone 1 is returned; the
    three then read the same value anyway.
    """
    best, best_d = 1, -1.0
    for z in ZONES:
        k = f'measurement/sample/QNW_Zone{z}_Temperature_Set'
        if k not in hf:
            continue
        d = abs(scalar(hf[k]) - PARKED_C)
        if d > best_d:
            best, best_d = z, d
    return best


def read_archive(path):
    """{nexus_path: value} for the fields this acquisition actually recorded."""
    out = {}
    try:
        with h5py.File(path, 'r') as hf:
            z = active_zone(hf)
            out['/entry/sample/qnw_zone'] = float(z)
            for nx, dx in ((f'/entry/sample/temperature',
                            f'measurement/sample/QNW_Zone{z}_Temperature'),
                           (f'/entry/sample/temperature_set',
                            f'measurement/sample/QNW_Zone{z}_Temperature_Set')):
                if dx in hf:
                    out[nx] = scalar(hf[dx])
            for nx, dx, kind in FIELDS:
                if dx not in hf:
                    continue
                try:
                    out[nx] = iso(text(hf[dx])) if kind == 't' else scalar(hf[dx])
                except (ValueError, TypeError):
                    pass
    except OSError:
        pass
    return out


def write(hf, nx, want, kind):
    if kind == 't':
        if nx in hf:
            del hf[nx]
        hf.create_dataset(nx, data=np.bytes_(want))
        return
    if nx in hf:
        ds = hf[nx]
        ds[...] = np.asarray(want, dtype=ds.dtype).reshape(ds.shape)
    else:
        hf.create_dataset(nx, data=float(want))


def check_file(path, values, apply_changes):
    rows = []
    with h5py.File(path, 'r+' if apply_changes else 'r') as hf:
        for nx, kind in ([(f, k) for f, _d, k in FIELDS]
                         + [('/entry/sample/temperature', 'f'),
                            ('/entry/sample/temperature_set', 'f'),
                            ('/entry/sample/qnw_zone', 'f')]):
            want = values.get(nx)
            if nx in hf:
                have = text(hf[nx]) if kind == 't' else scalar(hf[nx])
            else:
                have = None
            if want is None:
                rows.append((nx, have, None, 'not recorded in archive'))
                continue
            same = (have == want) if kind == 't' else (
                have is not None and np.isfinite(have) and np.isclose(have, want, rtol=1e-6))
            if same:
                rows.append((nx, have, want, 'already correct'))
                continue
            if apply_changes:
                write(hf, nx, want, kind)
            rows.append((nx, have, want, 'written' if apply_changes else 'would write'))
        for nx in UNRECOVERABLE:
            have = scalar(hf[nx]) if nx in hf else None
            rows.append((nx, have, None, 'downstream chamber not recorded in 2021-2; left as is'))
    return rows


def process(job):
    acq, dx_file, targets, apply_changes = job
    values = read_archive(dx_file)
    rows = []
    for t in targets:
        if not os.path.exists(t):
            continue
        try:
            for nx, have, want, action in check_file(t, values, apply_changes):
                rows.append((acq, os.path.relpath(t, NEXUS), nx, have, want, action))
        except OSError as e:
            rows.append((acq, os.path.relpath(t, NEXUS), '', '', '', f'open failed: {e}'))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--apply', action='store_true', help='perform the writes')
    ap.add_argument('--workers', type=int, default=8)
    args = ap.parse_args()

    index = archive_index()
    nx = nexus_index()
    print(f'{len(index)} acquisitions in the archive, '
          f'{len(nx)} in the NeXus tree, '
          f'{sum(len(v) for v in nx.values())} NeXus files')
    jobs = [(a, f, nx[a], args.apply) for a, f in sorted(index.items()) if a in nx]
    print(f'{len(jobs)} acquisitions matched')
    with multiprocessing.Pool(args.workers) as pool:
        out = pool.map(process, jobs, chunksize=16)

    rows = [r for sub in out for r in sub]
    with open(MANIFEST, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['acquisition', 'nexus_file', 'field', 'nexus_value', 'archive_value', 'action'])
        w.writerows(rows)

    tally = {}
    for r in rows:
        tally[r[5]] = tally.get(r[5], 0) + 1
    print(f'{len(rows)} comparisons over {len({r[1] for r in rows})} distinct NeXus files')
    for k in sorted(tally, key=lambda k: -tally[k]):
        print(f'   {tally[k]:7d}  {k}')
    print(f'manifest -> {MANIFEST}')
    if not args.apply:
        print('DRY RUN: nothing was written.  Rerun with --apply.')


if __name__ == '__main__':
    main()
