"""Restore into the NeXus tree the measured values the 2025 rewrite dropped.

BEAMLINE-ONLY.  This script has to see /gdata, which is mounted on the 10.x
hosts (amber) and not on the 164.x analysis machines, so it cannot run from a
clone of the repository.  It is kept here for provenance: it documents exactly
how the data the figures are built from came to be what they are.

WHY IT EXISTS
-------------
The 2022-1 beamtime was recorded in the Data Exchange schema.  In 2025 the
metadata were rewritten into ``babnigg202203_nexus`` so the current
``boost_corr_bin``, which reads NeXus only, could reprocess the 2022 frames.
The rewrite carried the detector data across intact but wrote constants in
place of four quantities that were actually measured:

    NeXus destination                                       value after rewrite
    /entry/instrument/incident_beam/incident_beam_intensity      133.49012809232
    /entry/instrument/incident_beam/transmitted_beam_intensity            0.0001
    /entry/instrument/incident_beam/ring_current                             0.0
    /entry/start_time, /entry/end_time                       the 2025 rewrite date

The first two are the upstream and downstream ion chambers.  abs_xsec.py turns
them into the incident flux and the sample transmission, so with the constants
in place every absolute cross section is wrong: the implied transmission is
8.6e-07 instead of 0.35 and the implied flux is negative.

The originals survive in the 2022 archive, which the NeXus rewrite never
touched, under measurement/instrument/source_begin and source_end.  This script
copies them back.

WHAT IT WRITES
--------------
Only the derived ``babnigg202203_nexus`` tree.  The 2022 archive is opened
read-only and is never modified.  Within the NeXus tree both levels are
written: the per-acquisition ``*_metadata.hdf``, so that a future reprocessing
inherits the measured values instead of the constants, and the
``reprocess_results/*_results.hdf`` that average_ranges.py reads.

A field is written only when it still holds the exact constant the 2025 rewrite
left.  A field already holding the archive value is skipped, so the script is
idempotent; anything else is left alone and reported, so a value somebody has
already corrected by hand cannot be clobbered.

Every restored value is written to a manifest CSV, which is committed to the
repository as the record of what this script changed.

Usage
-----
    python restore_dx_metadata.py              # dry run; writes the manifest only
    python restore_dx_metadata.py --apply      # perform the writes
    python restore_dx_metadata.py --limit 50   # first 50 acquisitions, for a smoke test
"""

import os
import re
import csv
import sys
import glob
import argparse
import multiprocessing

from datetime import datetime

import numpy as np
import h5py

ARCHIVE = '/gdata/s8id-dmdtn/2022-1/babnigg202203'
NEXUS   = '/gdata/s8id-dmdtn/2022-1/babnigg202203_nexus'
RESULTS = os.path.join(NEXUS, 'reprocess_results')

MANIFEST = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        'dx_restore_manifest.csv')

# Scalars: NeXus path, Data Exchange path, and the constant the 2025 rewrite
# left behind.  Only that exact constant is overwritten.
SCALARS = (
    ('/entry/instrument/incident_beam/incident_beam_intensity',
     'measurement/instrument/source_end/I0Monitor',           133.49012809232),
    ('/entry/instrument/incident_beam/transmitted_beam_intensity',
     'measurement/instrument/source_end/TransmissionMonitor',          0.0001),
    ('/entry/instrument/incident_beam/ring_current',
     'measurement/instrument/source_end/current',                         0.0),
)

# Timestamps: NeXus path, Data Exchange path.  The rewrite put its own run date
# here, so the guard is "the stored text is a 2025 date" rather than a constant.
TIMES = (
    ('/entry/start_time', 'measurement/instrument/source_begin/datetime'),
    ('/entry/end_time',   'measurement/instrument/source_end/datetime'),
)

DX_TIME_FORMAT  = '%a %b %d %H:%M:%S %Y'      # 'Fri Mar 04 23:31:26 2022'
OUT_TIME_FORMAT = '%Y-%m-%d %H:%M:%S'         # what nexus_read.py parses
STALE_TIME_RE   = re.compile(r'^20(2[3-9]|[3-9]\d)-')   # 2023 or later


def scalar(obj):
    """First element of a Data Exchange field, which may be (1,) or (1, 1)."""
    return float(np.asarray(obj[()]).ravel()[0])


def text(obj):
    """A Data Exchange string field as str, however it happens to be wrapped."""
    v = np.asarray(obj[()]).ravel()[0]
    return v.decode('utf-8') if isinstance(v, bytes) else str(v)


def archive_file(stem):
    """The 2022 Data Exchange metadata file of one acquisition, or None."""
    direct = os.path.join(ARCHIVE, stem, f'{stem}_0001-100000.hdf')
    if os.path.exists(direct):
        return direct
    found = sorted(glob.glob(os.path.join(ARCHIVE, stem, '*_0001-*.hdf')))
    return found[0] if found else None


def read_archive(path):
    """The measured values of one acquisition: {nexus_path: value}."""
    out = {}
    with h5py.File(path, 'r') as hf:
        for nx_path, dx_path, _const in SCALARS:
            if dx_path in hf:
                out[nx_path] = scalar(hf[dx_path])
        for nx_path, dx_path in TIMES:
            if dx_path not in hf:
                continue
            # a handful of acquisitions stored the string 'None' here
            try:
                stamp = datetime.strptime(text(hf[dx_path]).strip(), DX_TIME_FORMAT)
            except ValueError:
                continue
            out[nx_path] = stamp.strftime(OUT_TIME_FORMAT)
    return out


def restore(path, values, apply_changes, tally):
    """Write the archive values into one NeXus file.  Returns what was done."""
    size_before = os.path.getsize(path)
    mode = 'r+' if apply_changes else 'r'
    actions = {}
    with h5py.File(path, mode) as hf:
        for nx_path, _dx_path, const in SCALARS:
            if nx_path not in hf or nx_path not in values:
                actions[nx_path] = 'absent'
                continue
            current = scalar(hf[nx_path])
            if current == values[nx_path]:
                actions[nx_path] = 'already'
            elif current == const:
                actions[nx_path] = 'write'
                if apply_changes:
                    hf[nx_path][()] = values[nx_path]
            else:
                actions[nx_path] = 'unexpected'

        for nx_path, _dx_path in TIMES:
            if nx_path not in hf or nx_path not in values:
                actions[nx_path] = 'absent'
                continue
            current = text(hf[nx_path]).strip()
            if current == values[nx_path]:
                actions[nx_path] = 'already'
            elif STALE_TIME_RE.match(current):
                actions[nx_path] = 'write'
                if apply_changes:
                    hf[nx_path][()] = values[nx_path]
            else:
                actions[nx_path] = 'unexpected'

    # The scalars are fixed-width, so a file whose timestamps were already
    # correct must not have changed size.  The timestamps are variable-length
    # strings and do move the file's global heap, so only check when no
    # timestamp was written.
    wrote_time = any(actions.get(p) == 'write' for p, _ in TIMES)
    if apply_changes and not wrote_time:
        size_after = os.path.getsize(path)
        if size_after != size_before:
            raise RuntimeError(f'{path}: size changed {size_before} -> {size_after}')

    for state in actions.values():
        tally[state] = tally.get(state, 0) + 1
    return actions


def acquisitions():
    """Every acquisition reachable from either level of the NeXus tree."""
    stems = {f[:-len('_results.hdf')]
             for f in os.listdir(RESULTS)
             if f.endswith('_results.hdf')
             and not f.startswith(('Average', 'Avg', 'Stack'))}
    stems |= {d for d in os.listdir(NEXUS)
              if d != 'reprocess_results' and os.path.isdir(os.path.join(NEXUS, d))}
    return sorted(stems)


def targets(stem):
    """The NeXus files of one acquisition that exist, tagged by level."""
    out = []
    res = os.path.join(RESULTS, f'{stem}_results.hdf')
    if os.path.exists(res):
        out.append(('results', res))
    meta = os.path.join(NEXUS, stem, f'{stem}_metadata.hdf')
    if os.path.exists(meta):
        out.append(('metadata', meta))
    return out


def process_one(stem):
    """Restore one acquisition.  Top level so a process pool can call it."""
    src = archive_file(stem)
    if src is None:
        return stem, None, {}, [], []
    values = read_archive(src)
    missing = [p for p, _d, _c in SCALARS if p not in values] + \
              [p for p, _d in TIMES if p not in values]
    tally, unexpected = {}, []
    for _level, path in targets(stem):
        actions = restore(path, values, _APPLY, tally)
        if 'unexpected' in actions.values():
            unexpected.append((path, actions))
    return stem, values, tally, unexpected, missing


# Set once in each worker by the pool initialiser; the writes are independent
# per acquisition, so the only shared state is this flag.
_APPLY = False


def _init(apply_changes):
    global _APPLY
    _APPLY = apply_changes


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--apply', action='store_true',
                    help='perform the writes (default is a dry run)')
    ap.add_argument('--limit', type=int, default=0,
                    help='only the first N acquisitions')
    ap.add_argument('--jobs', type=int, default=16,
                    help='worker processes; the work is NFS latency, not CPU')
    args = ap.parse_args()

    stems = acquisitions()
    if args.limit:
        stems = stems[:args.limit]
    what = 'APPLY' if args.apply else 'DRY RUN'
    print(f'{what}: {len(stems)} acquisitions'
          f' on {args.jobs} workers')
    print(f'  archive (read-only) : {ARCHIVE}')
    print(f'  nexus   (written)   : {NEXUS}')

    tally = {}
    no_archive, incomplete, unexpected = [], [], []
    rows = []
    with multiprocessing.Pool(args.jobs, _init, (args.apply,)) as pool:
        for i, out in enumerate(pool.imap_unordered(process_one, stems, chunksize=16), 1):
            stem, values, counts, odd, missing = out
            if i % 2000 == 0:
                print(f'  ... {i}/{len(stems)}', flush=True)
            if values is None:
                no_archive.append(stem)
                continue
            if missing:
                incomplete.append((stem, missing))
            unexpected.extend(odd)
            for state, n in counts.items():
                tally[state] = tally.get(state, 0) + n
            rows.append([stem,
                         values.get(SCALARS[0][0], ''),
                         values.get(SCALARS[1][0], ''),
                         values.get(SCALARS[2][0], ''),
                         values.get(TIMES[0][0], ''),
                         values.get(TIMES[1][0], '')])
    rows.sort()

    with open(MANIFEST, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['acquisition', 'incident_beam_intensity',
                    'transmitted_beam_intensity', 'ring_current',
                    'start_time', 'end_time'])
        w.writerows(rows)

    print(f'\nfields by outcome: {dict(sorted(tally.items()))}')
    print(f'acquisitions with no archive file : {len(no_archive)}')
    print(f'acquisitions missing some field   : {len(incomplete)}')
    print(f'files with an unexpected value    : {len(unexpected)}')
    for stem in no_archive[:5]:
        print(f'   no archive: {stem}')
    for stem, missing in incomplete[:5]:
        print(f'   incomplete: {stem} -> {missing}')
    for path, actions in unexpected[:5]:
        print(f'   unexpected: {path} -> {actions}')
    print(f'manifest -> {MANIFEST} ({len(rows)} rows)')
    if not args.apply:
        print('\nnothing was written; rerun with --apply')


if __name__ == '__main__':
    main()
