"""Reduce the 2021-2 chain-length ramps the way average_ranges.py reduces 2022-1.

BEAMLINE-ONLY.  Reads /gdata, which mounts on the 10.x hosts (adamite, amber)
and not on the 164.x analysis machines.  Everything downstream reads only what
this writes into ramps/data/, so plot_ramps.py runs anywhere.

THE DATA
--------
Three chain lengths, each taken up from 20 C and back down, on the QNW stage,
21 July 2021:

    (VPAVG)10   B039 up to 45.2 C    B040 down
    (VPAVG)20   T214 up to 41.4 C    T215 down
    (VPAVG)30   H041 up to 33.1 C    H042 down

The ceiling falls with chain length, which is the point of the comparison: the
same ramp protocol reaches a different part of each phase diagram.

THE SAMPLE IS NOT ALWAYS IN ZONE 1.  The stage holds three cells and parks the
idle zones at 6 C.  (VPAVG)10 and (VPAVG)20 ramp in zone 1, (VPAVG)30 in zone 3.
restore_dx_metadata_2021-2.py wrote the ACTIVE zone's reading to
/entry/sample/temperature and the zone index to /entry/sample/qnw_zone, so this
script reads those rather than choosing a zone itself; it falls back to picking
the zone furthest from the 6 C park if an older file lacks them.

GROUPING IS BY TEMPERATURE, NOT BY COUNT.  The 2021-3 reducer blocks on
acquisition count because its runs are isothermal holds and the x axis is
elapsed time.  Here the x axis is temperature and the three samples are to be
compared against each other, so equal-count blocks would put the samples on
different temperature grids.  Fixed BIN_C bins put every sample on one grid.
On a constant-rate ramp the two are the same thing anyway, up to the ends.

WHAT IS WRITTEN
---------------
1. TRACE: one row per acquisition with frame, time, temperature, set point,
   zone, incident monitor and the mean I(Q) over 0.004-0.008 A^-1.  No
   averaging, so an onset can be placed in temperature.
2. AVERAGES: one per temperature bin per run, through average_ranges.process_group
   so the two outlier cuts are exactly the ones every other group in either
   paper receives.

NO ABSOLUTE SCALE HERE.  2021-2 never recorded the downstream ion chamber, so
the transmitted-beam field is the placeholder the NeXus rewrite left and must
not be divided by.  Conversion is the plotting script's job, through
abs_xsec.abs_xsec_coef_2021_2(), which supplies the transmission from 2022-1.

    python reduce_ramps.py
"""

import csv
import os
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'SAXPCS_8id'))
sys.path.insert(0, os.path.join(HERE, '..'))
import average_ranges as ar                                       # noqa: E402

PREFIX = '/gdata/s8id-dmdtn/2021-2/babnigg202107_2_nexus/reprocess_results'
OUT_DIR = os.path.join(HERE, 'data')
TRACE_CSV = os.path.join(OUT_DIR, 'trace_2021-2.csv')

# run -> (sample label, direction)
RUNS = {
    'H041_VPAVG_H05_2_500C10p': ('(VPAVG)30', 'up'),
    'H042_VPAVG_H05_2_200C10p': ('(VPAVG)30', 'down'),
    'T214_H04_2nd_500C10p':     ('(VPAVG)20', 'up'),
    'T215_H04_2nd_200C10p':     ('(VPAVG)20', 'down'),
    'B039_VPAVG_H02_500C10p':   ('(VPAVG)10', 'up'),
    'B040_VPAVG_H02_200C10p':   ('(VPAVG)10', 'down'),
}

BIN_C = 2.0        # C per averaged group
BIN_LO = 20.0      # C; every ramp starts here
MIN_ACQ = 20       # a bin with fewer acquisitions than this is not averaged
PARKED_C = 6.0

ar.prefix = PREFIX
ar.out_dir = OUT_DIR

TEMP = '/entry/sample/temperature'
TEMP_SET = '/entry/sample/temperature_set'
ZONE = '/entry/sample/qnw_zone'


def scalar(hf, path):
    return float(np.asarray(hf[path][()]).ravel()[0])


def temperature(hf):
    """(T, T_set, zone) for one acquisition, preferring the repaired fields."""
    if TEMP in hf:
        z = int(scalar(hf, ZONE)) if ZONE in hf else 0
        tset = scalar(hf, TEMP_SET) if TEMP_SET in hf else float('nan')
        return scalar(hf, TEMP), tset, z
    # older file: choose the zone whose set point is furthest from the park
    best, best_d = 1, -1.0
    for z in (1, 2, 3):
        k = f'/entry/sample/qnw{z}_temperature_set'
        if k in hf and abs(scalar(hf, k) - PARKED_C) > best_d:
            best, best_d = z, abs(scalar(hf, k) - PARKED_C)
    return (scalar(hf, f'/entry/sample/qnw{best}_temperature'),
            scalar(hf, f'/entry/sample/qnw{best}_temperature_set'), best)


def trace(run):
    """One row per acquisition, in frame order."""
    files = ar.group_files(run)
    if not files:
        return []
    q = ar.static_q(files[0])
    band = (q >= ar.SPIKE_BAND[0]) & (q < ar.SPIKE_BAND[1])
    rows = []
    for f in files:
        with h5py.File(f, 'r') as hf:
            T, tset, z = temperature(hf)
            saxs = np.asarray(hf[ar.SAXS_PATH][()]).reshape(-1)
            rows.append([run, ar.extract_frame(f), ar.get_start_time(f),
                         T, tset, z, scalar(hf, ar.INCIDENT_PATH),
                         float(np.nanmean(saxs[band]))])
    return rows


def bins(rows):
    """Contiguous frame ranges, one per BIN_C of temperature.

    A ramp is monotonic in temperature, so a temperature bin is a contiguous
    run of frames and process_group can average it directly.  Frames are taken
    from the first to the last acquisition falling in the bin, which closes the
    one-frame gaps that noise on the RTD would otherwise punch.
    """
    frames = np.array([r[1] for r in rows])
    T = np.array([r[3] for r in rows])
    out = []
    edges = np.arange(BIN_LO, np.nanmax(T) + BIN_C, BIN_C)
    for lo in edges:
        sel = (T >= lo) & (T < lo + BIN_C)
        if sel.sum() < MIN_ACQ:
            continue
        f = frames[sel]
        out.append((lo + BIN_C / 2, int(f.min()), int(f.max()), int(sel.sum())))
    return out


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    all_rows = []
    for run, (label, direction) in RUNS.items():
        rows = trace(run)
        if not rows:
            print(f'{run}: NO RESULT FILES under {PREFIX}')
            continue
        all_rows += rows
        bb = bins(rows)
        T = [r[3] for r in rows]
        print(f'{run}  {label} {direction}: {len(rows)} acquisitions, '
              f'T {min(T):.1f}-{max(T):.1f} C, zone {rows[0][5]}, {len(bb)} bins')
        for Tc, f0, f1, n in bb:
            print(f'    {Tc:5.1f} C  frames {f0:5d}-{f1:5d}  n={n:4d}')
        ar.process_group(run, [(f0, f1) for _Tc, f0, f1, _n in bb])

    with open(TRACE_CSV, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['run', 'frame', 'start_time', 'T_C', 'T_set_C', 'qnw_zone',
                    'incident', 'I_band'])
        w.writerows(all_rows)
    print(f'\ntrace -> {TRACE_CSV}  ({len(all_rows)} rows)')
    print(f'averages -> {OUT_DIR}')


if __name__ == '__main__':
    main()
