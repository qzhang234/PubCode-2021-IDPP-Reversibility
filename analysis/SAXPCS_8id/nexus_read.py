"""Where things live inside an 8-ID-I averaged HDF file, and how to read them.

The averaged NeXus files written by `average_ranges.py` all have the same
layout, and four of the figure scripts here (`saxpcs.py`, `saxs_evolution.py`,
`g2_grid_SI.py`, `thermal_cycle.py`) need the same handful of quantities out of
them: when the acquisition started, the time-averaged I(Q), and the correlation
functions.  Each script used to carry its own copy of these readers, so a reader
comparing two scripts had to check four near-identical blocks to be sure they
really were identical.  They live here once instead, and every script imports
them, so all four figures provably read the files the same way.

Nothing here modifies a file: every caller opens the HDF read-only, and these
functions only fetch.
"""

import os
import re
from datetime import datetime

import numpy as np

# --- where each quantity lives inside the file ---
START_TIME_PATH = '/entry/start_time'
TIME_FORMAT     = '%Y-%m-%d %H:%M:%S'
FRAME_TIME_PATH = '/entry/instrument/detector_1/frame_time'
DELAY_PATH      = '/xpcs/multitau/delay_list'
G2_PATH         = '/xpcs/multitau/normalized_g2'
G2_ERR_PATH     = '/xpcs/multitau/normalized_g2_err'
DYN_Q_PATH      = '/xpcs/qmap/dynamic_v_list_dim0'
SAXS_PATH       = '/xpcs/temporal_mean/scattering_1d'
STATIC_MAP_PATH = '/xpcs/qmap/static_index_mapping'
STATIC_Q_PATH   = '/xpcs/qmap/static_v_list_dim0'
STATIC_PHI_PATH = '/xpcs/qmap/static_v_list_dim1'

# 'Average_B0147_..._801_1000_results.hdf' -> header 'B0147', frames 801-1000
_NAME_RE = re.compile(r'Average_([A-Za-z]\d+)_.*?_(\d+)_(\d+)_results')


def parse_name(fname):
    """Pull ('B0147', first_frame, last_frame) out of an averaged file's name.

    Returns (None, -1, -1) if the name does not match, so a caller can sort a
    directory listing by frame range without crashing on a stray file.
    """
    m = _NAME_RE.search(os.path.basename(fname))
    if m is None:
        return None, -1, -1
    return m.group(1), int(m.group(2)), int(m.group(3))


def read_start_time(hf):
    """Acquisition start time as a datetime, for computing elapsed times.

    The stored value is a text timestamp, but different writers wrap it
    differently -- as a bare string, as bytes, or as a one-element array -- so
    it is unwrapped before parsing.
    """
    raw = hf[START_TIME_PATH][()]
    if isinstance(raw, np.ndarray):
        raw = raw.reshape(-1)[0]
    if isinstance(raw, bytes):
        raw = raw.decode('utf-8')
    return datetime.strptime(str(raw).strip(), TIME_FORMAT)


def read_saxs_iq(hf, phi_average=True):
    """Time-averaged intensity, collapsed from the (Q, phi) map down to (Q, I).

    The file stores one intensity per (Q, phi) bin, flattened into a single
    list, and `static_index_mapping` says which bin each entry belongs to.  The
    bins run phi-fastest, so integer-dividing that index by the number of phi
    sectors recovers the Q index of every entry.  With phi_average=True the
    entries sharing a Q index are averaged -- the azimuthal average used in
    every figure; otherwise only the first phi sector is kept.
    """
    intensity = np.asarray(hf[SAXS_PATH][()]).reshape(-1)
    idx_map = hf[STATIC_MAP_PATH][()]
    q_list = hf[STATIC_Q_PATH][()]
    n_phi = hf[STATIC_PHI_PATH].shape[0]
    q_idx = idx_map // n_phi
    uq = np.unique(q_idx)

    inten = []
    for qi in uq:
        entries = intensity[q_idx == qi]          # every phi sector at this Q
        if phi_average and n_phi > 1:
            inten.append(np.nanmean(entries))
        else:
            inten.append(entries[0])
    return q_list[uq], np.array(inten)


def read_g2(hf):
    """Correlation functions: (tau in seconds, g2, g2_err, Q values).

    `delay_list` is stored in FRAMES, so it is multiplied by the frame time to
    get seconds.  g2 and g2_err are 2-D, indexed (delay, q bin).
    """
    # frame_time is a single number, stored as a tiny array in some files and as
    # a plain scalar in others; .item() pulls the number out of either.
    t0 = hf[FRAME_TIME_PATH][()]
    if isinstance(t0, np.ndarray):
        t0 = t0.item()
    tau = hf[DELAY_PATH][()] * t0
    if tau.ndim > 1:
        tau = tau[:, 0]                           # some files store it as a column
    return tau, hf[G2_PATH][()], hf[G2_ERR_PATH][()], hf[DYN_Q_PATH][()]
