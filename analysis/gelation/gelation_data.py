"""Read the reduced ramp data in data/ and put I(Q) on the absolute scale.

Every gelation figure script imports this, so all of them read the files, pick
the buffer and calibrate the intensity the same way.  Nothing here writes.

ABSOLUTE SCALE.  Each section carries the mean ion-chamber readings of the
acquisitions averaged into it, so each curve is scaled by its own flux, the
way ../SAXPCS_8id/abs_xsec.py prescribes:

    I_abs(Q) = coef_sample * I_sample(Q) - coef_buffer * I_buffer(Q)

2022-1 recorded both chambers, so the coefficient uses the measured
transmission (abs_xsec.coef_from_monitors).  2021-2 recorded only the
incident chamber, so its transmission is the 2022-1 value of the same
construct in the same cell (abs_xsec.coef_from_monitors_2021_2).  The water
blank D029 is given that same transmission, so for 2021-2 the subtraction
weight is the ratio of the incident fluxes alone.  Lengths are in mm, so
I_abs is in mm^-1.

The Overleaf version of these figures used one fixed pair of coefficients,
6.93e4 and 7.62e4, for both beamtimes.  They had been computed for a 1 mm path
through the cell, which is 2 mm thick, and the axes were labelled cm^-1.  For
2022-1 the per-section coefficients here are within 5 % of half those values;
the 2021-2 ones follow from that beamtime's own incident flux.
"""

import os
import sys

import h5py

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')
sys.path.insert(0, os.path.join(HERE, '..', 'SAXPCS_8id'))
from abs_xsec import coef_from_monitors, coef_from_monitors_2021_2   # noqa: E402

# chain length -> (label, heating key, cooling key, buffer key)
SAMPLES = {
    10: (r'(VPAVG)$_{10}$', 'B039_n10', 'B040_n10', 'D029'),
    20: (r'(VPAVG)$_{20}$', 'T214_n20', 'T215_n20', 'D029'),
    30: (r'(VPAVG)$_{30}$', 'E0110_n10', 'E0111_n10', 'D0138'),
}

# The Overleaf figures dropped the last 15 Q bins, and so do these.  They lie in
# the detector corners: by 0.032 A^-1 only one or two of the eight phi sectors
# still hold pixels, a few hundred in all, falling to ten in the last bin.
Q_MAX = 0.032          # A^-1


def load(key):
    """Every array and attribute of data/<key>.h5, as a dict."""
    out = {}
    with h5py.File(os.path.join(DATA, f'{key}.h5'), 'r') as h:
        out.update(h.attrs)
        for name in ('tau', 'q_static', 'q_dynamic', 'g2', 'g2_err', 'saxs_1d'):
            out[name] = h[name][()]
        for grp in ('section', 'acquisition'):
            for name, ds in h[grp].items():
                out[f'{grp}/{name}'] = ds[()]
    return out


def coefficients(d):
    """Absolute cross-section coefficient of every section of one dataset."""
    up = d['section/incident']
    if d['beamtime'] == '2022-1':
        return coef_from_monitors(up, d['section/transmitted'])
    return coef_from_monitors_2021_2(up)


def absolute_iq(key, buffer_key):
    """(Q, I_abs) with I_abs of shape (sections, Q), buffer subtracted, mm^-1."""
    d, b = load(key), load(buffer_key)
    sample = coefficients(d)[:, None] * d['saxs_1d']
    buffer = coefficients(b)[0] * b['saxs_1d'][0]
    return d['q_static'], sample - buffer[None, :]
