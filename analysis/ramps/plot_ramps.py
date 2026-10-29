"""Figure 3: structure and dynamics of three chain lengths through a thermal ramp.

Reads only ramps/data/, written by reduce_ramps.py, so this runs anywhere; no
/gdata, no beamline account.

THE COMPARISON.  One protocol, three chain lengths, up and back down:

    (VPAVG)10   B039 up, B040 down    20 -> 45.2 C
    (VPAVG)20   T214 up, T215 down    20 -> 41.4 C
    (VPAVG)30   H041 up, H042 down    20 -> 33.1 C

The ceiling falls with chain length because the stage ran one protocol per
cell, so each sample reaches a different part of its own phase diagram.  The
rows are therefore read down a column -- same sample, up then down -- and
across a row only at temperatures all three reached.

LAYOUT.  Three columns, one per chain length.  Top row I(Q) on the way up,
bottom row g2 at the lowest usable Q, colour running with temperature on one
scale shared by every panel so a colour means the same temperature throughout.
The down-ramp is drawn in the same panels with open markers, which is what
makes reversibility or its absence visible in one look: for (VPAVG)30 the two
directions lie on top of each other, and for (VPAVG)10 they do not.

ABSOLUTE SCALE.  Through abs_xsec_coef_2021_2(), not abs_xsec_coef(): 2021-2
never recorded the downstream ion chamber, so the transmission comes from
2022-1 and the scale carries about 6% of systematic as a common factor.  It
moves no peak, no exponent and no ratio within a profile, and it is the same
factor for all three samples, so it cannot touch the comparison this figure
is read for.  See the comment at TRANSMISSION_2021_2.
"""

import glob
import os
import re
import sys

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'SAXPCS_8id'))
sys.path.insert(0, os.path.join(HERE, '..'))
from common.acs_style import (DOUBLE_COL, MS, MEW, LW_THIN, LW_DATA,      # noqa: E402
                              apply_style, add_minor_grid, label_panels,
                              q_log_ticks, save_fig)
from abs_xsec import abs_xsec_coef_2021_2                                 # noqa: E402
from nexus_read import read_g2                                            # noqa: E402

DATA = os.path.join(HERE, 'data')
TRACE = os.path.join(DATA, 'trace_2021-2.csv')

SAMPLES = [('(VPAVG)$_{10}$', 'B039_VPAVG_H02_500C10p', 'B040_VPAVG_H02_200C10p'),
           ('(VPAVG)$_{20}$', 'T214_H04_2nd_500C10p',   'T215_H04_2nd_200C10p'),
           ('(VPAVG)$_{30}$', 'H041_VPAVG_H05_2_500C10p', 'H042_VPAVG_H05_2_200C10p')]

T_LO, T_HI = 20.0, 46.0          # C, the colour scale shared by every panel
CMAP = plt.cm.coolwarm
QI_G2 = 0                        # lowest dynamic bin, the largest length scale
FIG_SIZE = (DOUBLE_COL, 4.4)


def mean_temperature(run, f0, f1):
    """Mean recorded temperature over the frames an average was built from."""
    T = []
    with open(TRACE) as fh:
        next(fh)
        for line in fh:
            p = line.rstrip('\n').split(',')
            if p[0] == run and f0 <= int(p[1]) <= f1:
                T.append(float(p[3]))
    return float(np.mean(T)) if T else np.nan


def averages(run):
    """(T, file) for every averaged bin of a run, in temperature order."""
    out = []
    for fp in sorted(glob.glob(os.path.join(DATA, f'Average_{run}_*_results.hdf'))):
        m = re.search(r'_(\d{5})_(\d{5})_results\.hdf$', fp)
        if not m:
            continue
        f0, f1 = int(m.group(1)), int(m.group(2))
        out.append((mean_temperature(run, f0, f1), fp))
    return sorted(out)


def read_iq(fp):
    """Absolute, azimuthally averaged I(Q) for one averaged file."""
    with h5py.File(fp, 'r') as hf:
        inten = np.asarray(hf['/xpcs/temporal_mean/scattering_1d'][()], float).reshape(-1)
        idx = hf['/xpcs/qmap/static_index_mapping'][()]
        qlist = hf['/xpcs/qmap/static_v_list_dim0'][()]
        nphi = hf['/xpcs/qmap/static_v_list_dim1'].shape[0]
        coef = abs_xsec_coef_2021_2(hf)
    qidx = idx // nphi
    uq = np.unique(qidx)
    return qlist[uq], np.array([np.nanmean(inten[qidx == i]) for i in uq]) * coef


apply_style()
norm = Normalize(vmin=T_LO, vmax=T_HI)
fig, axes = plt.subplots(2, 3, figsize=FIG_SIZE)

for col, (label, run_up, run_dn) in enumerate(SAMPLES):
    ax_i, ax_g = axes[0, col], axes[1, col]
    for run, filled in ((run_up, True), (run_dn, False)):
        for T, fp in averages(run):
            if not np.isfinite(T):
                continue
            c = CMAP(norm(T))
            q, I = read_iq(fp)
            ok = (q > 0.0035) & (q < 0.03) & np.isfinite(I) & (I > 0)
            ax_i.plot(q[ok], I[ok], '-' if filled else '--', color=c, lw=LW_DATA)
            with h5py.File(fp, 'r') as hf:
                tau, g2, g2e, _qd = read_g2(hf)
            u = (tau > 0) & np.isfinite(g2[:, QI_G2])
            ax_g.plot(tau[u], g2[u, QI_G2], 'o' if filled else 's', color=c,
                      ms=MS * 0.6, mfc=c if filled else 'none', mew=MEW, lw=0)
    ax_i.set_xscale('log'); ax_i.set_yscale('log')
    q_log_ticks(ax_i)
    ax_i.set_title(label, fontsize=plt.rcParams['font.size'])
    ax_i.set_xlabel(r'$Q$ ($\AA^{-1}$)')
    ax_g.set_xscale('log')
    ax_g.set_xlabel(r'$\Delta t$ (s)')
    for a in (ax_i, ax_g):
        add_minor_grid(a)
axes[0, 0].set_ylabel(r'$I(Q)$ (mm$^{-1}$)')
axes[1, 0].set_ylabel(r'$g_2$')

# one colour bar for the whole figure: a colour means one temperature everywhere
fig.tight_layout(rect=(0, 0.08, 1, 1))
cax = fig.add_axes([0.12, 0.045, 0.76, 0.022])
cb = fig.colorbar(ScalarMappable(norm=norm, cmap=CMAP), cax=cax, orientation='horizontal')
cb.set_label(r'Temperature ($^\circ$C);  solid / filled = heating,'
             r'  dashed / open = cooling', labelpad=2)
cb.outline.set_linewidth(LW_THIN)
label_panels(list(axes.ravel()))
save_fig(fig, 'Figure3_Chain_Length_Ramps.pdf')
