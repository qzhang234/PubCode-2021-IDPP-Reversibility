"""Figures S11 and S12: no structural or dynamical change within an exposure.

Figure S11, after the Overleaf Figure S10.  I(Q) of the twelve consecutive
8,333-frame subsets of one 2 s acquisition, for the first (20.0 C) and the
last (45.2 C) acquisition of the (VPAVG)10 heating ramp, B039.  The subsets
are those the correlator stores, multiplied back by their mean intensity and
put on the photons-per-pixel-per-frame scale of the whole acquisition (see
reduce_gelation.segments); data/B039_segments.h5.

Figure S12, after Old/H02_RadDam_Check/H02_RampUp_XPCS_Plot_Combined.ipynb.
g2 at the lowest Q bin for the ten heating sections of B039, correlated from
all 100,000 frames of each acquisition (circles), from the first 50,000
(squares) and from the last 50,000 (triangles); data/B039_n10.h5 and
B039_part{1,2}_n10.h5.  Two changes from the original: each section is drawn
once (the original redrew every curve once per acquisition in it), and points
below g2 = 0.98 are clipped by the axis instead of being dropped.
"""

import os
import sys

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..'))
from common.acs_style import (DOUBLE_COL, SINGLE_COL, MS, LW_THIN,          # noqa: E402
                              apply_style, add_minor_grid, label_panels,
                              q_log_ticks, save_fig)
from gelation_data import DATA, Q_MAX, load                                 # noqa: E402

QI = 0
apply_style()

# --- Figure S11 ------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(DOUBLE_COL, 2.9))
with h5py.File(os.path.join(DATA, 'B039_segments.h5'), 'r') as h:
    for ax, name in zip(axes, sorted(h)):
        g = h[name]
        q, sub = g['q'][()], g['saxs_1d_subsets'][()]
        ok = q <= Q_MAX
        spread = np.nanstd(sub[:, ok], axis=0) / np.nanmean(sub[:, ok], axis=0)
        lowq = q[ok] < 0.01
        print(f"{name}: {g.attrs['temperature']:.2f} C, {len(sub)} subsets; subset-to-subset "
              f"spread of I(Q) {np.median(spread[lowq]):.1%} median below 0.01 A^-1, "
              f"{np.median(spread):.1%} over the whole range")
        for k, curve in enumerate(sub):
            ax.plot(q[ok], curve[ok], '-', color=plt.cm.viridis(k / (len(sub) - 1)),
                    lw=LW_THIN, label=f'{k + 1}')
        ax.set_xscale('log')
        ax.set_yscale('log')
        q_log_ticks(ax)
        ax.set_title(f"first acquisition, {g.attrs['temperature']:.1f} °C" if name.endswith('00001')
                     else f"last acquisition, {g.attrs['temperature']:.1f} °C",
                     fontsize=plt.rcParams['font.size'])
        ax.set_xlabel(r'$Q$ ($\AA^{-1}$)')
        add_minor_grid(ax)
axes[0].set_ylabel('$I(Q)$ (photons pixel$^{-1}$ frame$^{-1}$)')
axes[0].legend(title='subset', ncol=3, loc='lower left', handlelength=1.2,
               columnspacing=0.8)
label_panels(axes)
fig.tight_layout(pad=0.4, w_pad=1.5)
save_fig(fig, 'FigureS11_SAXS_Stability.pdf')

# --- Figure S12 ------------------------------------------------------------------
runs = [('B039_n10', 'o', 'all 100,000 frames'), ('B039_part1_n10', 's', 'first 50,000'),
        ('B039_part2_n10', '^', 'last 50,000')]
data = [load(k) for k, _, _ in runs]
T = data[0]['section/T_mean']
fig = plt.figure(figsize=(SINGLE_COL, 3.75))
gs = fig.add_gridspec(2, 1, height_ratios=[3.2, 0.95], hspace=0.36)
ax, axl = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
axl.axis('off')
for i in range(len(T)):
    c = plt.cm.jet(i / len(T))
    for d, (_, marker, _) in zip(data, runs):
        ax.errorbar(d['tau'], d['g2'][i, :, QI], yerr=d['g2_err'][i, :, QI], fmt=marker,
                    color=c, mfc='none', ms=MS * 0.8, mew=LW_THIN, elinewidth=LW_THIN)
# how far the half-exposure correlations sit from the full one, delay by delay
for d, (k, _, _) in zip(data[1:], runs[1:]):
    n = d['g2'].shape[1]
    dev = np.abs(d['g2'][:, :, QI] - data[0]['g2'][:, :n, QI])
    sig = np.hypot(d['g2_err'][:, :, QI], data[0]['g2_err'][:, :n, QI])
    late = d['tau'] > 1e-3
    print(f'{k}: median |g2 - g2(full)| = {np.median(dev[:, late]):.4f} for delays > 1 ms '
          f'({np.median(dev[:, late] / sig[:, late]):.2f} sigma)')
ax.set_xscale('log')
ax.set_ylim(0.98, 1.24)
ax.set_xlabel(r'$\Delta t$ (s)')
ax.set_ylabel(r'$g_2(Q = 0.0038\ \AA^{-1},\ \Delta t)$')
add_minor_grid(ax)
temps = [Line2D([], [], ls='none', marker='o', color=plt.cm.jet(i / len(T)), mfc='none',
                ms=MS, mew=LW_THIN, label=f'{T[i]:.1f} °C') for i in range(len(T))]
shapes = [Line2D([], [], ls='none', marker=m, color='k', mfc='none', ms=MS, mew=LW_THIN,
                 label=lab) for _, m, lab in runs]
leg = axl.legend(handles=temps, ncol=5, loc='upper center', handlelength=0.9,
                 columnspacing=0.6, handletextpad=0.3)
axl.add_artist(leg)
axl.legend(handles=shapes, ncol=3, loc='lower center', handlelength=0.9, columnspacing=0.6,
           handletextpad=0.3)
fig.subplots_adjust(left=0.17, right=0.97, top=0.98, bottom=0.02)
save_fig(fig, 'FigureS12_XPCS_Stability.pdf')
