"""Figure S9: 2x2 grid of XPCS g2(tau) for the 4 non-primary q bins.

This figure is an EXPANSION of Figure 3b onto the four higher q bins: Figure 3b
shows q index 0, this shows q indices 1-4.  It therefore has to show the same
model, fitted the same way -- so the curves here come from exactly the same
global fit that produces Figure 3c and Figure S8, imported from xpcs_fit.py.

For each elapsed time the five q bins are fitted SIMULTANEOUSLY with the two
stretching exponents shared across q and the contrast fixed at the measured
instrumental value; the curve drawn in each panel is that global solution
evaluated at that panel's q.  (Previously this script carried its own private
copy of the model which fitted each q bin independently with both exponents
frozen at 0.5, so the caption's description of the fits did not match the
curves; the parameters shown here are now the ones the paper reports.)
"""

import glob
import os
import sys

import numpy as np
import h5py
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from common.acs_style import (DOUBLE_COL, MS, MEW, LW_THIN, LW_DATA,
                              apply_style, add_minor_grid, label_panels, save_fig)
from matplotlib.lines import Line2D
from xpcs_fit import double_exp, fit_g2_joint

FIG_SIZE = (DOUBLE_COL, 5.0)         # 2x2 panels + the shared key row beneath
LEGEND_H = 0.05                      # height fraction reserved for that key
G2_YLIM = (0.98, 1.20)               # trims 2 of 1220 points above and 4 below,
                                     # plus 4 error-bar caps, all at the noisiest
                                     # shortest delays

# --- PARAMETERS ---
data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
XPCS_HEADER = 'B0147'
N_LAST = 5
# XPCS elapsed-time colours: same scale as saxpcs.py -- the N_LAST times are
# mapped to evenly-spaced positions of this colormap for maximum separation.
XPCS_CMAP = plt.cm.plasma
CMAP_LO, CMAP_HI = 0.10, 0.88

# The model, the measured contrast and the global fit come from xpcs_fit.py --
# the same implementation Figure 3 and Figure S8 use.
fit_q_indices = [0, 1, 2, 3, 4]      # all five bins enter the global fit ...
grid_q_indices = [1, 2, 3, 4]        # ... these four are the ones plotted here

# Field locations and readers are shared with the other 8-ID figure scripts.
from nexus_read import parse_name, read_start_time, read_g2


# --- DISCOVER FILES ---
# every range average in data/, including the thermal-cycle groups that belong to
# Figure S4; the loop below keeps only XPCS_HEADER
file_paths = sorted(glob.glob(os.path.join(data_dir, 'Average_*.hdf')))
b0147, start_times = [], {}
for fp in file_paths:
    if parse_name(fp)[0] != XPCS_HEADER:
        continue
    with h5py.File(fp, 'r') as hf:
        start_times[fp] = read_start_time(hf)
    b0147.append(fp)
b0147.sort(key=lambda p: parse_name(p)[1])

xpcs_files = b0147[-N_LAST:]                    # last N only (NOT the 1-200 file)
t_ref = start_times[b0147[0]]                   # time origin = first B0147 file
elapsed = {fp: (start_times[fp] - t_ref).total_seconds() for fp in xpcs_files}

# Elapsed-time colour, matching saxpcs.py: N_LAST times mapped to evenly spaced
# colormap positions (sorted by time) so consecutive times are easy to tell apart.
_xpcs_sorted = sorted(xpcs_files, key=lambda p: elapsed[p])
_xpcs_pos = np.linspace(CMAP_LO, CMAP_HI, len(_xpcs_sorted))
xpcs_color = {fp: XPCS_CMAP(p) for fp, p in zip(_xpcs_sorted, _xpcs_pos)}

# Pre-load g2 once per file.
data = {}
for fp in xpcs_files:
    with h5py.File(fp, 'r') as hf:
        data[fp] = read_g2(hf)

# --- FIGURE ---
apply_style()
fig, axes = plt.subplots(2, 2, figsize=FIG_SIZE, sharex=True, sharey=True)
label_panels(axes.flat)

# One global fit per elapsed time, over all five q bins at once -- identical to
# saxpcs.py, so the curves below are the fits the paper reports.
# The exponents are shared across the whole waiting-time series, exactly as in
# saxpcs.py, so these curves are the fits the paper reports.
_joint, _per_time = fit_g2_joint([data[fp][:3] for fp in xpcs_files],
                                 [data[fp][3] for fp in xpcs_files],
                                 fit_q_indices)
fits = dict(zip(xpcs_files, _per_time))
print(f"  joint fit: chi2/dof = {_joint['red_chi2']:.3f}, "
      f"shared p_fast = {_joint['p1']:.3f}")

for ax, q_idx in zip(axes.flat, grid_q_indices):
    for fp in xpcs_files:
        tau, g2, g2_err, q_vals = data[fp]
        color = xpcs_color[fp]
        m = tau > 0
        ax.errorbar(tau[m], g2[m, q_idx], yerr=g2_err[m, q_idx], fmt='o', color=color,
                    mfc='none', markersize=MS, mew=LW_THIN, capsize=1.5,
                    elinewidth=LW_THIN, capthick=LW_THIN, alpha=0.85, linestyle='none')
        r = fits[fp]
        if r is not None and q_idx in r['per_q']:
            pq = r['per_q'][q_idx]
            t_fit = np.logspace(np.log10(tau[m].min()), np.log10(tau[m].max()), 200)
            ax.plot(t_fit, double_exp(t_fit, pq['tau_fast'], pq['f'], pq['tau_slow'],
                                      r['p1'], r['p2']), color=color, lw=LW_DATA)
    ax.set_xscale('log')
    ax.set_ylim(*G2_YLIM)
    ax.set_yticks([1.0, 1.05, 1.10, 1.15])
    ax.set_title(rf'$Q = {q_vals[q_idx]:.5f}\ \AA^{{-1}}$')
    add_minor_grid(ax)

# elapsed-time key in the artwork (ACS prefers this over explaining colours in
# the caption).  The same five colours serve all four panels, so -- as in
# Figure 3 -- the key belongs to the figure and sits in one row along the
# bottom; out there it costs no panel any headroom, which is what lets the y
# axis stop at 1.20.
t_handles = [Line2D([], [], marker='o', ls='none', mfc='none', mew=MEW,
                    mec=xpcs_color[fp], markersize=MS, label=f'{elapsed[fp]:.0f} s')
             for fp in sorted(xpcs_files, key=lambda p: elapsed[p])]
leg = fig.legend(handles=t_handles, loc='lower center', bbox_to_anchor=(0.5, 0.0),
                 ncol=len(t_handles), handlelength=1.2, handletextpad=0.4,
                 columnspacing=1.2, borderaxespad=0.0, borderpad=0.35)
leg.get_frame().set_linewidth(LW_THIN)

for ax in axes[-1, :]:
    ax.set_xlabel(r'Delay Time, $\Delta t$ (s)')
for ax in axes[:, 0]:
    ax.set_ylabel('$g_2$')

# sharex/sharey strip the inner tick labels, so the panels can sit close
fig.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.6, rect=(0, LEGEND_H, 1, 1))
save_fig(fig, 'FigureS9_g2_Grid.pdf')
plt.show()
