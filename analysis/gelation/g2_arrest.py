"""Figure 3: g2 of (VPAVG)10 through the heating ramp, at the lowest Q.

A port of Old/H02/H02_RampUp_XPCS_Plot.ipynb, which made Figure 3 of the
Overleaf version.  Reads data/B039_n10.h5: the 2780 acquisitions of the heating
run B039, 20.0-45.2 C at 0.10 C/min, in ten consecutive sections of 278, each
averaged after the cross-correlation outlier cut.

Changes from the Overleaf figure: the key to the colours, one temperature per
section, is in the artwork instead of the caption, and the Q of the bin is
stated as it is, 0.0038 A^-1 (the original label read 0.03 nm^-1).
"""

import os
import sys

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..'))
from common.acs_style import (SINGLE_COL, MS, LW_THIN,                     # noqa: E402
                              apply_style, add_minor_grid, save_fig)
from gelation_data import load                                            # noqa: E402

QI = 0                      # lowest dynamic Q bin

d = load('B039_n10')
tau, g2, g2e = d['tau'], d['g2'][:, :, QI], d['g2_err'][:, :, QI]
T, Tlo, Thi = d['section/T_mean'], d['section/T_min'], d['section/T_max']
print(f"Q = {d['q_dynamic'][QI]:.5f} A^-1; sections of {d['section/n_files'][0]} acquisitions, "
      f"{d['section/n_kept'].min()}-{d['section/n_kept'].max()} kept")
for i in range(len(T)):
    print(f'  {Tlo[i]:5.2f}-{Thi[i]:5.2f} C (mean {T[i]:5.2f}): g2 = {g2[i, 0]:.3f} at '
          f'{tau[0] * 1e6:.0f} us, {g2[i, -1]:.3f} at {tau[-1]:.2f} s')

apply_style()
fig = plt.figure(figsize=(SINGLE_COL, 3.45))
gs = fig.add_gridspec(2, 1, height_ratios=[3.2, 0.62], hspace=0.38)
ax, axl = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
axl.axis('off')
for i in range(len(T)):
    c = plt.cm.jet(i / len(T))
    ax.errorbar(tau, g2[i], yerr=g2e[i], fmt='o', color=c, mfc='none', ms=MS * 0.85,
                mew=LW_THIN, elinewidth=LW_THIN, label=f'{T[i]:.1f} °C')
ax.set_xscale('log')
ax.set_ylim(0.97, 1.26)
ax.set_xlabel(r'$\Delta t$ (s)')
ax.set_ylabel(r'$g_2(Q = 0.0038\ \AA^{-1},\ \Delta t)$')
add_minor_grid(ax)
h = [Line2D([], [], ls='none', marker='o', color=plt.cm.jet(i / len(T)), mfc='none',
            ms=MS, mew=LW_THIN, label=f'{T[i]:.1f} °C') for i in range(len(T))]
axl.legend(handles=h, ncol=5, loc='center', handlelength=1.0, columnspacing=0.6,
           handletextpad=0.3, title='heating, section mean temperature')
fig.subplots_adjust(left=0.17, right=0.97, top=0.98, bottom=0.03)
save_fig(fig, 'Figure3_Dynamic_Arrest.pdf')
