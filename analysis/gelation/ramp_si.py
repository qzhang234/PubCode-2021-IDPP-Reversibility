"""Figures S2-S4: one heating and cooling ramp per chain length.

A port of Old/H02_SI_RampUpDown, Old/H04_SI and Old/H06_SI (the H0x_Plot
notebooks), which made Figures S2-S4 of the Overleaf version.  Each figure has
the heating ramp on the left and the cooling ramp that followed on the right:

    (a), (b)  temperature of every acquisition against elapsed time
    (c), (d)  I(Q) of each section, buffer subtracted, absolute scale
    (e), (f)  g2 at the lowest Q bin

    Figure S2  (VPAVG)10  B039 / B040, 10 sections each   2021-2
    Figure S3  (VPAVG)20  T214 / T215, 20 sections each   2021-2
    Figure S4  (VPAVG)30  E0110 / E0111, 10 sections each 2022-1

Colour is the section's mean temperature on one scale per figure, so a colour
means the same temperature in both columns.

g2 CURVES ARE OMITTED where five or more of their 61 points lie outside
1.00-1.25, the rule the original applied to Figure S4.  Here it is applied to
all three figures, which removes the seven coldest heating sections of
(VPAVG)30, whose coherent intensity at this Q is too weak for a correlation,
its two coldest cooling sections, and one cooling section of (VPAVG)10 at
26.3 C.  In that one, 48 consecutive acquisitions (B040 02076 onward) carry a
flat g2 of 1.7-1.8 at the lowest Q bin only; the cross-correlation cut compares
shapes, not levels, and keeps them, so the section average sits at 1.29.  The
original hid that curve by stopping the axis at 1.23.  The script prints what
it omits.

Other changes from the Overleaf figures: the time axes are the recorded start
times (5 s per acquisition in 2021-2, 6 s in 2022-1; the original assumed 6 s
throughout), I is in mm^-1 (see gelation_data), Q is in A^-1, and a colour bar
keys the temperatures.
"""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..'))
from common.acs_style import (DOUBLE_COL, MS, LW_THIN,                    # noqa: E402
                              apply_style, add_minor_grid, label_panels,
                              q_log_ticks, save_fig)
from gelation_data import Q_MAX, SAMPLES, absolute_iq, load              # noqa: E402

QI = 0
G2_RANGE, G2_MAX_OUT = (1.0, 1.25), 5
FIGURES = {10: 'FigureS2_Ramp_VPAVG10.pdf', 20: 'FigureS3_Ramp_VPAVG20.pdf',
           30: 'FigureS4_Ramp_VPAVG30.pdf'}


def column(axes, key, buffer_key, norm):
    """Draw one ramp into a column of three axes; returns the g2 curves it omitted."""
    d = load(key)
    q, I = absolute_iq(key, buffer_key)
    T = d['section/T_mean']
    colors = plt.cm.jet(norm(T))
    sec = d['acquisition/section']
    axes[0].scatter(d['acquisition/elapsed_s'] / 1000, d['acquisition/temperature'],
                    c=colors[sec], s=0.6, lw=0, rasterized=True)
    show = np.arange(len(q)) % 3 == 0
    omitted = []
    for i in range(len(T)):
        ok = show & (q <= Q_MAX) & (I[i] > 0)
        axes[1].plot(q[ok], I[i][ok], 'o-', color=colors[i], lw=LW_THIN, ms=MS * 0.45,
                     mew=LW_THIN, mfc='none')
        g, e = d['g2'][i, :, QI], d['g2_err'][i, :, QI]
        if np.sum((g < G2_RANGE[0]) | (g > G2_RANGE[1])) >= G2_MAX_OUT:
            omitted.append(T[i])
            continue
        axes[2].errorbar(d['tau'], g, yerr=e, fmt='o', color=colors[i], ms=MS * 0.45,
                         mew=LW_THIN, mfc='none', elinewidth=LW_THIN)
    return omitted


def make(n):
    label, up, down, buf = SAMPLES[n]
    Ts = np.concatenate([load(k)['section/T_mean'] for k in (up, down)])
    norm = Normalize(vmin=Ts.min(), vmax=Ts.max())

    apply_style()
    fig, ax = plt.subplots(3, 2, figsize=(DOUBLE_COL, 6.3), sharey='row',
                           gridspec_kw=dict(hspace=0.5, wspace=0.06))
    for c, (key, title) in enumerate(((up, 'heating'), (down, 'cooling'))):
        omitted = column(ax[:, c], key, buf, norm)
        ax[0, c].set_title(f'{label}, {title}', fontsize=plt.rcParams['font.size'])
        print(f'{label} {title} ({key}): g2 omitted at '
              + (', '.join(f'{t:.1f} C' for t in omitted) or 'none'))
    for c in range(2):
        ax[0, c].set_xlabel('Elapsed time (1000 s)')
        ax[1, c].set_xscale('log')
        ax[1, c].set_yscale('log')
        q_log_ticks(ax[1, c])
        ax[1, c].set_xlabel(r'$Q$ ($\AA^{-1}$)')
        ax[2, c].set_xscale('log')
        ax[2, c].set_xlabel(r'$\Delta t$ (s)')
        for r in range(3):
            add_minor_grid(ax[r, c])
    ax[0, 0].set_ylabel('Temperature (°C)')
    ax[1, 0].set_ylabel(r'$I(Q)$ (mm$^{-1}$)')
    ax[2, 0].set_ylabel(r'$g_2(Q = 0.0038\ \AA^{-1})$')
    ax[2, 0].set_ylim(0.97, 1.26)

    fig.subplots_adjust(left=0.085, right=0.985, top=0.96, bottom=0.16)
    cax = fig.add_axes([0.25, 0.062, 0.5, 0.016])
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=plt.cm.jet), cax=cax,
                      orientation='horizontal')
    cb.set_label('Section mean temperature (°C)', labelpad=2)
    cb.outline.set_linewidth(LW_THIN)
    # the right column shares its y axis and has no tick labels to clear
    for r in range(3):
        label_panels([ax[r, 0]], dx=-30, start=2 * r)
        label_panels([ax[r, 1]], dx=-10, start=2 * r + 1)
    save_fig(fig, FIGURES[n])


for n in FIGURES:
    make(n)
