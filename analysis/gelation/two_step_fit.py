"""Figure 4 and its table: the two-step decay of g2 for three chain lengths.

A port of Old/g2_Summary/g2_uniform_mid_check.ipynb, which made Figure 4 of the
Overleaf version.  One section of each heating ramp near 31.2 C, from the finer
binning the original used for this comparison (data/, see reduce_gelation.py):

    (VPAVG)10   B039_n20_skip70, section 8    31.10 C, 136 acquisitions
    (VPAVG)20   T214_n20, section 10           31.25 C, 119 acquisitions
    (VPAVG)30   E0110_n20, section 17          31.24 C,  60 acquisitions

g2 at the lowest Q bin is turned into g1^2 = (g2 - 1) / beta with the contrast
beta the original set for each beamtime from the short-delay level of g2
(0.19 and 0.195 for the two 2021-2 samples, 0.135 for 2022-1), and fitted with

    g1^2 = [a1 exp(-(t/tau1)^p1) + (1 - a1) exp(-(t/tau2)^p2)]^2

over every delay, with the original starting values and bounds.

ONE CHANGE.  The original weighted all three fits with the (VPAVG)10 error
bars, a copy-paste slip; here each sample is fitted with its own.  The script
prints the table both ways, and the first reproduces the Overleaf table.

Writes Figure4_Two_Step_Decay.pdf and Table_Two_Step_Fit.tex, which main.tex
inputs below the figure so that the numbers in the text cannot drift from it.
"""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from scipy.optimize import curve_fit

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..'))
from common.acs_style import (SINGLE_COL, MS, MEW, LW_THIN, LW_DATA,      # noqa: E402
                              apply_style, add_minor_grid, save_fig)
from gelation_data import load                                           # noqa: E402

QI = 0
# label, dataset, section, contrast, upper bound on p1 and p2, colour, marker
SAMPLES = [(r'(VPAVG)$_{10}$', 'B039_n20_skip70', 8, 0.19, 2, 'r', 'o'),
           (r'(VPAVG)$_{20}$', 'T214_n20', 10, 0.195, 3, 'b', 's'),
           (r'(VPAVG)$_{30}$', 'E0110_n20', 17, 0.135, 3, 'k', 'D')]
P0 = [0.5, 1e-3, 1, 1, 1]


def two_step(t, a1, tau1, tau2, p1, p2):
    return (a1 * np.exp(-(t / tau1) ** p1) + (1 - a1) * np.exp(-(t / tau2) ** p2)) ** 2


def fit(t, y, sigma, p_max):
    bounds = ([0, 1e-6, 1e-6, 0, 0], [1, 1e3, 1e3, p_max, p_max])
    p, cov = curve_fit(two_step, t, y, p0=P0, bounds=bounds, sigma=sigma)
    return p, np.sqrt(np.diag(cov))


def row(label, p, e):
    return (f'{label:16s} a1={p[0]:.3f}+/-{e[0]:.3f}  tau1={p[1] * 1e3:.1f}+/-{e[1] * 1e3:.1f} ms  '
            f'tau2={p[2]:.1f}+/-{e[2]:.1f} s  p1={p[3]:.2f}+/-{e[3]:.2f}  p2={p[4]:.2f}+/-{e[4]:.2f}')


data = []
for label, key, sec, beta, pmax, color, marker in SAMPLES:
    d = load(key)
    y = (d['g2'][sec, :, QI] - 1) / beta
    s = d['g2_err'][sec, :, QI] / beta
    data.append(dict(label=label, t=d['tau'], y=y, s=s, pmax=pmax, color=color,
                     marker=marker, T=d['section/T_mean'][sec],
                     n=d['section/n_kept'][sec], q=d['q_dynamic'][QI]))

print('As the original (every fit weighted with the (VPAVG)10 errors):')
for x in data:
    print('  ' + row(x['label'], *fit(x['t'], x['y'], data[0]['s'], x['pmax'])))
print('Each fit weighted with its own errors (used below):')
for x in data:
    x['p'], x['e'] = fit(x['t'], x['y'], x['s'], x['pmax'])
    print(f"  {row(x['label'], x['p'], x['e'])}   [{x['T']:.2f} C, {x['n']} acq., Q = {x['q']:.5f}]")

# --- figure -------------------------------------------------------------------
apply_style()
fig, ax = plt.subplots(figsize=(SINGLE_COL, 2.9))
tt = np.logspace(np.log10(2e-5), np.log10(2), 200)
handles = []
for x in data:
    ax.errorbar(x['t'], x['y'], yerr=x['s'], fmt=x['marker'] + '-', color=x['color'],
                mfc='none', ms=MS * 0.8, mew=MEW, lw=LW_THIN, elinewidth=LW_THIN)
    ax.plot(tt, two_step(tt, *x['p']), '-', color=x['color'], lw=LW_DATA)
    handles.append(Line2D([], [], color=x['color'], marker=x['marker'], mfc='none',
                          ms=MS, mew=MEW, lw=LW_THIN,
                          label=f"{x['label']} at {x['T']:.1f} °C"))
ax.set_xscale('log')
ax.set_xlim(1.3e-5, 3)
ax.set_ylim(0, 1.18)
ax.set_xlabel(r'$\Delta t$ (s)')
ax.set_ylabel(r'$g_1^2(Q = 0.0038\ \AA^{-1},\ \Delta t)$')
add_minor_grid(ax)
ax.legend(handles=handles, loc='lower left')
fig.tight_layout(pad=0.4)
save_fig(fig, 'Figure4_Two_Step_Decay.pdf')

# --- table ----------------------------------------------------------------------
lines = [r'% Written by analysis/gelation/two_step_fit.py; do not edit by hand.',
         r'\begin{tabular}{lccccc}', r'\toprule',
         r' & $a_1$ & $\tau_1$ (ms) & $\tau_2$ (s) & $p_1$ & $p_2$ \\', r'\midrule']
for x in data:
    p, e = x['p'], x['e']
    name = x['label'].replace('$_{', '$_{').replace('(VPAVG)', r'(VPAVG)')
    lines.append(f"{name} & ${p[0]:.3f} \\pm {e[0]:.3f}$ & ${p[1] * 1e3:.1f} \\pm {e[1] * 1e3:.1f}$ & "
                 f"${p[2]:.1f} \\pm {e[2]:.1f}$ & ${p[3]:.2f} \\pm {e[3]:.2f}$ & "
                 f"${p[4]:.2f} \\pm {e[4]:.2f}$ \\\\")
lines += [r'\bottomrule', r'\end{tabular}']
with open(os.path.join(HERE, 'Table_Two_Step_Fit.tex'), 'w') as fh:
    fh.write('\n'.join(lines) + '\n')
print('wrote Table_Two_Step_Fit.tex')
