"""Figure S5: reversibility, cross-beta-like peaks and Guinier radius at 12-ID-B.

(a) The merged SAXS+WAXS profiles of two (VPAVG)30 aliquots at 10 and 30 C.
    The Measurement aliquot had been through more than ten heating/cooling
    cycles, the Reference aliquot had stayed at 4 C.  The dashed line is the
    power law fitted to the Reference at 30 C over 0.012-0.040 A^-1; the inset
    shows the WAXS region of that profile with the two Gaussian-plus-background
    fits of saxs_waxs_fits.py.
(b) Guinier plot of the Reference aliquot at 10 C, with the weighted fit of
    guinier.py.

Reads only reduced_data/, written by Read_12ID_SAWAXS.py, so it runs from a
clone.  Panel (a) is the single-column figure of Plot_12ID.py (Figure 1 of the
28 Oct version) drawn as one panel of this one, with the same fits.  The
Overleaf version of panel (b) was a screenshot of the beamline software's fit
(R_g = 22.754 +/- 0.418 A); this panel redoes the fit from the committed data
with the same Q range, which gives 22.65 +/- 0.50 A.
"""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FixedFormatter, FixedLocator

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..'))
from common.acs_style import (DOUBLE_COL, MS, MS_DENSE, MEW, LW_THIN, LW_DATA,  # noqa: E402
                              apply_style, add_minor_grid, label_panels, save_fig)
from guinier import guinier_fit, load                                          # noqa: E402
from saxs_waxs_fits import POWER_LAW_RANGE, gauss_bg, power_law, waxs_peaks    # noqa: E402

REF30 = 'Merged_Reference_30C.csv'
# The two aliquots overlap at both temperatures, which is the result, so the
# Reference is drawn as a line through the Measurement markers and neither
# hides the other.
CURVES = [('Merged_Reference_10C.csv',   'r', None, 'Reference, 10 °C'),
          ('Merged_Measurement_10C.csv', 'b', 's',  'Measurement, 10 °C'),
          (REF30,                        'k', None, 'Reference, 30 °C'),
          ('Merged_Measurement_30C.csv', 'g', 'D',  'Measurement, 30 °C')]

fit = guinier_fit()
print(f"Guinier: R_g = {fit['R_g']:.2f} +/- {fit['R_g_err']:.2f} A, "
      f"I_0 = {fit['I_0']:.3f} +/- {fit['I_0_err']:.3f}, "
      f"Q_max R_g = {fit['q_max'] * fit['R_g']:.3f}, reduced chi2 = {fit['red_chi2']:.3f}")
for name in (REF30, 'Merged_Measurement_30C.csv'):
    slope, _ = power_law(name)
    pk = waxs_peaks(name)
    print(f"{name}: power law Q^{slope:.2f}; peaks "
          + ', '.join(f"{x['center']:.3f} A^-1 (d = {x['d']:.2f} A, FWHM {x['fwhm']:.3f}, "
                      f"2pi/FWHM = {2 * np.pi / x['fwhm']:.0f} A)" for x in pk))

apply_style()
fig, (ax, axg) = plt.subplots(1, 2, figsize=(DOUBLE_COL, 3.2))

# (a) colour and marker shape both encode the condition
for name, color, marker, label in CURVES:
    Q, I, _ = load(name)
    ok = I > 0
    if marker is None:
        ax.plot(Q[ok], I[ok], '-', color=color, lw=LW_DATA, label=label, zorder=3)
    else:
        ax.plot(Q[ok], I[ok], ls='none', marker=marker, color=color, mfc='none',
                ms=MS_DENSE, mew=LW_THIN, label=label)
slope, intercept = power_law(REF30)
q_line = np.linspace(*POWER_LAW_RANGE, 50)
ax.plot(q_line, 10 ** intercept * q_line ** slope, 'r--', lw=LW_DATA, zorder=5)
ax.text(POWER_LAW_RANGE[0], 10 ** intercept * POWER_LAW_RANGE[0] ** slope * 1.8,
        rf'~$Q^{{{slope:.1f}}}$', ha='left', va='bottom')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(3.5e-3, 2.4)
ax.set_ylim(7e-3, 5e3)
ax.set_xlabel(r'$Q$ ($\AA^{-1}$)')
ax.set_ylabel(r'$I(Q)$ (arbs.)')
add_minor_grid(ax)
ax.legend(loc='lower left', markerscale=1.6, handlelength=1.2, handletextpad=0.3)

# inset: the WAXS region of the Reference at 30 C, with the two peak fits
ins = ax.inset_axes([0.50, 0.52, 0.48, 0.46])
Q, I, _ = load(REF30)
waxs = Q > 0.22
ins.plot(Q[waxs], I[waxs], 'o', color='k', mfc='none', ms=MS_DENSE, mew=LW_THIN)
for x in waxs_peaks(REF30):
    ins.plot(x['q'], gauss_bg(x['q'], *x['p']), 'r-', lw=LW_DATA)
    ins.axvline(x['center'], ymax=0.78, color='r', ls='--', lw=LW_THIN)
    ins.text(x['center'], 0.97, f"{x['center']:.2f} " + r'$\AA^{-1}$',
             transform=ins.get_xaxis_transform(), color='r', ha='center', va='top')
ins.set_xlim(0.2, 1.09 * Q[waxs].max())
ins.set_ylim(I[waxs].min() * 0.96, I[waxs].max() * 1.18)
ins.set_xlabel(r'$Q$ ($\AA^{-1}$)', labelpad=1.5)
ins.set_ylabel(r'$I(Q)$', labelpad=1.5)
ins.xaxis.set_major_locator(FixedLocator([0.5, 1.0, 1.5]))
ins.yaxis.set_major_locator(FixedLocator([0.02, 0.03, 0.04]))
ins.yaxis.set_major_formatter(FixedFormatter(['0.02', '0.03', '0.04']))
ins.tick_params(pad=1.5)
add_minor_grid(ins)

# (b) ln I against Q^2 a little beyond the fitted range, with the fit over it
Q, I, I_err = load()
show = (Q <= 0.075) & (I > 0)
axg.errorbar(Q[show] ** 2, np.log(I[show]), yerr=I_err[show] / I[show], fmt='o',
             color='r', mfc='none', ms=MS, mew=MEW, elinewidth=LW_THIN,
             capsize=1.5, capthick=LW_THIN, label='Reference, 10 °C')
x = np.linspace(fit['q_min'], fit['q_max'], 50) ** 2
axg.plot(x, fit['slope'] * x + fit['intercept'], 'k-', lw=LW_DATA, label='Guinier fit')
axg.text(0.97, 0.95,
         f"$R_g$ = {fit['R_g']:.2f} $\\pm$ {fit['R_g_err']:.2f} $\\AA$\n"
         f"$I_0$ = {fit['I_0']:.3f} $\\pm$ {fit['I_0_err']:.3f}\n"
         f"$Q_{{\\max}}R_g$ = {fit['q_max'] * fit['R_g']:.3f}",
         transform=axg.transAxes, va='top', ha='right',
         bbox=dict(boxstyle='square,pad=0.35', fc='white', ec='black', lw=LW_THIN))
axg.set_xlabel(r'$Q^2$ ($\AA^{-2}$)')
axg.set_ylabel(r'$\ln I(Q)$ (arbs.)')
axg.ticklabel_format(style='sci', axis='x', scilimits=(0, 0), useMathText=True)
add_minor_grid(axg)
axg.legend(loc='lower left')

label_panels([ax, axg])
fig.tight_layout(pad=0.4, w_pad=1.6)
save_fig(fig, 'FigureS5_12ID_Reversibility.pdf')
