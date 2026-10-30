"""Figure S3 of the 28 Oct version: Guinier analysis of the 10 C Reference aliquot.

The present manuscript draws the same fit as Figure S5b (reversibility_12id.py).

Weighted straight-line fit of ln I(q) against q^2 over 0.020 <= q <= 0.055
A^-1, giving the apparent radius of gyration of the unassembled ELP.  The data
are the merged SAXS+WAXS profile written by Read_12ID_SAWAXS.py into
reduced_data/, so this script needs nothing from the beamline storage.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from common.acs_style import (SINGLE_COL, MS, MEW, LW_THIN, LW_DATA,
                              apply_style, add_minor_grid, save_fig)

# 1. READ DATA
# script-relative, so the script runs correctly from any directory
data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'reduced_data')
df = pd.read_csv(os.path.join(data_dir, 'Merged_Reference_10C.csv'))

# Q is strictly in A^-1
Q = df['Q(A^-1)'].values
I = df['I(Q)'].values
# Intensity uncertainty (propagated in Read_12ID_SAWAXS.py); enables the same
# weighted Guinier fit the 12-ID beamline software performs.
I_err = df['I_err(Q)'].values if 'I_err(Q)' in df.columns else None

# 2-3. WEIGHTED GUINIER FIT over 0.020-0.055 A^-1, the region of the
#      unassembled ELP (~23 A).  The fit itself lives in guinier.py, shared with
#      reversibility_12id.py and gelation/saxs_scaling.py.
from guinier import guinier_fit                              # noqa: E402
fit = guinier_fit()
q_min, q_max = fit['q_min'], fit['q_max']
slope, intercept = fit['slope'], fit['intercept']
red_chi2 = fit['red_chi2']
R_g, R_g_err = fit['R_g'], fit['R_g_err']
I_0, I_0_err = fit['I_0'], fit['I_0_err']
x = Q[(Q >= q_min) & (Q <= q_max) & (I > 0)] ** 2

print("--- Guinier Fit Results (weighted) ---")
print(f"I_0: {I_0:.4f} ± {I_0_err:.4f} (arbs., the units of the merged profile)")
print(f"R_g: {R_g:.3f} ± {R_g_err:.3f} Å")
print(f"Max Q*R_g: {q_max * R_g:.4f} (Target < 1.3)")
print(f"reduced chi^2: {red_chi2:.3f}  (N = {len(x)} points)")

# 4. EXECUTE PLOT
apply_style()

fig, ax = plt.subplots(figsize=(SINGLE_COL, 3.0))

# Plot slightly wider raw data range for context, with propagated error bars
mask_plot = (Q >= 0.005) & (Q <= 0.08) & (I > 0)
xp = Q[mask_plot]**2
yp = np.log(I[mask_plot])
if I_err is not None:
    yperr = I_err[mask_plot] / I[mask_plot]
    ax.errorbar(xp, yp, yerr=yperr, fmt='o', mfc='none', mec='r', ecolor='r',
                markersize=MS, mew=MEW, capsize=1.5, elinewidth=LW_THIN,
                capthick=LW_THIN, ls='none', label='Reference 10 °C')
else:
    ax.plot(xp, yp, 'ro', fillstyle='none', markersize=MS, mew=MEW, label='Reference 10 °C')

# Plot fit line
fit_line = slope * x + intercept
ax.plot(x, fit_line, 'k-', linewidth=LW_DATA, label='Guinier fit')

# Fit-result annotation (mirrors the 12-ID beamline readout).  It sits in the
# upper right, where the data has already decayed away.
txt = (f'$I_0$ = {I_0:.3f} ± {I_0_err:.3f}\n'
       f'$R_g$ = {R_g:.2f} ± {R_g_err:.2f} $\\AA$\n'
       f'$Q_{{max}}$·$R_g$ = {q_max * R_g:.3f}')
ax.text(0.97, 0.95, txt, transform=ax.transAxes, va='top', ha='right',
        # square frame in black, matching the legend box in the same panel
        # (rcParam legend.fancybox is False everywhere in this paper)
        bbox=dict(boxstyle='square,pad=0.35', fc='white', ec='black', lw=LW_THIN))

# Formatting
add_minor_grid(ax)
ax.set_xlabel(r'$Q^2$ ($\AA^{-2}$)')
ax.set_ylabel(r'$\ln[I(Q)]$ (arbs.)')
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0), useMathText=True)
# centre right: the lower-left corner is where the longest low-Q error bars
# reach, so a key there would sit on top of them
ax.legend(loc='center right')

fig.tight_layout(pad=0.4)
save_fig(fig, 'FigureS3_Guinier.pdf')
plt.show()