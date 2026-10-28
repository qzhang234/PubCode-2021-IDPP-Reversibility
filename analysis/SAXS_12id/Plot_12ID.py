"""Figure 2: temperature-dependent SAXS/WAXS of (VPAVG)30, with a WAXS inset.

Single-column ACS figure (3.33 in).  Layout is unchanged from the original --
one log-log I(Q) panel with the WAXS peak region inset at the top right -- but
everything is drawn through common/acs_style.py, and three things were moved so
that no text sits on top of data at 8 pt:

  * the main y-axis floor is set just under the data and the four-entry key is
    made narrow enough to fit in the clear lower-left corner (see the note at
    set_ylim);
  * the inset is given headroom above its highest point, and the two peak
    labels sit in that clear strip instead of over the rising edge at 0.22 A^-1
    (the real-space spacings they used to carry are already in the caption);
  * each condition gets its own marker shape, so the four curves are still
    separable without colour (ACS: "avoid relying on color alone").

None of the fitted ranges, models or numbers changed.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter
from scipy.optimize import curve_fit
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from common.acs_style import (SINGLE_COL, MS_DENSE, LW_THIN, LW_DATA,
                              apply_style, add_minor_grid, save_fig)

# 1. READ GENERATED CSV FILES
# script-relative, so the script runs correctly from any directory
data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'reduced_data')
df_ref10 = pd.read_csv(os.path.join(data_dir, 'Merged_Reference_10C.csv'))
df_ref30 = pd.read_csv(os.path.join(data_dir, 'Merged_Reference_30C.csv'))
df_meas10 = pd.read_csv(os.path.join(data_dir, 'Merged_Measurement_10C.csv'))
df_meas30 = pd.read_csv(os.path.join(data_dir, 'Merged_Measurement_30C.csv'))

def gauss_bg(x, a, x0, sigma, m, c):
    """A Gaussian peak of height a, centre x0 and width sigma, on the sloping
    background m*x + c that the wing of the neighbouring scattering leaves under
    each WAXS peak.

    curve_fit returns the fitted parameters in this same order, so the peak
    centre is element [1] and its width is element [2] of the result below.
    """
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2)) + m * x + c

# 2. EXECUTE PLOT
apply_style()

fig, ax = plt.subplots(figsize=(SINGLE_COL, 3.15))

# Main plot.  Colour AND marker shape encode the condition.
CURVES = [(df_ref10,  'r', 'o', 'Reference, 10 °C'),
          (df_ref30,  'k', 's', 'Reference, 30 °C'),
          (df_meas10, 'g', '^', 'Measurement, 10 °C'),
          (df_meas30, 'b', 'D', 'Measurement, 30 °C')]
for df, color, marker, label in CURVES:
    # MS_DENSE + a 0.5 pt edge: ~500 points per curve over four overlapping
    # curves, so the default marker merges into a solid band.
    ax.plot(df['Q(A^-1)'].values, df['I(Q)'].values, color=color, marker=marker,
            fillstyle='none', markersize=MS_DENSE, mew=LW_THIN,
            linestyle='none', label=label)

# Formatting main plot
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$Q$ ($\AA^{-1}$)')
ax.set_ylabel(r'$I(Q)$ (arbs.)')
add_minor_grid(ax)
# The key sits in the empty bottom-left corner, so the floor only has to clear
# the key itself, not a whole spare decade.  What caps how high the key can ride
# is its WIDTH: the lowest curve is at I ~ 1.6 for Q < 0.02 but has fallen to
# ~0.14 by Q = 0.15, so a wide box has to stay low.  Short handles keep the box
# narrow (its right edge lands near Q = 0.11, where the ceiling is still ~0.3),
# which lets the floor come up to just under the data and closes the gap between
# the key and the curves.  Data then fills ~95 % of the axes instead of ~78 %.
ax.set_ylim(7e-3, 5e3)
ax.set_xlim(3.5e-3, 2.4)
ax.legend(loc='lower left', borderaxespad=0.4,
          handlelength=0.8, handletextpad=0.25)

# --- ADD Q-SCALING FOR REFERENCE 30C ---
q_ref30_main = df_ref30['Q(A^-1)'].values
i_ref30_main = df_ref30['I(Q)'].values

# Define the low-Q fitting range for 3D domain scaling
q_min_scale, q_max_scale = 0.012, 0.04
mask_scale = (q_ref30_main > q_min_scale) & (q_ref30_main < q_max_scale)

# Fit the linear relationship in log-log space
log_q = np.log10(q_ref30_main[mask_scale])
log_i = np.log10(i_ref30_main[mask_scale])
slope, intercept = np.polyfit(log_q, log_i, 1)

# Generate points for the trendline directly through the data
q_line = np.linspace(q_min_scale, q_max_scale, 50)
i_line = (10**intercept) * (q_line**slope)

# Plot the scaling line.  Red, not black: the fit runs straight through the
# 30 C data it was fitted to, and in black it was indistinguishable from the
# markers at print size.  Red is already this figure's colour for a fit -- the
# inset uses it for the two Gaussian peak fits.
ax.plot(q_line, i_line, 'r--', lw=LW_DATA, zorder=5)

# Add the scaling label at the top-left start of the fit line (offset vertically)
q_start = q_min_scale
i_start = (10**intercept) * (q_start**slope) * 1.8
ax.text(q_start, i_start, rf'~$Q^{{{slope:.1f}}}$', color='k', ha='left', va='bottom')
# ---------------------------------------

# 3. ADD INSET FIGURE (Top Right)
ax_ins = ax.inset_axes([0.50, 0.52, 0.48, 0.46])

q_ref30 = df_ref30['Q(A^-1)'].values
i_ref30 = df_ref30['I(Q)'].values

# Plot inset data
mask = q_ref30 > 0.22
ax_ins.plot(q_ref30[mask], i_ref30[mask], color='k', fillstyle='none', marker='o',
            markersize=MS_DENSE, mew=LW_THIN, linestyle='none')

print("\n--- Gaussian Fit Results for Reference 30C WAXS Peaks ---")

# --- FIT PEAK 1 (~0.74 A^-1) ---
m1 = (q_ref30 > 0.48) & (q_ref30 < 1.00)
q1, i1 = q_ref30[m1], i_ref30[m1]
p1_opt, _ = curve_fit(gauss_bg, q1, i1, p0=[np.ptp(i1), 0.74, 0.15, 0, np.min(i1)])

print(f"\nPeak 1 (Inter-sheet spacing):")
print(f"  Position (Q0): {p1_opt[1]:.4f} A^-1  (d = {2*np.pi/p1_opt[1]:.2f} A)")
print(f"  Width (Sigma): {p1_opt[2]:.4f} A^-1")
# FWHM = 2*sqrt(2*ln 2)*sigma = 2.355*sigma for a Gaussian
print(f"  Width (FWHM):  {2.355 * p1_opt[2]:.4f} A^-1")

ax_ins.plot(q1, gauss_bg(q1, *p1_opt), 'r-', lw=LW_DATA)

# --- FIT PEAK 2 (~1.44 A^-1) ---
m2 = (q_ref30 > 1.18) & (q_ref30 < 1.70)
q2, i2 = q_ref30[m2], i_ref30[m2]
p2_opt, _ = curve_fit(gauss_bg, q2, i2, p0=[np.ptp(i2), 1.44, 0.15, 0, np.min(i2)])

print(f"\nPeak 2 (Inter-strand hydrogen bonding):")
print(f"  Position (Q0): {p2_opt[1]:.4f} A^-1  (d = {2*np.pi/p2_opt[1]:.2f} A)")
print(f"  Width (Sigma): {p2_opt[2]:.4f} A^-1")
print(f"  Width (FWHM):  {2.355 * p2_opt[2]:.4f} A^-1\n")

ax_ins.plot(q2, gauss_bg(q2, *p2_opt), 'r-', lw=LW_DATA)

# Format inset.  A narrow strip of headroom above the highest point carries the
# two peak labels; each label is centred on its peak.  The dashed marker line
# stops short of that strip (ymax is an axes fraction) rather than running the
# full height, so the line points at the peak without striking through the label
# it belongs to.
i_ins = i_ref30[mask]
# a little margin past the last point: keeps it off the right spine, and
# gives the 1.38 label room to stay inside the frame
ax_ins.set_xlim(0.2, 1.09 * max(q_ref30[mask]))
ax_ins.set_ylim(i_ins.min() * 0.96, i_ins.max() * 1.18)
for p_opt in (p1_opt, p2_opt):
    ax_ins.axvline(p_opt[1], ymax=0.78, color='r', linestyle='--', lw=LW_THIN)
    ax_ins.text(p_opt[1], 0.97, f'{p_opt[1]:.2f} ' + r'$\AA^{-1}$',
                transform=ax_ins.get_xaxis_transform(), color='r',
                ha='center', va='top')
ax_ins.set_xlabel(r'$Q$ ($\AA^{-1}$)', labelpad=1.5)
ax_ins.set_ylabel(r'$I(Q)$', labelpad=1.5)
ax_ins.xaxis.set_major_locator(FixedLocator([0.5, 1.0, 1.5]))
ax_ins.yaxis.set_major_locator(FixedLocator([0.02, 0.03, 0.04]))
ax_ins.yaxis.set_major_formatter(FixedFormatter(['0.02', '0.03', '0.04']))
ax_ins.tick_params(pad=1.5)
add_minor_grid(ax_ins)

fig.tight_layout(pad=0.4)
save_fig(fig, 'Figure2_SAXS_WAXS.pdf')
plt.show()
