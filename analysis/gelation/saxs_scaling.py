"""Figure 2: SAXS of (VPAVG)30 through the heating ramp, and the length xi(T).

A port of the four pieces that made Figure 2 of the Overleaf version
(Old/H06: SAXS_Global_Scale, SAXS_Scale_Fit, SAXS_Global_Fit, saxs_tools.py).
Reads data/E0110_n10.h5 and data/D0138.h5 through gelation_data, and the
12-ID-B Guinier radius through ../SAXS_12id/guinier.py.  Writes the figure and
prints every number the text quotes.

THE ANALYSIS, as in the original
1. I(Q, T) of the ten heating sections, buffer subtracted, on the absolute
   scale, positive points with Q <= Q_MAX.
2. Self-similar scaling.  Each curve is moved in log-log by q_scale along Q and
   i_scale along I onto the previous, already moved, curve, minimising the mean
   squared difference of log10 I over their overlap.  The search is a grid:
   q_scale in steps of 0.01 inside a window per curve, i_scale in steps of
   0.001 from 0.001 to 1.5.  The windows are the original ones; the script
   stops if a minimum lands on a window edge.  The two hottest curves do not
   collapse and stay out of the master curve.
3. eta.  The master curve, normalised by its maximum y, is linearised as
   ln(1/y - 1) = eta ln x + c and fitted by least squares.
4. The generalised Lorentzian I = I0 / [1 + (xi Q)^eta], with that eta fixed,
   is fitted to every curve but the hottest (dashed lines).
5. xi(T) = R_g q_scale(T), with R_g the 12-ID-B Guinier radius of the
   dispersed chains, and xi = a (T0 - T)^gamma over the eight collapsing
   curves.  For each trial T0 = max(T) + delta, delta in [0.1, 3] C, ln xi is
   fitted against ln(T0 - T) by least squares; T0 minimises the residual, and
   gamma and its standard error come from the fit at that T0.

WHAT CHANGED FROM THE OVERLEAF FIGURE
* I is in mm^-1 through the per-section calibration of gelation_data, which
  halves it relative to the original's fixed coefficients (computed for a 1 mm
  path; the cell is 2 mm) and corrects the unit, which was labelled cm^-1.  A
  common factor; steps 2-5 do not see it.
* Q is in A^-1 like every other figure here.
* The inset abscissa is Q xi with both in A; the original multiplied Q in A^-1
  by R_g in nm.  Its ordinate is divided by the lowest-Q intensity of the
  reference curve, computed rather than typed in.
* Panel (b) is labelled T0 - T, which is what is plotted, and shows all eight
  points of the fit; the original read T - T0 and its axis stopped at 1 C,
  hiding the 29.6 C point.
* The dotted line on the hottest curve is a least-squares fit, in log I, of
  the generalised Lorentzian with eta free to Q >= 0.008 A^-1.  The original
  drew it with hand-chosen parameters, eta = 2.8.
* R_g = 22.65 +/- 0.50 A is fitted from the committed 12-ID-B data rather than
  typed in from the beamline software (22.8 +/- 0.4 A).  It scales xi and
  moves neither T0 nor gamma.
* Step 5 is unweighted.  The original weighted both fits by a q_scale
  'uncertainty' that was the chi^2 GRADIENT at the grid minimum, a number of
  arbitrary size, and gave the reference curve 1e-6, which pinned both fits to
  it.  Fed the Overleaf inputs, this script's step 2 returns the Overleaf
  q_scale values exactly, and the original weighting returns T0 = 29.96 C and
  gamma = -0.55 +/- 0.02.  But the same weighting gives -0.63 once the four
  corner bins beyond Q_MAX are dropped, and -0.92 with the calibrated buffer
  weight, so that number belongs to the weights rather than to the data.
  Unweighted, every one of those inputs gives T0 = 30.5 C and
  gamma = -0.80 +/- 0.04.
* The error bars in (b) are the R_g uncertainty alone, common to every point,
  so they do not enter the fit.  The q_scale grid step, 0.01, moves xi by
  0.5 % at most.
"""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FixedLocator, FixedFormatter, NullFormatter
from scipy.optimize import curve_fit, minimize

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..'))
sys.path.insert(0, os.path.join(HERE, '..', 'SAXS_12id'))
from common.acs_style import (SINGLE_COL, MS, MEW, LW_THIN, LW_DATA,      # noqa: E402
                              apply_style, add_minor_grid, label_panels,
                              q_log_ticks, save_fig)
from gelation_data import Q_MAX, absolute_iq, load                       # noqa: E402
from guinier import guinier_fit                                          # noqa: E402

# q_scale search window for curve i (the original's, by curve index)
WINDOWS = {1: (0.5, 2.0), 2: (0.5, 2.0), 3: (0.5, 2.0), 4: (0.5, 2.0),
           5: (1.9, 4.2), 6: (1.9, 4.2), 7: (5.5, 7.0), 8: (9.5, 11.0),
           9: (11.0, 13.5)}
I_SCALES = np.arange(0.001, 1.5, 0.001)
N_MASTER = 8                     # curves 0-7 collapse; 8 and 9 do not
HOT_Q_MIN = 0.008                # A^-1, lower end of the eta fit on the hottest curve
MARKERS = ['o', 's', 'D', 'v', '^', '<', '>', 'p', '*', 'h']


def chi2_grid(q_scales, ref_x, ref_y, x, y):
    """chi^2 of every (q_scale, i_scale) pair, as saxs_tools.chi_squared_direct.

    For a given q_scale the overlap and the interpolated reference do not
    depend on i_scale, so the i_scale axis is done in one array operation.  The
    arithmetic is the original's, term for term.
    """
    rx, ry, tx, ty = np.log10(ref_x), np.log10(ref_y), np.log10(x), np.log10(y)
    out = np.full((len(q_scales), len(I_SCALES)), np.inf)
    for n, qs in enumerate(q_scales):
        sx = tx + np.log10(qs)
        lo, hi = max(sx.min(), rx.min()), min(sx.max(), rx.max())
        ov = (sx > lo) & (sx < hi)
        ref_i = np.interp(sx[ov], rx, ry)
        sy = ty[ov][None, :] + np.log10(I_SCALES)[:, None]
        out[n] = np.mean((sy - ref_i[None, :]) ** 2, axis=1)
    return out


def iterative_scaling(q, curves):
    """Step 2: [q_scale, i_scale, chi2] per curve, and the moved curves."""
    results = [[1.0, 1.0, 0.0]]
    moved = [(q, curves[0])]
    for i in range(1, len(curves)):
        lo, hi = WINDOWS[i]
        q_scales = np.arange(lo, hi, 0.01)
        grid = chi2_grid(q_scales, *moved[-1], q, curves[i])
        n, m = np.unravel_index(np.argmin(grid), grid.shape)
        if n in (0, len(q_scales) - 1) or m in (0, len(I_SCALES) - 1):
            # the two hottest curves never collapse, and in the original the
            # last one already sat on the i_scale floor; only the master curve
            # has to be free of the grid's edges
            if i < N_MASTER:
                sys.exit(f'curve {i}: minimum on the edge of its search window')
            print(f'  (curve {i} is on the edge of its search window; it is not '
                  f'part of the master curve)')
        qs, iscale = q_scales[n], I_SCALES[m]
        results.append([qs, iscale, grid[n, m]])
        moved.append((q * qs, curves[i] * iscale))
    return np.array(results), moved


def lorentzian(q, i0, xi, eta):
    return i0 / (1 + (xi * q) ** eta)


def power_law(T, xi):
    """Step 5: (T0, gamma, gamma error, ln a) of ln xi = ln a + gamma ln(T0 - T).

    T0 = max(T) + delta with delta in [0.1, 3] C, chosen to minimise the
    residual sum of squares of the straight-line fit in log-log.
    """
    def rss(delta):
        x = np.log(np.max(T) + delta[0] - T)
        c = np.polyfit(x, np.log(xi), 1)
        return np.sum((np.log(xi) - np.polyval(c, x)) ** 2)

    res = minimize(rss, [1.0], bounds=[(0.1, 3)])
    T0 = np.max(T) + res.x[0]
    x, y = np.log(T0 - T), np.log(xi)
    A = np.vstack([x, np.ones_like(x)]).T
    coef = np.linalg.lstsq(A, y, rcond=None)[0]
    s2 = np.sum((y - A @ coef) ** 2) / (len(y) - 2)
    cov = s2 * np.linalg.inv(A.T @ A)
    return T0, coef[0], np.sqrt(cov[0, 0]), coef[1]


# --- data --------------------------------------------------------------------
d = load('E0110_n10')
T = d['section/T_mean']
q_all, I_all = absolute_iq('E0110_n10', 'D0138')
keep = (q_all <= Q_MAX) & np.all(I_all > 0, axis=0)
q, I = q_all[keep], I_all[:, keep]
print(f'{I.shape[0]} heating curves, {T.min():.1f}-{T.max():.1f} C, '
      f'{keep.sum()} Q points {q.min():.4f}-{q.max():.4f} A^-1')

# --- step 2: self-similar scaling ----------------------------------------------
res, moved = iterative_scaling(q, I)
for i, (qs, iscale, c2) in enumerate(res):
    print(f'  curve {i} ({T[i]:.1f} C): q_scale={qs:.2f}  i_scale={iscale:.3f}  chi2={c2:.2e}')

# --- step 3: eta from the master curve ----------------------------------------
mx = np.concatenate([m[0] for m in moved[:N_MASTER]])
my = np.concatenate([m[1] for m in moved[:N_MASTER]])
yn = np.clip(my / my.max(), 1e-10, 1 - 1e-10)
eta = np.polyfit(np.log(mx), np.log(1 / yn - 1), 1)[0]
print(f'eta (master curve) = {eta:.4f}')

# --- step 4: fixed-eta Lorentzian per curve ---------------------------------
fits = []
for i in range(len(I) - 1):
    p, _ = curve_fit(lambda x, i0, xi: lorentzian(x, i0, xi, eta), q, I[i],
                     p0=[I[i].max(), 1], bounds=([0, 0], [np.inf, np.inf]))
    fits.append(p)

# the hottest curve at Q >= HOT_Q_MIN, eta free, least squares in log I
hot = q >= HOT_Q_MIN
p_hot, c_hot = curve_fit(lambda x, a, b, c: np.log(lorentzian(x, a, b, c)), q[hot],
                         np.log(I[-1][hot]), p0=[I[-1][hot][0] * 2, 200, 2.8],
                         bounds=([0, 0, 0], [np.inf, np.inf, 10]))
eta_hot, eta_hot_err = p_hot[2], np.sqrt(c_hot[2, 2])
print(f'eta at {T[-1]:.1f} C, Q >= {HOT_Q_MIN} A^-1: {eta_hot:.2f} +/- {eta_hot_err:.2f}')

# --- step 5: xi(T) and its power law -----------------------------------------
g = guinier_fit()
Rg, Rg_err = g['R_g'], g['R_g_err']                         # A
xs = res[:N_MASTER, 0]
Tm = T[:N_MASTER]
T0, gamma, gamma_err, ln_a = power_law(Tm, xs)
dT = T0 - Tm
xi_nm = Rg * xs / 10
xi_err_nm = Rg_err / Rg * xi_nm
print(f'R_g = {Rg:.2f} +/- {Rg_err:.2f} A;  T0 = {T0:.2f} C;  gamma = {gamma:.2f} +/- {gamma_err:.2f}')
print('  xi(T) nm: ' + ', '.join(f'{t:.1f} C {x:.2f}' for t, x in zip(Tm, xi_nm)))

# --- figure ----------------------------------------------------------------------
apply_style()
fig = plt.figure(figsize=(SINGLE_COL, 5.2))
gs = fig.add_gridspec(3, 1, height_ratios=[3.1, 1.25, 0.62], hspace=0.42)
ax, axb, axl = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[2])
axl.axis('off')
colors = [plt.cm.jet(i / len(I)) for i in range(len(I))]

for i in range(len(I)):
    ax.plot(q[::2], I[i][::2], MARKERS[i], color=colors[i], mfc='none', ms=MS * 0.8,
            mew=LW_THIN, label=f'{T[i]:.1f} °C')
for i, p in enumerate(fits):
    ax.plot(q, lorentzian(q, *p, eta), '--', color='k', lw=LW_THIN)
ax.plot(q[hot], lorentzian(q[hot], *p_hot), ':', color='k', lw=LW_DATA * 1.5)
ax.set_xscale('log')
ax.set_yscale('log')
q_log_ticks(ax)
ax.set_xlabel(r'$Q$ ($\AA^{-1}$)')
ax.set_ylabel(r'$I(Q)$ (mm$^{-1}$)')
add_minor_grid(ax)
# headroom above the hottest curve, where the inset sits clear of the data
ax.set_ylim(0.02, 400)

# inset: the master curve in reduced units
ins = ax.inset_axes([0.60, 0.635, 0.37, 0.335])
i_ref = I[0][0]
for i in range(N_MASTER):
    ins.plot(moved[i][0] * Rg, moved[i][1] / i_ref, 'o', color=colors[i], ms=1.0, mew=0)
ins.set_xscale('log')
ins.set_yscale('log')
ins.set_xlabel(r'$Q\,\xi(T)$', labelpad=1)
ins.set_ylabel(r'$I/I_{0}(T)$', labelpad=1)
ins.tick_params(pad=1.5)
add_minor_grid(ins)
ins.patch.set_alpha(1.0)

# (b) xi against T0 - T
for i in range(N_MASTER):
    axb.errorbar(dT[i], xi_nm[i], yerr=xi_err_nm[i], fmt=MARKERS[i], color=colors[i],
                 mfc='none', ms=MS, mew=MEW, elinewidth=LW_THIN, capsize=1.5,
                 capthick=LW_THIN)
xx = np.logspace(np.log10(dT.min() * 0.8), np.log10(dT.max() * 1.2), 100)
axb.plot(xx, Rg / 10 * np.exp(ln_a) * xx ** gamma, '--', color='k', lw=LW_THIN)
axb.set_xscale('log')
axb.set_yscale('log')
axb.set_xlim(0.2, 20)
axb.set_ylim(1.5, 30)
axb.xaxis.set_major_locator(FixedLocator([0.3, 1, 3, 10]))
axb.xaxis.set_major_formatter(FixedFormatter(['0.3', '1', '3', '10']))
axb.yaxis.set_major_locator(FixedLocator([2, 5, 10, 20]))
axb.yaxis.set_major_formatter(FixedFormatter(['2', '5', '10', '20']))
axb.yaxis.set_minor_formatter(NullFormatter())
axb.set_xlabel(r'$T_0 - T$ (°C)')
axb.set_ylabel(r'$\xi(T)$ (nm)')
add_minor_grid(axb)
axb.text(0.97, 0.92, f'$T_0$ = {T0:.2f} °C\n$\\gamma$ = {gamma:.2f} $\\pm$ {gamma_err:.2f}',
         transform=axb.transAxes, ha='right', va='top',
         bbox=dict(boxstyle='square,pad=0.3', fc='white', ec='black', lw=LW_THIN))

h, lab = ax.get_legend_handles_labels()
axl.legend(h, lab, ncol=5, loc='center', handlelength=1.0, columnspacing=0.6,
           handletextpad=0.3, markerscale=1.2)
label_panels([ax, axb], dx=-34)
fig.subplots_adjust(left=0.17, right=0.97, top=0.97, bottom=0.03)
save_fig(fig, 'Figure2_SAXS_Scaling.pdf')
