"""Combined SAXS + XPCS analysis, reading the averaged HDF files in data/.

This one script writes three of the paper's figures.

Figure 3, three panels:
  (a) SAXS I(q) for B0146, the first B0147 file (frames 1-200) and the last five
      (frames 801-1313), with the D0138 buffer subtracted.
  (b) g2(tau) at one q for those same files, with the fits.
  (c) fitted fast fraction f against elapsed time, one series per q bin.
Figure S8, 2x2, the fit parameters: (a) the shared exponents p1 and p2 against
elapsed time, (b) tau_fast against Q, (c) the power-law exponents gamma_fast and
gamma_slow from fitting tau = A * Q**gamma, and (d) tau_slow against Q.
Figure S6, 2 panels, documents the absolute-scale calibration: the ion-chamber
to photon fit and the air transmission.  See abs_xsec.py.

The SAXS panel is on an absolute scale, with the coefficient computed PER FILE
from that file's own range-averaged ion-chamber readings, so a drift in incident
flux across the series is handled correctly.  The axis carries the units, mm^-1.

Colour encodes elapsed time and means the same thing in all three panels of
Figure 3, so ONE key serves the whole figure and sits in a row along the bottom
rather than inside any panel; that way no panel needs headroom for it.  Marker
SHAPE encodes the q bin, and its key stays inside panel (c).  B0146 is a 6 C
reference taken before the isothermal run, so it is not one of the elapsed times
and is keyed separately in panel (a), beside its own curve.

Geometry follows common/acs_style.py: double column, 7.0 in, 8 pt Arial.

FIT MODEL.  A double stretched exponential in the Siegert form,

    g2 = contrast * ( f e^-(tau/tau_fast)^p1 + (1-f) e^-(tau/tau_slow)^p2 )^2 + 1

with the contrast fixed at beta = 0.13139, measured on a static reference by
contrast_calibration.py (Figure S7), and the baseline fixed at 1.  At each
elapsed time all five q bins are fitted SIMULTANEOUSLY: p1 and p2 are shared
across q, while tau_fast, f and tau_slow are free in each bin.  The shared
exponents are the physically motivated constraint that relaxation shape is a
property of the sample state at a given age rather than of the q bin, and they
stabilise the slow mode, whose tau lies beyond the 1.6 s delay window and so
cannot be pinned bin by bin.

UNCERTAINTIES.  The fit minimises error-weighted residuals using the g2_err
stored in the averaged files (absolute_sigma convention).  The covariance is the
pseudo-inverse of J^T J at the solution, and the 1-sigma errors are the square
roots of its diagonal; see the note in xpcs_fit.py for why a pseudo-inverse.
Sharing p1 and p2 propagates their uncertainty into f, so the f errors here are
larger than a per-q fixed-exponent fit would report.  That difference is the
exponent systematic, folded in honestly.
"""

import glob
import os
import sys

import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import (FixedLocator, FixedFormatter, NullFormatter,
                               FuncFormatter, LogLocator)

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from common.acs_style import (DOUBLE_COL, MS, MS_SPARSE, MEW, LW_THIN, LW_DATA,
                              apply_style, add_minor_grid, label_panels,
                              q_log_ticks, save_fig)

# --- PARAMETERS ---
data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

# --- ACS FIGURE GEOMETRY (see common/acs_style.py) ---
# All three have two or more panels, so all three are double-column figures at
# the 7.0 in (504 pt) ACS maximum.  Depths are far inside the 9.167 in cap.
FIG1_SIZE = (DOUBLE_COL, 3.15)        # 3 panels + the shared key row beneath them
FIG1_LEGEND_H = 0.075                # fraction of the height reserved at the
                                     # bottom for that shared key
FIG2_SIZE = (DOUBLE_COL, 2.6)        # 1x3 fit-parameter panels
FIG3_SIZE = (DOUBLE_COL, 2.9)        # 2 calibration panels
# Panel (a) carries no key of its own except the one-line 6 C reference in the
# bottom-left corner, so its limits sit just clear of the data: a little over a
# third of a decade below, half a decade above.
SAXS_YLO_PAD, SAXS_YHI_PAD = 0.45, 3.0

# --- DATASET SELECTION ---
BACKGROUND_HEADER = 'D0138'          # buffer, subtracted as background in the SAXS panel
XPCS_HEADER = 'B0147'                # header used for the g2 / fit panels
# The absolute cross-section coefficients coef_sam / coef_buf are no longer
# hard-coded: they are computed per file from that file's own (range-averaged)
# ion-chamber readings by abs_xsec_coef(), defined in the ABSOLUTE SCATTERING
# CROSS-SECTION section below.
N_LAST = 5                           # number of last (highest-frame) files
target_q_idx = 0                     # q bin shown in the g2 panel
fit_q_indices = [0, 1, 2, 3, 4]      # q bins fitted for the f-vs-time panel
Q_MARKERS = ['o', 's', '^', 'D', 'v']   # one marker per fit q bin
CHI2_MAX = 10.0                      # skip fits worse than this reduced chi^2.
                                     # A guard only: the five fitted files give
                                     # 0.78-1.89, so it never fires.  The liquid
                                     # t=0 file, which it was written for, is not
                                     # in the fit set at all (chi^2 ~ 2e3 if it
                                     # were: it does not obey the arrested model)
PHI_AVERAGE = True                   # azimuthally average I(q) over phi sectors
# XPCS elapsed-time colours: the N_LAST times are mapped to evenly-spaced
# positions of this colormap (maximises separation so times are easy to tell
# apart), independent of the SAXS-only SI figure.  The 0 s and 6 C reference
# datasets get their own distinct colours (below), outside this scale.
XPCS_CMAP = plt.cm.plasma
CMAP_LO, CMAP_HI = 0.10, 0.88        # colormap span used for the N_LAST times
COLOR_0S = 'black'                   # first B0147 file (elapsed = 0 s)
COLOR_6C = '#1f77b4'                 # B0146 (6 C reference, before isothermal)

# --- FIT MODEL ---
# The model, the measured contrast and the global fit live in xpcs_fit.py so
# that Figure 3b, Figure S8 and Figure S9 provably share one implementation.
from xpcs_fit import (double_exp, fit_g2_joint,   # noqa: E402
                      mean_tau)


# ============================================================================
# ABSOLUTE SCATTERING CROSS-SECTION
# ============================================================================
# The calibration constants, the IC->photon linear fit, the air transmission
# and the per-file coefficient function abs_xsec_coef() all live in the shared
# module abs_xsec.py.  saxpcs.py and saxs_evolution.py both import from it so
# the absolute scale is computed identically in every figure.  See abs_xsec.py
# for the full derivation, its provenance and the unit notes.
from abs_xsec import (                                     # noqa: E402
    CAL_A, CAL_B, CAL_UPIC, CAL_PHOTONS, CAL_CROP,
    AIR_TRANS_SERIES, AIR_TRANSMISSION, AIR_TRANS_STD,
    abs_xsec_coef, calibration_summary,
)

print(calibration_summary())


# --- READING THE AVERAGED HDF FILES ---
# The field locations and the four readers are shared with the other three
# 8-ID figure scripts, so all four provably read the files the same way.
from nexus_read import (parse_name, read_start_time,   # noqa: E402
                        read_saxs_iq, read_g2)


def fit_powerlaw(Q, tau, tau_err):
    """Weighted log-log fit of tau = A * Q**gamma at one elapsed time.

    Minimises the error-weighted residual of ln(tau) vs ln(Q) (a weighted
    linear regression), with sigma_ln(tau) = tau_err / tau (propagated from
    the tau uncertainty). Returns (gamma, gamma_err, A), with gamma_err from
    the weighted-least-squares covariance (absolute_sigma convention) and A
    the fitted prefactor (for overlaying the fit line).
    """
    x = np.log(np.asarray(Q, dtype=float))
    y = np.log(np.asarray(tau, dtype=float))
    sigma_y = np.asarray(tau_err, dtype=float) / np.asarray(tau, dtype=float)
    w = 1.0 / sigma_y**2
    X = np.vstack([x, np.ones_like(x)]).T
    cov = np.linalg.inv(X.T @ (w[:, None] * X))
    gamma, intercept = cov @ (X.T @ (w * y))
    gamma_err = np.sqrt(cov[0, 0])
    return gamma, gamma_err, np.exp(intercept)


# --- DISCOVER FILES ---
# every range average in data/, including the thermal-cycle groups that belong to
# Figure S4; the headers wanted here are selected out of by_header below
file_paths = sorted(glob.glob(os.path.join(data_dir, 'Average_*.hdf')))
assert file_paths, f'no HDF files found in {data_dir}'
by_header, start_times = {}, {}
for fp in file_paths:
    header = parse_name(fp)[0]
    with h5py.File(fp, 'r') as hf:
        start_times[fp] = read_start_time(hf)
    by_header.setdefault(header, []).append(fp)
for h in by_header:
    by_header[h].sort(key=lambda p: parse_name(p)[1])

b0147 = by_header[XPCS_HEADER]
first_file = b0147[0]                          # first B0147 file (frames 1-200)
xpcs_files = b0147[-N_LAST:]                    # g2 / fit: last N only (NOT 1-200)
saxs_files = [first_file] + xpcs_files          # SAXS panel: first + last N
t_ref = start_times[first_file]                 # time origin = first B0147 file


def elapsed(fp):
    return (start_times[fp] - t_ref).total_seconds()


# XPCS elapsed-time colours: map the N_LAST files (sorted by time) to evenly
# spaced colormap positions so consecutive times are maximally distinguishable.
# The same mapping is reused in every XPCS panel (and in g2_grid_SI.py).
_xpcs_sorted = sorted(xpcs_files, key=elapsed)
_xpcs_pos = np.linspace(CMAP_LO, CMAP_HI, len(_xpcs_sorted))
_xpcs_color = {fp: XPCS_CMAP(p) for fp, p in zip(_xpcs_sorted, _xpcs_pos)}


def ecolor(fp):
    """Colour for an XPCS elapsed time; 0 s and the 6 C ref use distinct colours."""
    if fp == first_file:
        return COLOR_0S
    return _xpcs_color[fp]


# elapsed-time value -> colour (xpcs_files only; used where fp isn't in scope)
time_color = {elapsed(fp): ecolor(fp) for fp in xpcs_files}


print('SAXS files:', [f'{parse_name(f)[1]}-{parse_name(f)[2]} ({elapsed(f):.0f} s)' for f in saxs_files])
print('XPCS files:', [f'{parse_name(f)[1]}-{parse_name(f)[2]} ({elapsed(f):.0f} s)' for f in xpcs_files])

# --- GLOBAL FITS (one per elapsed time; p1, p2 shared across q) ---
# Fit each XPCS file once and reuse the result in every panel.
g2_data = {}   # fp -> (tau, g2, g2_err, q_vals)
for fp in xpcs_files:
    with h5py.File(fp, 'r') as hf:
        g2_data[fp] = read_g2(hf)

# One fit of the whole series: p_fast is shared across every elapsed time so
# that relaxation times ARE comparable between them (see fit_g2_joint); p_slow
# stays free at each time.  tau_fast is NOT free in each q bin -- Equation 2
# ties it to 1/(D0 q^2), so the q^-2 scaling is imposed, not fitted.  What
# justifies imposing it is the free fit, whose q exponents are -2.08 to -2.37.
joint, _per_time = fit_g2_joint([g2_data[fp][:3] for fp in xpcs_files],
                                [g2_data[fp][3] for fp in xpcs_files],
                                fit_q_indices)
fits = dict(zip(xpcs_files, _per_time))
MEAN_TAU_FACTOR = mean_tau(1.0, joint['p1'])   # <tau>/tau_SE, a constant here
print(f"  joint fit: {joint['n_par']} par, chi2/dof = {joint['red_chi2']:.3f}, "
      f"shared p_fast = {joint['p1']:.3f}+/-{joint['p1_err']:.3f}, "
      f"<tau>/tau_SE = {MEAN_TAU_FACTOR:.3f}")
for fp in xpcs_files:
    r = fits[fp]
    if r is not None:
        fmt = ' '.join(f'f{qi}={r["per_q"][qi]["f"]:.3f}' for qi in r['per_q'])
        print(f'  fit {elapsed(fp):5.0f} s: chi2={r["red_chi2"]:.2f} '
              f'p_slow={r["p2"]:.3f}+/-{r["p2_err"]:.3f}  {fmt}')

# --- FIGURE ---
apply_style()
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=FIG1_SIZE)
label_panels((ax1, ax2, ax3))

# ============================================================
# PANEL 1: SAXS I(q), background-subtracted (B0146 + first + last N B0147)
# ============================================================
# Background buffer curve AND its own absolute coefficient (coef_buf).  Both are
# taken from the D0138 buffer file; coef_buf multiplies the buffer term in every
# subtraction below.
bg_I = None
coef_buf = None
if BACKGROUND_HEADER in by_header:
    with h5py.File(by_header[BACKGROUND_HEADER][0], 'r') as hf:
        _, bg_I = read_saxs_iq(hf, PHI_AVERAGE)
        coef_buf = abs_xsec_coef(hf)
    print(f'coef_buf ({BACKGROUND_HEADER}) = {coef_buf:.3e}')

time_handles = []                         # the six elapsed times -> the key row
ref_handles = []                          # the 6 C reference -> its own small key
saxs_I_lo, saxs_I_hi = np.inf, 0.0        # data range, for the panel-(a) headroom
for fp in by_header.get('B0146', []):                        # 6 C reference (distinct)
    with h5py.File(fp, 'r') as hf:
        q, I = read_saxs_iq(hf, PHI_AVERAGE)
        coef_sam = abs_xsec_coef(hf)                          # this file's own coefficient
    print(f'coef_sam (B0146 {parse_name(fp)[1]}-{parse_name(fp)[2]}) = {coef_sam:.3e}')
    # Put the profile on the absolute scale, then take the buffer off.  The
    # buffer already carries its own coefficient, so both terms are absolute.
    I = coef_sam * I
    if bg_I is not None:
        I = I - coef_buf * bg_I
    pos = I > 0
    saxs_I_lo, saxs_I_hi = min(saxs_I_lo, I[pos].min()), max(saxs_I_hi, I[pos].max())
    ax1.plot(q[pos], I[pos], color=COLOR_6C, marker='s', ls='none', ms=MS,
             mfc='none', mew=MEW)
    ref_handles.append(Line2D([], [], color=COLOR_6C, marker='s', ls='none',
                              mfc='none', mew=MEW, ms=MS,
                              label='6 °C ref'))

for fp in saxs_files:                                        # first (0 s) + last N B0147
    with h5py.File(fp, 'r') as hf:
        q, I = read_saxs_iq(hf, PHI_AVERAGE)
        coef_sam = abs_xsec_coef(hf)                          # this file's own coefficient
    print(f'coef_sam (B0147 {parse_name(fp)[1]}-{parse_name(fp)[2]}) = {coef_sam:.3e}')
    # Put the profile on the absolute scale, then take the buffer off.  The
    # buffer already carries its own coefficient, so both terms are absolute.
    I = coef_sam * I
    if bg_I is not None:
        I = I - coef_buf * bg_I
    pos = I > 0
    saxs_I_lo, saxs_I_hi = min(saxs_I_lo, I[pos].min()), max(saxs_I_hi, I[pos].max())
    color = ecolor(fp)
    ax1.plot(q[pos], I[pos], color=color, marker='o', ls='none', ms=MS,
             mfc='none', mew=MEW)
    time_handles.append(Line2D([], [], color=color, marker='o', ls='none',
                               mfc='none', mew=MEW, ms=MS,
                               label=f'{elapsed(fp):.0f} s'))

ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlabel(r'$Q$ ($\AA^{-1}$)')
ax1.set_ylabel(r'$I(Q)$ (mm$^{-1}$)')
add_minor_grid(ax1)
# Three named Q ticks instead of the single 10^-2 the automatic log locator
# finds in this ~1-decade range; shared with Figures S4b and S8 so the three
# absolute-scale profiles can be compared tick for tick.  The limits sit well
# outside the data so no point sits on the frame.
ax1.set_xlim(2.6e-3, 4.4e-2)
q_log_ticks(ax1)

# The elapsed-time key lives at the bottom of the FIGURE (below), so panel (a)
# needs no headroom for it.  Only the one-line 6 C reference key sits inside,
# in the bottom-left corner next to the curve it labels: that corner is empty
# because I(Q) falls to the right.
ax1.set_ylim(saxs_I_lo * SAXS_YLO_PAD, saxs_I_hi * SAXS_YHI_PAD)
leg_r = ax1.legend(handles=ref_handles, loc='lower left', borderaxespad=0.3,
                   handlelength=1.2, handletextpad=0.3)
leg_r.get_frame().set_linewidth(LW_THIN)

# ============================================================
# PANEL 2: g2(tau) at target q, with fits
# ============================================================
q_val_target = None
for fp in xpcs_files:
    tau, g2, g2_err, q_vals = g2_data[fp]
    q_val_target = q_vals[target_q_idx]
    color = ecolor(fp)
    m = tau > 0
    # Data drawn at the THIN end of the allowed range (0.5 pt = 0.18 mm edges)
    # and below the fits in z, so the fit line reads as a line and not as one
    # more marker outline of the same colour.
    ax2.errorbar(tau[m], g2[m, target_q_idx], yerr=g2_err[m, target_q_idx], fmt='o',
                 color=color, mfc='none', markersize=MS, mew=LW_THIN,
                 capsize=1.5, elinewidth=LW_THIN, capthick=LW_THIN,
                 alpha=0.75, linestyle='none', zorder=2)
    r = fits[fp]
    if r is not None and target_q_idx in r['per_q'] and r['red_chi2'] < CHI2_MAX:
        pq = r['per_q'][target_q_idx]
        tt = tau[m]
        t_fit = np.logspace(np.log10(tt.min()), np.log10(tt.max()), 200)
        g2_fit = double_exp(t_fit, pq['tau_fast'], pq['f'], pq['tau_slow'],
                            r['p1'], r['p2'])
        # LW_DATA (1.0 pt) against 0.6 pt marker outlines is enough contrast on
        # its own -- deliberately NO white casing under the line, which would
        # blank out the very data points the fit is meant to be judged against.
        ax2.plot(t_fit, g2_fit, color=color, lw=LW_DATA, solid_capstyle='round',
                 zorder=3)

ax2.set_xscale('log')
ax2.set_xlabel(r'Delay Time, $\Delta t$ (s)')
ax2.set_ylabel(r'$g_2$')
ax2.set_ylim(1.0, 1.18)
# two-decimal ticks: '1.000'-style labels are wide enough at 8 pt to squeeze
# the three panels together, and the extra digit carries no information here.
ax2.yaxis.set_major_locator(FixedLocator([1.00, 1.05, 1.10, 1.15]))
ax2.yaxis.set_major_formatter(FixedFormatter(['1.00', '1.05', '1.10', '1.15']))
add_minor_grid(ax2)
# tau covers ~5 decades; label every other one so the exponents stay legible
# and well separated, keeping a tick (and grid line) on every decade.  The
# limits leave a clear margin either side of the shortest / longest delay.
ax2.set_xlim(9e-6, 5.0)
ax2.xaxis.set_major_locator(LogLocator(base=10.0, numticks=12))
def every_other_decade(x, _pos):
    """Label 1e-4, 1e-2, 1e0 and leave the odd decades blank, to avoid crowding."""
    power = int(round(np.log10(x)))
    if power % 2 != 0:
        return ''
    return rf'$10^{{{power}}}$'


ax2.xaxis.set_major_formatter(FuncFormatter(every_other_decade))

# ============================================================
# PANEL 3: fitted fast fraction f vs elapsed time (colour=time, marker=q)
# ============================================================
fit_rows = []
q_val_of = {}                                    # q_idx -> Q value (for legends)
for fp in xpcs_files:
    r = fits[fp]
    if r is None or r['red_chi2'] >= CHI2_MAX:   # drop unreliable global fits
        continue
    _, _, _, q_vals = g2_data[fp]
    for q_idx, pq in r['per_q'].items():
        q_val_of[q_idx] = q_vals[q_idx]
        fit_rows.append({'q_index': q_idx, 'q_val': q_vals[q_idx],
                         'elapsed': elapsed(fp), 'f': pq['f'], 'f_err': pq['f_err']})

for q_idx, marker in zip(fit_q_indices, Q_MARKERS):
    rows = sorted([r for r in fit_rows if r['q_index'] == q_idx], key=lambda r: r['elapsed'])
    if not rows:
        continue
    xs = [r['elapsed'] for r in rows]
    ys = [r['f'] for r in rows]
    es = [r['f_err'] for r in rows]
    ax3.plot(xs, ys, '-', color='0.75', lw=LW_THIN, zorder=1)        # per-q guide line
    for x, y, e in zip(xs, ys, es):
        ax3.errorbar(x, y, yerr=e, marker=marker, color=time_color[x],
                     markersize=MS_SPARSE, capsize=1.5, mfc='none', mew=MEW,
                     elinewidth=LW_THIN, capthick=LW_THIN, zorder=2)

# Headroom above the highest fitted point, for the q key.  Three columns make
# that key two rows deep instead of five, so 0.18 of the span is enough where a
# single column needed 0.42 -- the data now fills most of the panel.
_f_hi = max(r['f'] + r['f_err'] for r in fit_rows)
_f_lo = min(r['f'] - r['f_err'] for r in fit_rows)
_span = _f_hi - _f_lo
ax3.set_ylim(_f_lo - 0.08 * _span, _f_hi + 0.18 * _span)
# widen the elapsed-time axis so the first and last groups are not on the frame
_t_lo = min(r['elapsed'] for r in fit_rows)
_t_hi = max(r['elapsed'] for r in fit_rows)
_t_pad = 0.16 * (_t_hi - _t_lo)
ax3.set_xlim(_t_lo - _t_pad, _t_hi + _t_pad)

# legend: marker shape -> q bin (open black markers).  The Q values are carried
# by a legend TITLE with the 10^-3 factor pulled out, so each entry is a short
# mantissa instead of '$Q = 0.00376 \AA^{-1}$'.  Three columns keep the box
# shallow, so it clears the descending 5040 s points on little headroom.
q_handles = [Line2D([], [], marker=mk, ls='none', mfc='none', mec='k', mew=MEW,
                    markersize=MS_SPARSE, label=f'{1e3 * q_val_of[qi]:.2f}')
             for qi, mk in zip(fit_q_indices, Q_MARKERS) if qi in q_val_of]
leg3 = ax3.legend(handles=q_handles, loc='upper right', ncol=3,
                  borderaxespad=0.3, title=r'$Q$ ($10^{-3}\ \AA^{-1}$)',
                  handlelength=0.8, handletextpad=0.25, columnspacing=0.5,
                  labelspacing=0.2, borderpad=0.3)
leg3.get_frame().set_linewidth(LW_THIN)

ax3.set_xlabel('Elapsed Time (s)')
ax3.set_ylabel('Fast Fraction, $f$')
add_minor_grid(ax3)

# --- right axis: mean collective diffusion coefficient of the fast mode ---
# <D> = 1/(<tau> Q^2), averaged over the q bins where the fast mode is measured
# (the same det test panel (b) of Figure S8 uses).  Plotted beside f because the
# comparison is the point: f falls monotonically to zero while <D> moves by
# about a factor of two with no systematic trend.  The two endpoint times carry
# large errors -- at 5040 s the two-step decay is only just resolved, at 7863 s
# f has nearly vanished -- and those error bars are left visible rather than the
# points dropped.  Filled black squares against the open coloured f markers.
ax3r = ax3.twinx()
# D_eff is a fitted parameter of Equation 2, not a quantity derived afterwards:
# the fast mode is tied to tau_fast(Q) = 1/(D0 Q^2) in the model, and D_eff is
# D0 rescaled by p/Gamma(1/p) so that it belongs with the MEAN relaxation time.
# Its error bar therefore comes from the covariance matrix.
d_rows = [(elapsed(fp), fits[fp]['D_eff'], fits[fp]['D_eff_err'])
          for fp in xpcs_files
          if fits[fp] is not None and fits[fp]['red_chi2'] < CHI2_MAX]
d_rows.sort()
d_rows.sort()
xd = [r[0] for r in d_rows]
yd = [r[1] / 1e7 for r in d_rows]
ed = [r[2] / 1e7 for r in d_rows]
ax3r.errorbar(xd, yd, yerr=ed, marker='*', ls='none', color='k',
              markersize=MS_SPARSE + 2.5, capsize=1.5, mfc='none', mew=MEW,
              elinewidth=LW_THIN, capthick=LW_THIN, zorder=4)
ax3r.set_ylabel(r'$D_{\mathrm{eff}}$ ($10^{7}\ \AA^{2}\,$s$^{-1}$)')
ax3r.set_ylim(0, max(y + e for y, e in zip(yd, ed)) * 1.30)
print('  D_eff (1e7 A^2/s): ' + '  '.join(f'{x:.0f}s {y:.2f}+/-{e:.2f}'
                                        for x, y, e in zip(xd, yd, ed)))

# --- shared elapsed-time key, one row along the bottom of the figure ---
# Colour carries the same meaning in all three panels, so this key belongs to
# the figure, not to any single panel.  Out here it costs no panel any y range.
leg_fig = fig.legend(handles=time_handles, loc='lower center',
                     bbox_to_anchor=(0.5, 0.0), ncol=len(time_handles),
                     handlelength=1.2, handletextpad=0.4, columnspacing=1.2,
                     borderaxespad=0.0, borderpad=0.35)
leg_fig.get_frame().set_linewidth(LW_THIN)

# Reserve the bottom strip for that key; w_pad keeps the y labels of panels (b)
# and (c) clear of the neighbouring panel's tick labels.
fig.tight_layout(pad=0.4, w_pad=1.1, rect=(0, FIG1_LEGEND_H, 1, 1))
# Figure 3 of the manuscript is now the \SI{28.50}{\degreeCelsius} loading of the
# 2021-3 series, written by temperature_series_2021-3.py.  This panel set is kept
# for the 2022-1 analysis but is written under its own name so the two scripts
# never write the same file.
save_fig(fig, 'Figure_2022-1_Isothermal_SAXPCS.pdf')

# ============================================================
# FIGURE S8 (1x3): (a) p_slow vs elapsed time;  (b) <tau_slow> vs Q;
#                  (c) gamma_slow vs elapsed time.
# The fast mode gets no panel: Equation 2 ties it to tau_fast = 1/(D0 q^2), so
# a plot against q would be that line and nothing else, and the one number it
# carries, D_eff, is in Figure 3c.  What is left is the slow mode, free at
# every elapsed time.
# All XPCS colours use the same elapsed-time scale as Figure 3.
# ============================================================
fig2, (axps, axs, axg) = plt.subplots(1, 3, figsize=FIG2_SIZE)
label_panels((axps, axs, axg))   # reading order: a b c

# The fast mode has no panel of its own here.  Equation 2 ties it to
# tau_fast(Q) = 1/(D0 Q^2), so a plot of it against Q would be that line and
# nothing else; the one number it contains, D_eff, is in Figure 3c.  p_fast is
# likewise a single fitted value for the series and is quoted in the caption.
# What remains are the slow-mode parameters, which are free at every elapsed
# time.

# --- (a) slow-mode stretching exponent vs elapsed time ---
p2_rows = sorted((elapsed(fp), fits[fp]['p2'], fits[fp]['p2_err'], fp)
                 for fp in xpcs_files
                 if fits[fp] is not None and fits[fp]['red_chi2'] < CHI2_MAX)
axps.plot([r[0] for r in p2_rows], [r[1] for r in p2_rows], '-',
          color='0.75', lw=LW_THIN, zorder=1)
for x, y, e, fp in p2_rows:
    axps.errorbar(x, y, yerr=e, marker='s', color=ecolor(fp),
                  markersize=MS_SPARSE, capsize=1.5, mfc='none', mew=MEW,
                  elinewidth=LW_THIN, capthick=LW_THIN, zorder=2)
axps.set_xlabel('Elapsed Time (s)')
axps.set_ylabel(r'$p_{\mathrm{slow}}$')
add_minor_grid(axps)
_tp = 0.16 * (p2_rows[-1][0] - p2_rows[0][0])
axps.set_xlim(p2_rows[0][0] - _tp, p2_rows[-1][0] + _tp)

# --- (b) mean slow relaxation time vs Q, one curve per elapsed time; each is
# fitted to a power law tau = A Q**gamma for panel (c). ---
gamma_rows = []
for fp in xpcs_files:
    r = fits[fp]
    if r is None or r['red_chi2'] >= CHI2_MAX:
        continue
    _, _, _, q_vals = g2_data[fp]
    qs = sorted(r['per_q'])
    Q  = np.array([q_vals[qi] for qi in qs])
    # Mean relaxation time, not the bare Kohlrausch tau: <tau> = Gamma(1/p) tau / p.
    ts     = mean_tau(np.array([r['per_q'][qi]['tau_slow'] for qi in qs]), r['p2'])
    ts_err = mean_tau(np.array([r['per_q'][qi]['tau_slow_err'] for qi in qs]), r['p2'])
    color  = ecolor(fp)
    axs.errorbar(Q, ts, yerr=ts_err, marker='s', ls='none', color=color,
                 markersize=MS_SPARSE, capsize=1.5, mfc='none', mew=MEW,
                 elinewidth=LW_THIN, capthick=LW_THIN, zorder=2)
    if len(Q) >= 3:
        gs, gs_err, gs_A = fit_powerlaw(Q, ts, ts_err)
        Q_line = np.linspace(Q.min(), Q.max(), 50)
        axs.plot(Q_line, gs_A * Q_line**gs, '--', color=color, lw=LW_DATA, zorder=1)
        gamma_rows.append({'elapsed': elapsed(fp), 'fp': fp,
                           'gamma_slow': gs, 'gamma_slow_err': gs_err})
        print(f'  gamma_slow {elapsed(fp):5.0f} s: {gs:.3f}+/-{gs_err:.3f}')

# Q ticks: label the mantissa (4..8); the x10^-3 factor is folded into the
# shared x-axis label on the bottom row instead of a separate corner text.
_q_major = [4e-3, 5e-3, 6e-3, 7e-3, 8e-3]
# integer multiples of 10^-3 only; see the note on Figure 3a above
_q_minor = [9e-3]
for a in (axs,):
    a.set_xscale('log')
    a.set_yscale('log')
    a.set_xlim(3.25e-3, 9.6e-3)          # margin either side of Q = 3.76-8.27
    add_minor_grid(a)                              # sets minorticks_on + grids
    a.xaxis.set_major_locator(FixedLocator(_q_major))
    a.xaxis.set_minor_locator(FixedLocator(_q_minor))
    a.xaxis.set_major_formatter(FixedFormatter(['4', '5', '6', '7', '8']))
    a.xaxis.set_minor_formatter(NullFormatter())
axs.set_ylabel(r'$\langle\tau_{\mathrm{slow}}\rangle$ (s)')
axs.set_xlabel(r'$Q$ ($\times 10^{-3}\ \AA^{-1}$)')

# --- (c) power-law scaling exponent gamma_slow vs elapsed time ---
gamma_rows.sort(key=lambda r: r['elapsed'])
xg = [r['elapsed'] for r in gamma_rows]
# gamma_fast is imposed at -2 in the refit that yields <D>, so it is not a
# fitted quantity and is not plotted; the free-fit values that justify imposing
# it (-2.08 to -2.37) are quoted in the text.  The reference line marks -2.
axg.plot(xg, [r['gamma_slow'] for r in gamma_rows], '-', color='0.75', lw=LW_THIN, zorder=1)
for r in gamma_rows:
    axg.errorbar(r['elapsed'], r['gamma_slow'], yerr=r['gamma_slow_err'], marker='s',
                 color=ecolor(r['fp']), markersize=MS_SPARSE, capsize=1.5, mfc='none', mew=MEW,
                 elinewidth=LW_THIN, capthick=LW_THIN, zorder=2)
# Bottom right, in two columns.  NOT the lower left, which holds the 5040 s
# gamma_slow point (-3.78 +/- 0.37), the lowest in the panel; a box there hides
# it completely.  On the right the lowest thing is the 7863 s gamma_fast error
# bar, which stops at -3.33.
axg.set_xlabel('Elapsed Time (s)')
axg.set_ylabel(r'Slow-mode exponent ($\gamma_{\mathrm{slow}}$)')
add_minor_grid(axg)

# Widen the elapsed-time axis so the first and last points are not sitting on
# the frame (matplotlib's 5 % default is too tight for five widely-spaced
# groups).
_xe = [r['elapsed'] for r in gamma_rows]
_t_pad2 = 0.16 * (max(_xe) - min(_xe))
axg.set_xlim(min(_xe) - _t_pad2, max(_xe) + _t_pad2)

# matplotlib's default 5 % y margin leaves the extreme error-bar caps ~1 mm off
# the frame; widen it for some breathing room.
for a in (axs, axg):
    a.set_ymargin(0.11)
    a.autoscale_view()

fig2.tight_layout(pad=0.4, w_pad=1.4, h_pad=0.8)
# Renamed off the SI numbering: 'S8' now belongs to the 28.00 C hold of the
# 2021-3 series.  These 2022-1 fit parameters are in neither document; the
# script is kept because the absolute-scale calibration below is, and that
# calibration is shared by both beamtimes.
save_fig(fig2, 'Figure_2022-1_Fit_Parameters.pdf')

# ============================================================
# FIGURE S3 (2 panels): ion-chamber -> photon calibration
#   (left)  air transmission across the calibration measurements (mean ~0.868)
#   (right) photon counts vs upstream IC, with the linear fit CAL_A/CAL_B
# The same CAL_A, CAL_B and AIR_TRANSMISSION feed abs_xsec_coef(), so this
# figure documents exactly the calibration the SAXS panels are scaled by.
# ============================================================
fig3, (axt, axc) = plt.subplots(1, 2, figsize=FIG3_SIZE)
label_panels((axt, axc))

# --- left: air transmission (mean +/- standard deviation reported in legend) ---
axt.plot(np.arange(len(AIR_TRANS_SERIES)), AIR_TRANS_SERIES, 'ko:', ms=MS, mew=MEW, lw=LW_THIN)
axt.axhline(AIR_TRANSMISSION, color='b', ls='-', lw=LW_DATA,
            label=f'mean = {AIR_TRANSMISSION:.4f} ± {AIR_TRANS_STD:.4f}')
axt.set_xlabel('Measurement #')
axt.set_ylabel('Transmission of Air')
# bottom-left: the transmission falls to the right, so that corner is empty
axt.legend(loc='lower left', borderaxespad=0.4)
add_minor_grid(axt)

# --- right: ion chamber -> photon linear calibration ---
_up = CAL_UPIC[CAL_CROP:]
axc.plot(_up, CAL_PHOTONS[CAL_CROP:], 'ko', ms=MS, mew=MEW, label='measured')
_uline = np.linspace(_up.min(), _up.max(), 100)
axc.plot(_uline, CAL_A * _uline + CAL_B, 'b-', lw=LW_DATA,
         label=f'fit: {CAL_A:.2e}·Up_IC {CAL_B:+.2e}')
axc.set_xlabel('Upstream Ion Chamber')
axc.set_ylabel('Number of Photons')
axc.legend(loc='upper left')
add_minor_grid(axc)

fig3.tight_layout(pad=0.4, w_pad=1.4)
save_fig(fig3, 'FigureS5_Calibration.pdf')

plt.show()
