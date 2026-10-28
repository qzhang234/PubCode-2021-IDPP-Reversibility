"""Shared XPCS two-mode fit: the model, the measured contrast, and the global fit.

Figure 3b,c and panels (b)-(g) of Figures S7-S10 must all show the SAME model
applied the SAME way -- (d)-(g) are an expansion of (b) onto the four higher q
bins -- so the model and the fitting routine live here and every figure script
imports them.  (Before this module, g2_grid_SI.py carried its own copy that
fitted each q bin independently with the stretching exponents frozen at 0.5, so
that grid was showing a different model from the one its caption described.)

CONTRAST
--------
beta is the instrumental speckle contrast: a property of the beam coherence and
the detector pixel size, not of the sample.  It is therefore measured once, on a
static reference, and held fixed for every sample fit.

The value below was measured on the 10 nm nano-porous glass standard (Doraglas
S10-10-1200-50, 10 mm x 1.2 mm, 10 nm pores), dataset F0145, at 6 C.  The value
comes from the group average average_ranges.py writes into data/ -- 41 of the 50
acquisitions, the survivors of the same two outlier cuts every other group in
the paper gets -- and it is measured in the SAME five q bins the sample is
fitted in, which is where it is used.  A constant is fitted to the averaged g2
of each, as a static sample requires; the five are flat (reduced chi^2 0.71 to
1.03) and return 0.12948, 0.13124, 0.13291, 0.13145 and 0.13186, whose mean is

    g2 = 1 + beta,   beta = 0.13139 +/- 0.00056,

the error being the standard error of that mean.  The strongest-counted flat bin
anywhere, bin 16 at q = 0.02067 A^-1, gives 0.13031, within one standard error,
so nothing hinges on the choice: moving between the two shifts every fitted fast
fraction by less than 0.02.
"""

import h5py
import numpy as np
from scipy.optimize import least_squares
from scipy.special import gamma as gamma_fn

from nexus_read import read_g2

# --- measured instrumental contrast (see the module docstring) ---
# This is the 2022-1 value.  The 2021-3 scripts measure their own beta with
# measure_beta() below and overwrite CONTRAST before fitting anything; a script
# that forgets to do so silently fits the 2021-3 data at the 2022-1 contrast,
# which is 6.8% low and drags every fitted f and the scaling temperature with
# it.  Measuring and overwriting is one line, and it belongs beside the
# constant it replaces, which is why measure_beta lives here rather than in
# each figure script.
CONTRAST = 0.13139
BASELINE = 1.0

# --- measuring beta on a static reference -----------------------------------
# NOT measured over the five fitted bins.  The 2021-3 static reference is the
# NIST glassy carbon standard, whose surface scratches add parasitic scattering
# at low Q.  That scattering is static, so it dilutes the speckle and depresses
# the APPARENT contrast exactly where the paper fits.  The measurement shows it
# plainly: beta climbs from 0.1226 at 0.0036 A^-1 to a plateau of 0.140-0.143
# by 0.0086 A^-1, a 12% rise.  The 2022-1 reference, a different sample, rises
# by only 2% over the same interval, so the effect belongs to the standard and
# not to the beamline.
#
# Contrast is a property of the optics and the detector, near enough constant
# in Q over this range, so it is taken from the clean plateau and applied to
# every bin.  BETA_BAND stops at 0.0165 A^-1, the last bin before beta leaves
# the plateau: it dips to 0.126 at 0.0175-0.0185 and reads 0.1429 at 0.0214,
# the bin holding the stuck pixel, against 0.1399 either side.  The dip was
# once attributed here to a detector module gap; the XSPA-500k has none, so no
# cause is assumed and the band is set on the measurement alone.  A single
# scalar rather than a per-bin value: the bin-to-bin scatter inside the plateau
# is instrumental noise, and letting it into f would inject that noise into the
# very ordering the paper reports.
BETA_BAND = (0.0085, 0.0165)          # A^-1, the flat part of the standard


def flat_level(tau, g, e):
    """Inverse-variance weighted mean of g2 over the measurable delays, and its
    standard error.  Every delay carries its own g2 uncertainty, so each enters
    with weight 1/e^2 and the noisy short delays count for little."""
    ok = (tau > 0) & (tau < 2) & np.isfinite(g) & (e > 0)
    w = 1.0 / e[ok] ** 2
    return float(np.sum(w * g[ok]) / np.sum(w)), float(np.sqrt(1.0 / np.sum(w)))


def measure_beta(path, label, fitted_bins=(0, 1, 2, 3, 4)):
    """beta from the flat correlation level of a static reference, over
    BETA_BAND.  A static sample cannot relax, so each of its correlation
    functions is flat in delay and sits at 1 + beta."""
    with h5py.File(path, 'r') as hf:
        tau, g2, g2e, q_std = read_g2(hf)
    lv = [flat_level(tau, g2[:, i], g2e[:, i]) for i in range(g2.shape[1])]
    b = np.array([x[0] - 1.0 for x in lv])
    e = np.array([x[1] for x in lv])
    sel = (q_std >= BETA_BAND[0]) & (q_std <= BETA_BAND[1]) & np.isfinite(b)
    w = 1.0 / e[sel] ** 2
    beta = float(np.sum(w * b[sel]) / np.sum(w))
    print(f'beta ({label}) = {beta:.5f} from {sel.sum()} bins over '
          f'{BETA_BAND[0]}-{BETA_BAND[1]} A^-1 (spread {b[sel].min():.4f}-{b[sel].max():.4f})')
    print('   over the five FITTED bins it would read '
          + ', '.join(f'{b[i]:.4f}' for i in fitted_bins) + '  <- depressed by the standard')
    return beta

# per-q parameter bounds / start (tau_fast, f, tau_slow) and shared (p1, p2)
PQ_P0    = [1e-3, 0.5, 100.0]
PQ_LO    = [1e-6, 0.0, 1.0]
PQ_HI    = [10.0, 1.0, 10000.0]
P_EXP_P0 = [0.5, 0.5]
P_EXP_LO = [0.2, 0.2]
P_EXP_HI = [3.0, 3.0]

# Mean collective diffusion coefficient, used only by fit_g2_diffusive.  <D> is
# of order 1e7 A^2/s while every other fitted parameter is of order 1, so it is
# fitted in units of D_SCALE.  Without that rescaling its Jacobian column is
# ~1e-8 smaller than the rest, J^T J is numerically singular in that direction,
# and the pseudo-inverse returns a zero uncertainty for it.
D_SCALE = 1e7
D_LO, D_HI = 1e-5, 1e3          # in units of D_SCALE


def double_exp(tau, tau_fast, f, tau_slow, p1, p2):
    """Two-mode g2: a fraction f of the scattering relaxes fast, 1-f slowly.

    Each mode is a stretched exponential (Kohlrausch-Williams-Watts): p = 1 is
    a plain exponential, p < 1 stretches the decay, p > 1 compresses it.  The
    two are added as FIELDS and then squared, because g2 is the intensity
    correlation and the field correlation is what decays (Siegert relation).
    """
    decay_fast = f * np.exp(-(tau / tau_fast)**p1)
    decay_slow = (1 - f) * np.exp(-(tau / tau_slow)**p2)
    return CONTRAST * (decay_fast + decay_slow)**2 + BASELINE


def bin_params(p, i, offset=2):
    """The three parameters belonging to q bin number i.

    The fit works on one long list of numbers, p.  The first `offset` entries
    are shared: the exponents p1 and p2 when they are free, nothing when they
    are held fixed.  After those come three numbers for each q bin, always in
    the order tau_fast, f, tau_slow.  So with the default offset bin 0 starts
    at position 2, bin 1 at position 5, and bin i at position 2 + 3 * i.
    """
    start = offset + 3 * i
    return p[start], p[start + 1], p[start + 2]


def fit_g2_global(tau, g2, g2_err, q_indices, p1_fixed=None):
    """Global fit of several q bins for one elapsed time, sharing p1 and p2.

With p1_fixed the FAST stretching exponent is held at the value fit_g2_joint
    returns for the whole series, so that tau_fast is comparable between elapsed
    times; the parameter vector then drops its first entry.  p_slow stays free
    at every elapsed time -- tau_slow lies beyond the acquisition window and is
    not used for scaling, so there is nothing to gain by constraining it and the
    fit is better without.

    Parameter vector = [p1, p2, (tau_fast, f, tau_slow) x nq].  Minimises the
    error-weighted residual (model - g2) / g2_err over all q simultaneously,
    using the g2_err stored in the file directly (absolute_sigma convention).
    Parameter 1-sigma errors are sqrt(diag(inv(J^T J))).

    Returns a dict:
      {'p1','p1_err','p2','p2_err','red_chi2',
       'per_q': {q_idx: {'tau_fast','tau_fast_err','f','f_err',
                         'tau_slow','tau_slow_err'}}}
    or None if no q bin has usable data.
    """
    data = []
    for qi in q_indices:
        # A delay point is usable when the delay is positive, g2 is a real
        # number, and its uncertainty is a real number above zero, since the
        # fit divides by that uncertainty.
        usable = (tau > 0) & np.isfinite(g2[:, qi]) & (g2_err[:, qi] > 0)
        if usable.sum() >= 5:
            data.append((qi, tau[usable], g2[usable, qi], g2_err[usable, qi]))
    nq = len(data)
    if nq == 0:
        return None

    # Fixing p_fast removes one entry from the head of the parameter vector, so
    # every per-q block shifts down by one and the offset is 1 instead of 2.
    off = 1 if p1_fixed is not None else 2

    def residual(p):
        p1, p2 = (p1_fixed, p[0]) if p1_fixed is not None else (p[0], p[1])
        parts = []
        for i, (qi, tv, gv, ev) in enumerate(data):
            tf, f, ts = bin_params(p, i, offset=off)
            parts.append((double_exp(tv, tf, f, ts, p1, p2) - gv) / ev)
        return np.concatenate(parts)

    # "PQ_P0 * nq" repeats that list of three numbers once per q bin, which
    # builds the starting guess for the whole parameter vector in one line.
    # head = [p_slow] when p_fast is fixed, [p_fast, p_slow] when neither is
    x0 = (P_EXP_P0[1:] if p1_fixed is not None else list(P_EXP_P0)) + PQ_P0 * nq
    lo = (P_EXP_LO[1:] if p1_fixed is not None else list(P_EXP_LO)) + PQ_LO * nq
    hi = (P_EXP_HI[1:] if p1_fixed is not None else list(P_EXP_HI)) + PQ_HI * nq
    res = least_squares(residual, x0, bounds=(lo, hi), max_nfev=40000)

    ndof = max(len(res.fun) - len(res.x), 1)
    red_chi2 = float(np.sum(res.fun**2) / ndof)
    # Parameter uncertainties from the Gauss-Newton approximation to the
    # Hessian of the error-weighted residuals.
    #
    # Pseudo-inverse rather than plain inverse, as a safeguard.  If an amplitude
    # f ever settles exactly on its bound of 0, its tau_fast stops affecting the
    # model, J^T J becomes singular, and inv() returns nan for EVERY parameter,
    # including the well-determined ones.  pinv drops only the useless direction
    # and leaves the rest with honest uncertainties.  This is not hypothetical:
    # f refines to exactly 0 in the lowest q bin at t_w = 7863 s, where the fast
    # decay has gone, so that one direction really is useless and pinv is what
    # keeps the other 60 parameters' errors finite.
    jtj = res.jac.T @ res.jac
    cov = np.linalg.pinv(jtj, rcond=1e-12)
    perr = np.sqrt(np.abs(np.diag(cov)))

    if p1_fixed is not None:
        p1v, p1e = p1_fixed, 0.0
        p2v, p2e = res.x[0], perr[0]
    else:
        p1v, p1e = res.x[0], perr[0]
        p2v, p2e = res.x[1], perr[1]
    out = {'p1': p1v, 'p1_err': p1e, 'p2': p2v, 'p2_err': p2e,
           'red_chi2': red_chi2, 'n_par': len(res.x), 'per_q': {}}
    for i, (qi, tv, gv, ev) in enumerate(data):
        tf, f, ts = bin_params(res.x, i, offset=off)
        tf_err, f_err, ts_err = bin_params(perr, i, offset=off)
        out['per_q'][qi] = {'tau_fast': tf, 'tau_fast_err': tf_err,
                            'f': f, 'f_err': f_err,
                            'tau_slow': ts, 'tau_slow_err': ts_err}
    return out


def mean_tau(tau_se, p):
    """Mean relaxation time of a stretched exponential.

    A Kohlrausch decay exp[-(t/tau_se)^p] does not relax on the time tau_se;
    its mean is

        <tau> = tau_se * Gamma(1/p) / p,

    which is the quantity to compare when p is not the same everywhere (Guo et
    al., Phys. Rev. Lett. 109, 055901 (2012), Eq. 1 and the text below it).
    With p held global across the series this is a constant rescaling, so the
    q dependence is unchanged and only the absolute time shifts; the reason to
    apply it is that <tau> is the physically meaningful relaxation time and is
    what the XPCS literature quotes.
    """
    return tau_se * gamma_fn(1.0 / p) / p


def eff_diffusion(D0, p):
    """Effective diffusion coefficient from the Kohlrausch one.

    Equation 2 parameterises the fast mode as tau_fast(Q) = 1/(D0 Q^2), so D0 is
    built from the Kohlrausch time.  A stretched exponential does not relax on
    that time: its mean is <tau> = Gamma(1/p) tau / p (Guo et al., Phys. Rev.
    Lett. 109, 055901 (2012)).  The effective coefficient, the one that belongs
    with a mean relaxation time and the one the XPCS literature quotes, is

        D_eff = 1/(<tau_fast> Q^2) = [p / Gamma(1/p)] D0.

    With p shared across the waiting-time series this is a constant rescaling,
    but it is what keeps D_eff comparable between fits if p ever differs.
    """
    return D0 * p / gamma_fn(1.0 / p)


def fit_g2_diffusive(tau, g2, g2_err, q_indices, q_vals, p1_fixed=None):
    """Two-mode fit with the fast mode constrained to diffuse: Equation 2.

    tau_fast is not free in each q bin.  It is tied to a single coefficient,

        tau_fast(Q) = 1 / (D0 Q^2),

    which is the model the paper reports.  Diffusion is the first thing to try,
    and it describes the data.  A stretched g2 is still diffusive provided p
    does not depend on Q and tau goes as Q^-2: the first holds by construction
    here, and the second is what an unconstrained fit returns, with Q exponents
    of -2.08 to -2.37 at the elapsed times where the fast mode is resolved.

    p1_fixed holds p_fast at the value fit_g2_joint shares across the series.
    p_slow is always free: tau_slow lies beyond the acquisition window and is
    never used for scaling, so constraining it buys nothing.

    Parameter vector = [p_fast (omitted when fixed), p_slow, D0] followed by
    (f, tau_slow) for each q bin.  D0 is fitted in units of D_SCALE because it
    is of order 1e7 A^2/s while every other parameter is of order 1; without
    that rescaling its Jacobian column is small enough that the pseudo-inverse
    returns a zero uncertainty for it.

    Returns a fit_g2_global-style dict with 'D0', 'D0_err', 'D_eff' and
    'D_eff_err' added, or None if no q bin has usable data.
    """
    data = []
    for qi in q_indices:
        usable = (tau > 0) & np.isfinite(g2[:, qi]) & (g2_err[:, qi] > 0)
        if usable.sum() >= 5:
            data.append((qi, tau[usable], g2[usable, qi], g2_err[usable, qi]))
    nq = len(data)
    if nq == 0:
        return None
    qs = np.array([q_vals[qi] for qi, _, _, _ in data])

    # Seed from the unconstrained fit, whose tau_fast gives a starting D0.
    free = fit_g2_global(tau, g2, g2_err, q_indices, p1_fixed=p1_fixed)
    if free is None:
        return None
    d_seed = float(np.median([1.0 / (free['per_q'][qi]['tau_fast'] * q ** 2)
                              for (qi, _, _, _), q in zip(data, qs)])) / D_SCALE

    # Head of the parameter vector: p_fast only when it is free, then p_slow and
    # D0.  HEAD is where the per-bin (f, tau_slow) pairs begin.
    HEAD = 2 if p1_fixed is not None else 3

    def unpack(p):
        pf = p1_fixed if p1_fixed is not None else p[0]
        return pf, p[HEAD - 2], p[HEAD - 1] * D_SCALE        # p_fast, p_slow, D0

    def pair(p, i):
        return p[HEAD + 2 * i], p[HEAD + 2 * i + 1]           # f, tau_slow

    def residual(p):
        pf, ps, D0 = unpack(p)
        parts = []
        for i, (qi, tv, gv, ev) in enumerate(data):
            f, ts = pair(p, i)
            parts.append((double_exp(tv, 1.0 / (D0 * qs[i] ** 2), f, ts, pf, ps)
                          - gv) / ev)
        return np.concatenate(parts)

    x0 = ([] if p1_fixed is not None else [free['p1']]) + [free['p2'], d_seed]
    lo = ([] if p1_fixed is not None else [P_EXP_LO[0]]) + [P_EXP_LO[1], D_LO]
    hi = ([] if p1_fixed is not None else [P_EXP_HI[0]]) + [P_EXP_HI[1], D_HI]
    for qi, _, _, _ in data:
        pq = free['per_q'][qi]
        x0 += [pq['f'], pq['tau_slow']]
        lo += [PQ_LO[1], PQ_LO[2]]
        hi += [PQ_HI[1], PQ_HI[2]]
    x0 = [min(max(v, l), h) for v, l, h in zip(x0, lo, hi)]
    res = least_squares(residual, x0, bounds=(lo, hi), max_nfev=120000)

    ndof = max(len(res.fun) - len(res.x), 1)
    cov = np.linalg.pinv(res.jac.T @ res.jac, rcond=1e-12)
    perr = np.sqrt(np.abs(np.diag(cov)))
    pf, ps, D0 = unpack(res.x)
    D0_err = perr[HEAD - 1] * D_SCALE
    scale = pf / gamma_fn(1.0 / pf)
    out = {'p1': pf, 'p1_err': 0.0 if p1_fixed is not None else perr[0],
           'p2': ps, 'p2_err': perr[HEAD - 2],
           'D0': D0, 'D0_err': D0_err,
           'D_eff': scale * D0, 'D_eff_err': scale * D0_err,
           'red_chi2': float(np.sum(res.fun**2) / ndof),
           'n_par': len(res.x), 'per_q': {}}
    for i, (qi, tv, gv, ev) in enumerate(data):
        tf = 1.0 / (D0 * qs[i] ** 2)
        f, ts = pair(res.x, i)
        f_err, ts_err = pair(perr, i)
        out['per_q'][qi] = {'tau_fast': tf, 'tau_fast_err': tf * D0_err / D0,
                            'f': f, 'f_err': f_err,
                            'tau_slow': ts, 'tau_slow_err': ts_err}
    return out


def fit_g2_joint(datasets, q_vals_list, q_indices):
    """Fit the whole waiting-time series with Equation 2, sharing p_fast.

    Each elapsed time gets its own D0, p_slow, and (f, tau_slow) per q bin;
    p_fast is one number for the series.  Sharing it costs nothing -- reduced
    chi^2 is no worse than letting it float, for four fewer parameters -- and
    the chi^2 profile is flat to within 1 % between p_fast 0.55 and 0.70, so
    the data cannot resolve five separate values.  Letting a poorly determined
    p float also lets it trade against tau_fast, which is what makes relaxation
    times incomparable between elapsed times in the first place.

    Returns (joint, per_time) with joint = {'p1','p1_err','red_chi2','n_par'}.
    """
    seeds = [fit_g2_global(t, g, e, q_indices) for t, g, e in datasets]
    if not any(s is not None for s in seeds):
        return None, None
    p1_0 = float(np.mean([s['p1'] for s in seeds if s is not None]))

    blocks, qss = [], []
    for (tau, g2, g2_err), q_vals in zip(datasets, q_vals_list):
        d = []
        for qi in q_indices:
            usable = (tau > 0) & np.isfinite(g2[:, qi]) & (g2_err[:, qi] > 0)
            if usable.sum() >= 5:
                d.append((qi, tau[usable], g2[usable, qi], g2_err[usable, qi]))
        blocks.append(d)
        qss.append(np.array([q_vals[qi] for qi, _, _, _ in d]))

    # [p_fast] then, per elapsed time, [p_slow, D0, (f, tau_slow) x nq].
    starts, x0, lo, hi = [], [p1_0], [P_EXP_LO[0]], [P_EXP_HI[0]]
    for d, qs, s in zip(blocks, qss, seeds):
        starts.append(len(x0))
        d0 = float(np.median([1.0 / (s['per_q'][qi]['tau_fast'] * q ** 2)
                              for (qi, _, _, _), q in zip(d, qs)])) / D_SCALE
        x0 += [s['p2'], d0]; lo += [P_EXP_LO[1], D_LO]; hi += [P_EXP_HI[1], D_HI]
        for qi, _, _, _ in d:
            pq = s['per_q'][qi]
            x0 += [pq['f'], pq['tau_slow']]
            lo += [PQ_LO[1], PQ_LO[2]]
            hi += [PQ_HI[1], PQ_HI[2]]
    x0 = [min(max(v, l), h) for v, l, h in zip(x0, lo, hi)]

    def residual(p):
        parts = []
        for d, qs, st in zip(blocks, qss, starts):
            ps, D0 = p[st], p[st + 1] * D_SCALE
            for i, (qi, tv, gv, ev) in enumerate(d):
                f, ts = p[st + 2 + 2 * i], p[st + 3 + 2 * i]
                parts.append((double_exp(tv, 1.0 / (D0 * qs[i] ** 2), f, ts,
                                         p[0], ps) - gv) / ev)
        return np.concatenate(parts)

    res = least_squares(residual, x0, bounds=(lo, hi), max_nfev=300000)
    ndof = max(len(res.fun) - len(res.x), 1)
    cov = np.linalg.pinv(res.jac.T @ res.jac, rcond=1e-12)
    perr = np.sqrt(np.abs(np.diag(cov)))
    joint = {'p1': res.x[0], 'p1_err': perr[0],
             'red_chi2': float(np.sum(res.fun**2) / ndof), 'n_par': len(res.x)}
    per_time = [fit_g2_diffusive(t, g, e, q_indices, qv, p1_fixed=res.x[0])
                for (t, g, e), qv in zip(datasets, q_vals_list)]
    return joint, per_time
