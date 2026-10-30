"""Fits to the 12-ID-B profiles at 30 C: the low-Q power law and the two WAXS peaks.

The fits live here once because two scripts draw them: Plot_12ID.py (Figure 1
of the 28 Oct version) and reversibility_12id.py (Figure S5a).

* Power law: a straight line in log10 I against log10 Q over
  0.012 < Q < 0.040 A^-1.
* WAXS peaks: a Gaussian on a sloping background, a exp[-(Q - Q0)^2 / 2 s^2]
  + m Q + c, fitted independently over 0.48 < Q < 1.00 and 1.18 < Q < 1.70
  A^-1.  The background is the wing of the neighbouring scattering under each
  peak.  FWHM = 2.355 s.
"""

import numpy as np
from scipy.optimize import curve_fit

from guinier import load

POWER_LAW_RANGE = (0.012, 0.040)                      # A^-1
PEAK_WINDOWS = [((0.48, 1.00), 0.74), ((1.18, 1.70), 1.44)]   # window, starting Q0


def gauss_bg(x, a, x0, sigma, m, c):
    """Gaussian of height a, centre x0 and width sigma on the background m x + c."""
    return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2)) + m * x + c


def power_law(name):
    """(slope, log10 prefactor) of the power law fitted to one merged profile."""
    Q, I, _ = load(name)
    m = (Q > POWER_LAW_RANGE[0]) & (Q < POWER_LAW_RANGE[1])
    slope, intercept = np.polyfit(np.log10(Q[m]), np.log10(I[m]), 1)
    return slope, intercept


def waxs_peaks(name):
    """One dict per peak: fitted parameters, the Q it was fitted over, centre,
    FWHM and real-space spacing 2 pi / Q0 (A)."""
    Q, I, _ = load(name)
    out = []
    for (lo, hi), x0 in PEAK_WINDOWS:
        m = (Q > lo) & (Q < hi)
        p, _ = curve_fit(gauss_bg, Q[m], I[m], p0=[np.ptp(I[m]), x0, 0.15, 0, np.min(I[m])])
        out.append(dict(p=p, q=Q[m], center=p[1], fwhm=2.355 * p[2], d=2 * np.pi / p[1]))
    return out
