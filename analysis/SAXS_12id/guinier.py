"""Weighted Guinier fit of the 10 C Reference aliquot from 12-ID-B.

The fit lives here once because three scripts need it: Guinier_Plot.py and
reversibility_12id.py draw it, and gelation/saxs_scaling.py takes the fitted
R_g as the length of the dispersed state from which xi(T) is scaled.

ln I = ln I_0 - (R_g^2 / 3) Q^2 over 0.020 <= Q <= 0.055 A^-1, each point
weighted by 1 / sigma_lnI^2 with sigma_lnI = sigma_I / I, the uncertainty that
Read_12ID_SAWAXS.py propagated into reduced_data/.  This is the weighting the
12-ID beamline software uses; an unweighted fit ignores the measured errors and
reports an artificially small uncertainty.
"""

import os

import numpy as np
import pandas as pd

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'reduced_data')
Q_MIN, Q_MAX = 0.020, 0.055        # A^-1


def load(name='Merged_Reference_10C.csv'):
    """(Q, I, I_err) of one merged profile; Q in A^-1, I in arbitrary units."""
    df = pd.read_csv(os.path.join(DATA, name))
    return df['Q(A^-1)'].values, df['I(Q)'].values, df['I_err(Q)'].values


def guinier_fit(q_min=Q_MIN, q_max=Q_MAX):
    """Fit the Reference profile at 10 C; returns a dict of results.

    R_g is in A.  Uncertainties are from the covariance of the weighted fit
    (absolute-sigma convention), propagated through sqrt(-3 slope) and exp().
    """
    Q, I, I_err = load()
    m = (Q >= q_min) & (Q <= q_max) & (I > 0) & (I_err > 0)
    x, y = Q[m] ** 2, np.log(I[m])
    sigma = I_err[m] / I[m]
    w = 1.0 / sigma ** 2
    X = np.vstack([x, np.ones_like(x)]).T
    cov = np.linalg.inv(X.T @ (w[:, None] * X))
    slope, intercept = cov @ (X.T @ (w * y))
    R_g = np.sqrt(-3 * slope)
    resid = (y - (slope * x + intercept)) / sigma
    return dict(R_g=R_g, R_g_err=(3.0 / (2.0 * R_g)) * np.sqrt(cov[0, 0]),
                I_0=np.exp(intercept), I_0_err=np.exp(intercept) * np.sqrt(cov[1, 1]),
                slope=slope, intercept=intercept, q_min=q_min, q_max=q_max,
                red_chi2=np.sum(resid ** 2) / (len(x) - 2), n=len(x))
