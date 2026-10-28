"""Figure S7: measurement of the instrumental speckle contrast beta.

beta is a property of the beamline -- the transverse coherence of the incident
beam and the speckle-to-pixel size ratio -- not of the sample, so it is measured
once on a STATIC reference and then held fixed in every g2 fit.  A static sample
has no decorrelation, so its correlation function is flat at

    g2(tau) = 1 + beta.

Reference: nano-porous glass, 10 nm nominal pore diameter, 10 x 10 mm, 1.2 mm
thick (Doraglas S10-10-1200-50).  Dataset F0145, 6 C, 50 repeat acquisitions of
100,000 frames each, taken at 19:49 on 2022-03-06 -- about two hours before the
isothermal series of Figure 3 (21:57), in the same beamline configuration.  This
script reads the group average that average_ranges.py writes into data/, on the
same terms as every other group in the paper: the two outlier cuts leave 41 of
the 50 acquisitions, and those 41 are averaged.  Nothing here touches the
beamline storage.

Panel (a): the averaged g2 at q = 0.00376 A^-1, the lowest of the five bins the
sample is analysed in, with its measured pointwise uncertainties and the
straight-line fit in log(delay) that a static sample requires.  It is flat to
0.3 sigma with a reduced chi^2 of 0.71.

Panel (b): the same constant fitted bin by bin, giving beta(q).  The ADOPTED
value is the unweighted mean over the five bins the sample is analysed in,

    beta = 0.13139 +/- 0.00056,

the error being the standard error of that mean; the five are 0.12948, 0.13124,
0.13291, 0.13145 and 0.13186 at q = 0.00376, 0.00489, 0.00602, 0.00714 and
0.00827 A^-1, and every one of them is flat (reduced chi^2 <= 1.03).  Measuring
beta in the same q range it is applied in removes the only assumption the older
choice needed, that the plateau in panel (b) extends down from q = 0.02067 A^-1
to the sample range.  It does: the strongest-counted flat bin, bin 16 at
q = 0.02067 A^-1, returns 0.13031, within one standard error of the adopted
mean, and the bin-to-bin spread of beta falls as 1/sqrt(N) with the number of
acquisitions averaged, so it is counting noise rather than a real q dependence.
The two choices differ by 0.8 %, and moving between them shifts every fitted
fast fraction by less than 0.02.

Four of the 27 bins are excluded from panel (b) and flagged in the print-out
(bins 7, 12, 15 and 22): their averaged correlation function is not flat,
reduced chi^2 > 5 against a constant.  All four lie at q >= 0.0105 A^-1, above
the range used for the sample.

Panel (b) is flat at about 0.132 out to q ~ 0.021 A^-1 and then falls off,
reaching 0.107 in the outermost bin at q = 0.0331 A^-1.  That falloff is real
and instrumental, not a defect of the standard: the 0.03 % bandwidth passed by
the Si(111) monochromator at 10.94 keV gives a finite LONGITUDINAL coherence
length, and the path-length difference across the scattering volume grows with
scattering angle, so contrast is progressively lost as q rises.  All five bins
used for the sample sit well below that onset.
"""

import os
import sys

import numpy as np
import h5py
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from common.acs_style import (DOUBLE_COL, MS, MS_DENSE, MEW, LW_THIN, LW_DATA,
                              apply_style, add_minor_grid, label_panels, save_fig)

# The group average written by average_ranges.py.  The BadpixRm tag records
# that this group was correlated against a qmap with three hot detector pixels
# masked; see GROUP_SUFFIX there.
GLASS = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data',
                     'Average_F0145_10nm_Glass_006C_att00_Rq0_00001_00050'
                     '_BadpixRm_results.hdf')
SAMPLE_Q = [0, 1, 2, 3, 4]   # 0-based; the five bins saxpcs.py fits the sample in
                             # (fit_q_indices there).  beta is the mean over these.
Q_INDEX = SAMPLE_Q[0]        # the bin panel (a) shows: q = 0.00376 A^-1, the
                             # lowest of the five and the worst-counted of them
CHI2_FLAT_MAX = 5.0          # a q bin whose averaged g2 is this far from flat is
                             # an artefact, not a contrast measurement
FIG_SIZE = (DOUBLE_COL, 2.9)


def load(path):
    """Delay times, q bins, the averaged g2 and g2_err, and how many averaged."""
    with h5py.File(path, 'r') as hf:
        # frame_time is stored as a tiny array, shape (1,) or (1, 1) depending
        # on the file, so flatten to 1-D and take the first entry.
        ft = float(np.asarray(hf['/entry/instrument/detector_1/frame_time'][()]).reshape(-1)[0])
        t = np.asarray(hf['/xpcs/multitau/delay_list'][()])
        if t.ndim > 1:
            t = t[:, 0]                  # some files store the delays as a column
        tau = t * ft
        q = hf['/xpcs/qmap/dynamic_v_list_dim0'][()]
        G = hf['/xpcs/multitau/normalized_g2'][()]
        E = hf['/xpcs/multitau/normalized_g2_err'][()]
        n = hf['/xpcs/average/file_list'].shape[0]
    return tau, np.asarray(q), np.asarray(G), np.asarray(E), n


def flat_level(tau, g, e):
    """Weighted mean level of a flat g2, its error, and the reduced chi^2 of the
    flat hypothesis."""
    ok = (tau > 0) & (tau < 2) & np.isfinite(g) & (e > 0)
    weight = 1.0 / e[ok]**2              # least squares weights each point by 1/sigma^2
    level = np.sum(weight * g[ok]) / np.sum(weight)
    level_err = 1.0 / np.sqrt(np.sum(weight))
    chi2 = np.sum(weight * (g[ok] - level)**2) / (ok.sum() - 1)
    return level, level_err, chi2


tau, q, G, E, n_avg = load(GLASS)
print(f'{n_avg} acquisitions of the 10 nm glass standard averaged, {len(q)} q bins')

# --- the working bin, and the straight-line fit ---
gbar, ebar = G[:, Q_INDEX], E[:, Q_INDEX]
ok = (tau > 0) & (tau < 2) & np.isfinite(gbar) & (ebar > 0)
x, y, w = np.log10(tau[ok]), gbar[ok], 1.0 / ebar[ok]**2
# Weighted straight-line fit, y = slope * x + intercept, done with the normal
# equations.  Each row of X is one data point, [x, 1].  Writing the weights as a
# column, w[:, None], multiplies each ROW of X by that point's weight.
X = np.vstack([x, np.ones_like(x)]).T
cov = np.linalg.inv(X.T @ (w[:, None] * X))
slope, intercept = cov @ (X.T @ (w * y))
level, _, chi2 = flat_level(tau, gbar, ebar)
print(f'q bin {Q_INDEX+1} (q = {q[Q_INDEX]:.5f} A^-1), shown in panel (a):')
print(f'  straight-line slope = {slope:+.2e} +/- {np.sqrt(cov[0,0]):.2e} '
      f'({abs(slope)/np.sqrt(cov[0,0]):.1f} sigma from flat), reduced chi^2 = {chi2:.2f}')

# --- beta(q) over every bin ---
bq, bb, bad = [], [], []             # usable bin numbers, their beta, rejected bins
for qi in range(len(q)):
    lvl, _, c2 = flat_level(tau, G[:, qi], E[:, qi])
    if c2 > CHI2_FLAT_MAX:
        bad.append(qi)
    else:
        bq.append(qi)
        bb.append(lvl - 1.0)
bq = np.array(bq)
bb = np.array(bb)
print(f'  beta(q): {len(bq)} usable bins, {len(bad)} rejected as non-flat '
      f'(bins {[i+1 for i in bad]})')

# --- the adopted value: the mean over the bins the sample is analysed in ---
# Unweighted, because the per-bin fit errors (~1e-5) are far smaller than the
# bin-to-bin spread and so do not describe it: the multi-tau delay points share
# frames and are strongly correlated, which makes the fit error optimistic.  The
# spread of the five is the honest measure, and its standard error is quoted.
bsample = []
csample = []
for i in SAMPLE_Q:
    lvl, _, c2 = flat_level(tau, G[:, i], E[:, i])
    bsample.append(lvl - 1.0)
    csample.append(c2)
bsample = np.array(bsample)
beta = bsample.mean()
beta_err = bsample.std(ddof=1) / np.sqrt(len(bsample))
levels = ', '.join(f'{v:.5f}' for v in bsample)
chi2s = ', '.join(f'{v:.2f}' for v in csample)
print(f'  bins used for the sample: {levels}  (chi2 {chi2s})')
print(f'  ==> beta = {beta:.5f} +/- {beta_err:.5f} '
      f'(mean of the five, standard error of that mean)')
print(f'  for comparison, bin 16 at q = {q[15]:.5f} A^-1 gives '
      f'{flat_level(tau, G[:, 15], E[:, 15])[0] - 1:.5f}')

# --- figure ---
apply_style()
fig, (axg, axq) = plt.subplots(1, 2, figsize=FIG_SIZE)
label_panels((axg, axq))

m = tau > 0
axg.errorbar(tau[m], gbar[m], yerr=ebar[m], fmt='o', color='k', mfc='none',
             ms=MS_DENSE, mew=LW_THIN, capsize=1.5, elinewidth=LW_THIN,
             capthick=LW_THIN, zorder=2)
axg.plot(tau[m], intercept + slope * np.log10(tau[m]), 'r-', lw=LW_DATA, zorder=3)
axg.set_xscale('log')
axg.set_xlabel(r'Delay Time, $\Delta t$ (s)')
axg.set_ylabel('$g_2$')
axg.set_ylim(1.10, 1.16)
axg.set_title(rf'$Q = {q[Q_INDEX]:.5f}\ \AA^{{-1}}$')
add_minor_grid(axg)
axg.text(0.03, 0.06, rf'$\beta(Q) = {level - 1:.4f}$', transform=axg.transAxes,
         ha='left', va='bottom', color='r')

# Filled markers for the five bins the adopted value is the mean of, open for
# the rest, so the reader can see which points the red line averages.
inb = np.isin(bq, SAMPLE_Q)
axq.plot(q[bq[~inb]], bb[~inb], 'ko', mfc='none', ms=MS, mew=MEW)
axq.plot(q[bq[inb]], bb[inb], 'ko', mfc='k', ms=MS, mew=MEW)
axq.axhline(beta, color='r', ls='-', lw=LW_DATA)
axq.axvspan(0.00376, 0.00827, color='0.88', zorder=0)
# Left-aligned ON the left edge of the shaded band, not centred on it: the
# label is wider than the band, so centring pushed it 2.9 pt past the y axis.
axq.text(0.00376, 0.1465, 'range used\nfor the sample', ha='left', va='top', color='0.35')
axq.text(0.985, 0.06, rf'$\beta = {beta:.4f}$', transform=axq.transAxes,
         ha='right', va='bottom', color='r')
axq.set_xlabel(r'$Q$ ($\AA^{-1}$)')
axq.set_ylabel(r'Contrast, $\beta$')
axq.set_ylim(0.10, 0.15)
add_minor_grid(axq)

fig.tight_layout(pad=0.4, w_pad=1.4)
save_fig(fig, 'FigureS7_Contrast.pdf')
plt.show()
