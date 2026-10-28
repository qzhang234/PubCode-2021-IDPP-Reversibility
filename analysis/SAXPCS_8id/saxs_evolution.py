"""Figure S5: full SAXS 1d evolution for all B0147 files (D0138 subtracted).

Every B0147 averaged file is plotted, coloured by elapsed time (time origin =
first B0147 file, frames 1-200).  B0146 (6 C reference, before the 30 C
isothermal) is drawn in the same blue square as the 6 C reference in Figure 3a,
and the y axis carries the same label, so the two figures key and label that
dataset identically.

Elapsed time is keyed by a COLOURBAR rather than by an 11-entry legend box.  The
curves are stacked in y across the whole q range, so there is no corner of the
axes a box that size can sit in without covering data; the colourbar carries the
same information outside the frame.  The colormap is plasma, the same
elapsed-time colour language used in the main XPCS figure (and, unlike the turbo
scale used before, it is legible with colour-vision deficiency) -- but note the
two figures normalise over different time spans, so read the value off this bar
rather than matching hues between figures.
"""

import glob
import os
import sys

import h5py
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

# Absolute scattering cross-section calibration lives in the shared module so
# this figure uses exactly the same constants and per-file coefficient as the
# main saxpcs.py figure.  See abs_xsec.py for the full derivation.
from abs_xsec import abs_xsec_coef, calibration_summary

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from common.acs_style import (SINGLE_COL, MS, MEW, LW_THIN,
                              apply_style, add_minor_grid, q_log_ticks, save_fig)

# --- PARAMETERS ---
data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
BACKGROUND_HEADER = 'D0138'
# Absolute scattering cross-section: convert the raw SAXS 1d to an absolute
# differential cross section (mm^-1) as
#     I_abs(q) = coef_sam * I_sample(q) - coef_buf * I_buffer(q)
# coef_sam / coef_buf are NOT hard-coded any more: each is computed by
# abs_xsec_coef() from that file's own range-averaged ion-chamber readings, so a
# drift in incident flux across the time series is handled per file.  coef_buf
# is the D0138 buffer's own coefficient; coef_sam is recomputed for every curve.
PHI_AVERAGE = True
CMAP = plt.cm.plasma
COLOR_6C = '#1f77b4'                 # identical to COLOR_6C in saxpcs.py

# Field locations and readers are shared with the other 8-ID figure scripts.
from nexus_read import parse_name, read_start_time, read_saxs_iq


# --- DISCOVER FILES ---
# every range average in data/, including the thermal-cycle groups that belong to
# Figure S4; only XPCS_HEADER and the 6 C reference are kept below
file_paths = sorted(glob.glob(os.path.join(data_dir, 'Average_*.hdf')))
by_header, start_times = {}, {}
for fp in file_paths:
    header = parse_name(fp)[0]
    with h5py.File(fp, 'r') as hf:
        start_times[fp] = read_start_time(hf)
    by_header.setdefault(header, []).append(fp)
for h in by_header:
    by_header[h].sort(key=lambda p: parse_name(p)[1])

b0147 = by_header['B0147']
t_ref = start_times[b0147[0]]
elapsed = {fp: (start_times[fp] - t_ref).total_seconds() for fp in b0147}
norm = Normalize(vmin=0, vmax=max(elapsed.values()))

print(calibration_summary())

# Buffer curve AND its own absolute coefficient (coef_buf), both from D0138.
bg_I = None
coef_buf = None
if BACKGROUND_HEADER in by_header:
    with h5py.File(by_header[BACKGROUND_HEADER][0], 'r') as hf:
        _, bg_I = read_saxs_iq(hf, PHI_AVERAGE)
        coef_buf = abs_xsec_coef(hf)
    print(f'coef_buf ({BACKGROUND_HEADER}) = {coef_buf:.3e}')

# --- FIGURE ---
apply_style()
fig, ax = plt.subplots(figsize=(SINGLE_COL, 3.3))

for fp in by_header.get('B0146', []):
    with h5py.File(fp, 'r') as hf:
        q, I = read_saxs_iq(hf, PHI_AVERAGE)
        coef_sam = abs_xsec_coef(hf)                       # this file's own coefficient
    # Put the profile on the absolute scale, then take the buffer off.  The
    # buffer already carries its own coefficient, so both terms are absolute.
    I = coef_sam * I
    if bg_I is not None:
        I = I - coef_buf * bg_I
    pos = I > 0
    # same colour and marker as the 6 C reference in Figure 3a, so the two
    # figures key that dataset identically
    ax.plot(q[pos], I[pos], color=COLOR_6C, marker='s', ls='none', ms=MS,
            mfc='none', mew=MEW, label='6 °C ref')

for fp in b0147:
    with h5py.File(fp, 'r') as hf:
        q, I = read_saxs_iq(hf, PHI_AVERAGE)
        coef_sam = abs_xsec_coef(hf)                       # this file's own coefficient
    # Put the profile on the absolute scale, then take the buffer off.  The
    # buffer already carries its own coefficient, so both terms are absolute.
    I = coef_sam * I
    if bg_I is not None:
        I = I - coef_buf * bg_I
    pos = I > 0
    ax.plot(q[pos], I[pos], color=CMAP(norm(elapsed[fp])), marker='o', ls='none',
            ms=MS, mfc='none', mew=MEW)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$Q$ ($\AA^{-1}$)')
q_log_ticks(ax)   # same three named Q ticks as Figure 3a
# same axis label as Figure 3a: the quantity is an absolute differential cross
# section, and the units are on the axis (see the note in saxpcs.py)
ax.set_ylabel(r'$I(Q)$ (mm$^{-1}$)')
add_minor_grid(ax)
# only the 6 C reference needs a legend entry; elapsed time is the colourbar
ax.legend(loc='lower left', borderaxespad=0.4)

cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=CMAP), ax=ax, pad=0.03)
cb.set_label('Elapsed time (s)')
cb.outline.set_linewidth(LW_THIN)
cb.ax.tick_params(width=LW_THIN)
# a tick on the bar at every acquisition, so the reader can see that the times
# are unevenly spaced and where the Figure 3 subset (5040-7863 s) sits
cb.ax.hlines(sorted(elapsed.values()), 0, 1, colors='w', lw=LW_THIN)

fig.tight_layout(pad=0.4)
save_fig(fig, 'FigureS5_SAXS_Evolution.pdf')
plt.show()
