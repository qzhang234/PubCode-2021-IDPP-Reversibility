"""Figure S10: the ion-chamber calibration behind every absolute intensity here.

(a) PIN-diode photon count against the upstream ion-chamber readout across the
    attenuation series, with the straight line whose slope and intercept are
    CAL_A and CAL_B;
(b) transmission of the air path between the two ion chambers across the same
    series, and its mean.

READS NO DATA FILES.  Everything plotted is a constant in abs_xsec.py, where
the eight calibration points are recorded and CAL_A, CAL_B and
AIR_TRANSMISSION are derived from them, the same numbers that put every
absolute I(Q) of this paper on its scale.  The panel order is that of the
Overleaf Figure S9; the companion repository draws the same two panels
(PubCode-2021-IDPP, analysis/SAXPCS_8id/calibration_figure.py).
"""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..'))
from common.acs_style import (DOUBLE_COL, MS, MEW, LW_THIN, LW_DATA,      # noqa: E402
                              apply_style, add_minor_grid, label_panels, save_fig)
from abs_xsec import (AIR_TRANS_SERIES, AIR_TRANSMISSION, AIR_TRANS_STD,  # noqa: E402
                      CAL_UPIC, CAL_PHOTONS, CAL_CROP, CAL_A, CAL_B,
                      calibration_summary)

print(calibration_summary())

apply_style()
fig, (axc, axt) = plt.subplots(1, 2, figsize=(DOUBLE_COL, 2.9))

# (a) the unattenuated point lies off the line the others define and is left
# out of the fit (CAL_CROP), so it is left out here too
up = CAL_UPIC[CAL_CROP:]
axc.plot(up, CAL_PHOTONS[CAL_CROP:] / 1e9, 'ko', ms=MS, mew=MEW, label='measured')
line = np.linspace(0, up.max() * 1.05, 100)
axc.plot(line, (CAL_A * line + CAL_B) / 1e9, 'b-', lw=LW_DATA,
         label=f'fit: {CAL_A:.3g} $\\times$ IC {CAL_B:+.3g}')
axc.set_xlabel('Upstream ion chamber (counts)')
axc.set_ylabel(r'Photons ($10^9$)')
add_minor_grid(axc)
axc.legend(loc='upper left')

# (b) air-path transmission over the same seven points
axt.plot(np.arange(1, len(AIR_TRANS_SERIES) + 1), AIR_TRANS_SERIES, 'ko:', ms=MS, mew=MEW,
         lw=LW_THIN)
axt.axhline(AIR_TRANSMISSION, color='b', ls='-', lw=LW_DATA,
            label=f'mean {AIR_TRANSMISSION:.4f} $\\pm$ {AIR_TRANS_STD:.4f}')
axt.set_xlabel('Measurement')
axt.set_ylabel('Transmission of air')
add_minor_grid(axt)
axt.legend(loc='lower left')

label_panels([axc, axt])
fig.tight_layout(pad=0.4, w_pad=1.6)
save_fig(fig, 'FigureS10_IC_Calibration.pdf')
