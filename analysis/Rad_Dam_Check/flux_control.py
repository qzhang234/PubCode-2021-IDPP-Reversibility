"""Figure S13: SAXS and XPCS of diluted (VPAVG)30 at four incident fluxes.

After the Overleaf Figure S12, and the same figure the companion repository
draws (PubCode-2021-IDPP, analysis/Rad_Dam_Check/g2_SAXPCS_Rad_Cali.py).  A
20 mg/mL (VPAVG)30 sample, measured after the APS Upgrade at four attenuator
settings, E0171-E0174.  Reads data/flux_control.h5, which reduce_flux_control.py
writes from the beamline's correlation results on /gdata; how those
differ from the correlation behind the Overleaf figure is in its docstring.

The legend gives the photon flux on the sample, recovered for each measurement
from its upstream ion-chamber reading, because the nominal attenuator settings
(a7, a9, a11, a14) are not monotonic in flux; the Overleaf legend listed the
nominal attenuation ratios.

    IC readout (entry/instrument/incident_beam/incident_beam_intensity, x 1e-9)
        --(Pind_calibration.csv, fit through the origin)--> calibrated PIN current [A]
        --(Flux_Cal.py formula)---------------------------> photon flux [ph/s]
"""

import os
import sys

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
from common.acs_style import (DOUBLE_COL, MS, MEW, LW_THIN,                  # noqa: E402
                              apply_style, add_minor_grid, label_panels,
                              q_log_ticks, save_fig)

DATA = os.path.join(HERE, 'data', 'flux_control.h5')
COLORS = {'a7': 'C0', 'a9': 'C1', 'a11': 'C2', 'a14': 'C3'}
MARKERS = {'a7': 'o', 'a9': 's', 'a11': '^', 'a14': 'D'}
IC_SCALE = 1e-9                           # ion-chamber readout -> A

# Flux_Cal.py, without its attenuation-ratio factor: the IC reading here already
# corresponds to the attenuated beam on the sample.
SPEC_RES = 0.1882                         # Si PIN-diode responsivity [A/W]
PHOTON_ENERGY = 12.4e3                    # eV
EV2J = 1.6e-19                            # J per eV
LOSS_FACTOR = 0.94 * 0.94 * 0.735         # two mirrors, BPM and CVD diamond window


def pind_to_flux(current):
    """Calibrated PIN-diode current [A] -> photon flux on the sample [ph/s]."""
    return current / SPEC_RES / (PHOTON_ENERGY * EV2J) / LOSS_FACTOR


def ic_to_pind_slope():
    """Least-squares slope, through the origin, of calibrated PIN against upstream IC."""
    df = pd.read_csv(os.path.join(HERE, 'Pind_calibration.csv'))
    ic = pd.to_numeric(df['Up IC (120 nA, 1e9 scale)'], errors='coerce')
    pind = pd.to_numeric(df['Calibrated pind'], errors='coerce')
    ok = ic.notna() & pind.notna() & (ic > 0) & (pind > 0)
    x, y = ic[ok].values, pind[ok].values
    return np.sum(x * y) / np.sum(x * x)


slope = ic_to_pind_slope()
data = {}
with h5py.File(DATA, 'r') as h:
    for key in COLORS:
        g = h[key]
        ic = g.attrs['ic_readout'] * IC_SCALE
        data[key] = dict(flux=pind_to_flux(slope * ic), q=g['q_static'][()], I=g['saxs_1d'][()],
                         tau=g['tau'][()], g2=g['g2'][:, 0], g2_err=g['g2_err'][:, 0],
                         q_dyn=g['q_dynamic'][0])
        print(f"{key:>4}: IC {ic:.3e} A -> flux {data[key]['flux']:.2e} ph/s")


def flux_label(f):
    e = int(np.floor(np.log10(f)))
    return rf'${f / 10 ** e:.1f}\times10^{{{e}}}$ ph/s'


apply_style()
fig, axs = plt.subplots(1, 2, figsize=(DOUBLE_COL, 2.9))
for key in sorted(data, key=lambda k: data[k]['flux']):       # key in order of flux
    d, c, m = data[key], COLORS[key], MARKERS[key]
    axs[0].plot(d['q'], d['I'], ls='none', marker=m, mfc='none', mec=c, mew=MEW, ms=MS,
                label=flux_label(d['flux']))
    axs[1].errorbar(d['tau'], d['g2'], d['g2_err'], ls='none', marker=m, mfc='none',
                    mec=c, color=c, mew=MEW, ms=MS, elinewidth=LW_THIN, capsize=1.5,
                    capthick=LW_THIN, label=flux_label(d['flux']))
axs[0].set_xscale('log')
axs[0].set_yscale('log')
q_log_ticks(axs[0])
axs[0].set_xlabel(r'$Q$ ($\AA^{-1}$)')
axs[0].set_ylabel('$I(Q)$ (photons pixel$^{-1}$ frame$^{-1}$)')
add_minor_grid(axs[0])
axs[0].legend(title='Flux on sample', loc='lower left')
axs[1].set_xscale('log')
axs[1].set_xlabel(r'$\Delta t$ (s)')
axs[1].set_ylabel(rf"$g_2(Q = {data['a7']['q_dyn']:.4f}\ \AA^{{-1}},\ \Delta t)$")
add_minor_grid(axs[1])
axs[1].legend(title='Flux on sample')
label_panels(axs)
fig.tight_layout(pad=0.4, w_pad=1.4)
save_fig(fig, 'FigureS13_Flux_Control.pdf')
