"""Absolute scattering cross-section calibration (shared helper).

saxpcs.py (Figure S5), thermal_cycle.py (Figure S4b) and
temperature_series_2021-3.py (Figure 3a and Figures S7-S10a) all import from
here, so they use identical constants and one per-file coefficient function.

PROVENANCE: worked out originally in two beamline notebooks, IC_pind4_convert
and Abs_Scatt_Cross, which are not part of this repository.  The code below is
the authoritative implementation and reproduces their numbers.

GOAL.  Turn the azimuthally averaged I(Q), in photons per pixel per frame, into
an absolute differential cross section with the buffer subtracted:

    I_abs(Q) = coef_sam * I_sample(Q) - coef_buf * I_buffer(Q)

Following Sheyfer et al., Phys. Rev. Lett. 125, 125504 (2020), Eq. S1, the
coefficient of one measurement is

    coef = num_frames / (t_exp * F * T * f * dOmega)

    F      incident flux              [photons / s]
    T      transmission               [dimensionless]
    t_exp  total exposure time        [s] = frame_time * num_frames
    f      sample thickness
    dOmega solid angle of one pixel   = (pixel / distance)^2

F and T come from the upstream and downstream ion chambers stored per run at
    /entry/instrument/incident_beam/incident_beam_intensity     (Up_IC)
    /entry/instrument/incident_beam/transmitted_beam_intensity  (Dn_IC)
average_ranges.py averages both over each frame range, so every averaged file
carries the mean Up_IC and Dn_IC of the frames that went into it, and each
curve is scaled by ITS OWN flux and transmission.

STEP 1, ion chamber to photon flux.  With no sample in the beam, Up_IC was
recorded at eight attenuation ratios alongside the downstream PIN diode.  A
straight line photons = CAL_A * Up_IC + CAL_B is fitted to seven of the eight;
the unattenuated point is dropped (CAL_CROP = 1) because it lies off the line
the others define.  Keeping it leaves residuals of -2 to -148 % across the set,
dropping it leaves about 1 % over the middle five.  A real run then has
F = CAL_A * (Up_IC / t_exp) + CAL_B.

STEP 2, air transmission.  With no sample, (Dn_IC - dark) / (Up_IC - dark) is
the transmission of the air path between the chambers, flat at about 0.868.  So
the sample transmission is T = Dn_IC / Up_IC / AIR_TRANSMISSION.

STEP 3, coefficient.  abs_xsec_coef() evaluates coef for one open file from its
range-averaged monitors and the constants below.

THICKNESS IS 2 mm, the cell used in this experiment, and every number in the
paper rests on it; coef_buf = 3.75e4 for the D0138 buffer.  Setting it to 1 mm
reproduces the original notebook, which had been written for a 1 mm cell, to
2 %, which is how this implementation was checked.  Thickness enters only as
1 / f, so the two differ by exactly a factor of two.

UNITS.  Lengths are in millimetres and dOmega is a pure ratio, so coef * I(Q)
comes out in mm^-1 sr^-1 and that is what every absolute-scale figure plots,
with no further scaling.  The literature more often quotes cm^-1, a factor of
ten, which is deliberately not applied anywhere here.
"""

import numpy as np


# --- fixed beamline constants (from the experiment log / conversation notes) ---
DET_DIST_MM      = 7800.0    # sample-to-detector distance, 7.8 m
PIXEL_MM         = 0.076     # detector pixel size, 76 um (same in x and y)
SAMPLE_THICKNESS = 2.0       # sample / capillary thickness, 2 mm
NUM_FRAMES       = 100000    # frames per measurement
FRAME_TIME       = 20e-6     # acquisition (frame) time, 20 us
T_EXP            = FRAME_TIME * NUM_FRAMES          # total exposure time [s]
DELTA_OMEGA      = (PIXEL_MM / DET_DIST_MM) ** 2    # pixel solid angle

# --- IC -> photon calibration data (no sample; only the air gap attenuates) ---
UPIC_DARK   = 83.6           # upstream IC dark reading
DNIC_DARK   = 204.4          # downstream IC dark reading
# Upstream IC, downstream IC and PIN-diode photon counts at 8 attenuation ratios
# (dark already subtracted from the photon counts).
CAL_UPIC    = np.array([257609, 125204, 55586.5, 27036.8, 12296.8, 6028.8, 2714.4, 1361.1])
CAL_DNIC    = np.array([223888, 108873, 48413.4, 23613.7, 10807.5, 5367.4, 2485.5, 1310.9])
CAL_PHOTONS = np.array([1.49e10, 8.12e9, 3.75e9, 1.97e9, 1.06e9, 6.71e8, 4.66e8, 3.83e8]) - 3.04e8
CAL_CROP    = 1              # drop the first (unattenuated) point before fitting

# STEP 1: linear fit  photons = CAL_A * Up_IC + CAL_B
CAL_A, CAL_B = np.polyfit(CAL_UPIC[CAL_CROP:], CAL_PHOTONS[CAL_CROP:], 1)

# STEP 2: air transmission = mean of (Dn_IC - dark) / (Up_IC - dark).  The
# spread of the per-measurement values (AIR_TRANS_STD) is a fair error-from-mean
# estimate and is reported in the calibration figure legend.
AIR_TRANS_SERIES = (CAL_DNIC[CAL_CROP:] - DNIC_DARK) / (CAL_UPIC[CAL_CROP:] - UPIC_DARK)
AIR_TRANSMISSION = float(np.mean(AIR_TRANS_SERIES))
AIR_TRANS_STD    = float(np.std(AIR_TRANS_SERIES, ddof=1))   # sample standard deviation

# HDF fields holding the (range-averaged) ion-chamber monitors.
INCIDENT_PATH    = '/entry/instrument/incident_beam/incident_beam_intensity'
TRANSMITTED_PATH = '/entry/instrument/incident_beam/transmitted_beam_intensity'


def abs_xsec_coef(hf):
    """Absolute cross-section coefficient for one open XPCS/SAXS file.

    Reads the range-averaged upstream (Up_IC) and downstream (Dn_IC)
    ion-chamber intensities and combines them with the beamline constants and
    the IC->photon calibration (see the module docstring):

        F    = CAL_A * (Up_IC / T_EXP) + CAL_B          incident flux [photons/s]
        T    = Dn_IC / Up_IC / AIR_TRANSMISSION         transmission  [-]
        coef = NUM_FRAMES / (T_EXP * F * T * SAMPLE_THICKNESS * DELTA_OMEGA)

    Returns the scalar coefficient that multiplies this file's I(Q).
    """
    # These are stored as tiny arrays whose shape varies between files -- (1,)
    # in some, (1, 1) in others -- so flatten to 1-D and take the first entry.
    up_ic = float(np.asarray(hf[INCIDENT_PATH][()]).ravel()[0])
    dn_ic = float(np.asarray(hf[TRANSMITTED_PATH][()]).ravel()[0])
    flux = CAL_A * (up_ic / T_EXP) + CAL_B
    transmission = dn_ic / up_ic / AIR_TRANSMISSION
    return NUM_FRAMES / (T_EXP * flux * transmission * SAMPLE_THICKNESS * DELTA_OMEGA)


# --- 2021-2: incident chamber recorded, downstream chamber not ---------------
# The 2021-2 archive records the upstream chamber as source_end/I0Monitor and
# never recorded the downstream one (see restore_dx_metadata_2021-2.py), so the
# incident flux follows from the ordinary calibration but the transmission has
# no Dn_IC to come from.  It is taken from 2022-1 instead, on three grounds the
# beamtimes share: CAL_A and CAL_B describe the chambers and their electronics
# rather than the run; the cell is the same 2 mm aluminium block behind the same
# 0.127 mm polycarbonate, so the window and path absorption are unchanged; and
# the sample is the same construct at the same nominal concentration in the same
# buffer, so its own absorption barely moves.
#
# The number is the mean over the two 2022-1 ELP loadings, which is also the
# honest statement of its uncertainty: the S1 cycling loading transmits 0.3497
# (0.3461-0.3526 over 21 averaged files) and the isothermal loading 0.3877
# (0.3859-0.3896 over 10), so two preparations of the same material in the same
# cell differ by 11%.  A third preparation is bracketed by them and not pinned.
#
# What that costs: the absolute scale of any 2021-2 intensity carries about 6%
# of systematic from this choice alone, ON TOP of the calibration's own.  It is
# a SCALE FACTOR.  It moves no peak, no exponent and no ratio within a profile,
# and it is common to every 2021-2 curve, so it cannot move the comparison
# between chain lengths or between temperatures, which is what the ramp figures
# are read for.  Quote 2021-2 intensities to one significant figure.
TRANSMISSION_2021_2 = 0.369
TRANSMISSION_2021_2_RANGE = (0.3461, 0.3896)


def abs_xsec_coef_2021_2(hf, transmission=TRANSMISSION_2021_2):
    """Absolute cross-section coefficient for a 2021-2 file.

    Identical to abs_xsec_coef except that the transmission is supplied rather
    than measured, because 2021-2 has no downstream ion chamber.  Everything
    else -- the IC-to-photon calibration, the thickness, the pixel solid angle,
    the frame count -- is the same constant it is for every other beamtime here.
    """
    up_ic = float(np.asarray(hf[INCIDENT_PATH][()]).ravel()[0])
    if not np.isfinite(up_ic) or up_ic <= 0:
        raise ValueError('2021-2 incident chamber missing; run '
                         'restore_dx_metadata_2021-2.py --apply first')
    flux = CAL_A * (up_ic / T_EXP) + CAL_B
    return NUM_FRAMES / (T_EXP * flux * transmission * SAMPLE_THICKNESS * DELTA_OMEGA)


def calibration_summary():
    """One-line-per-item summary string for logging at import in the scripts."""
    return (f'IC->photon calibration : photons = {CAL_A:.3e} * Up_IC {CAL_B:+.3e}\n'
            f'air transmission        : {AIR_TRANSMISSION:.4f} +/- {AIR_TRANS_STD:.4f}')
