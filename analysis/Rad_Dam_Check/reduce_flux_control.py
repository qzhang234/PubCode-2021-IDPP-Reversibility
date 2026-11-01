"""Reduce the four flux-control measurements behind Figure S13 into data/.

READS /gdata, SO IT RUNS ON ADAMITE (10.54.116.17; see CLAUDE.md).  Its input is
the beamline's own correlation of each measurement, the Multitau results of
experiment qzhang202503 (8-ID-I, 2025-1), copied from the experiment directory
to /gdata/s8id-dmdtn/2025-1/qzhang202503/Multitau on 1 Nov 2026.  Nothing is read from
/home/8-id-i.  The raw frames were archived off /gdata by DM
(/gdata/dm/8ID/8IDI/2025-1/qzhang202503/DM_ARCHIVE_README.txt), so these
results cannot be re-correlated.

flux_control.py reads only what this writes, so `make figures` never runs it.

THE DATA
--------
A 20 mg/mL (VPAVG)30 sample at four attenuator settings, 100,000 frames each:

    key   acquisition
    a7    E0174_H06-c6c5-4_a0007_f100000_r00001
    a9    E0173_H06-c6c5-4_a0009_f100000_r00001
    a11   E0172_H06-c6c5-4_a0011_f100000_r00001
    a14   E0171_H06-c6c5-4_a0014_f100000_r00001

THE BEAMLINE MASK.  The correlator ran against rigaku3m_qmap_S360_D36_Lin_refined.h5
(36 dynamic and 360 static Q bins in one phi sector), whose mask leaves in a
few bad pixels.  Two things follow, and this script corrects both:

  * HOT PIXELS.  Nine to twelve pixels per measurement, mostly the same ones in
    all four, count far above their Q bin; the worst, (1086, 1822), records
    about 35,000 photons against a bin mean of about 14 and alone raises I(Q)
    at 0.0114 A^-1 by a factor of 2.1-3.3.  I(Q) is therefore re-averaged from the
    time-averaged image the correlator stores
    (/xpcs/temporal_mean/scattering_2d), over the same static bins, leaving out
    every pixel with more than HOT_FACTOR times its bin's mean count and at least
    HOT_MIN photons.  In bins without such a pixel this reproduces the stored
    I(Q) to float precision, which the script checks.
  * THE ONE-FRAME DELAY.  At the shortest delay, 20 us, g2 at the lowest Q
    reads 1.5-1.8 where every later delay continues smoothly from about 1.2.
    The correlator stores no per-pixel G2, so the pixels responsible cannot be
    removed; the first delay is dropped.  No hot pixel lies in that dynamic
    bin, which holds 119 pixels more than under the q-map below.

THE OVERLEAF FIGURE used a later correlation of the same frames against
rigaku_s360_d36_lin_MC.hdf, whose mask excludes every one of the hot pixels and
whose g2 runs smoothly through the first delay; its results reached this
repository from Archive/New/ and are not on /gdata.  That q-map
puts the bins on 75 um pixels instead of 76 um, so its Q values are 1.3 % lower
(the g2 bin of the figure was 0.0030 A^-1, here 0.0031).  Checked against it on
Amber, 1 Nov 2026, for all four measurements: the ion-chamber readings, frame
time and delays are identical; bin for bin, I(Q) agrees within 4 % (median
ratio 1.000), and g2 at the lowest Q within 0.023 at every delay from the
second on (median 0.003).  The flux in the legend and the conclusion drawn
from the figure are unchanged.

    python reduce_flux_control.py
"""

import hashlib
import os

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SOURCE = '/gdata/s8id-dmdtn/2025-1/qzhang202503/Multitau'
OUT = os.path.join(HERE, 'data', 'flux_control.h5')
FILES = {'a7':  'E0174_H06-c6c5-4_a0007_f100000_r00001_results.hdf',
         'a9':  'E0173_H06-c6c5-4_a0009_f100000_r00001_results.hdf',
         'a11': 'E0172_H06-c6c5-4_a0011_f100000_r00001_results.hdf',
         'a14': 'E0171_H06-c6c5-4_a0014_f100000_r00001_results.hdf'}
HOT_FACTOR, HOT_MIN = 20, 20      # hot: > 20 x the bin's mean count and >= 20 photons
FIRST_DELAY = 1                   # drop the one-frame delay

IC = '/entry/instrument/incident_beam/incident_beam_intensity'
G2, G2_ERR = '/xpcs/multitau/normalized_g2', '/xpcs/multitau/normalized_g2_err'
DELAY, FRAME_TIME = '/xpcs/multitau/delay_list', '/entry/instrument/detector_1/frame_time'
SAXS, SAXS_2D = '/xpcs/temporal_mean/scattering_1d', '/xpcs/temporal_mean/scattering_2d'
INT_T = '/xpcs/spatial_mean/intensity_vs_time'
S_MAP, S_Q, D_Q = ('/xpcs/qmap/static_roi_map', '/xpcs/qmap/static_v_list_dim0',
                   '/xpcs/qmap/dynamic_v_list_dim0')
QMAP = '/xpcs/multitau/config/qmap'


def md5(path):
    with open(path, 'rb') as f:
        return hashlib.md5(f.read()).hexdigest()


def bin_mean(values, label, use, nbins):
    """Mean of values over the pixels in use, for static bins 1..nbins."""
    s = np.bincount(label[use], values[use], nbins + 1)
    n = np.bincount(label[use], minlength=nbins + 1)
    return s[1:] / np.maximum(n[1:], 1)


def hot_pixels(counts, label, nbins):
    """Pixels with more than HOT_FACTOR times their bin's mean count and at least
    HOT_MIN photons.  Repeated with the flagged pixels left out of the means, so
    that one very hot pixel does not shelter another in its bin."""
    hot = np.zeros(label.shape, dtype=bool)
    while True:
        mean = np.r_[0, bin_mean(counts, label, (label > 0) & ~hot, nbins)]
        new = (label > 0) & (counts > np.maximum(HOT_FACTOR * mean[label], HOT_MIN))
        if np.array_equal(new, hot):
            return hot
        hot = new


def reduce_one(path):
    with h5py.File(path, 'r') as hf:
        label = hf[S_MAP][()].astype(np.int64)            # static bin + 1, 0 = masked
        image = np.squeeze(hf[SAXS_2D][()]).astype(np.float64)
        stored = np.squeeze(hf[SAXS][()])
        nframes = hf[INT_T].shape[1]
        nbins = stored.size
        hot = hot_pixels(image * nframes, label, nbins)
        saxs = bin_mean(image, label, (label > 0) & ~hot, nbins)
        clean = np.setdiff1d(np.arange(nbins), label[hot] - 1)
        assert np.allclose(saxs[clean], stored[clean], rtol=1e-5), 'image does not reproduce I(Q)'
        t0 = float(np.ravel(hf[FRAME_TIME][()])[0])
        data = dict(tau=(np.ravel(hf[DELAY][()]) * t0)[FIRST_DELAY:],
                    q_static=hf[S_Q][()], q_dynamic=hf[D_Q][()], saxs_1d=saxs,
                    g2=hf[G2][()][FIRST_DELAY:], g2_err=hf[G2_ERR][()][FIRST_DELAY:],
                    hot_pixels=np.argwhere(hot))
        attrs = dict(qmap=hf[QMAP][()].decode(), ic_readout=float(np.ravel(hf[IC][()])[0]),
                     frame_time=t0, n_frames=nframes,
                     max_saxs_change=float(np.max(np.abs(stored / saxs - 1))))
    return data, attrs


def main():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with h5py.File(OUT, 'w') as h:
        h.attrs.update(source=SOURCE, written_by='reduce_flux_control.py',
                       hot_pixel_rule=f'> {HOT_FACTOR} x bin mean count and >= {HOT_MIN} photons',
                       delays_dropped=FIRST_DELAY)
        for key, name in FILES.items():
            path = os.path.join(SOURCE, name)
            data, attrs = reduce_one(path)
            g = h.create_group(key)
            g.attrs.update(acquisition=name.replace('_results.hdf', ''), source=path,
                           md5=md5(path), **attrs)
            for k, v in data.items():
                g[k] = v
            for k, v in (('tau', 's'), ('q_static', 'A^-1'), ('q_dynamic', 'A^-1'),
                         ('saxs_1d', 'photons/pixel/frame')):
                g[k].attrs['units'] = v
            print(f"{key:>4}: {g.attrs['acquisition']}  md5 {g.attrs['md5']}\n"
                  f"      IC readout {attrs['ic_readout']:.4g}, {len(data['hot_pixels'])} hot pixels, "
                  f"which moved I(Q) by up to {attrs['max_saxs_change']:.0%}")
    print(f'  -> {os.path.relpath(OUT, HERE)}')


if __name__ == '__main__':
    main()
