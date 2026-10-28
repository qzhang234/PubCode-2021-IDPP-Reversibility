"""Helpers for reading and averaging 8-ID-I XPCS cluster-result files.

Only two functions here are on the path from the reduced data in this
repository to a figure -- ``outlier_removal`` and ``average_datasets``, both
called by SAXPCS_8id/average_ranges.py.  They are pure array code and need
nothing but numpy.

Everything else is the beamline-side tooling that produced the raw results in
the first place: ``read_keys_from_files``, ``read_keys_from_files_parallel``,
``read_temperature_from_files``, ``process_group`` and
``process_group_by_range``.  Those walk a directory of per-acquisition result
files on the beamline storage, and the first two also import ``pyxpcsviewer``
lazily, which is deliberately NOT in environment.yml.  They are kept for
provenance and for reprocessing at the beamline; they cannot run from a clone
of this repository alone, and no figure script calls them.  None of them
carries a default data directory -- the caller must pass the ``prefix`` that
applies to the beamtime being reduced.
"""

import numpy as np
import os
import h5py
import glob
import matplotlib.pyplot as plt
from multiprocessing import Pool
import logging



logger = logging.getLogger(__name__)


def split(arr, n_segments: int):
    """split arr into n parts

    Parameters
    ----------
    arr : list
        input part
    n_segments : int 
        number of segments

    Returns
    -------
    list
        a list of n_segments lists
    """
    k, m = divmod(len(arr), n_segments)
    return list(arr[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] 
                for i in range(n_segments))


def read_keys_from_files_parallel(flist_list, num_cores=24):
    """
    Read keys from a list of files and process them in parallel
    """  
    with Pool(num_cores) as p:
        result = p.map(read_keys_from_files, flist_list)
        
    result = [x for x in result if x is not None]
    return result

    
def read_keys_from_files(flist, keys=('g2', 'g2_err', 'saxs_1d')):
    """Read the named datasets out of each beamline result file.

    BEAMLINE-ONLY: needs pyxpcsviewer and the raw per-acquisition files, so it
    cannot run from a clone of this repository.

    Parameters
    ----------
    flist : list of str
        Paths to per-acquisition XPCS result files.
    keys : tuple of str
        Which datasets to pull out of each file.

    Returns
    -------
    dict
        One key per entry in ``keys``; each value is an array with one row per
        file, in the order ``flist`` was given.  A 2-D ``saxs_1d`` is collapsed
        to one curve by averaging over its rows, ignoring NaNs.
    """
    from pyxpcsviewer import XpcsFile as XF

    def get_field(xf_obj, key):
        if key == 'saxs_1d':
            x = xf_obj.saxs_1d['Iq']
            assert x.ndim in [1, 2]
            if x.ndim == 2:
                try:
                    x = np.nanmean(x, axis=0)
                except:
                    print(xf_obj.saxs_1d['Iq'])
                    raise
                return x
        else:
            return getattr(xf_obj, key)

    data_dict = {}
    for key in keys:
        data_dict[key] = []

    for f in flist:
        try:
            dset = XF(f, fields=keys)
            for k in keys:
                data_dict[k].append(get_field(dset, k))
        except Exception:
            logger.info(f'failed to read file {f}, skip this file')
            pass
    
    for k in keys:
        data_dict[k] = np.array(data_dict[k])
        
    return data_dict


def average_datasets(flist=None, data_dict=None, mask=None):
    """Average the acquisitions a mask keeps, one averaged curve per dataset.

    Give it either an already-stacked ``data_dict`` (what average_ranges.py
    does) or a file list to read first.  Each array in ``data_dict`` has one row
    per acquisition; ``mask`` is the keep-mask from the outlier cuts.

    'g2' and 'saxs_1d' are averaged with nanmean.  'g2_err' is NOT: averaging N
    measurements whose individual uncertainties are sigma_i gives the mean an
    uncertainty of sqrt(sum sigma_i^2) / N, which is what the else-branch below
    computes.  N is counted per point (the number of non-NaN entries), so a
    point that is missing from some acquisitions still gets the right divisor.
    """
    if data_dict is None and flist is not None:
        keys=('g2', 'g2_err', 'saxs_1d')
        data_dict = read_keys_from_files(flist)
    elif data_dict is not None:
        keys = data_dict.keys()
    else:
        print('must provide flist or data_dict')
        raise

    average_dict = {}

    if mask is None:
        num_valid_files = data_dict['g2'].shape[0]
        mask = np.ones_like(num_valid_files, dtype=bool)
    else:
        num_valid_files = np.sum(mask)

    for k in keys:
        if k != 'g2_err':
            average_dict[k] = np.nanmean(data_dict[k][mask], axis=0)
        else:
            average_dict[k] = np.sqrt(np.nansum(data_dict[k][mask] ** 2, axis=0))
            average_dict[k] /= np.sum(~np.isnan(data_dict[k][mask]), axis=0)
        
    return average_dict


def apply_cross_corr_threshold(x0, percentile=5, style='linear', debug_fig_ax=None,
                               label='debug'):
    # Work on a copy: this routine rewrites NaNs and non-positive values (e.g.
    # x[x <= 0] = 1 below) purely to make the correlation metric well-defined.
    # If we mutated x0 in place it would corrupt data_dict, and those injected
    # 1.0's would later be averaged into saxs_1d as huge single-q spikes.
    x = np.array(x0, dtype=np.float64)
    x[np.isnan(x)] = 0
    if style == 'log':
        x[np.isnan(x)] = 1
        x[x <= 0] = 1
        x = np.log10(x)

    # Similarity between every pair of acquisitions: the cosine similarity of
    # their curves, i.e. the dot product divided by the two vector lengths.  It
    # is 1 for two curves of identical shape and drops as the shapes diverge;
    # dividing by the lengths is what makes it blind to overall level.
    num = x.shape[0]
    result = np.zeros((num, num), dtype=np.float64)
    for n in range(num):
        for m in range(n + 1, num):
            val = np.sum(x[n] * x[m]) / np.sqrt(np.sum(x[n] ** 2) * np.sum(x[m] ** 2))
            result[n, m] = val
            result[m, n] = val

    result_1d = np.sum(result, axis=1)
    low_cutoff = np.percentile(result_1d, percentile)
    if debug_fig_ax is not None:
        title = f'{label}_{style=}_{percentile=}'
        debug_fig_ax.plot(result_1d, 'o')
        debug_fig_ax.hlines(low_cutoff, xmin=-1, xmax=len(result_1d))
        debug_fig_ax.set_title(title)

    mask = result_1d >= low_cutoff
    return mask


def outlier_removal(data_dict, label='testrun', percentile=5, plot_debug=False):
    """Flag acquisitions whose curves do not look like the rest of the group.

    For each dataset in ``data_dict`` (g2, g2_err, saxs_1d), every acquisition
    is compared with every other by cross-correlation, giving each acquisition
    a similarity score against the group.  The lowest ``percentile`` of those
    scores is rejected.  An acquisition has to pass on ALL datasets to be kept.

    This cut sees the SHAPE of a curve, not its level: an acquisition that is
    anomalously bright in a few q bins but normal elsewhere stays close to the
    group mean by this measure and survives.  That is what the second,
    level-based cut in SAXPCS_8id/average_ranges.py is for.

    Parameters
    ----------
    data_dict : dict of str -> array
        One entry per dataset; each array has one row per acquisition.
    label : str
        Name used in the log line and in the debug figure filename.
    percentile : float
        Percentage of acquisitions rejected per dataset (5 means the least
        similar 5 %).
    plot_debug : bool
        Write ``debug/debug_outlier_<label>.png`` showing the scores.

    Returns
    -------
    numpy.ndarray of bool
        One entry per acquisition, True to keep.
    """
    mask_all = np.ones(len(data_dict['g2']), dtype=bool)
    num_features = len(data_dict.keys()) + 1

    if plot_debug:
        fig, ax = plt.subplots(num_features, 1, figsize=(4, 2.4 * num_features),
                               sharex=True)
    else:
        ax = [None for _ in range(num_features)]

    for n, key in enumerate(data_dict.keys()):
        if key == 'saxs_1d':
            style = 'log'
        else:
            style = 'linear'
        mask = apply_cross_corr_threshold(data_dict[key], debug_fig_ax=ax[n],
                                          style=style,
                                          label=key, percentile=percentile)

        mask_all = np.logical_and(mask_all, mask)
    logger.info(f'{label=}: remove {np.sum(mask_all==False)} datasets out of {len(mask_all)}')

    if plot_debug:
        ax[-1].plot(mask_all, 'o')
        ax[-1].set_title('combined_axis')
        plt.tight_layout()
    
        if not os.path.isdir('debug'):
            os.mkdir('debug')
        plt.savefig(f'debug/debug_outlier_{label}.png', dpi=300)
        plt.close(fig)

    return mask_all


def average_datasets_without_outlier(args):
    """
    Average datasets without outlier.
    """
    mask = outlier_removal(args[0], label=args[1], percentile=5)
    avg_dict = average_datasets(data_dict=args[0], mask=mask)
    return avg_dict


def average_datasets_without_outlier_parallel(args, num_cores=24):
    with Pool(num_cores) as p:
        result = p.map(average_datasets_without_outlier, args)
    return result


def get_temperature(fname, zone_idx='auto'):
    if zone_idx == 'auto':
        zone_idx = (ord(os.path.basename(fname)[0]) - ord('A')) // 3 + 1
        
    assert 1 <= zone_idx <= 3, 'zone_idx must be in [1, 2, 3]'
    
    # key = f'/measurement/sample/QNW_Zone{zone_idx}_Temperature'
    key = f"/entry/sample/qnw{zone_idx}_temperature"
        
    try:
        with h5py.File(fname) as f:
            val = float(f[key][()][0])
    except Exception:
        val = np.nan
        logger.info(f'Failed to get temperature for file {fname}, return nan')
    return val


def read_temperature_from_files(fnames, zone_idx=1):
    return np.array([get_temperature(f, zone_idx) for f in fnames])


def _read_xpcs_hdf(fname):
    """Read tau, g2, g2_err, and q_vals from a single XPCS HDF file using h5py."""
    with h5py.File(fname, 'r') as hf:
        t0 = float(np.asarray(hf['/entry/instrument/detector_1/frame_time'][()]).flat[0])
        delay_list = hf['/xpcs/multitau/delay_list'][()]
        if delay_list.ndim > 1:
            delay_list = delay_list[:, 0]   # some files store the delays as a column
        tau = delay_list * t0
        g2 = hf['/xpcs/multitau/normalized_g2'][()]
        g2_err = hf['/xpcs/multitau/normalized_g2_err'][()]
        q_vals = hf['/xpcs/qmap/dynamic_v_list_dim0'][()]
    return tau, g2, g2_err, q_vals


def process_group(group,
                  prefix,
                  num_sections=10,
                  zone_idx=1,
                  num_cores=24,
                  skip_first_files=0,
                  skip_last_files=0):
    
    
    # read file list in the folder
    flist = glob.glob(os.path.join(prefix, f'{group}*.hdf'))
    flist.sort()
    assert len(flist) > 0, f'no dataset is found in {prefix}'

    logger.info(f'total number of files in {group}  is {len(flist)}')

    logger.info(f'{skip_first_files=}, {skip_last_files=}')

    skip_end = len(flist) - skip_last_files
    flist = flist[skip_first_files: skip_end]
    assert len(flist) > 0, 'no dataset after the skip filter'
    
    # get the temperature
    temperature_list = read_temperature_from_files(flist, zone_idx=zone_idx)

    from pyxpcsviewer import XpcsFile as XF

    if num_sections <= 0:
        keys=('g2', 'g2_err', 'saxs_1d')
        data_dict = read_keys_from_files(flist, keys=keys)
        axf = XF(flist[0])
        t_el = axf.t_el
        ql_dyn = axf.dqlist
        ql_sta = axf.sqlist
        data_dict['temperature'] = temperature_list

        return data_dict, t_el, ql_dyn, ql_sta

    flist_sections = split(flist, num_sections)
    idx_sections = split(np.arange(len(flist)), num_sections)
    print('index\t T-min(C)\t T-max(C)\t T-mean(C)\t points')
    for n in range(num_sections):
        temp = temperature_list[idx_sections[n]]
        print(f'{n=:02d}\t {round(np.min(temp), 4):8}\t {round(np.max(temp), 4):8}\t {round(np.mean(temp), 4):8}\t {len(idx_sections[n])}')

    # read main field for averag
    data_dict_all = read_keys_from_files_parallel(flist_sections, num_cores=num_cores)

    # do outlier removal and average
    args = [(data_dict_all[n], f'{group}_section_{n:02d}') for n in range(num_sections)]
    avg_all = average_datasets_without_outlier_parallel(args, num_cores=num_cores)

    for n, avg_dict in enumerate(avg_all):
        avg_dict['temperature'] = temperature_list[idx_sections[n]]
        avg_dict['temperature_x'] = idx_sections[n]

    # get additonal field for plotting
    axf = XF(flist[0])
    t_el = axf.t_el
    ql_dyn = axf.dqlist
    ql_sta = axf.sqlist
    return avg_all, t_el, ql_dyn, ql_sta


def process_group_by_range(group,
                            prefix,
                            file_ranges=None,
                            zone_idx=1,
                            num_cores=24):
    """Like process_group, but selects files by explicit frame-number ranges instead of
    splitting into equal sections.  Uses h5py directly — no pyxpcsviewer required.

    Parameters
    ----------
    group : str
        Dataset group label used as glob prefix (e.g. 'B0147').
    prefix : str
        Directory containing the cluster-result HDF files.
    file_ranges : list of (int, int)
        Each tuple ``(start, end)`` selects files whose embedded frame number
        satisfies ``start <= frame_number <= end`` (inclusive).
        E.g. ``[(951, 1050), (1051, 1150), (1151, 1250), (1251, 1313)]``.
    zone_idx : int or 'auto'
        Accepted for signature compatibility with ``process_group``; this
        function reads no temperature and ignores it.

    Returns
    -------
    avg_all : list of dict
        One averaged data dict per range with keys 'g2' and 'g2_err'.
    t_el : ndarray
        Delay times in seconds, read from the first file in the first range.
    ql_dyn : ndarray
        Dynamic q values from the first file.
    ql_sta : None
        Placeholder kept for API compatibility with process_group.
    """
    import re

    if file_ranges is None:
        raise ValueError('file_ranges must be a list of (start, end) tuples')

    flist_all = glob.glob(os.path.join(prefix, f'{group}*.hdf'))
    flist_all.sort()
    assert len(flist_all) > 0, f'no dataset found in {prefix} for group {group}'

    logger.info(f'total number of files in {group} is {len(flist_all)}')

    # Extract the frame number from filenames like B0147_..._00950_results.hdf
    _frame_re = re.compile(r'_(\d+)_results')

    def _extract_frame(fname):
        m = _frame_re.search(os.path.basename(fname))
        if m is None:
            return -1
        return int(m.group(1))

    frame_numbers = np.array([_extract_frame(f) for f in flist_all])

    avg_all = []
    t_el = None
    ql_dyn = None

    print('range\t\t files')
    for start, end in file_ranges:
        sel = (frame_numbers >= start) & (frame_numbers <= end)
        section_files = [f for f, m in zip(flist_all, sel) if m]
        assert len(section_files) > 0, f'no files found for range ({start}, {end})'
        print(f'{start}-{end}\t {len(section_files)} files')

        g2_stack, g2_err_stack = [], []
        for f in section_files:
            try:
                tau, g2, g2_err, q_vals = _read_xpcs_hdf(f)
                if t_el is None:
                    t_el = tau
                    ql_dyn = q_vals
                g2_stack.append(g2)
                g2_err_stack.append(g2_err)
            except Exception as exc:
                logger.info(f'failed to read {f}: {exc}')

        g2_arr     = np.array(g2_stack)
        g2_err_arr = np.array(g2_err_stack)
        avg_g2     = np.nanmean(g2_arr, axis=0)
        count      = np.sum(~np.isnan(g2_err_arr), axis=0)
        avg_g2_err = np.sqrt(np.nansum(g2_err_arr ** 2, axis=0)) / np.maximum(count, 1)

        avg_all.append({'g2': avg_g2, 'g2_err': avg_g2_err})

    return avg_all, t_el, ql_dyn, None


# This module is a library only; it has no __main__ entry point.  The reduction
# that this repository actually performs is SAXPCS_8id/average_ranges.py, which
# imports outlier_removal and average_datasets from here.  To use the
# beamline-side helpers interactively, import them and pass the prefix of the
# beamtime you are reducing, e.g.::
#
#     from common.utils import process_group
#     avg, t_el, ql_dyn, ql_sta = process_group('B0147', prefix=<beamtime dir>)
    exit(0)