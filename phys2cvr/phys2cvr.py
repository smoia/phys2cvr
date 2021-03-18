#!/usr/bin/env python3

import datetime
import logging
import os
import sys

import numpy as np

import matplotlib.pyplot as plt
import nibabel as nib
import scipy.interpolate as spint
import scipy.stats as sct
from peakutils.io import load_physio
from scipy.signal import butter, filtfilt

from phys2cvr import _version
from phys2cvr.cli.run import _get_parser

from . import __version__

SET_DPI = 100
FIGSIZE = (18, 10)
LGR = logging.getLogger(__name__)
LGR.setLevel(logging.INFO)
EXT_1D = ['.txt', '.csv', '.tsv', '.1d', '.par', '.tsv.gz']
EXT_NIFTI = ['.nii', '.nii.gz']


def save_bash_call(outdir):
    """
    Save the bash call into file `p2d_call.sh`.

    Parameters
    ----------
    metric : function
        Metric function to retrieve arguments for
    metric_args : dict
        Dictionary containing all arguments for all functions requested by the
        user
    """
    arg_str = ' '.join(sys.argv[1:])
    call_str = f'phys2cvr {arg_str}'
    outdir = os.path.abspath(outdir)
    log_path = os.path.join(outdir, 'code', 'logs')
    os.makedirs(log_path, exist_ok=True)
    isotime = datetime.datetime.now().strftime('%Y-%m-%dT%H%M%S')
    f = open(os.path.join(log_path, f'p2c_call_{isotime}.sh'), "a")
    f.write(f'#!bin/bash \n{call_str}')
    f.close()


def check_ext(all_ext, fname):
    """
    Check which extension a file has.

    Parameters
    ----------
    all_ext: list
        All possibel extensions to check within
    fname: str
        The filename to check

    Returns
    -------
    has_ext: boolean
        True if the extension is found, false otherwise.
    """
    has_ext = False
    for ext in all_ext:
        if fname.endswith(ext):
            has_ext = True
            break

    return has_ext


def check_nifti_dim(fname, data, dim=4):
    """
    Remove extra dimensions.
    """
    if len(data.shape) < dim:
        raise Exception(f'{fname} does not seem to be a {dim}D file. '
                        f'Plase provide a {dim}D nifti file.')
    if len(data.shape) > dim:
        LGR.warning(f'{fname} has more than {dim} dimensions. Removing D > {dim}.')
        for ax in range(dim, len(data.shape)):
            data = np.delete(data, np.s_[1:], axis=ax)

    return np.squeeze(data)


def load_nifti_get_mask(fname, is_mask=False):
    """
    Load a nifti file and returns its data, its image, and a 3d mask.

    Parameters
    ----------
    fname: str
        The filename to read in
    dim: int
        The number of dimensions expected

    Returns
    -------
    data: np.ndarray
        Data from nifti file.
    mask: np.ndarray

    """
    img = nib.load(fname)
    data = img.get_fdata()

    if is_mask:
        data = check_nifti_dim(fname, data, dim=3)
        mask = (data < 0) + (data > 0)
    else:
        data = check_nifti_dim(fname, data)
        mask = np.squeeze(np.any(data, axis=-1))

    return data, mask, img


def create_hrf(freq=40):
    """
    Create a canonical haemodynamic response function sampled at the given frequency.

    Parameters
    ----------
    freq: float
        Sampling frequency of the haemodynamic response function.

    Returns
    -------
    hrf: np.ndarray
        Haemodynamic response function.

    """
    # Create HRF
    RT = 1/freq
    fMRI_T = 16
    p = [6, 16, 1, 1, 6, 0, 32]

    # Modelled hemodynamic response function - {mixture of Gammas}
    dt = RT / fMRI_T
    u = np.arange(0, p[6]/dt+1, 1) - p[5]/dt
    a1 = p[0] / p[2]
    b1 = 1 / p[3]
    a2 = p[1] / p[3]
    b2 = 1 / p[3]
    hrf = (sct.gamma.pdf(u*dt, a1, scale=b1) - sct.gamma.pdf(u*dt, a2, scale=b2)/p[4])/dt
    time_axis = np.arange(0, int(p[6]/RT+1), 1) * fMRI_T
    hrf = hrf[time_axis]
    min_hrf = 1e-9*min(hrf[hrf > 10*np.finfo(float).eps])

    if min_hrf < 10*np.finfo(float).eps:
        min_hrf = 10*np.finfo(float).eps

    hrf[hrf == 0] = min_hrf
    hrf = hrf/max(hrf)

    plt.figure()
    plt.plot(time_axis, hrf)

    return hrf


def filter_signal(data, tr, lowcut='', highcut=''):
    """
    Create a bandpass filter given a lowcut and a highcut, then filter data.

    Parameters
    ----------
    data: np.ndarray
        Data to filter (along last axis)
    tr: float
        TR of functional files
    lowcut: float
        Lower frequency in the bandpass
    highcut: float
        Higher frequency in the bandpass

    Returns
    -------
    filt_data: np.ndarray
        Input `data`, but filtered.

    """
    if not lowcut:
        lowcut = 0.02
    if not highcut:
        highcut = 0.04
    nyq = (1 / tr) / 2
    low = lowcut / nyq
    high = highcut / nyq
    a, b = butter(9, [low, high], btype='band')
    filt_data = filtfilt(a, b, data, axis=-1)
    return filt_data


def convolve_petco2(co2, pidx, freq, fname_co2, outdir):
    # Extract PETco2
    hrf = create_hrf(freq)
    co2_lenght = len(co2)
    nx = np.linspace(0, co2_lenght, co2_lenght)
    f = spint.interp1d(pidx, co2[pidx], fill_value='extrapolate')
    petco2 = f(nx)

    # Plot PETco2
    plt.figure(figsize=FIGSIZE, dpi=SET_DPI)
    plt.title('CO2 and PetCO2')
    plt.plot(co2, '-', petco2, '-')
    plt.savefig(f'{fname_co2}_petco2.png', dpi=SET_DPI)
    plt.close()

    # Demean and export
    petco2 = petco2 - petco2.mean()
    np.savetxt(f'{fname_co2}_petco2.1D', petco2, fmt='%.18f')

    # Convolve, and then rescale to have same amplitude (?)
    co2_conv = np.convolve(petco2, hrf)
    co2_conv = np.interp(co2_conv, (co2_conv.min(), co2_conv.max()),
                         (petco2.min(), petco2.max()))

    plt.figure(figsize=FIGSIZE, dpi=SET_DPI)
    plt.title('PetCO2 and convolved regressor (PetCO2hrf)')
    plt.plot(co2_conv, '-', petco2, '-')
    plt.savefig(f'{fname_co2}_co2_conv.png', dpi=SET_DPI)
    plt.close()

    np.savetxt(f'{fname_co2}_co2_conv.1D', co2_conv, fmt='%.18f')

    return co2_conv


def export_regressor(regr_x, co_shift, GM_x, GM_name, suffix='_co_regr'):
    f = spint.interp1d(regr_x, co_shift, fill_value='extrapolate')
    co_tr = f(GM_x)
    co_tr = co_tr - co_tr.mean()
    textname = GM_name + suffix + '.1D'
    np.savetxt(textname, co_tr, fmt='%.18f')


def get_regr(func_avg, petco2hrf, tr=1.5, freq=40, BH_len='', nBH='', ext='.1D'):
    sequence_tps = len(func_avg)

    regr_x = np.arange(0, ((sequence_tps-1) * tr + 1/freq), 1/freq)
    func_x = np.linspace(0, (sequence_tps - 1) * tr, sequence_tps)

    regr_len = len(regr_x)
    BH_len_upsampled = int(BH_len*freq)
    nBH = int(nBH)

    f = spint.interp1d(func_x, func_avg, fill_value='extrapolate')
    func_avg_upsampled = f(regr_x)

    # Preparing central breathhold and CO2 trace for Xcorr
    # CO2 trace should have the equivalent of
    # ten tr of bad data at the beginning of the file
    if BH_len:
        last_tp = BH_len_upsampled*(nBH-1)
    else:
        last_tp = -1

    func_avg_cut = func_avg_upsampled[BH_len_upsampled:last_tp]
    petco2hrf_cut = petco2hrf[BH_len_upsampled:]

    func_avg_cut_len = len(func_avg_cut)
    nrep = len(petco2hrf_cut) - func_avg_cut_len
    if BH_len and nrep > BH_len_upsampled:
        nrep = BH_len_upsampled

    func_avg_co_r = np.zeros(nrep)
    for i in range(0, nrep):
        func_avg_co_r[i] = np.corrcoef(func_avg_cut, petco2hrf[0+i:func_avg_cut_len+i].T)[1, 0]

    optshift = int(func_avg_co_r.argmax())
    textname = func_avg_name + '_optshift.1D'
    # #!#
    optshiftout = np.array((optshift/freq,0))
    np.savetxt(textname, optshiftout, fmt='%.4f')
    co_shift = petco2hrf[optshift:optshift+regr_len]

    # preparing for and exporting figures of shift
    time_axis = np.arange(0, nrep/freq, 1/freq)
    # #!# I should change to following line but not tested yet
    # time_axis = np.linspace(0, nrep*freq, nrep)
    if nrep < len(time_axis):
        time_axis = time_axis[:nrep]
    elif nrep > len(time_axis):
        time_axis = np.pad(time_axis, (0, int(nrep - len(time_axis))), 'linear_ramp')

    plt.figure(figsize=FIGSIZE, dpi=SET_DPI)
    plt.plot(time_axis, func_avg_co_r)
    plt.title('optshift')
    plt.savefig(func_avg_name + '_optshift.png', dpi=SET_DPI)

    plt.figure(figsize=FIGSIZE, dpi=SET_DPI)
    plt.plot(sct.zscore(co_shift), '-', sct.zscore(func_avg_upsampled), '-')
    plt.title('GM and shift')
    plt.savefig(func_avg_name + '_co_regr.png', dpi=SET_DPI)

    export_regressor(regr_x, co_shift, func_x, func_avg_name, '_co_regr')

    # Create folder
    func_avg_dir = func_avg_name + '_regr_shift'
    if not os.path.exists(func_avg_dir):
        os.makedirs(func_avg_dir)

    # Set num of fine shifts: 9 seconds is a bit more than physiologically feasible
    nrep = int(9 * freq)

    # Padding regressor for shift, and padding optshift too
    if nrep > optshift:
        left_pad = nrep - optshift
    else:
        left_pad = 0

    if (optshift + nrep + regr_len) > len(petco2hrf):
        right_pad = (optshift + nrep + regr_len) - len(petco2hrf)
    else:
        right_pad = 0

    co_padded = np.pad(petco2hrf, (int(left_pad), int(right_pad)), 'mean')
    optshift_padded = optshift + left_pad

    for i in range(-nrep, nrep):
        co_shift = co_padded[optshift_padded-i:optshift_padded-i+regr_len]
        suffix = '/shift_' + '%04d' % (i + nrep)
        export_regressor(regr_x, co_shift, GM_x, GM_dir, suffix)


def phys2cvr(fname_func, fname_co2='', fname_pidx='', fname_mask='', outdir='',
             freq='', tr='', trial_len='', n_trials='', highcut='', lowcut='',
             apply_filter='',
             do_regression=False, no_phys=False, lagged_regression=True,
             quiet=False, debug=False):
    """
    Run main workflow of phys2cvr.
    """
    # Add logger and suff
    outdir = os.path.abspath(outdir)
    os.makedirs(outdir, exist_ok=True)
    petco2log_path = outdir
    # petco2log_path = os.path.join(outdir, 'code', 'petco2log')
    # os.makedirs(petco2log_path, exist_ok=True)

    # Create logfile name
    basename = 'phys2cvr_'
    extension = 'tsv'
    isotime = datetime.datetime.now().strftime('%Y-%m-%dT%H%M%S')
    logname = os.path.join(petco2log_path, (basename + isotime + '.' + extension))

    # Set logging format
    log_formatter = logging.Formatter(
        '%(asctime)s\t%(name)-12s\t%(levelname)-8s\t%(message)s',
        datefmt='%Y-%m-%dT%H:%M:%S')

    # Set up logging file and open it for writing
    log_handler = logging.FileHandler(logname)
    log_handler.setFormatter(log_formatter)
    sh = logging.StreamHandler()

    if quiet:
        logging.basicConfig(level=logging.WARNING,
                            handlers=[log_handler, sh], format='%(levelname)-10s %(message)s')
    elif debug:
        logging.basicConfig(level=logging.DEBUG,
                            handlers=[log_handler, sh], format='%(levelname)-10s %(message)s')
    else:
        logging.basicConfig(level=logging.INFO,
                            handlers=[log_handler, sh], format='%(levelname)-10s %(message)s')

    version_number = _version.get_versions()['version']
    LGR.info(f'Currently running phys2cvr version {version_number}')
    LGR.info(f'Input file is {fname_func}')

    # Save call.sh
    arg_str = ' '.join(sys.argv[1:])
    call_str = f'phys2bids {arg_str}'
    f = open(os.path.join(petco2log_path, 'call.sh'), "a")
    f.write(f'#!bin/bash \n{call_str}')
    f.close()

    # Check func type and read it
    func_is_1d = check_ext(EXT_1D, fname_func)
    func_is_nifti = check_ext(EXT_NIFTI, fname_func)

    if func_is_1d:
        if tr:
            func_avg = np.genfromtxt(fname_func)
            if apply_filter:
                LGR.info('Applying butterworth filter to {fname_func}')
                func_avg = filter_signal(func_avg, tr, lowcut, highcut)
        else:
            raise Exception('Provided functional signal, but no TR specified! '
                            'Rerun specifying the TR')
    elif func_is_nifti:
        func, mask, img = load_nifti_get_mask(fname_func)
        if len(func.shape) < 4:
            raise Exception(f'Provided functional file {fname_func} is not a 4D file!')
        # Read TR or declare its overwriting
        if tr:
            LGR.warning(f'Forcing functional TR to be {tr} seconds')
        else:
            tr = img.header['pixdim'][4]

        # Read mask if provided
        if fname_mask:
            _, mask, img = load_nifti_get_mask(fname_mask)
            if func.shape[:3] != mask.shape:
                raise Exception(f'{fname_mask} and {fname_func} have different sizes!')
        else:
            LGR.warning(f'No mask specified, using any voxel different from 0 in {fname_func}')

        if apply_filter:
            LGR.info('Obtaining filtered average of {fname_func}')
            func_filt = filter_signal(func, tr, lowcut, highcut)
            func_avg = func_filt[mask].mean(axis=-1)
        else:
            func_avg = func[mask].mean(axis=-1)

    else:
        raise Exception(f'{fname_func} file type is not supported yet.')

    if fname_co2 == '' and not no_phys:
        raise Exception('The pipeline for no physiological files was not selected, and '
                        'no file for physiological regressor was found. Please rerun with '
                        'one option or the other.')
    elif no_phys:
        LGR.info(f'Computing "CVR" maps using {fname_func} only')
        if func_is_1d:
            LGR.warning('Using an average signal only, solution might be unoptimal.')

            if not apply_filter:
                LGR.warning('No filter applied to the input average! You know '
                            'what you are doing, right?')

        petco2hrf = func_avg

    else:
        co2_is_phys = check_ext('.phys', fname_co2)
        co2_is_1d = check_ext(EXT_1D, fname_co2)

        if co2_is_1d:
            if fname_pidx:
                pidx = np.genfromtxt(fname_pidx)
            else:
                raise Exception(f'{fname_co2} file is a text file, but no '
                                'file containing its peaks was provided. '
                                ' Please provide peak file!')

            if not freq:
                raise Exception(f'{fname_co2} file is a text file, but no '
                                'frequency was specified. Please provide peak '
                                ' file!')

            co2 = np.genfromtxt(fname_co2)
        elif co2_is_phys:
            # Read a phys file!
            phys = load_physio

            if freq:
                LGR.warning(f'Forcing CO2 frequency to be {freq} Hz')
        else:
            raise Exception(f'{fname_co2} file type is not supported yet.')

        petco2hrf = convolve_petco2(co2, pidx, freq, fname_co2, outdir)

    os.makedirs('regr', exist_ok=True)

    get_regr(func_avg, petco2hrf, tr, freq, trial_len, n_trials)

    # Add internal regression if required!
    if func_is_nifti and do_regression:
        print('Running regression!')
        pass


def _main(argv=None):
    options = _get_parser().parse_args(argv)

    save_bash_call(options.outdir)

    phys2cvr(**vars(options))


if __name__ == '__main__':
    _main(sys.argv[1:])


if __name__ == '__main__':
    _main()