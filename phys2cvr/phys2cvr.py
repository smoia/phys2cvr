#!/usr/bin/env python3

import datetime
import logging
import os
import sys
from copy import deepcopy

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as spint
import scipy.stats as sct
from peakdet.io import load_physio
from scipy.signal import butter, filtfilt

from phys2cvr import io, _version
from phys2cvr.cli.run import _get_parser


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


def convolve_petco2(co2, pidx, freq, outname):
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
    plt.savefig(f'{outname}_petco2.png', dpi=SET_DPI)
    plt.close()

    # Demean and export
    petco2 = petco2 - petco2.mean()
    np.savetxt(f'{outname}_petco2.1D', petco2, fmt='%.18f')

    # Convolve, and then rescale to have same amplitude (?)
    co2_conv = np.convolve(petco2, hrf)
    co2_conv = np.interp(co2_conv, (co2_conv.min(), co2_conv.max()),
                         (petco2.min(), petco2.max()))

    plt.figure(figsize=FIGSIZE, dpi=SET_DPI)
    plt.title('PetCO2 and convolved regressor (PetCO2hrf)')
    plt.plot(co2_conv, '-', petco2, '-')
    plt.savefig(f'{outname}_petco2hrf.png', dpi=SET_DPI)
    plt.close()

    np.savetxt(f'{outname}_petco2hrf.1D', co2_conv, fmt='%.18f')

    return co2_conv


def x_corr(func, co2, lastrep, firstrep=0, offset=0):
    if len(func) + offset > len(co2):
        raise ValueError(f'The specified offset of {offset} is too high to '
                         f'compare func of length {len(func)} with co2 of '
                         f'length {len(co2)}')
    if firstrep + offset < 0:
        firstrep = -offset
    if lastrep + offset + len(func) > len(co2):
        lastrep = len(co2) - offset - len(func)

    xcorr = np.empty(lastrep+firstrep)
    for i in range(firstrep, lastrep):
        xcorr[i] = np.corrcoef(func, co2[0+i+offset:len(func)+i+offset].T)[1, 0]

    return xcorr.max(), (xcorr.argmax() + firstrep + offset), xcorr


def get_regr(func_avg, petco2hrf, tr, freq, outname, maxlag=9, trial_len='',
             n_trials='', no_pad=False, ext='.1D', lagged_regression=True):
    # Setting up some variables
    first_tp = 0
    last_tp = -1

    if trial_len and n_trials:
        # If both are specified, disregard two extreme _trial from matching.
        LGR.info(f'Specified {n_trials} trials lasting {trial_len} seconds')
        if n_trials > 2:
            LGR.info('Ignoring first trial to improve first bulk shift estimation')
            first_tp = int(trial_len*freq)
        else:
            LGR.info('Using all trials for bulk shift estimation')
        if n_trials > 3:
            LGR.info('Ignoring last trial to improve first bulk shift estimation')
            last_tp = first_tp*(n_trials-1)

    elif trial_len and not n_trials:
        LGR.warning('The length of trial was specified, but the number of '
                    'trials was not. Using all trials for bulk shift estimation')
    elif not trial_len and n_trials:
        LGR.warning('The number of trials was specified, but the length of '
                    'trial was not. Using all trials for bulk shift estimation')
    else:
        LGR.info('Using all trials for bulk shift estimation.')

    # Upsample functional signal
    func_len = len(func_avg)
    regr_x = np.arange(0, ((func_len-1) * tr + 1/freq), 1/freq)
    func_x = np.linspace(0, (func_len - 1) * tr, func_len)
    f = spint.interp1d(func_x, func_avg, fill_value='extrapolate')
    func_upsampled = f(regr_x)
    len_upd = len(func_upsampled)

    # Preparing breathhold and CO2 trace for Xcorr
    func_cut = func_upsampled[first_tp:last_tp]
    petco2hrf_cut = petco2hrf[first_tp:]

    nrep = len(petco2hrf_cut) - len(func_cut)

    _, optshift, xcorr = x_corr(func_cut, petco2hrf, nrep)
    LGR.info(f'First cross correlation estimated bulk shift at {optshift/freq} seconds')

    if trial_len and n_trials and n_trials > 2:
        LGR.info('Running second bulk shift estimation')
        if len(func_upsampled) + nrep > len(petco2hrf):
            pad = len(func_upsampled) + nrep - len(petco2hrf)
            petco2hrf_padded = np.pad(petco2hrf, pad, mode='mean')
        else:
            petco2hrf_padded = petco2hrf

        _, optshift, xcorr = x_corr(func_upsampled, petco2hrf_padded, nrep, -nrep, optshift)
        LGR.info(f'Second cross correlation estimated bulk shift at {optshift/freq} seconds')

    # Export estimated optimal shift in seconds
    with open(f'{outname}_optshift.1D', 'w') as f:
        print(f'{(optshift/freq):.4f}', file=f)

    petco2hrf_shift = petco2hrf[optshift:optshift+len_upd]

    # preparing for and exporting figures of shift
    time_axis = np.arange(0, nrep/freq, 1/freq)

    if nrep < len(time_axis):
        time_axis = time_axis[:nrep]
    elif nrep > len(time_axis):
        time_axis = np.pad(time_axis, (0, int(nrep - len(time_axis))), 'linear_ramp')

    plt.figure(figsize=FIGSIZE, dpi=SET_DPI)
    plt.plot(time_axis, xcorr)
    plt.title('optshift')
    plt.savefig(f'{outname}_optshift.png', dpi=SET_DPI)

    plt.figure(figsize=FIGSIZE, dpi=SET_DPI)
    plt.plot(sct.zscore(petco2hrf_shift), '-', sct.zscore(func_upsampled), '-')
    plt.title('GM and shift')
    plt.savefig(f'{outname}_petco2hrf.png', dpi=SET_DPI)

    petco2hrf_demean = io.export_regressor(regr_x, petco2hrf_shift, func_x, outname, 'petco2hrf', ext)

    if lagged_regression:
        outprefix = os.path.join(os.path.split(outname)[0], 'regr', os.path.split(outname)[1])
        os.makedirs(os.path.join(os.path.split(outname)[0], 'regr'), exist_ok=True)

        # Set num of fine shifts: 9 seconds is a bit more than physiologically feasible
        nrep = int(maxlag * freq)

        # Padding regressor for shift, and padding optshift too
        if (optshift - nrep) < 0:
            lpad = nrep - optshift
        else:
            lpad = 0

        if (optshift + nrep + len_upd) > len(petco2hrf):
            rpad = (optshift + nrep + len_upd) - len(petco2hrf)
        else:
            rpad = 0

        petco2hrf_padded = np.pad(petco2hrf, (int(lpad), int(rpad)), 'mean')

        for i in range(-nrep, nrep):
            petco2hrf_shift = petco2hrf_padded[optshift+lpad-i:optshift+lpad-i+len_upd]
            io.export_regressor(regr_x, petco2hrf_shift, func_x, outprefix, f'_{(i + nrep):04g}', ext)

    return petco2hrf_demean


def phys2cvr(fname_func, fname_co2='', fname_pidx='', fname_mask='', outdir='',
             freq='', tr='', trial_len='', n_trials='', highcut='', lowcut='',
             apply_filter='', no_pad='', no_phys=False, do_regression=False,
             lagged_regression=True, maxlag=9, denoise_matrix='',
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

    # Check func type and read it
    func_is_1d = io.check_ext(EXT_1D, fname_func)
    func_is_nifti = io.check_ext(EXT_NIFTI, fname_func)

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
        func, dmask, img = io.load_nifti_get_mask(fname_func)
        if len(func.shape) < 4:
            raise Exception(f'Provided functional file {fname_func} is not a 4D file!')
        # Read TR or declare its overwriting
        if tr:
            LGR.warning(f'Forcing TR to be {tr} seconds')
        else:
            tr = img.header['pixdim'][4]

        # Read mask if provided
        if fname_mask:
            _, mask, _ = io.load_nifti_get_mask(fname_mask)
            if func.shape[:3] != mask.shape:
                raise Exception(f'{fname_mask} and {fname_func} have different sizes!')
            mask = mask * dmask
        else:
            mask = dmask
            LGR.warning(f'No mask specified, using any voxel different from 0 in {fname_func}')

        if apply_filter:
            LGR.info('Obtaining filtered average of {fname_func}')
            func_filt = filter_signal(func, tr, lowcut, highcut)
            func_avg = func_filt[mask].mean(axis=-1)
        else:
            func_avg = func[mask].mean(axis=-1)

    else:
        raise Exception(f'{fname_func} file type is not supported yet, or '
                        'the extension was not specified.')

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
        co2_is_phys = io.check_ext('.phys', fname_co2)
        co2_is_1d = io.check_ext(EXT_1D, fname_co2)

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
            phys = load_physio(fname_co2)

            co2 = phys.data
            pidx = phys.peaks
            if freq:
                LGR.warning(f'Forcing CO2 frequency to be {freq} Hz')
            else:
                freq = phys.fs
        else:
            raise Exception(f'{fname_co2} file type is not supported yet, or '
                            'the extension was not specified.')

        # Set output file & path - calling splitext twice cause .gz
        basename_co2 = os.path.splitext(os.path.splitext(os.path.basename(fname_co2))[0])[0]
        outname = os.join(outdir, basename_co2)
        petco2hrf = convolve_petco2(co2, pidx, freq, outname)

    regr = get_regr(func_avg, petco2hrf, tr, freq, outname, maxlag, trial_len,
                    n_trials, no_pad, '.1D', lagged_regression)

    # Add internal regression if required!
    if func_is_nifti and do_regression:
        print('Running regression!')
        pass

        outfuncname = os.path.splitext(os.path.splitext(fname_func)[0])[0]
        newdim = deepcopy(img.header['dim'])
        newdim[0], newdim[4] = 3, 1
        oimg = deepcopy(img)
        oimg.header['dim'] = newdim

        # Read in eventual denoising factors
        if denoise_matrix:
            mat_conf = np.genfromtxt(denoise_matrix)
            if not np.any(np.all(mat_conf == 1, axis=1)):
                mat_conf = np.stack(mat_conf, np.ones_like(regr))
        else:
            mat_conf = np.ones_like(regr)

        mat = np.stack([regr, mat_conf])

        # #!# func has to be SPC!

        model = OLSModel(mat).fit(func)

        beta = model.predicted()

        # #!# beta is not cvr!

        io.export_nifti(beta, img, f'{outfuncname}_cvr_simple')

        if lagged_regression:
            nrep = int(maxlag * freq * 2)
            outprefix = os.path.join(os.path.split(outname)[0], 'regr', os.path.split(outname)[1])
            r_square = np.empty((func.shape[0], func.shape[1], func.shape[2], nrep))
            beta_all = np.empty((func.shape[0], func.shape[1], func.shape[2], nrep))

            for i in range(nrep):
                regr = np.genfromtxt(f'{outprefix}_{(i + nrep):04g}')

                mat = np.stack([regr, mat_conf])

                model = OLSModel(mat).fit(func)

                beta_all[:, :, :, i] = model.predicted()
                r_square[:, :, :, i] = model.r_square()

            lag_idx = np.argmax(r_square, axis=-1)
            lag = lag_idx / freq
            beta = beta_all[:, :, :, lag_idx]

            io.export_nifti(beta, oimg, f'{outfuncname}_cvr')
            io.export_nifti(lag, oimg, f'{outfuncname}_lag')

    elif do_regression:
        LGR.warning('The input file is not a nifti volume. At the moment, '
                    'regression is not supported for other formats.')

    LGR.info('phys2cvr finished! Enjoy your outputs!')


def _main(argv=None):
    options = _get_parser().parse_args(argv)

    save_bash_call(options.outdir)

    phys2cvr(**vars(options))


if __name__ == '__main__':
    _main(sys.argv[1:])
