#!/usr/bin/env python3

import datetime
import logging
import os
import sys
from copy import deepcopy

import numpy as np
from peakdet.io import load_physio

from phys2cvr import io, signal, stats, _version
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
    log_path = os.path.join(outdir, 'logs')
    os.makedirs(log_path, exist_ok=True)
    isotime = datetime.datetime.now().strftime('%Y-%m-%dT%H%M%S')
    f = open(os.path.join(log_path, f'p2c_call_{isotime}.sh'), 'a')
    f.write(f'#!bin/bash \n{call_str}')
    f.close()


def phys2cvr(fname_func, fname_co2='', fname_pidx='', fname_mask='', outdir='',
             freq='', tr='', trial_len='', n_trials='', highcut='', lowcut='',
             apply_filter=False, no_pad=False, do_regression=False,
             lagged_regression=True, maxlag=9, d_lag='', l_degree=0, denoise_matrix=[],
             scale_factor='', lag_map='', regr_dir='', skip_conv=False,
             quiet=False, debug=False):
    """
    Run main workflow of phys2cvr.
    """
    # Add logger and suff
    if outdir:
        outdir = os.path.abspath(outdir)
    else:
        outdir = os.path.join(os.path.split(fname_func)[0], 'phys2cvr')
    outdir = os.path.abspath(outdir)
    petco2log_path = os.path.join(outdir, 'logs')
    os.makedirs(petco2log_path, exist_ok=True)

    # Create logfile name
    basename = 'phys2cvr_'
    extension = 'tsv'
    isotime = datetime.datetime.now().strftime('%Y-%m-%dT%H%M%S')
    logname = os.path.join(petco2log_path, f'{basename}{isotime}.{extension}')

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
                func_avg = signal.filter_signal(func_avg, tr, lowcut, highcut)
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
            _, mask, _ = io.load_nifti_get_mask(fname_mask, is_mask=True)
            if func.shape[:3] != mask.shape:
                raise Exception(f'{fname_mask} and {fname_func} have different sizes!')
            mask = mask * dmask
        else:
            mask = dmask
            LGR.warning(f'No mask specified, using any voxel different from 0 in {fname_func}')

        if apply_filter:
            LGR.info(f'Obtaining filtered average of {fname_func}')
            func_filt = signal.filter_signal(func, tr, lowcut, highcut)
            func_avg = func_filt[mask].mean(axis=0)
        else:
            func_avg = func[mask].mean(axis=0)

    else:
        raise Exception(f'{fname_func} file type is not supported yet, or '
                        'the extension was not specified.')

    if fname_co2 == '':
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
            elif not skip_conv:
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
            phys = load_physio(fname_co2, allow_pickle=True)

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
        outname = os.path.join(outdir, basename_co2)

        # Unless user asks to skip this step, convolve the end tidal signal.
        if skip_conv or not fname_co2:
            petco2hrf = co2
        else:
            petco2hrf = signal.convolve_petco2(co2, pidx, freq, outname)

    # If a regressor directory is not specified, compute the regressors.
    if not regr_dir:
        regr, regr_shifts = stats.get_regr(func_avg, petco2hrf, tr, freq, outname,
                                           maxlag, trial_len, n_trials, no_pad,
                                           '.1D', lagged_regression)
    elif do_regression:
        try:
            regr = np.genfromtxt(f'{outname}_petco2hrf.1D')
        except:
            regr, _ = stats.get_regr(func_avg, petco2hrf, tr, freq, outname, maxlag,
                                     trial_len, n_trials, no_pad, '.1D')

    # Run internal regression if required and possible!
    if func_is_nifti and do_regression:
        LGR.info('Running regression!')

        # Change dimensions in image header before export
        LGR.info('Prepare output image')
        outfuncname = os.path.splitext(os.path.splitext(fname_func)[0])[0]
        outfuncname = os.path.join(outdir, outfuncname)
        newdim = deepcopy(img.header['dim'])
        newdim[0], newdim[4] = 3, 1
        oimg = deepcopy(img)
        oimg.header['dim'] = newdim

        # Compute signal percentage change of functional data
        m = func.mean(axis=-1)[..., np.newaxis]
        func = (func - m) / m

        # Start computing the polynomial regressor (at least average)
        LGR.info(f'Compute Legendre polynomials up to order {l_degree}')
        mat_conf = stats.get_legendre(l_degree, regr.size)

        # Read in eventual denoising factors
        if denoise_matrix:
            for matrix in denoise_matrix:
                LGR.info(f'Read confounding factor from {matrix}')
                conf = np.genfromtxt(matrix)
                mat_conf = np.hstack([mat_conf, conf])

        LGR.info('Compute simple CVR estimation (bulk shift only)')
        beta, tstat, _ = stats.regression(func, dmask, regr, mat_conf)

        LGR.info('Export bulk shift results')
        if not scale_factor:
            LGR.warning('Remember: CVR might not be in %BOLD/mmHg!')
        else:
            beta = beta * scale_factor
        # Scale beta by scale factor while exporting (useful to transform V in mmHg)
        io.export_nifti(beta, oimg, f'{outfuncname}_cvr_simple')
        io.export_nifti(tstat, oimg, f'{outfuncname}_tstat_simple')

        if lagged_regression:
            LGR.info(f'Running lagged CVR estimation with max lag = {maxlag}!'
                     '(might take a while...)')

            nrep = int(maxlag * freq * 2)
            if d_lag:
                step = int(d_lag * freq)
            else:
                step = 1

            if regr_dir:
                outprefix = os.path.join(regr_dir, os.path.split(outname)[1])

            # If user specified a lag map, use that one to regress things
            if lag_map:
                lag, _, _ = io.load_nifti_get_mask(lag_map)
                if func.shape[:3] != lag.shape:
                    raise Exception(f'{lag_map} and {fname_func} have different sizes!')

                # Read d_lag from file (or try to)
                lag_list = np.unique(lag)
                if not d_lag:
                    d_lag = np.unique(lag_list[1:] - lag_list[:-1])
                    if d_lag.size > 1:
                        raise Exception(f'phys2cvr found different delta lags in {lag_map}')
                    else:
                        d_lag = float(d_lag)
                        LGR.warning(f'phys2cvr detected a delta lag of {d_lag} seconds')
                if d_lag:
                    LGR.warning(f'Forcing delta lag to be {d_lag}')

                lag = lag * dmask

                lag_idx = (lag + maxlag) * freq / step

                lag_idx_list = np.unique[lag_idx]

                # Prepare empty matrices
                beta = np.empty_like(lag)
                tstat = np.empty_like(lag)

                for i in lag_idx_list:
                    LGR.info(f'Perform L-GLM number {i+1} of {nrep // step}')
                    try:
                        regr = regr_shifts[:, (i*step)]
                    except NameError:
                        regr = np.genfromtxt(f'{outprefix}_{(i*step):04g}')

                    (beta[lag_idx == i],
                     tstat[lag_idx == i],
                     _) = stats.regression(func[lag_idx == i], [lag_idx == i],
                                           regr, mat_conf)

            else:
                # Prepare empty matrices
                r_square = np.empty(list(func.shape[:3]) + [nrep // step])
                beta_all = np.empty(list(func.shape[:3]) + [nrep // step])
                tstat_all = np.empty(list(func.shape[:3]) + [nrep // step])

                for n, i in enumerate(range(0, nrep, step)):
                    LGR.info(f'Perform L-GLM number {n+1} of {nrep // step}')
                    try:
                        regr = regr_shifts[:, i]
                    except NameError:
                        regr = np.genfromtxt(f'{outprefix}_{i:04g}')

                    (beta_all[:, :, :, n],
                     tstat_all[:, :, :, n],
                     r_square[:, :, :, n]) = stats.regression(func, dmask, regr,
                                                              mat_conf)

                # Find the right lag for CVR estimation
                lag_idx = np.argmax(r_square, axis=-1)
                lag = (lag_idx * step) / freq - (dmask * maxlag)
                # Express lag map relative to median of the mask
                lag_rel = lag - (dmask * np.median(lag[mask]))

                # Run through indexes to pick the right value
                lag_idx_list = np.unique(lag_idx)
                beta = np.empty_like(lag)
                tstat = np.empty_like(lag)
                for i in lag_idx_list:
                    beta[lag_idx == i] = beta_all[:, :, :, i][lag_idx == i]
                    tstat[lag_idx == i] = tstat_all[:, :, :, i][lag_idx == i]

            LGR.info('Export fine shift results')
            if not scale_factor:
                LGR.warning('Remember: CVR might not be in %BOLD/mmHg!')
            else:
                beta = beta * scale_factor

            io.export_nifti(beta, oimg, f'{outfuncname}_cvr')
            io.export_nifti(tstat, oimg, f'{outfuncname}_tstat')
            if not lag_map:
                io.export_nifti(lag, oimg, f'{outfuncname}_lag')
                io.export_nifti(lag_rel, oimg, f'{outfuncname}_lag_mkrel')

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
