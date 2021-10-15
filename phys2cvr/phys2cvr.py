#!/usr/bin/env python3
"""
Main file for phys2cvr.

Attributes
----------
LGR :
    Logger
"""

import datetime
import logging
import os
import sys
from copy import deepcopy

import numpy as np
from peakdet.io import load_physio

from phys2cvr import io, signal, stats, _version
from phys2cvr.cli.run import _get_parser, _check_opt_conf
from phys2cvr.io import EXT_NIFTI, EXT_1D


LGR = logging.getLogger(__name__)
LGR.setLevel(logging.INFO)


def save_bash_call(fname, outdir):
    """
    Save the bash call into file `p2d_call.sh`.
    
    Parameters
    ----------
    outdir : str or path, optional
        output directory
    """
    arg_str = ' '.join(sys.argv[1:])
    call_str = f'phys2cvr {arg_str}'
    outdir = os.path.abspath(outdir)
    log_path = os.path.join(outdir, 'logs')
    os.makedirs(log_path, exist_ok=True)
    isotime = datetime.datetime.now().strftime('%Y-%m-%dT%H%M%S')
    fname = os.path.basename(fname).split('.')[0]
    f = open(os.path.join(log_path, f'p2c_call_{fname}_{isotime}.sh'), 'a')
    f.write(f'#!bin/bash \n{call_str}')
    f.close()


def phys2cvr(fname_func, fname_co2=None, fname_pidx=None, fname_mask=None,
             outdir=None, freq=None, tr=None, trial_len=None, n_trials=None,
             highcut=None, lowcut=None, apply_filter=False,
             run_regression=False, lagged_regression=True, lag_max=None,
             lag_step=None, l_degree=0, denoise_matrix=[], scale_factor=None,
             lag_map=None, regr_dir=None, run_conv=True, quiet=False, debug=False):
    """
    Run main workflow of phys2cvr.
    
    Parameters
    ----------
    fname_func : str or path
        Filename of the functional input (nifti or txt)
    fname_co2 : str or path, optional
        Filename of the CO2 (physiological regressor) timeseries.
        Can be either peakdet's output or a txt file.
        If not declared, phys2cvr will consider the average temporal value 
        from the input.
        Default: empty
    fname_pidx : str or path, optional
        Filename of the CO2 (physiological regressor) timeseries' PEAKS.
        Required if CO2 file is a txt AND the convolution step is not skipped.
        If not declared AND the convolution step is not skipped, raises an exception.
        Default: empty
    fname_mask : str or path, optional
        Filename of the mask for nifti volume.
        If declared, phys2cvr will run only on these voxels.
        If not, phys2cvr will estimate a mask from the functional input.
        Ignored if input is a txt file.
        Default: empty
    outdir : str or path, optional
        Output directory
        Default: the directory where `fname_func` is.
    freq : str, int, or float, optional
        Sample frequency of the CO2 regressor. Required if CO2 input is TXT file.
        If declared with peakdet file, it will overwrite the file frequency. 
    tr : str, int, or float, optional
        TR of the timeseries. Required if input is TXT file.
        If declared with nifti file, it will overwrite the file TR. 
    trial_len : str or int, optional
        Length of each single trial for tasks that have more than one
        (E.g. BreathHold, CO2 challenges, ...)
        Used to improve cross correlation estimation.
        Default: None
    n_trials : str or int, optional
        Number of trials in the task.
        Default: None
    highcut : str, int, or float, optional
        High frequency limit for filter.
        Required if applying a filter.
        Default: empty
    lowcut : str, int, or float, optional
        Low frequency limit for filter.
        Required if applying a filter.
        Default: empty
    apply_filter : bool, optional
        Apply a Butterworth filter to the functional input.
        Default: False
    run_regression : bool, optional
        Also run the regression step within phys2cvr.
        By default, phys2cvr will only estimate the regressor(s) of interest.
        Default: False
    lagged_regression : bool, optional
        Estimates regressors to run a lagged regression approach.
        If `run_regression` is True, also run the lagged regression.
        Can be turned off.
        Default: True
    lag_max : int or float, optional
        Limits (both positive and negative) of the temporal area to explore,
        expressed in seconds (i.e. ±9 seconds).
        Default: None 
    lag_step : int or float, optional
        Step of the lag to take into account in seconds.
        Default: None
    l_degree : int, optional
        Only used if performing the regression step.
        Highest order of the Legendre polynomial to add to the denoising matrix.
        phys2cvr will add all polynomials up to the specified order 
        (e.g. if user specifies 3, orders 0, 1, 2, and 3 will be added). 
        Default is 0, which will model only the mean of the timeseries.
    denoise_matrix : list of np.array(s), optional
        Add one or multiple denoising matrices to the regression model.
        Ignored if not performing the regression step.
    scale_factor : str, int, or float, optional
        A scale factor to apply to the CVR map before exporting it.
        Useful when using inputs recorded/stored in Volts that have a meaningful 
        unit of measurement otherwise, e.g. (CO2 traces in mmHg).
        V->mmHg: CO2[mmHg] = (Patm-Pvap)[mmHg] * 10*CO2[V]/100[V]
        Default: None
    lag_map : str or path, optional
        Filename of a lag map to get lags from.
        Ignored if not running a lagged-GLM regression.
        Default: None
    regr_dir : str, optional
        Directory containing pre-generated lagged regressors, useful
        to (re-)run a GLM analysis.
        Default: None
    run_conv : bool, optional
        Run the convolution of the physiological trace.
        Can be turned off
        Default: True
    quiet : bool, optional
        Return to screen only warnings and errors.
        Default: False
    debug : bool, optional
        Return to screen more output.
        Default: False
    
    Raises
    ------
    Exception
        - If functional timeseries is lacking TR and the latter was not specified.
        - If functional nifti file is not at least 4D.
        - If mask was specified but it has different dimensions than the 
            functional nifti file.
        - If a file type is not supported yet.
        - If physiological file is a txt file and no peaks were provided.
        - If physiological file is lacking frequency and the latter was not specified.
        - If a lag map was specified but it has different dimensions than the 
            functional nifti file.
        - If a lag map was specified, lag_step was not, and the lag map seems 
            to have different lag_steps inside.
    """
    # If lagged regression is selected, make sure run_regression is true.
    if lagged_regression:
        run_regression = True
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

    # Check that all input values have right type
    tr = io.if_declared_force_type(tr, 'float', 'tr')
    freq = io.if_declared_force_type(freq, 'float', 'freq')
    trial_len = io.if_declared_force_type(trial_len, 'int', 'trial_len')
    n_trials = io.if_declared_force_type(n_trials, 'int', 'n_trials')
    highcut = io.if_declared_force_type(highcut, 'float', 'highcut')
    lowcut = io.if_declared_force_type(lowcut, 'float', 'lowcut')
    lag_max = io.if_declared_force_type(lag_max, 'float', 'lag_max')
    lag_step = io.if_declared_force_type(lag_step, 'float', 'lag_step')
    l_degree = io.if_declared_force_type(l_degree, 'int', 'l_degree')
    scale_factor = io.if_declared_force_type(scale_factor, 'float', 'scale_factor')

    if func_is_1d:
        if tr:
            func_avg = np.genfromtxt(fname_func)
            if apply_filter:
                LGR.info('Applying butterworth filter to {fname_func}')
                func_avg = signal.filter_signal(func_avg, tr, lowcut, highcut)
        else:
            raise NameError('Provided functional signal, but no TR specified! '
                            'Rerun specifying the TR')
    elif func_is_nifti:
        func, dmask, img = io.load_nifti_get_mask(fname_func)
        if len(func.shape) < 4:
            raise ValueError(f'Provided functional file {fname_func} is not a 4D file!')
        # Read TR or declare its overwriting
        if tr:
            LGR.warning(f'Forcing TR to be {tr} seconds')
        else:
            tr = img.header['pixdim'][4]

        # Read mask if provided
        if fname_mask:
            _, mask, _ = io.load_nifti_get_mask(fname_mask, is_mask=True)
            if func.shape[:3] != mask.shape:
                raise ValueError(f'{fname_mask} and {fname_func} have different sizes!')
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
        raise NotImplementedError(f'{fname_func} file type is not supported yet, or '
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
            elif run_conv:
                raise NameError(f'{fname_co2} file is a text file, but no '
                                'file containing its peaks was provided. '
                                ' Please provide peak file!')

            if not freq:
                raise NameError(f'{fname_co2} file is a text file, but no '
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
            raise NotImplementedError(f'{fname_co2} file type is not supported yet, or '
                                      'the extension was not specified.')

        # Set output file & path - calling splitext twice cause .gz
        basename_co2 = os.path.splitext(os.path.splitext(os.path.basename(fname_co2))[0])[0]
        outname = os.path.join(outdir, basename_co2)

        # Unless user asks to skip this step, convolve the end tidal signal.
        if not run_conv or not fname_co2:
            petco2hrf = co2
        else:
            petco2hrf = signal.convolve_petco2(co2, pidx, freq, outname)

    # If a regressor directory is not specified, compute the regressors.
    if not regr_dir:
        regr, regr_shifts = stats.get_regr(func_avg, petco2hrf, tr, freq, outname,
                                           lag_max, trial_len, n_trials,
                                           '.1D', lagged_regression)
    elif run_regression:
        try:
            regr = np.genfromtxt(f'{outname}_petco2hrf.1D')
        except IOError:
            LGR.warning(f'Regressor {outname}_petco2hrf.1D not found. '
                        'Estimating it.')
            regr, regr_shifts = stats.get_regr(func_avg, petco2hrf, tr, freq,
                                               outname, lag_max, trial_len,
                                               n_trials, '.1D')

    # Run internal regression if required and possible!
    if func_is_nifti and run_regression:
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
            beta = beta * float(scale_factor)
        # Scale beta by scale factor while exporting (useful to transform V in mmHg)
        io.export_nifti(beta, oimg, f'{outfuncname}_cvr_simple')
        io.export_nifti(tstat, oimg, f'{outfuncname}_tstat_simple')

        if lagged_regression and regr_shifts is not None and ((lag_max and lag_step) or lag_map):
            if lag_max:
                LGR.info(f'Running lagged CVR estimation with max lag = {lag_max}! '
                         '(might take a while...)')
            elif lag_map is not None:
                LGR.info(f'Running lagged CVR estimation with lag map {lag_map}! '
                         '(might take a while...)')

            nrep = int(lag_max * freq * 2)
            if lag_step:
                step = int(lag_step * freq)
            else:
                step = 1

            if regr_dir:
                outprefix = os.path.join(regr_dir, os.path.split(outname)[1])

            # If user specified a lag map, use that one to regress things
            if lag_map:
                lag, _, _ = io.load_nifti_get_mask(lag_map)
                if func.shape[:3] != lag.shape:
                    raise ValueError(f'{lag_map} and {fname_func} have different sizes!')

                # Read lag_step and lag_max from file (or try to)
                lag_list = np.unique(lag)
                if not lag_step:
                    lag_step = np.unique(lag_list[1:] - lag_list[:-1])
                    if lag_step.size > 1:
                        raise ValueError(f'phys2cvr found different delta lags in {lag_map}')
                    else:
                        LGR.warning(f'phys2cvr detected a delta lag of {lag_step} seconds')
                else:
                    LGR.warning(f'Forcing delta lag to be {lag_step}')

                if not lag_max:
                    lag_max = abs(np.min(lag_list))
                    LGR.warning(f'phys2cvr detected a max lag of {lag_max} seconds')
                else:
                    LGR.warning(f'Forcing max lag to be {lag_max}')

                lag = lag * dmask

                lag_idx = (lag + lag_max) * freq / step

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
                lag = (lag_idx * step) / freq - (dmask * lag_max)
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

    elif run_regression:
        LGR.warning('The input file is not a nifti volume. At the moment, '
                    'regression is not supported for other formats.')

    LGR.info('phys2cvr finished! Enjoy your outputs!')


def _main(argv=None):
    options = _get_parser().parse_args(argv)

    options = _check_opt_conf(options)

    save_bash_call(options.fname_func, options.outdir)

    phys2cvr(**vars(options))


if __name__ == '__main__':
    _main(sys.argv[1:])


"""
Copyright 2021, Stefano Moia.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""
