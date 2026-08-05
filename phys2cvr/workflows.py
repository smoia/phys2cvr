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

import nibabel as nib
import numpy as np

from phys2cvr import _version, io, signal, stats, utils
from phys2cvr.cli.run import _check_opt_conf, _get_parser
from phys2cvr.regressors import (
    compute_petco2hrf,
    create_legendre,
    create_physio_regressor,
    select_lag_avoid_boundary
    
)

LGR = logging.getLogger(__name__)
LGR.setLevel(logging.INFO)


def phys2cvr(
    fname_func,
    fname_co2=None,
    fname_pidx=None,
    fname_roi=None,
    fname_mask=None,
    outdir=None,
    freq=None,
    tr=None,
    trial_len=None,
    n_trials=None,
    abs_xcorr=False,
    skip_xcorr=False,
    highcut=0.04,
    lowcut=0.02,
    butter_order=9,
    apply_filter=False,
    run_regression=False,
    lagged_regression=True,
    r2model='full',
    lag_max=None,
    lag_min=None,
    lag_step=None,
    starting_lag_max=None,
    lag_increment=None,
    legacy=False,
    l_degree=0,
    denoise_matrix_file=None,
    orthogonalised_matrix_file=None,
    extra_matrix_file=None,
    scale_factor=None,
    lag_map=None,
    regr_dir=None,
    comp_endtidal=True,
    response_function='hrf',
    quiet=False,
    debug=False,
):
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
    fname_roi : str or path, optional
        Filename of the roi in a nifti volume.
        If declared, phys2cvr will use these voxels .
        If not, phys2cvr will use a mask, either the declared one or estimated
        from the functional input.
        Ignored if input is a txt file.
        Default: empty
    fname_mask : str or path, optional
        Filename of the mask in a nifti volume.
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
    abs_xcorr : bool, optional
        If True, the cross correlation will consider max(abs(xcorr)).
        If False, the cross correlation will consider max(xcorr).
        Default: False
    skip_xcorr : bool, optional
        If True, skip the cross correlation step.
        Default: False
    highcut : str, int, or float, optional
        High frequency limit for filter.
        Required if applying a filter.
        Default: 0.02
    lowcut : str, int, or float, optional
        Low frequency limit for filter.
        Required if applying a filter.
        Default: 0.04
    butter_order : int, optional
        Butterworth filter order.
        Default: 9
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
    r2model : {`full`, `partial`, `intercept`, `adj_full`, `adj_partial`, `adj_intercept`), optional
        R^2 model the regression should return (and hence used for lag selection).
        Potentially invariant if no orthogonalisation is introduced,
        will change results with orthogonalisations.
        See `stats.ols` help for more details.
        Default: `full`
    lag_max : int, float, or None, optional
        Upper limit of the temporal area to explore, expressed in seconds.
        Caution: this is not a pythonic range, but a real range, i.e. the upper limit is included.
        Default: None
    lag_min : int, float, or None, optional
        Lower limit of the temporal area to explore, expressed in seconds.
        If set to None, and lag_max is not None and is positive, lag_min defaults to -lag_max (symmetric range).
        Default: None
    lag_step : int, float, or None, optional
        Step of the lag to take into account in seconds.
        Default: None
    starting_lag_max : int, float, or None, optional
        Initial upper limit of the temporal area to explore, expressed in seconds. 
        If this value is specified, it is used as the initial upper limit of the search range for 
        identifying the optimal lag. If the optimal lag of a voxel
        is within 1 lag step of a boundary (either lag_min or starting-lag-max), the maximum lag for 
        that voxel will be iteratively increased by the amount specified with --lag-increment until the optimal lag
        is not at within 1 lag step of a boundary or --lag-max is reached, whichever comes first.
        Default: None
    lag_increment : int, float, or None, optional
        Step size (in seconds) to increase the maximum lag when the optimal lag is within 1 lag step
        of the boundary.  This value needs to be specified if starting_lag_max is not None.'
        Default: None
    legacy : bool, optional
        If True, use pythonic ranges when creating the regressors, i.e. exclude
        the upper range (e.g. [-9, +9) ).
        Default: False
    l_degree : int, optional
        Only used if performing the regression step.
        Highest order of the Legendre polynomial to add to the denoising matrix.
        phys2cvr will add all polynomials up to the specified order
        (e.g. if user specifies 3, orders 0, 1, 2, and 3 will be added).
        Default is 0, which will model only the mean of the timeseries.
    denoise_matrix_file : None, list of str(s) or path(s), optional
        Add one or multiple denoising matrices to the regression model.
        Ignored if not performing the regression step.
        Default is None.
    orthogonalised_matrix_file : None, list of str(s) or path(s), optional
        Add one or multiple denoising matrices to the regression model,
        AFTER orthogonalising them w.r.t. the task, the denoise matrix,
        and the extra matrix.
        Ignored if not performing the regression step.
        Default is None.
    extra_matrix_file : None, list of str(s) or path(s), optional
        Add one or multiple extra matrices to use in the orthogonalisation step.
        These matrices will not be added to the final regression model.
        Ignored if not performing the regression step.
        Default is None.
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
    comp_endtidal : bool, optional
        Compute the end-tidal interpolation of the CO2 signal.
    response_function : {`hrf`, `rrf`, `crf`, `None`}, None, str, path, or 1D array-like, optional
        Name of the response function to be used in the convolution of the regressor of
        interest. Default is `hrf`.
        If None, skips the convolution step.
        If file, loads it.
        For `rrf` and `crf`, `phys2denoise` must be installed (see extra installs).
        See `signal.compute_petco2hrf` for details.
    quiet : bool, optional
        Return to screen only warnings and errors.
        Default: False
    debug : bool, optional
        Return to screen more output.
        Default: False

    Raises
    ------
    ValueError
        - If the order of Legendre Polynomials is < 0.
        - If a wrong R^2 model was specified.
        - If functional nifti file is not at least 4D.
        - If mask was specified but it has different dimensions than the
            functional nifti file.
        - If ROI was specified but it has different dimensions than the
            functional nifti file.
        - If a lag map was specified but it has different dimensions than the
            functional nifti file.
        - If a lag map was specified, lag_step was not, and the lag map seems
            to have different lag_steps inside.
        - If the maximum lag is 0 or negative and no minimum lag is provided
        - If a minimum lag is provided and no maximum lag is provided
        - If the minimum lag is greater than or equal to the maximum lag
        - If a starting_lag_max is provided and no lag_increment is specified
        - If a starting_lag_max is provided and it is greater than or equal to lag_max
        - If a lag map is provided AND starting_lag_max is specified. (The lag map option does
        not work with the iterative lag search)
    NotImplementedError
        - If a file type is not supported yet.
    NameError
        - If functional timeseries is lacking TR and the latter was not specified.
        - If physiological file is a txt file and no peaks were provided.
        - If physiological file is lacking frequency and the latter was not specified.
    """
    # If lagged regression is selected, make sure run_regression is true.
    if extra_matrix_file is None:
        extra_matrix_file = []
    if orthogonalised_matrix_file is None:
        orthogonalised_matrix_file = []
    if denoise_matrix_file is None:
        denoise_matrix_file = []
    if lagged_regression:
        run_regression = True
    # Add logger and suff
    if outdir is None:
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
        datefmt='%Y-%m-%dT%H:%M:%S',
    )

    # Set up logging file and open it for writing
    log_handler = logging.FileHandler(logname)
    log_handler.setFormatter(log_formatter)
    sh = logging.StreamHandler()

    if quiet:
        logging.basicConfig(
            level=logging.WARNING,
            handlers=[log_handler, sh],
            format='%(levelname)-10s %(message)s',
        )
    elif debug:
        logging.basicConfig(
            level=logging.DEBUG,
            handlers=[log_handler, sh],
            format='%(levelname)-10s %(message)s',
        )
    else:
        logging.basicConfig(
            level=logging.INFO,
            handlers=[log_handler, sh],
            format='%(levelname)-10s %(message)s',
        )

    version_number = _version.get_versions()['version']
    LGR.info(f'Currently running phys2cvr version {version_number}')
    LGR.info(f'Input file is {fname_func}')

    # Check if func is 1d, save its extension in any case, and read it
    func_is_1d, _, fname_ext = utils.check_ext(io.EXT_ARRAY, fname_func, remove=True)

    # Check that all input values have right type
    tr = utils.if_declared_force_type(tr, 'float', 'tr')
    freq = utils.if_declared_force_type(freq, 'float', 'freq')
    trial_len = utils.if_declared_force_type(trial_len, 'int', 'trial_len')
    n_trials = utils.if_declared_force_type(n_trials, 'int', 'n_trials')
    highcut = utils.if_declared_force_type(highcut, 'float', 'highcut')
    lowcut = utils.if_declared_force_type(lowcut, 'float', 'lowcut')
    lag_max = utils.if_declared_force_type(lag_max, 'float', 'lag_max')
    lag_min = utils.if_declared_force_type(lag_min, 'float', 'lag_min')
    lag_step = utils.if_declared_force_type(lag_step, 'float', 'lag_step')
    l_degree = utils.if_declared_force_type(l_degree, 'int', 'l_degree')
    if lag_max is not None and lag_min is None:
        if lag_max > 0:
            lag_min = -lag_max
        else:
            raise ValueError(
                'Given maximum lag is 0 or negative, but no minimum lag was provided. Halting execution.'
            )
    if lag_max is None and lag_min is not None:
        raise ValueError(
            'A minimum lag was provided without providing a maximum lag. Please rerun providing both or none.'
        )
    if lag_max is not None and lag_min >= lag_max:
        raise ValueError(
            f'Invalid lag range: lag_min ({lag_min}) >= lag_max ({lag_max}). Please provide a range where lag_min < lag_max.'
        )
    if starting_lag_max is not None and lag_increment is None:
        raise ValueError(
            f"You provided a value for starting_lag_max but did not specify lag_increment, which determines "
            "the step size for increasing the maximum lag. If you are going to specify starting_lag_max, "
            "please also specify lag_increment!"
        ) 
    if starting_lag_max>=lag_max:
        raise ValueError(
            f'If you are going to specify a starting_lag_max, it should be less than lag_max.'
        ) 
    if lag_map is not None and starting_lag_max is not None:
        raise ValueError(
            f'The option to provide a lag map to calculate CVR is not supported with an iterative lag maximum.'
        )    
    if l_degree < 0:
        raise ValueError(
            'The specified order of the Legendre polynomials must be >= 0.'
        )
    scale_factor = utils.if_declared_force_type(scale_factor, 'float', 'scale_factor')
    if r2model not in stats.R2MODEL:
        raise ValueError(
            f'R^2 model {r2model} not supported. Supported models are {stats.R2MODEL}'
        )

    LGR.info('Load functional data')
    if func_is_1d:
        if tr:
            func_avg = io.load_array(fname_func)
            LGR.info(f'Loading {fname_func}')
            if apply_filter:
                LGR.info('Applying butterworth filter to {fname_func}')
                func_avg = signal.filter_signal(
                    func_avg, tr, lowcut, highcut, butter_order
                )
        else:
            raise NameError(
                'Provided functional signal, but no TR specified! '
                'Rerun specifying the TR'
            )
    else:
        _, _, fname_ext = utils.check_ext(io.EXT_NIMG, fname_func, remove=True)
        try:
            func, dmask, img = io.load_nifti_get_mask(fname_func, dim=4)
        except nib.filebasedimages.ImageFileError:
            raise NotImplementedError(
                f'{fname_func} file type is not supported yet, or '
                'the extension was not specified.'
            )

        if len(func.shape) < 4:
            raise ValueError(f'Provided functional file {fname_func} is not a 4D file!')
        # Read TR or declare its overwriting
        if tr:
            LGR.warning(f'Forcing TR to be {tr} seconds')
        else:
            tr = img.header['pixdim'][4]

        # Read mask (and mask func) if provided
        if fname_mask:
            LGR.info('Load mask to restrict operations on functional data')
            _, mask, _ = io.load_nifti_get_mask(fname_mask, is_mask=True)
            if func.shape[:3] != mask.shape:
                raise ValueError(f'{fname_mask} and {fname_func} have different sizes!')
            mask = mask * dmask
            LGR.info(
                f'Masking {os.path.basename(fname_func)} using {os.path.basename(fname_mask)}'
            )
            func = func * mask[..., np.newaxis]
            roiref = os.path.basename(fname_mask)
        else:
            mask = dmask
            LGR.warning(
                f'No mask specified, using any voxel different from 0 in '
                f'{os.path.basename(fname_func)}'
            )
            roiref = os.path.basename(fname_func)

        # Read roi if provided
        if fname_roi:
            LGR.info('Load ROI to obtain a reference from functional data')
            _, roi, _ = io.load_nifti_get_mask(fname_roi, is_mask=True)
            if func.shape[:3] != roi.shape:
                raise ValueError(f'{fname_roi} and {fname_func} have different sizes!')
            roi = roi * mask
            roiref = os.path.basename(fname_roi)
        else:
            roi = mask
            LGR.warning(
                f'No ROI specified, using any voxel different from 0 in {roiref}'
            )

        LGR.info(f'Obtaining average signal in {roiref}')
        func_avg = func[roi].mean(axis=0)

        if apply_filter:
            LGR.info(f'Obtaining filtered average signal in {roiref}')
            func_avg = signal.filter_signal(func_avg, tr, lowcut, highcut, butter_order)

    LGR.info('Load physiological data')
    if fname_co2 is None:
        LGR.info(f'Computing "CVR" (approximation) maps using {fname_func} only')
        if func_is_1d:
            LGR.warning('Using an average signal only, solution might be unoptimal.')

            if apply_filter is None:
                LGR.warning(
                    'No filter applied to the input average! You know '
                    'what you are doing, right?'
                )

        # Get the SPC of the average rather than the average of the SPC
        # The former is more robust to intrinsic data noise than the latter
        petco2hrf = signal.spc(func_avg)

        # Reassign fname_co2 to fname_func for later use
        _, basename_co2, _ = utils.check_ext(
            io.EXT_ALL, f'avg_{os.path.basename(fname_func)}', scan=True, remove=True
        )

        outprefix = os.path.join(outdir, basename_co2)

        # If freq was declared, upsample the average GM to that.
        # Otherwise, set freq to inverse of TR.
        if freq is None:
            freq = 1 / tr
            LGR.info(f'No frequency declared, using 1/tr ({freq}Hz)')
        else:
            LGR.info(f'Resampling the average fMRI timeseries at {freq}Hz')
            upsamp_tps = int(np.round(petco2hrf.shape[-1] * tr * freq))
            petco2hrf = signal.resample_signal_samples(petco2hrf, upsamp_tps)
    else:
        co2_is_phys, _ = utils.check_ext('.phys', fname_co2)
        co2_is_1d, _ = utils.check_ext(io.EXT_ARRAY, fname_co2)

        if co2_is_1d:
            if fname_pidx:
                pidx = io.load_array(fname_pidx)
                pidx = pidx.astype(int)
            elif comp_endtidal:
                raise NameError(
                    f'{fname_co2} file is a text file, but no '
                    'file containing its peaks was provided. '
                    ' Please provide peak file!'
                )
            else:
                # pidx to None for compatibility
                pidx = None

            if freq is None:
                raise NameError(
                    f'{fname_co2} file is a text file, but no '
                    'frequency was specified. Please provide peak '
                    ' file!'
                )

            co2 = io.load_array(fname_co2)
        elif co2_is_phys:
            # Read a phys file!
            if freq:
                LGR.warning(f'Forcing CO2 frequency to be {freq} Hz')
                co2, pidx, _ = io.load_physio(fname_co2)
            else:
                co2, pidx, freq = io.load_physio(fname_co2)
        else:
            raise NotImplementedError(
                f'{fname_co2} file type is not supported yet, or '
                'the extension was not specified.'
            )

        _, basename_co2, _ = utils.check_ext(
            io.EXT_ALL, os.path.basename(fname_co2), scan=True, remove=True
        )

        outprefix = os.path.join(outdir, basename_co2)

        petco2hrf = compute_petco2hrf(
            co2, pidx, freq, outprefix, comp_endtidal, response_function
        )

    # If the user provided a lag map, read it and extract information from it
    if lag_map:
        LGR.info('Load lag map')
        lag, _, _ = io.load_nifti_get_mask(lag_map)
        if func.shape[:3] != lag.shape:
            raise ValueError(f'{lag_map} and {fname_func} have different sizes!')

        # Read lag_step and lag_max from file (or try to)
        lag = lag * mask

        lag_list = np.unique(lag[mask])

        if lag_step is None:
            # np.unique sorts results already
            lag_step = np.unique(np.round(lag_list[1:] - lag_list[:-1], 3))
            if lag_step.size > 1:
                # Check if extra steps found are multiple of the first within numerical tolerance
                if not np.isclose(np.mod(lag_step, lag_step[0]), 0).all():
                    raise ValueError(
                        f'phys2cvr found different delta lags in {lag_map}: {lag_step}'
                    )
            lag_step = lag_step[0]

            LGR.warning(f'phys2cvr detected a delta lag of {lag_step} seconds')
        else:
            LGR.warning(f'Forcing delta lag to be {lag_step}')

        if lag_max is None:
            lag_max = np.asarray(lag_list).max()
            lag_min = np.asarray(lag_list).min()
            LGR.warning(f'phys2cvr detected a lag range of [{lag_min}, {lag_max}]')
        else:
            LGR.warning(f'Forcing lag range to be [{lag_min}, {lag_max}]')

    # If a regressor dir is specified, try load the data,
    # If failing or otherwise, compute the regressors.
    if run_regression and regr_dir is not None:
        try:
            regr = io.load_array(
                os.path.join(regr_dir, '..', f'{outprefix}_petco2hrf.1D')
            )
            regr_shifts = io.load_array(os.path.join(regr_dir, 'co2_shifts.1D'))
        except OSError:
            LGR.warning(f'Regressor {outprefix}_petco2hrf.1D not found. Estimating it.')
            # Set regr_dir to None to activate next if statement.
            regr_dir = None

    if regr_dir is None:
        regr, regr_shifts = create_physio_regressor(
            func_avg,
            petco2hrf,
            tr,
            freq,
            outprefix,
            lag_max,
            lag_min,
            trial_len,
            n_trials,
            '.1D',
            lagged_regression,
            legacy,
            abs_xcorr,
            skip_xcorr,
        )

    # Run internal regression if required and possible!
    if not func_is_1d and run_regression:
        LGR.info('Running regression!')

        # Change dimensions in image header before export
        LGR.info('Prepare output image')
        _, fname_out_func, _ = utils.check_ext(
            fname_ext, os.path.basename(fname_func), remove=True
        )
        fname_out_func = os.path.join(outdir, fname_out_func)
        newdim = deepcopy(img.header['dim'])
        newdim[0], newdim[4] = 3, 1
        oimg = deepcopy(img)
        oimg.header['dim'] = newdim

        # Compute signal percentage change of functional data
        func = signal.spc(func)

        # Generate polynomial regressors (at least average) and assign them to denoise_matrix
        LGR.info(f'Compute Legendre polynomials up to order {l_degree}')
        denoise_matrix = create_legendre(l_degree, regr.size)

        # Read in eventual denoising factors
        if denoise_matrix_file:
            denoise_matrix_file = utils.if_declared_force_type(
                denoise_matrix_file, 'list', 'denoise_matrix_file'
            )
            for matrix in denoise_matrix_file:
                LGR.info(f'Read confounding factor from {matrix}')
                conf = io.load_array(matrix)
                denoise_matrix = np.hstack([denoise_matrix, conf])
        # Read in eventual extra factors
        if extra_matrix_file:
            extra_matrix_file = utils.if_declared_force_type(
                extra_matrix_file, 'list', 'extra_matrix_file'
            )
            matlist = []
            for matrix in extra_matrix_file:
                LGR.info(f'Read extra factor for orthogonalisation from {matrix}')
                matlist += [io.load_array(matrix)]
            extra_matrix = np.hstack(matlist)
        else:
            extra_matrix = None
        # Read in eventual orthogonalisable factors
        if orthogonalised_matrix_file:
            orthogonalised_matrix_file = utils.if_declared_force_type(
                orthogonalised_matrix_file, 'list', 'orthogonalised_matrix_file'
            )
            matlist = []
            for matrix in orthogonalised_matrix_file:
                LGR.info(f'Read confounding factor from {matrix}')
                matlist += [io.load_array(matrix)]
            orthogonalised_matrix = np.hstack(matlist)
        else:
            orthogonalised_matrix = None

        LGR.info('Compute simple CVR estimation (bulk shift only)')
        x1D = os.path.join(outdir, 'mat', 'mat_simple.1D')
        regr = utils.check_array_dim('regr', regr, shape='rectangle')
        beta, tstat, r_square = stats.regression(
            func,
            regr,
            denoise_matrix,
            orthogonalised_matrix,
            extra_matrix,
            mask,
            r2model,
            debug,
            x1D,
        )

        LGR.info('Export bulk shift results')
        if scale_factor is None:
            LGR.warning('Remember: CVR might not be in %BOLD/mmHg!')
        else:
            beta = beta / float(scale_factor)
        # Scale beta by scale factor while exporting (useful to transform V in mmHg)
        LGR.info('Export CVR and T-stat of simple regression')
        io.export_nifti(beta, oimg, f'{fname_out_func}_cvr_simple{fname_ext}')
        io.export_nifti(tstat, oimg, f'{fname_out_func}_tstat_simple{fname_ext}')

        if debug:
            LGR.debug('Export R^2 volume of simple regression')
            io.export_nifti(
                r_square, oimg, f'{fname_out_func}_r_square_simple{fname_ext}'
            )

        if lagged_regression and regr_shifts is not None and (lag_max and lag_step):
            # If user specified a lag map, run regression based on it (see "Load lag map")
            if lag_map is not None:
                LGR.info(
                    f'Running lagged CVR estimation with lag map {lag_map}! '
                    '(might take a while...)'
                )

                step = int(lag_step * freq)

                lag_idx = np.round((lag - lag_min) * freq / step).astype(int)
                lag_idx_list = np.unique(lag_idx)

                # Prepare empty matrices
                beta = np.empty_like(lag, dtype='float32')
                tstat = np.empty_like(lag, dtype='float32')

                for n, i in enumerate(lag_idx_list):
                    LGR.info(f'Perform L-GLM number {n + 1} of {len(lag_idx_list)}')
                    regr = regr_shifts[(i * step), :, np.newaxis]

                    x1D = os.path.join(outdir, 'mat', f'mat_{i:04g}.1D')

                    (beta[lag_idx == i], tstat[lag_idx == i], _) = stats.regression(
                        func[lag_idx == i],
                        regr,
                        denoise_matrix,
                        orthogonalised_matrix,
                        extra_matrix,
                        mask[lag_idx == i],
                        r2model,
                        debug,
                        x1D,
                    )

            else:
                LGR.info(
                    f'Running lagged CVR estimation with lag range = [{lag_min}, {lag_max}]! '
                    '(might take a while...)'
                )

                if legacy:
                    nrep_neg = int(abs(lag_min) * freq)
                    nrep_pos = int(abs(lag_max) * freq)
                    nrep = nrep_neg + nrep_pos
                else:
                    nrep_neg = int(abs(lag_min) * freq)
                    nrep_pos = int(abs(lag_max) * freq)
                    nrep = nrep_neg + nrep_pos + 1
                # Check the number of repetitions first
                if lag_step:
                    step = int(lag_step * freq)
                else:
                    step = 1
                lag_range = list(range(0, nrep, step))
                # Prepare empty matrices
                r_square_all = np.zeros(
                    list(func.shape[:3]) + [len(lag_range)], dtype='float32'
                )
                beta_all = np.zeros(
                    list(func.shape[:3]) + [len(lag_range)], dtype='float32'
                )
                tstat_all = np.zeros(
                    list(func.shape[:3]) + [len(lag_range)], dtype='float32'
                )

                for n, i in enumerate(lag_range):
                    LGR.info(f'Perform L-GLM number {n + 1} of {len(lag_range)}')
                    regr = regr_shifts[i, :, np.newaxis]

                    x1D = os.path.join(outdir, 'mat', f'mat_{i:04g}.1D')
                    (
                        beta_all[:, :, :, n],
                        tstat_all[:, :, :, n],
                        r_square_all[:, :, :, n],
                    ) = stats.regression(
                        func,
                        regr,
                        denoise_matrix,
                        orthogonalised_matrix,
                        extra_matrix,
                        mask,
                        r2model,
                        debug,
                        x1D,
                    )

                if debug:
                    LGR.debug('Export all betas, tstats, and R^2 volumes.')
                    newdim_all = deepcopy(img.header['dim'])
                    newdim_all[0], newdim_all[4] = 4, int(len(lag_range))
                    oimg_all = deepcopy(img)
                    oimg_all.header['dim'] = newdim_all
                    io.export_nifti(
                        r_square_all,
                        oimg_all,
                        f'{fname_out_func}_r_square_all{fname_ext}',
                    )
                    io.export_nifti(
                        tstat_all, oimg_all, f'{fname_out_func}_tstat_all{fname_ext}'
                    )
                    io.export_nifti(
                        beta_all, oimg_all, f'{fname_out_func}_beta_all{fname_ext}'
                    )

                # Find the right lag for CVR estimation
                if starting_lag_max is not None:
                    lag_idx = select_lag_avoid_boundary(
                        r_square_all=r_square_all,
                        mask=mask,
                        starting_max=starting_lag_max,
                        final_max=lag_max,
                        lag_min=lag_min,
                        lag_step=lag_step,
                        freq=freq,
                        expand_by=lag_increment,
                    )
                    lag = (lag_idx * step) / freq + (mask * lag_min)
                else:
                    lag_idx = np.argmax(r_square_all, axis=-1)
                    lag = (lag_idx * step) / freq + (mask * lag_min)
                # Express lag map relative to median of the roi
                lag_rel = lag - (mask * np.median(lag[roi]))

                # Run through indexes to pick the right value
                lag_idx_list = np.unique(lag_idx)
                beta = np.zeros_like(lag, dtype='float32')
                tstat = np.zeros_like(lag, dtype='float32')
                for i in lag_idx_list:
                    beta[lag_idx == i] = beta_all[:, :, :, i][lag_idx == i]
                    tstat[lag_idx == i] = tstat_all[:, :, :, i][lag_idx == i]

            LGR.info('Export fine shift results')
            if scale_factor is None:
                LGR.warning('Remember: CVR might not be in %BOLD/mmHg!')
            else:
                beta = beta / float(scale_factor)

            io.export_nifti(beta, oimg, f'{fname_out_func}_cvr{fname_ext}')
            io.export_nifti(tstat, oimg, f'{fname_out_func}_tstat{fname_ext}')
            if not lag_map:
                io.export_nifti(lag, oimg, f'{fname_out_func}_lag{fname_ext}')
                io.export_nifti(lag_rel, oimg, f'{fname_out_func}_lag_mkrel{fname_ext}')

    elif run_regression:
        LGR.warning(
            'The input file is not a nifti volume. At the moment, '
            'regression is not supported for other formats.'
        )

    LGR.info('phys2cvr finished! Enjoy your outputs!')


def _main(argv=None):
    options = _get_parser().parse_args(argv)

    options = _check_opt_conf(options)

    utils.save_bash_call(options.fname_func, options.outdir)

    phys2cvr(**vars(options))


if __name__ == '__main__':
    _main(sys.argv[1:])


"""
Copyright 2021-2026, Stefano Moia & phys2cvr contributors.

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
