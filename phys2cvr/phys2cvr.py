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

from phys2cvr import _version, io, signal, stats
from phys2cvr.cli.run import _check_opt_conf, _get_parser
from phys2cvr.io import EXT_1D, EXT_NIFTI

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
    arg_str = " ".join(sys.argv[1:])
    call_str = f"phys2cvr {arg_str}"
    if outdir:
        outdir = os.path.abspath(outdir)
    else:
        outdir = os.path.join(os.path.split(fname)[0], "phys2cvr")
    log_path = os.path.join(outdir, "logs")
    os.makedirs(log_path, exist_ok=True)
    isotime = datetime.datetime.now().strftime("%Y-%m-%dT%H%M%S")
    fname, _ = io.check_ext(".nii.gz", os.path.basename(fname), remove=True)
    f = open(os.path.join(log_path, f"p2c_call_{fname}_{isotime}.sh"), "a")
    f.write(f"#!bin/bash \n{call_str}")
    f.close()


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
    r2model="full",
    lag_max=None,
    lag_step=None,
    legacy=False,
    l_degree=0,
    denoise_matrix_file=[],
    orthogonalised_matrix_file=[],
    extra_matrix_file=[],
    scale_factor=None,
    lag_map=None,
    regr_dir=None,
    run_conv=True,
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
    r2model : 'full', 'partial' or 'intercept', optional
        Submit the model of the R^2 the regression should return (hence, which
        R^2 model the lag regression should be based on).
        Possible options are 'full', 'partial', 'intercept'.
        See `stats.regression` help to understand them.
        Default: 'full'
    lag_max : int or float, optional
        Limits (both positive and negative) of the temporal area to explore,
        expressed in seconds (e.g. ±9 seconds). Caution: this is not a pythonic
        range, but a real range, i.e. the upper limit is included (e.g. [-9, +9]).
        Default: None
    lag_step : int or float, optional
        Step of the lag to take into account in seconds.
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
    denoise_matrix_file : list of str(s) or path(s), optional
        Add one or multiple denoising matrices to the regression model.
        Ignored if not performing the regression step.
    orthogonalised_matrix_file : list of str(s) or path(s), optional
        Add one or multiple denoising matrices to the regression model,
        AFTER orthogonalising them w.r.t. the task, the denoise matrix,
        and the extra matrix.
        Ignored if not performing the regression step.
    extra_matrix_file : list of str(s) or path(s), optional
        Add one or multiple extra matrices to use in the orthogonalisation step.
        These matrices will not be added to the final regression model.
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
    NotImplementedError
        - If a file type is not supported yet.
    NameError
        - If functional timeseries is lacking TR and the latter was not specified.
        - If physiological file is a txt file and no peaks were provided.
        - If physiological file is lacking frequency and the latter was not specified.
    """
    # If lagged regression is selected, make sure run_regression is true.
    if lagged_regression:
        run_regression = True
   # Set the output directory 
   # Add logger and suff
    if outdir:
        outdir = os.path.abspath(outdir)
    else:
        # Default output directory is "phys2cvr" folder in the same directory as the input file
        outdir = os.path.join(os.path.split(fname_func)[0], "phys2cvr")
    outdir = os.path.abspath(outdir)

    # Create a directory for logs
    petco2log_path = os.path.join(outdir, "logs")
    os.makedirs(petco2log_path, exist_ok=True)

    # Create logfile name based on the current time
    basename = "phys2cvr_"
    extension = "tsv"
    isotime = datetime.datetime.now().strftime("%Y-%m-%dT%H%M%S")
    logname = os.path.join(petco2log_path, f"{basename}{isotime}.{extension}")

    # Set logging format
    log_formatter = logging.Formatter(
        "%(asctime)s\t%(name)-12s\t%(levelname)-8s\t%(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
    )

    # Set up logging file and open it for writing
    # Configure logging to write to both a file and the console
    log_handler = logging.FileHandler(logname)
    log_handler.setFormatter(log_formatter)
    sh = logging.StreamHandler()

    if quiet:
        # Log only warnings and errors if "quiet" mode is enabled
        logging.basicConfig(
            level=logging.WARNING,
            handlers=[log_handler, sh],
            format="%(levelname)-10s %(message)s",
        )
    elif debug:
        # Log detailed debug information if "debug" mode is enabled
        logging.basicConfig(
            level=logging.DEBUG,
            handlers=[log_handler, sh],
            format="%(levelname)-10s %(message)s",
        )
    else:
         # Default logging level is INFO
        logging.basicConfig(
            level=logging.INFO,
            handlers=[log_handler, sh],
            format="%(levelname)-10s %(message)s",
        )

    # Log the version of the software and the input file being processed
    version_number = _version.get_versions()["version"]
    LGR.info(f"Currently running phys2cvr version {version_number}")
    LGR.info(f"Input file is {fname_func}")

    # Check the type of the input file (1D text file or NIfTI image) and read it
    func_is_1d = io.check_ext(EXT_1D, fname_func)
    func_is_nifti = io.check_ext(EXT_NIFTI, fname_func)

    # Validate and convert input parameters to the correct types    
    # Check that all input values have right type
    tr = io.if_declared_force_type(tr, "float", "tr") # Repetition time (TR)
    freq = io.if_declared_force_type(freq, "float", "freq") #Sampling frequency
    trial_len = io.if_declared_force_type(trial_len, "int", "trial_len") # Trial lens 
    n_trials = io.if_declared_force_type(n_trials, "int", "n_trials") # Number of trials
    highcut = io.if_declared_force_type(highcut, "float", "highcut") #High cutoff frequency 
    lowcut = io.if_declared_force_type(lowcut, "float", "lowcut") #Low cutoff frequency
    lag_max = io.if_declared_force_type(lag_max, "float", "lag_max") #Maximum lag
    
    # Added by Cristina (line 321 to 341)
    # Handle symmetric and asymmetric lag ranges
    if isinstance(lag_max, (list, tuple)) and len(lag_max) == 2:
         # If lag_max is a tuple or list with two values, treat it as an asymmetric range
        lag_min= float(lag_max[0])  # First value is the minimum lag
        lag_max= float(lag_max[1]) # Second value is the maximum lag
        LGR.info(f"Using asymmetric lag range: [{lag_min}, {lag_max}]")
    
    elif isinstance(lag_max, (int, float)):
        # If lag_max is a single value, treat it as a symmetric range
        lag_min =  -float(lag_max)
        lag_max = float(lag_max)
        LGR.info(f"Using symmetric lag range: [{lag_min}, {lag_max}]")

    else:  
        # Raise an error if lag_max is invalid
       raise ValueError(
            "Invalid lag_max value. Using default symmetric range."
        )
    
    # Validate and convert additional input parameters
    lag_step = io.if_declared_force_type(lag_step, "float", "lag_step")
    l_degree = io.if_declared_force_type(l_degree, "int", "l_degree")
    if l_degree < 0:
        raise ValueError(
            "The specified order of the Legendre polynomials must be >= 0."
        )
    scale_factor = io.if_declared_force_type(scale_factor, "float", "scale_factor")

    # Validate the R^2 model
    if r2model not in stats.R2MODEL:
         # Raise an error if the R^2 model is not supported
        raise ValueError(
            f"R^2 model {r2model} not supported. Supported models "
            f"are {stats.R2MODEL}"
        )

    # Check if the input file is a 1D text file
    if func_is_1d:
        if tr:
            # Load the functional signal from the 1D text file
            func_avg = np.genfromtxt(fname_func)
            LGR.info(f"Loading {fname_func}")
             # Apply a Butterworth filter to the functional signal
            if apply_filter:
                LGR.info("Applying butterworth filter to {fname_func}")
                func_avg = signal.filter_signal(
                    func_avg, tr, lowcut, highcut, butter_order
                )
        else:
             # Raise an error if TR is not specified for a 1D input file
            raise NameError(
                "Provided functional signal, but no TR specified! "
                "Rerun specifying the TR"
            )
    elif func_is_nifti:
         # Load the functional data from a NIfTI file
        func, dmask, img = io.load_nifti_get_mask(fname_func, dim=4)
        if len(func.shape) < 4:
            # Raise an error if the NIfTI file is not 4D
            raise ValueError(f"Provided functional file {fname_func} is not a 4D file!")
        # Read TR or declare its overwriting
        if tr:
            LGR.warning(f"Forcing TR to be {tr} seconds")
        else:
            tr = img.header["pixdim"][4]

        # Read mask (and mask func) if provided
        if fname_mask:
            _, mask, _ = io.load_nifti_get_mask(fname_mask, is_mask=True)
            if func.shape[:3] != mask.shape:
                 # Raise an error if the mask dimensions do not match the functional data
                raise ValueError(f"{fname_mask} and {fname_func} have different sizes!")
            mask = mask * dmask # Combine the mask with the data mask
            LGR.info(
                f"Masking {os.path.basename(fname_func)} using {os.path.basename(fname_mask)}"
            )
            func = func * mask[..., np.newaxis] # Apply the mask to the functional data
            roiref = os.path.basename(fname_mask)
        else:
            # Use the data mask if no mask file is provided
            mask = dmask
            LGR.warning(
                f"No mask specified, using any voxel different from 0 in "
                f"{os.path.basename(fname_func)}"
            )
            roiref = os.path.basename(fname_func)

        # Read roi if provided
        if fname_roi:
            _, roi, _ = io.load_nifti_get_mask(fname_roi, is_mask=True)
            if func.shape[:3] != roi.shape:
                # Raise an error if the ROI dimensions do not match the functional data
                raise ValueError(f"{fname_roi} and {fname_func} have different sizes!")
            roi = roi * mask # Combine the ROI with the mask
            roiref = os.path.basename(fname_roi)
        else:
            # Use the mask as the ROI if no ROI file is provided
            roi = mask
            LGR.warning(
                f"No ROI specified, using any voxel different from 0 in " f"{roiref}"
            )

        if apply_filter:
             # Apply a Butterworth filter to the functional data and compute the average signal
            LGR.info(f"Obtaining filtered average signal in {roiref}")
            func_filt = signal.filter_signal(func, tr, lowcut, highcut, butter_order)
            func_avg = func_filt[roi].mean(axis=0)
        else:
            # Compute the average signal without filtering
            LGR.info(f"Obtaining average signal in {roiref}")
            func_avg = func[roi].mean(axis=0)

    else:
        # Raise an error if the input file type is not supported
        raise NotImplementedError(
            f"{fname_func} file type is not supported yet, or "
            "the extension was not specified."
        )

    if fname_co2 is None:
        # If no CO2 file is provided, compute CVR maps using only the functional data
        LGR.info(f'Computing "CVR" (approximation) maps using {fname_func} only')
        if func_is_1d:
             # Warn the user that using only the average signal might not be optimal
            LGR.warning("Using an average signal only, solution might be unoptimal.")

            if apply_filter is None:
                # Warn the user if no filter is applied to the input average
                LGR.warning(
                    "No filter applied to the input average! You know "
                    "what you are doing, right?"
                )
        # Compute the SPC (signal percentage change) of the average signal (rather than the average of the SPC)
        # The former is more robust to intrinsic data noise than the latter
        petco2hrf = signal.spc(func_avg)

         # Generate a base name for the output CO2 file (using functional file name as a base)
        # Reassign fname_co2 to fname_func for later use - calling splitext twice cause .gz
        basename_co2 = os.path.splitext(
            os.path.splitext(f"avg_{os.path.basename(fname_func)}")[0]
        )[0]
        outname = os.path.join(outdir, basename_co2)

        # If freq was declared, upsample the average GM to that.
        # Otherwise, set freq to inverse of TR.
        if freq is None:
            freq = 1 / tr
            LGR.info(f"No frequency declared, using 1/tr ({freq}Hz)")
        else:
            # Resample the average fMRI timeseries to the specified frequency
            LGR.info(f"Resampling the average fMRI timeseries at {freq}Hz")
            petco2hrf = signal.resample_signal(petco2hrf, 1 / tr, freq)
    else:
         # If a CO2 file is provided, check its type (.phys or 1D)
        co2_is_phys = io.check_ext(".phys", fname_co2)
        co2_is_1d = io.check_ext(EXT_1D, fname_co2)

        if co2_is_1d:
            if fname_pidx:
                # Load the peaks file if provided
                pidx = np.genfromtxt(fname_pidx)
                pidx = pidx.astype(int)
            elif run_conv:
                # Raise an error if peaks are required but not provided
                raise NameError(
                    f"{fname_co2} file is a text file, but no "
                    "file containing its peaks was provided. "
                    " Please provide peak file!"
                )

            if freq is None:
                 # Raise an error if frequency is required but not provided
                raise NameError(
                    f"{fname_co2} file is a text file, but no "
                    "frequency was specified. Please provide peak "
                    " file!"
                )
             # Load the CO2 timeseries from the 1D text file
            co2 = np.genfromtxt(fname_co2)
        elif co2_is_phys:
            # Read a phys file!
            phys = load_physio(fname_co2, allow_pickle=True)

            co2 = phys.data # Extract the CO2 timeseries
            pidx = phys.peaks # Extract the CO2 peaks
            if freq:
                 # Warn the user if the frequency is being overwritten
                LGR.warning(f"Forcing CO2 frequency to be {freq} Hz")
            else:
                freq = phys.fs
        else:
             # Raise an error if the CO2 file type is not supported
            raise NotImplementedError(
                f"{fname_co2} file type is not supported yet, or "
                "the extension was not specified."
            )
        
        # Generate a base name for the output CO2 file
        # Set output file & path - calling splitext twice cause .gz
        basename_co2 = os.path.splitext(
            os.path.splitext(os.path.basename(fname_co2))[0]
        )[0]
        outname = os.path.join(outdir, basename_co2)

        # Unless user asks to skip this step, convolve the end tidal signal.
         # If the user asks to skip convolution, use the raw CO2 signal
        if run_conv is False:
            petco2hrf = co2
        else:
            # Convolve the CO2 signal with the peaks to generate the HRF
            petco2hrf = signal.convolve_petco2(co2, pidx, freq, outname)

    # If a regressor directory is not specified, compute the regressors.
    if regr_dir is None:
        # Generate regressors using the functional average and CO2 HRF
        regr, regr_shifts = stats.get_regr(
            func_avg, 
            petco2hrf, 
            tr,
            freq,
            outname, # output file base name
            lag_max,
            trial_len,
            n_trials,
            ".1D",
            lagged_regression, # whether to compute lagged regressors
            legacy, # Use legacy mode for ranges
            abs_xcorr, # Use absolute cross-correlation
            skip_xcorr, # Skip cross-correlation step
        )
    elif run_regression:
         # If regressors are precomputed, try to load them
        try:
            regr = np.genfromtxt(f"{outname}_petco2hrf.1D")
        except IOError:
              # If the regressor file is not found, compute it
            LGR.warning(
                f"Regressor {outname}_petco2hrf.1D not found. " "Estimating it."
            )
            regr, regr_shifts = stats.get_regr(
                func_avg,
                petco2hrf,
                tr,
                freq,
                outname,
                lag_max,
                trial_len,
                n_trials,
                ".1D",
                lagged_regression,
                legacy,
                abs_xcorr,
            )

    # Run internal regression if required and possible!
    if func_is_nifti and run_regression:
        LGR.info("Running regression!")

        # Change dimensions in image header before export
        LGR.info("Prepare output image")
        fname_out_func, _ = io.check_ext(
            ".nii.gz", os.path.basename(fname_func), remove=True
        )
        fname_out_func = os.path.join(outdir, fname_out_func)
        newdim = deepcopy(img.header["dim"])
        newdim[0], newdim[4] = 3, 1  # Set the header to 3D with 1 timepoint
        oimg = deepcopy(img)
        oimg.header["dim"] = newdim

        # Compute signal percentage change of functional data
        func = signal.spc(func)

        # Generate polynomial regressors (at least average) and assign them to denoise_matrix for denosing
        LGR.info(f"Compute Legendre polynomials up to order {l_degree}")
        denoise_matrix = stats.get_legendre(l_degree, regr.size)

        # Read in eventual denoising factors
        # Read additional denoising matrices if provided
        if denoise_matrix_file:
            denoise_matrix_file = io.if_declared_force_type(
                denoise_matrix_file, "list", "denoise_matrix_file"
            )
            for matrix in denoise_matrix_file:
                LGR.info(f"Read confounding factor from {matrix}")
                conf = np.genfromtxt(matrix)
                denoise_matrix = np.hstack([denoise_matrix, conf])

        # Read in eventual extra factors
        # Read extra matrices for orthogonalization if provided
        if extra_matrix_file:
            extra_matrix_file = io.if_declared_force_type(
                extra_matrix_file, "list", "extra_matrix_file"
            )
            matlist = []
            for matrix in extra_matrix_file:
                LGR.info(f"Read extra factor for orthogonalisation from {matrix}")
                matlist += [np.genfromtxt(matrix)]
            extra_matrix = np.hstack(matlist)
        else:
            extra_matrix = None

        # Read in eventual orthogonalisable factors
        # Read orthogonalizable matries if provided
        if orthogonalised_matrix_file:
            orthogonalised_matrix_file = io.if_declared_force_type(
                orthogonalised_matrix_file, "list", "orthogonalised_matrix_file"
            )
            matlist = []
            for matrix in orthogonalised_matrix_file:
                LGR.info(f"Read confounding factor from {matrix}")
                matlist += [np.genfromtxt(matrix)]
            orthogonalised_matrix = np.hstack(matlist)
        else:
            orthogonalised_matrix = None

        # Perform simple CVR estimation (bulk shift only)
        LGR.info("Compute simple CVR estimation (bulk shift only)")
        x1D = os.path.join(outdir, "mat", "mat_simple.1D")
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

        # Export the results of the simple regression
        LGR.info("Export bulk shift results")
        if scale_factor is None:
            LGR.warning("Remember: CVR might not be in %BOLD/mmHg!")
        else:
            beta = beta / float(scale_factor) # Scale beta by the scale factor
        # Scale beta by scale factor while exporting (useful to transform V in mmHg)
        # Export CVR and T-stat maps
        LGR.info("Export CVR and T-stat of simple regression")
        io.export_nifti(beta, oimg, f"{fname_out_func}_cvr_simple")
        io.export_nifti(tstat, oimg, f"{fname_out_func}_tstat_simple")

        if debug:
            # Export R^2 volume if in debug mode
            LGR.debug("Export R^2 volume of simple regression")
            io.export_nifti(r_square, oimg, f"{fname_out_func}_r_square_simple")

        if (
            lagged_regression
            and regr_shifts is not None
            and ((lag_max and lag_step) or lag_map)
        ):
            # If lagged regression is enabled and regressors are available
            # Perform lagged CVR estimation
            if lag_max:
                 # Log the minimum and maximum lag being used 
                LGR.info(
                    f"Running lagged CVR estimation with lag range: [{lag_min}, {lag_max}]"
                    "(might take a while...)"
                )
            elif lag_map is not None:
                # Log the use of a lag map
                LGR.info(
                    f"Running lagged CVR estimation with lag map {lag_map}! "
                    "(might take a while...)"
                )
             # Determine the number of repetitions (nrep) based on lag_max and frequency    
            if legacy:
                # Legacy mode: use pythonic ranges (exclude the upper limit)
                nrep = int(lag_max * freq * 2)
            else:
                # Legacy mode: use pythonic ranges (exclude the upper limit)
                nrep = int(lag_max * freq * 2) + 1

            if regr_dir:
                # If a regressor directory is specified, set the output prefix
                outprefix = os.path.join(regr_dir, os.path.split(outname)[1])

            # If user specified a lag map, use that one to regress things
            if lag_map:
                lag, _, _ = io.load_nifti_get_mask(lag_map)
                if func.shape[:3] != lag.shape:
                    raise ValueError(
                        f"{lag_map} and {fname_func} have different sizes!"
                    )

                # Read lag_step and lag_max from file (or try to)
                # Apply the mask to the lag map
                lag = lag * mask

                # Extract unique lag values from the lag map
                lag_list = np.unique(lag)

                if lag_step is None:
                    # If lag_step is not provided, calculate it from the lag map
                    lag_step = np.unique(lag_list[1:] - lag_list[:-1])
                    if lag_step.size > 1:
                         # Raise an error if the lag map has inconsistent step sizes
                        raise ValueError(
                            f"phys2cvr found different delta lags in {lag_map}"
                        )
                    else:
                        # Log the detected lag step
                        LGR.warning(
                            f"phys2cvr detected a delta lag of {lag_step} seconds"
                        )
                else:
                     # Log the forced lag step
                    LGR.warning(f"Forcing delta lag to be {lag_step}")

                # Convert lag_step to samples
                step = int(lag_step * freq)

                if lag_max is None:
                     # If lag_max is not provided, calculate it from the lag map
                    lag_max = np.abs(lag_list).max()
                    LGR.warning(f"phys2cvr detected a max lag of {lag_max} seconds")
                else:
                    # Log the forced lag_max
                    LGR.warning(f"Forcing max lag to be {lag_max}")

                # Convert lag values to indices
                lag_idx = (lag + lag_max) * freq / step

                lag_idx_list = np.unique[lag_idx]

                # Initialize empty matrices for beta and t-stat values
                beta = np.empty_like(lag, dtype="float32")
                tstat = np.empty_like(lag, dtype="float32")

                for i in lag_idx_list:
                     # Perform lagged regression for each unique lag index
                    LGR.info(
                        f"Perform L-GLM for lag {lag_list[i]} ({i + 1} of "
                        f"{len(lag_idx_list)}"
                    )
                    try:
                         # Use precomputed regressors if available
                        regr = regr_shifts[:, (i * step)]
                    except NameError:
                         # Otherwise, load regressors from file
                        regr = np.genfromtxt(f"{outprefix}_{i:04g}")
                    
                    # Perform regression for the current lag
                    x1D = os.path.join(outdir, "mat", f"mat_{i:04g}.1D")
                    (beta[lag_idx == i], tstat[lag_idx == i], _) = stats.regression(
                        func[lag_idx == i],
                        regr,
                        denoise_matrix,
                        orthogonalised_matrix,
                        extra_matrix,
                        [lag_idx == i],
                        r2model,
                        debug,
                        x1D,
                    )

            else:
                # If no lag map is provided, calculate lagged regressors dynamically
                # Check the number of repetitions first
                if lag_step:
                     # Convert lag_step to samples
                    step = int(lag_step * freq)
                else:
                    step = 1
                #this was the lag range implemment by stefano as it is symmetric
                # lag_range = list(range(0, nrep, step))
                
                # Create Cristina's lag range taking into account assymetric lag and symmetric lag
                lag_range = np.arange(lag_min, lag_max + lag_step, lag_step)
                LGR.info(f"Generated lag range: {lag_range}")

                # Prepare empty matrices to store results for all lags
                r_square_all = np.empty(
                    list(func.shape[:3]) + [len(lag_range)], dtype="float32"
                ) # Matrix to store R^2 values for all lags
                beta_all = np.empty(
                    list(func.shape[:3]) + [len(lag_range)], dtype="float32"
                ) # Matrix to store beta values for all lags
                tstat_all = np.empty(
                    list(func.shape[:3]) + [len(lag_range)], dtype="float32"
                ) # Matrix to store t-stat values for all lags

                # Loop through each lag in the generated lag range
                for n, i in enumerate(lag_range):
                    # Log the progress of the lagged regression
                    LGR.info(f"Perform L-GLM number {n + 1} of {len(lag_range)}")
                    try:
                        # Use precomputed regressors if available
                        regr = regr_shifts[:, i]
                        LGR.debug(f"Using shift {i} from matrix in memory: {regr}")
                    except NameError:
                        # Otherwise, load regressors from file
                        regr = np.genfromtxt(f"{outprefix}_{i:04g}")
                        LGR.debug(f"Reading shift {i} from file {outprefix}_{i:04g}")

                    # Define the output path for the current lag's design matrix
                    x1D = os.path.join(outdir, "mat", f"mat_{i:04g}.1D")
                    (
                        beta_all[:, :, :, n], # Store beta values for the current lag
                        tstat_all[:, :, :, n], # Store t-stat values for the current lag
                        r_square_all[:, :, :, n], # Store R^2 values for the current lag
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

                # If debug mode is enabled, export all betas, t-stats, and R^2 volumes
                if debug:
                    LGR.debug("Export all betas, tstats, and R^2 volumes.")
                    newdim_all = deepcopy(img.header["dim"])
                    newdim_all[0], newdim_all[4] = 4, int(len(lag_range))
                    oimg_all = deepcopy(img)
                    oimg_all.header["dim"] = newdim_all
                    # Export all beta, t-stat, and R^2 volumes for all lags
                    io.export_nifti(
                        r_square_all, oimg_all, f"{fname_out_func}_r_square_all"
                    )
                    io.export_nifti(tstat_all, oimg_all, f"{fname_out_func}_tstat_all")
                    io.export_nifti(beta_all, oimg_all, f"{fname_out_func}_beta_all")

                # Find the optimal lag for CVR estimation based on R^2 values
                lag_idx = np.argmax(r_square_all, axis=-1) # Index of the maximum R^2 value along the lag dimension
                lag = (lag_idx * step) / freq - (mask * lag_max) # Convert lag index to time in seconds
                
                # Express lag map relative to median of the roi
                lag_rel = lag - (mask * np.median(lag[roi]))

                # Run through indexes to pick the right value
                lag_idx_list = np.unique(lag_idx) # Get unique lag indices
                beta = np.empty_like(lag, dtype="float32") # Initialize empty matrix for beta values
                tstat = np.empty_like(lag, dtype="float32") # Initialize empty matrix for t-stat values
                # Loop through each unique lag index
                for i in lag_idx_list:
                    # Assign beta and t-stat values for the current lag index
                    beta[lag_idx == i] = beta_all[:, :, :, i][lag_idx == i]
                    tstat[lag_idx == i] = tstat_all[:, :, :, i][lag_idx == i]
            
            # Export the final results for the optimal lag  
            LGR.info("Export fine shift results")
            if scale_factor is None:
                # Warn the user if no scale factor is provided
                LGR.warning("Remember: CVR might not be in %BOLD/mmHg!")
            else:
                # Scale beta values by the scale factor
                beta = beta / float(scale_factor)

            # Export the CVR and T-stat maps 
            io.export_nifti(beta, oimg, f"{fname_out_func}_cvr")
            io.export_nifti(tstat, oimg, f"{fname_out_func}_tstat")
            
            # If no lag map is provided, export the lag map and relative lag map
            if not lag_map:
                io.export_nifti(lag, oimg, f"{fname_out_func}_lag")
                io.export_nifti(lag_rel, oimg, f"{fname_out_func}_lag_mkrel")

    elif run_regression:
        # Warning if the input file is not a NIfTI volume and regression is requested
        LGR.warning(
            "The input file is not a nifti volume. At the moment, "
            "regression is not supported for other formats."
        )

    LGR.info("phys2cvr finished! Enjoy your outputs!")


def _main(argv=None):
    """
    Main entry point for the script. Parses command-line arguments,
    validates them, and runs the `phys2cvr` function.
    """
    # Parse command-line arguments
    options = _get_parser().parse_args(argv)

    # Validate and configure the parsed options
    options = _check_opt_conf(options)

    # Save the bash call to a log file for reproducibility
    save_bash_call(options.fname_func, options.outdir)

    # Run the main phys2cvr function with the parsed options
    phys2cvr(**vars(options))


if __name__ == "__main__":

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
