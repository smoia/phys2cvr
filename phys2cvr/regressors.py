#!/usr/bin/env python3
"""
Module that creates the regressors used in phys2cvr.

Attributes
----------
LGR
    Logger
"""

import logging
import os

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view as swv

from phys2cvr.io import export_regressor
from phys2cvr.signal import (
    convolve_signal,
    endtidal_interpolation,
    resample_signal_freqs,
)
from phys2cvr.stats import x_corr
from phys2cvr.viz import plot_two_timeseries, plot_xcorr

LGR = logging.getLogger(__name__)
LGR.setLevel(logging.INFO)


def create_legendre(degree, length):
    """
    Produce the Legendre polynomials of order `degree`.

    Parameters
    ----------
    degree : int
        Highest number of desired orders.
    length : int
        Length of the desired polynomials (number of samples).

    Returns
    -------
    legendre : np.ndarray
        A `degree`*`length` array which includes all the polynomials up to order `degree`.
    """

    def _bonnet(d, x):
        """Use Bonnet method to create Leg polys."""
        if d == 0:
            return np.ones_like(x)
        if d == 1:
            return x
        return ((2 * d - 1) * x * _bonnet(d - 1, x) - (d - 1) * _bonnet(d - 2, x)) / d

    x = np.linspace(-1, 1, length)
    legendre = np.empty((length, degree + 1), dtype='float32')
    for n in range(degree + 1):
        legendre[:, n] = _bonnet(n, x)
    return legendre


def compute_petco2hrf(
    co2, pidx, freq, outprefix, comp_endtidal=True, response_function='hfr', mode='full'
):
    """
    Create the PetCO2 trace from CO2 trace, then convolve it to obtain the PetCO2hrf.

    Parameters
    ----------
    co2 : np.ndarray
        CO2 (or physiological) regressor
    pidx : np.ndarray
        Indices of peaks
    freq : str, int, or float
        Sample frequency of the CO2 regressor
    outprefix : str
        Prefix of the output file (i.e., the regressor of interest).
    comp_endtidal : bool
        If True, interpolate input regressor
    response_function : {`hrf`, `rrf`, `crf`}, None, str, path, or 1D array-like , optional
        Name of the response function to be used in the convolution of the regressor of
        interest or path to a 1D file containing one. Default is `hrf`.
        See `.signal.convolve_signal` for more information.
    mode : {'full', 'valid', 'same'} str, optional
        Convolution mode, see numpy.convolve.

    Returns
    -------
    petco2hrf : np.ndarray
        Convolved PetCO2 trace.

    Raises
    ------
    NotImplementedError
        If the provided CO2 is not a 1D array.
    """
    co2 = co2.squeeze()
    if co2.ndim > 1:
        raise NotImplementedError(
            'Arrays with more than 2 dimensions are not supported.'
        )

    if comp_endtidal and pidx is not None:
        petco2 = endtidal_interpolation(co2, pidx, axis=-1)

        # Plot PetCO2 vs CO2
        plot_two_timeseries(
            petco2, co2, f'{outprefix}_co2_vs_petco2.png', 'PetCO2', 'CO2', freq
        )

        # Demean and export
        petco2 = petco2 - petco2.mean()
        np.savetxt(f'{outprefix}_petco2.1D', petco2, fmt='%.18f')
    elif comp_endtidal and pidx is None:
        raise ValueError(
            'End-tidal interpolation was requested, but peaks were not provided.'
        )
    else:
        LGR.info(
            'Skipping End Tidal interpolation of PetCO2 trace (if you provided raw CO2 '
            'data, then you probably should not be doing this)'
        )
        petco2 = co2 - co2.mean()

    # Plot convolved PetCO2 vs PetCO2
    if response_function not in [None, 'none', 'None', 'NONE']:
        petco2hrf = convolve_signal(petco2, freq, response_function, mode)
        plot_two_timeseries(
            petco2hrf,
            petco2,
            f'{outprefix}_petco2_vs_petco2hrf.png',
            'Convolved PetCO2',
            'PetCO2',
            freq,
        )
    else:
        LGR.info('Skipping convolution of PetCO2 trace')
        petco2hrf = petco2

    return petco2hrf


def compute_bulk_shift(
    func_upsampled,
    petco2hrf,
    freq,
    outprefix,
    trial_len=None,
    n_trials=None,
    abs_xcorr=False,
):
    """
    Compute (initial) bulk shift of regressor.

    Parameters
    ----------
    func_upsampled : np.ndarray
        Functional timeseries average upsampled at the frequency of the regressor of interest.
    petco2hrf : np.ndarray
        Regressor of interest
    freq : str, int, or float
        Sample frequency of petco2hrf
    outprefix : list or path
        Path to output directory for regressors.
    trial_len : str or int, optional
        Length of each single trial for tasks that have more than one
        (E.g. BreathHold, CO2 challenges, ...)
        Used to improve cross correlation estimation.
        Default: None
    n_trials : str or int, optional
        Number of trials in the task.
        Default: None
    abs_xcorr : bool, optional
        If True, the cross correlation will consider the maximum absolute
        correlation, i.e. if a negative correlation is higher than the highest
        positive, the negative correlation will be chosen instead.

    Returns
    -------
    optshift : int
        The index of optimal shifting computed via Xcorr
    """
    first_tp, n_shifts = 0, None

    if trial_len and n_trials:
        # If both are specified, disregard two extreme _trial from matching.
        LGR.info(f'Specified {n_trials} trials lasting {trial_len} seconds')
        if n_trials > 3:
            LGR.info('Ignoring first trial to improve bulk shift estimation')
            first_tp = int(trial_len * freq)
        else:
            LGR.info('Using all trials for bulk shift estimation')
        if n_trials > 4:
            LGR.info('Ignoring last trial to improve bulk shift estimation')
            n_shifts = first_tp * (n_trials - 2)
    elif trial_len and not n_trials:
        LGR.warning(
            'The length of trial was specified, but the number of '
            'trials was not. Using all available trials for bulk shift estimation'
        )
    elif not trial_len and n_trials:
        LGR.warning(
            'The number of trials was specified, but the length of '
            'trial was not. Using all available trials for bulk shift estimation'
        )
    else:
        LGR.info('Using all trials for bulk shift estimation.')

    # Preparing breathhold and CO2 trace for Xcorr
    func_cut = func_upsampled[first_tp:]
    _, optshift, xcorr = x_corr(
        func_cut, petco2hrf, n_shifts=n_shifts, offset=first_tp, abs_xcorr=abs_xcorr
    )

    LGR.info(f'Cross correlation estimated a bulk shift of {optshift / freq} seconds')
    # Export estimated optimal shift in seconds
    with open(f'{outprefix}_optshift.1D', 'w') as f:
        print(f'{(optshift / freq):.4f}', file=f)

    # Export xcorr figure
    plot_xcorr(xcorr, outprefix, freq)

    # This shouldn't happen, but still check
    if optshift + func_upsampled.shape[0] > len(petco2hrf):
        raise Exception(
            f'The identified optimal shift {optshift / freq} removes too many samples to '
            'continue.'
        )

    return optshift


def create_fine_shift_regressors(
    petco2hrf,
    optshift,
    lag_max,
    lag_min,
    freq,
    func_size,
    func_upsamp_size,
    outprefix,
    ext='.1D',
    legacy=False,
):
    """
    Compute fine shifts to further optimize shifts.

    Parameters
    ----------
    petco2hrf : np.ndarray
        Regressor of interest
    optshift : int
        The index shift computed by the Xcorr/bulk shift
    lag_max : int, float, or None, optional
        Upper limit of the temporal area to explore, expressed in seconds.
        Caution: this is not a pythonic range, but a real range, i.e. the upper limit is included.
        Default: None
    lag_min : int, float, or None, optional
        Lower limit of the temporal area to explore, expressed in seconds.
        If set to None, and lag_max is not None and is positive, lag_min defaults to -lag_max (symmetric range).
        Default: None
    freq : str, int, or float
        Sample frequency of petco2hrf
    func_size : int
        Total timepoints of functional timeseries
    func_upsamp_size : int
        Total timepoints of functional timeseries, resampled at `freq` frequency
    outprefix : list or path
        Path to output directory for regressors.
    ext : str, optional
        Extension to be used for the exported regressors.
    legacy : bool, optional
        If True, exclude the upper lag limit from the regression estimation.
        If True, the maximum number of regressors will be `(freq*lag_max*2)`

    Returns
    -------
    petco2hrf_lagged : np.ndarray
        The shifted versions of the regressor of interest.
    """
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

    outdir, base = os.path.split(outprefix)
    regr_dir = os.path.join(outdir, 'regr')
    os.makedirs(regr_dir, exist_ok=True)
    outprefix = os.path.join(regr_dir, base)

    neg_shifts = int(abs(lag_min) * freq)
    pos_shifts = int(abs(lag_max) * freq) if legacy else int(abs(lag_max) * freq) + 1

    # Padding regressor right for shifts if not enough timepoints
    # Padding regressor left for shifts and update optshift if less than neg_shifts.
    rpad = max(0, int(func_upsamp_size + optshift + neg_shifts - petco2hrf.shape[0]))
    lpad = max(0, int(pos_shifts - optshift + 1))

    petco2hrf = np.pad(petco2hrf, (int(lpad), int(rpad)), 'mean')

    # Create sliding window view into petco2hrf, -1 because of reversed indexing
    neg_idx = optshift - pos_shifts + lpad - 1
    pos_idx = optshift + neg_shifts + lpad - 1
    # select the right windows the other way round
    petco2hrf_lagged = swv(petco2hrf, func_upsamp_size)[pos_idx:neg_idx:-1].copy()

    petco2hrf_lagged = export_regressor(
        petco2hrf_lagged, func_size, outprefix, 'shifts', ext
    )
    return petco2hrf_lagged


def create_physio_regressor(
    func_avg,
    petco2hrf,
    tr,
    freq,
    outprefix,
    lag_max=None,
    lag_min=None,
    trial_len=None,
    n_trials=None,
    ext='.1D',
    lagged_regression=True,
    legacy=False,
    abs_xcorr=False,
    skip_xcorr=False,
):
    """
    Create regressor(s) of interest for nifti GLM.

    Parameters
    ----------
    func_avg : np.ndarray
        Average functional timeseries (1D)
    petco2hrf : np.ndarray
        Regressor of interest (e.g., CO2 regressor)
    tr : str, int, or float
        Repetition time (TR) of timeseries
    freq : str, int, or float
        Sample frequency of petco2hrf
    outprefix : list or path
        Path to output directory for computed regressors.
    lag_max : int, float, or None, optional
        Upper limit of the temporal area to explore, expressed in seconds.
        Caution: this is not a pythonic range, but a real range, i.e. the upper limit is included.
        Default: None
    lag_min : int, float, or None, optional
        Lower limit of the temporal area to explore, expressed in seconds.
        If set to None, and lag_max is not None and is positive, lag_min defaults to -lag_max (symmetric range).
        Default: None
    trial_len : str or int, optional
        Length of each individual trial for timeseries which include more than one trial
        (e.g., multiple BreathHold trials, trials within CO2 challenges, ...)
        Used to improve cross correlation estimation.
        Default: None
    n_trials : str or int, optional
        Number of trials within the timeseries.
        Default: None
    ext : str, optional
        Extension to be used for the exported regressors (e.g., .txt, .csv)
    lagged_regression : bool, optional
        Estimate regressors for each possible lag of `petco2hrf`.
    legacy : bool, optional
        If True, exclude the upper (positive) lag limit from the regression estimation,
        i.e., the maximum number of regressors will be `(freq*lag_max*2)`
        If False, the maximum number of regressors will be `(freq*lag_max*2)+1`
    abs_xcorr : bool, optional
        If True, the cross correlation will consider the maximum absolute
        correlation, i.e., if a negative correlation is stronger than the strongest
        positive, the negative correlation will be used.
    skip_xcorr : bool, optional
        If True, skip the cross correlation step.

    Returns
    -------
    petco2hrf_demean : np.ndarray
        The demeaned petco2hrf regressor, central in time (not shifted).
    petco2hrf_lagged : np.ndarray
        The other shifted versions of the regressor.
    """
    # Upsample functional signal
    func_upsampled = resample_signal_freqs(func_avg, 1 / tr, freq)

    if skip_xcorr:
        LGR.info('Skipping Bulk Shift Computation')
        optshift = 0
    else:
        optshift = compute_bulk_shift(
            func_upsampled, petco2hrf, freq, outprefix, trial_len, n_trials, abs_xcorr
        )

    petco2hrf_shift = petco2hrf[optshift : optshift + func_upsampled.shape[0]]

    # Plot (shifted) regressor vs average ROI signal.
    plot_two_timeseries(
        petco2hrf_shift,
        func_upsampled,
        f'{outprefix}_petco2hrf_vs_avgroi.png',
        'Optimally shifted regressor',
        'Average ROI signal',
        freq,
        zscore=True,
    )

    petco2hrf_demean = export_regressor(
        petco2hrf_shift, func_avg.shape[-1], outprefix, 'petco2hrf_simple', ext
    )

    # Initialise the shifts first.
    petco2hrf_lagged = None
    if lagged_regression and lag_max:
        petco2hrf_lagged = create_fine_shift_regressors(
            petco2hrf,
            optshift,
            lag_max,
            lag_min,
            freq,
            func_avg.shape[-1],
            func_upsampled.shape[-1],
            outprefix,
            ext,
            legacy,
        )
    elif lagged_regression and not lag_max:
        LGR.warning(
            'Lagged regressors requested but maximum lag not provided. Skipping.'
        )
    else:
        LGR.info('Skipping generation of lagged regressors.')

    return petco2hrf_demean, petco2hrf_lagged


def select_lag_avoid_boundary(
    r_square_all,
    mask,
    starting_max,
    final_max,
    lag_min,
    lag_step,
    freq,
    expand_by=2,
):
    """
    Expand the maximum shift only in voxels with lag values at the boundaries.

    Parameters
    ----------
    r_square_all : ndarray
        R^2 values with shape (..., N_lags)
    mask : ndarray
        Boolean or 0/1 mask of valid voxels
    starting_max : float
        Initial maximum temporal shift to explore, expressed in seconds.
        Caution: this is not a pythonic range, but a real range, i.e. the upper limit is included.
        Default: None
    final_max : float
        Hard upper limit of the temporal area to explore, expressed in seconds.
        Caution: this is not a pythonic range, but a real range, i.e. the upper limit is included.
        Default: None
    lag_min : int, float, or None, optional
        Lower limit of the temporal area to explore, expressed in seconds.
        If set to None, and lag_max is not None and is positive, lag_min defaults to -lag_max (symmetric range).
        Default: None
    lag_step : float
        Lag step in seconds
    freq : float
        Sampling frequency (Hz)
    expand_by : float, optional
        Amount (seconds) to expand the lag window each iteration

    Returns
    -------
    final_lag_idx : ndarray
        Selected lag indices per voxel.
    """

    LGR.info(
        f'Identifying optimal shifts using a starting maximum of {starting_max}s and a final maximum of {final_max}s. If needed, the starting maximum will be increased in increments of {expand_by}s.'
    )

    final_lag_idx = np.full(mask.shape, 0, dtype=int) #stores final lag indices
    active = mask.astype(bool) #make all in mask voxels active
    current_max = starting_max

    while True:
        # Convert lag (seconds) -> index
        max_idx = int((current_max - lag_min) / lag_step)
        boundary_idx = max_idx - 1

        # Argmax only within current window
        lag_idx = np.argmax(r_square_all[..., : max_idx + 1], axis=-1)
        
        # Identify voxels that hit the boundary
        #at_boundary = (lag_idx >= boundary_idx) & active
        at_boundary = ((lag_idx >= boundary_idx) | (lag_idx == 0) | (lag_idx == 1)) & active

        # Accept non-boundary voxels and freeze them
        frozen = active & ~at_boundary
        final_lag_idx[frozen] = lag_idx[frozen]
    
        # Remove frozen voxels from active
        active[frozen] = 0  # now only boundary voxels remain active
    
        # Keep only boundary-stuck voxels
        active = active.astype(bool)  # ensure it's boolean

        if not np.any(active):
            break

        if current_max + expand_by > final_max:
            # Give up on remaining voxels
            final_lag_idx[active] = lag_idx[active]
            break

        current_max += expand_by

    return final_lag_idx

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
