#!/usr/bin/env python3
"""
Statistical module for phys2cvr.

Attributes
----------
LGR :
    Logger
"""

import logging
import os

import matplotlib.pyplot as plt
import numpy as np
from numpy.lib.stride_tricks import sliding_window_view as swv
from scipy.stats import zscore

from phys2cvr import io
from phys2cvr.io import FIGSIZE, SET_DPI
from phys2cvr.signal import resample_signal
from phys2cvr.stats import x_corr

R2MODEL = ["full", "partial", "intercept", "adj_full", "adj_partial", "adj_intercept"]

LGR = logging.getLogger(__name__)
LGR.setLevel(logging.INFO)


def create_legendre(degree, length):
    """
    Producesthe Legendre polynomials of order `degree`.

    Parameters
    ----------
    degree : int
        Highest order desired.
    length : int
        Number of samples of the polynomials.

    Returns
    -------
    legendre : np.ndarray
        A `degree`*`length` array with all the polynomials up to order `degree`
    """

    def _bonnet(d, x):
        if d == 0:
            return np.ones_like(x)
        elif d == 1:
            return x
        else:
            return (
                (2 * d - 1) * x * _bonnet(d - 1, x) - (d - 1) * _bonnet(d - 2, x)
            ) / d

    x = np.linspace(-1, 1, length)
    legendre = np.empty([length, degree + 1], dtype="float32")
    for n in range(degree + 1):
        legendre[:, n] = _bonnet(n, x)
    return legendre


def bulk_shift(
    func_upsampled,
    petco2hrf,
    freq,
    outname,
    trial_len=None,
    n_trials=None,
    abs_xcorr=False,
):
    """
    Compute bulk shift of regressor.
    """
    first_tp, last_tp = 0, -1

    if trial_len and n_trials:
        # If both are specified, disregard two extreme _trial from matching.
        LGR.info(f"Specified {n_trials} trials lasting {trial_len} seconds")
        if n_trials > 2:
            LGR.info("Ignoring first trial to improve first bulk shift estimation")
            first_tp = int(trial_len * freq)
        else:
            LGR.info("Using all trials for bulk shift estimation")
        if n_trials > 3:
            LGR.info("Ignoring last trial to improve first bulk shift estimation")
            last_tp = first_tp * (n_trials - 1)

    elif trial_len and not n_trials:
        LGR.warning(
            "The length of trial was specified, but the number of "
            "trials was not. Using all trials for bulk shift estimation"
        )
    elif not trial_len and n_trials:
        LGR.warning(
            "The number of trials was specified, but the length of "
            "trial was not. Using all trials for bulk shift estimation"
        )
    else:
        LGR.info("Using all trials for bulk shift estimation.")

    # Preparing breathhold and CO2 trace for Xcorr
    func_cut = func_upsampled[first_tp:last_tp]
    petco2hrf_cut = petco2hrf[:last_tp]
    _, optshift, xcorr = x_corr(
        func_cut, petco2hrf_cut, n_shifts=None, offset=None, abs_xcorr=abs_xcorr
    )
    LGR.info(f"Cross correlation estimated bulk shift at {optshift/freq} seconds")
    # Export estimated optimal shift in seconds
    with open(f"{outname}_optshift.1D", "w") as f:
        print(f"{(optshift/freq):.4f}", file=f)

    # Preparing time axis for plots
    time_axis = np.linspace(0, (len(xcorr) - 1) / freq, len(xcorr))

    # Export xcorr figure
    plt.figure(figsize=FIGSIZE, dpi=SET_DPI)
    plt.plot(time_axis, xcorr)
    plt.plot(time_axis[optshift], xcorr[optshift], "x")
    plt.legend(["Cross correlation value", "Optimal detected shift"])
    plt.title("Cross correlation and optimal shift")
    plt.tight_layout()
    plt.savefig(f"{outname}_optshift.png", dpi=SET_DPI)
    plt.close()

    # This shouldn't happen, but still check
    if optshift + func_upsampled.shape[0] > len(petco2hrf):
        raise Exception(
            f"The found optimal shift {optshift/freq} removes too many samples to "
            "continue. This error should not be possible."
        )

    return optshift


def fine_shift(
    petco2hrf,
    optshift,
    lag_max,
    freq,
    func_size,
    func_upsamp_size,
    outname,
    ext=".1D",
    legacy=False,
):
    """
    Compute fine shifts.
    """
    # func_size = func_avg.shape[-1]
    # func_upsamp_size = len_upd

    outprefix = os.path.join(
        os.path.split(outname)[0], "regr", os.path.split(outname)[1]
    )
    os.makedirs(os.path.join(os.path.split(outname)[0], "regr"), exist_ok=True)

    # Set num of fine shifts
    neg_shifts = int(lag_max * freq)
    pos_shifts = neg_shifts if legacy is True else (neg_shifts + 1)

    # Padding regressor right for shifts if not enough timepoints
    rpad = (
        func_upsamp_size - (petco2hrf.shape[0] - (optshift + pos_shifts))
        if petco2hrf.shape[0] - (optshift + pos_shifts) < func_upsamp_size
        else 0
    )
    # Padding regressor left for shifts and update optshift if less than neg_shifts.
    lpad = neg_shifts - optshift if (optshift - neg_shifts) < 0 else 0
    optshift = neg_shifts if (optshift - neg_shifts) < 0 else optshift

    petco2hrf = np.pad(petco2hrf, (int(lpad), int(rpad)), "mean")

    # Create sliding window view into petco2hrf, then select those
    negp, posp = optshift - neg_shifts, optshift + pos_shifts
    petco2hrf_lagged = swv(petco2hrf, func_upsamp_size)[negp:posp].copy()

    petco2hrf_lagged = io.export_regressor(
        petco2hrf_lagged, func_size, outprefix, "", ext
    )

    return petco2hrf_lagged


def create_physio_regressor(
    func_avg,
    petco2hrf,
    tr,
    freq,
    outname,
    lag_max=None,
    trial_len=None,
    n_trials=None,
    ext=".1D",
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
        Functional timeseries (1D)
    petco2hrf : np.ndarray
        Regressor of interest
    tr : str, int, or float
        TR of timeseries
    freq : str, int, or float
        Sample frequency of petco2hrf
    outname : list or path
        Path to output directory for regressors.
    lag_max : int or float, optional
        Limits (both positive and negative) of the temporal area to explore,
        expressed in seconds.
        Default: 9 (i.e. ±9 seconds)
    trial_len : str or int, optional
        Length of each single trial for tasks that have more than one
        (E.g. BreathHold, CO2 challenges, ...)
        Used to improve cross correlation estimation.
        Default: None
    n_trials : str or int, optional
        Number of trials in the task.
        Default: None
    ext : str, optional
        Extension to be used for the exported regressors.
    lagged_regression : bool, optional
        Estimate regressors for each possible lag of `petco2hrf`.
        If True, the maximum number of regressors will be `(freq*lag_max*2)+1`
    legacy : bool, optional
        If True, exclude the upper lag limit from the regression estimation.
        If True, the maximum number of regressors will be `(freq*lag_max*2)`
    abs_xcorr : bool, optional
        If True, the cross correlation will consider the maximum absolute
        correlation, i.e. if a negative correlation is higher than the highest
        positive, the negative correlation will be chosen instead.
    skip_xcorr : bool, optional
        If True, skip the cross correlation step.

    Returns
    -------
    petco2hrf_demean : np.ndarray
        The central, demeaned petco2hrf regressor.
    petco2hrf_lagged : np.ndarray
        The other shifted versions of the regressor.
    """
    # Upsample functional signal
    upsamp_tps = int(np.round(func_avg.shape[-1] * tr * freq))
    func_upsampled = resample_signal(func_avg, upsamp_tps)

    if not skip_xcorr:
        optshift = bulk_shift(
            func_upsampled, petco2hrf, freq, outname, trial_len, n_trials, abs_xcorr
        )
    else:
        optshift = 0

    petco2hrf_shift = petco2hrf[optshift : optshift + func_upsampled.shape[0]]

    # Exporting figures of shift
    plt.figure(figsize=FIGSIZE, dpi=SET_DPI)
    plt.plot(zscore(petco2hrf_shift), "-", zscore(func_upsampled), "-")
    plt.title("Optimally shifted regressor and average ROI signal")
    plt.legend(["Optimally shifted regressor", "Average ROI signal"])
    plt.tight_layout()
    plt.savefig(f"{outname}_petco2hrf_vs_avgroi.png", dpi=SET_DPI)
    plt.close()

    petco2hrf_demean = io.export_regressor(
        petco2hrf_shift, func_avg.shape[-1], outname, "petco2hrf_simple", ext
    )

    # Initialise the shifts first.
    petco2hrf_lagged = None

    if lagged_regression and lag_max:
        petco2hrf_lagged = fine_shift(
            petco2hrf,
            optshift,
            lag_max,
            freq,
            func_avg.shape[-1],
            func_upsampled[0],
            outname,
            ext,
            legacy,
        )
    elif lagged_regression and lag_max is None:
        LGR.warning(
            "The generation of lagged regressors was requested, "
            "but the maximum lag was not specified. Skipping "
            "lagged regressor generation."
        )
    else:
        LGR.info("Skipping lag regressors generation.")

    return petco2hrf_demean, petco2hrf_lagged


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
