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

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as spint
import scipy.stats as sct

from phys2cvr import io
from phys2cvr.io import SET_DPI, FIGSIZE


LGR = logging.getLogger(__name__)
LGR.setLevel(logging.INFO)


def x_corr(func, co2, lastrep, firstrep=0, offset=0):
    """
    Cross correlation between `func` and `co2`.
    
    Parameters
    ----------
    func : np.ndarray
        Timeseries, must be SHORTER (or of equal length) than `co2`
    co2 : np.ndarray
        Second timeseries, can be LONGER than `func` 
    lastrep : int
        Last index to take into account in `func`
    firstrep : int, optional
        First index totake into account in `func`
    offset : int, optional
        Optional amount of offset desired for `func`
    
    Returns
    -------
    float :
        Highest correlation
    int :
        Index of higher correlation
    xcorr : np.ndarray
        Full Xcorr
    
    Raises
    ------
    ValueError
        If `offset` is too high for the functional file
    """
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


def get_regr(func_avg, petco2hrf, tr, freq, outname, lag_max=None,
             trial_len=None, n_trials=None, ext='.1D', lagged_regression=True):
    """
    Create regressor(s) of interest for nifti GLM.
    
    Parameters
    ----------
    func_avg : np.ndarray
        Functional timeseries
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
        If True, the maximum number of regressors will be `(freq*lag_max*2)-1`
    
    Returns
    -------
    petco2hrf_demean : np.ndarray
        The central, demeaned petco2hrf regressor.
    petco2hrf_shifts : np.ndarray
        The other shifted versions of the regressor.
    """
    # Setting up some variables
    first_tp, last_tp = 0, -1

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
    regr_t = np.arange(0, ((func_len-1) * tr + 1/freq), 1/freq)
    func_t = np.linspace(0, (func_len - 1) * tr, func_len)
    f = spint.interp1d(func_t, func_avg, fill_value='extrapolate')
    func_upsampled = f(regr_t)
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
    plt.close()

    plt.figure(figsize=FIGSIZE, dpi=SET_DPI)
    plt.plot(sct.zscore(petco2hrf_shift), '-', sct.zscore(func_upsampled), '-')
    plt.title('GM and shift')
    plt.savefig(f'{outname}_petco2hrf.png', dpi=SET_DPI)
    plt.close()

    petco2hrf_demean = io.export_regressor(regr_t, petco2hrf_shift, func_t, outname,
                                           'petco2hrf', ext)

    # Initialise the shifts first.
    petco2hrf_shift = None
    if lagged_regression and lag_max:
        outprefix = os.path.join(os.path.split(outname)[0], 'regr', os.path.split(outname)[1])
        os.makedirs(os.path.join(os.path.split(outname)[0], 'regr'), exist_ok=True)

        # Set num of fine shifts: 9 seconds is a bit more than physiologically feasible
        nrep = int(lag_max * freq)

        petco2hrf_shifts = np.empty([func_len, nrep*2])

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
            petco2hrf_lagged = petco2hrf_padded[optshift+lpad-i:optshift+lpad-i+len_upd]
            petco2hrf_shifts[:, i] = io.export_regressor(regr_t, petco2hrf_lagged,
                                                         func_t, outprefix,
                                                         f'{(i + nrep):04g}', ext)

    elif not lag_max:
        LGR.warning('The generation of lagged regressors was requested, '
                    'but the maximum lag was not specified. Skipping '
                    'lagged regressor generation.')

    return petco2hrf_demean, petco2hrf_shifts


def get_legendre(degree, length):
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
        if(d == 0):
            return np.ones_like(x)
        elif(d == 1):
            return x
        else:
            return ((2*d-1)*x*_bonnet(d-1, x)-(d-1)*_bonnet(d-2, x))/d

    x = np.linspace(-1, 1, length)
    legendre = np.empty([length, degree+1])
    for n in range(degree+1):
        legendre[:, n] = _bonnet(n, x)
    return legendre


def regression(data, mask, regr, mat_conf):
    """
    Estimate regression parameters.
    
    Parameters
    ----------
    data : np.ndarray
        Dependent variable of the model (i.e. Y), as a 4D volume.
    mask : np.ndarray
        A 3D mask to reduce the number of voxels to run the regression for.
    regr : np.ndarray
        Regressor of interest
    mat_conf : np.ndarray
        Confounding effects (regressors)
    
    Returns
    -------
    bout : np.ndarray
        Beta map
    tout : np.ndarray
        T-stats map
    rout : np.ndarray
        R^2 map
    
    Raises
    ------
    Exception
        If `mat_conf` and `regr` do not have at least one common dimension.
    """
    # Obtain data in 2d
    data_2d = data[mask]
    # Check that regr has "two" dimensions
    if len(regr.shape) < 2:
        regr = regr[..., np.newaxis]
    if regr.shape[0] != mat_conf.shape[0]:
        regr = regr.T
        if regr.shape[0] != mat_conf.shape[0]:
            raise ValueError('The provided confounding matrix does not match '
                             'the dimensionality of the PetCO2hrf regressor!')
    # Stack mat and solve least square with it demeaned
    mat = np.hstack([mat_conf, regr])
    betas = np.linalg.lstsq(mat-mat.mean(axis=0),
                            data_2d.T, rcond=None)[0]

    # compute t-values of betas (estimates)
    # first compute number of degrees of freedom
    dofs = data_2d.shape[1] - mat.shape[1]
    # compute sigma:
    # RSS = sum{[mdata - (X * betas)]^2}
    # sigma = RSS / Degrees_of_Freedom
    RSS = np.sum(np.power(data_2d.T - np.dot(mat, betas), 2), axis=0)
    sigma = (RSS / dofs)
    sigma = sigma[..., np.newaxis]
    # Copmute std of betas:
    # C = (mat^T * mat)_ii^(-1)
    # std(betas) = sqrt(sigma * C)
    C = np.diag(np.linalg.pinv(np.dot(mat.T, mat)))
    C = C[..., np.newaxis]
    std_betas = np.sqrt(np.dot(sigma, C.T))
    tstats = betas / std_betas.T

    # Compute R^2 coefficient of multiple determination!
    # R^2 = 1 - RSS/TSS, where TSS (total sum of square) is variance*samples
    # #!# Check the axis here
    TSS = data_2d.var(axis=1) * data_2d.shape[1]
    r_square = np.ones(data_2d.shape[0]) - (RSS / TSS)

    # Assign betas, Rsquare and tstats to new volume
    bout = mask * 1.0
    tout = mask * 1.0
    rout = mask * 1.0
    bout[mask] = betas[-1, :]
    tout[mask] = tstats[-1, :]
    rout[mask] = r_square
    return bout, tout, rout


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
