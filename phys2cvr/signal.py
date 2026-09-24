#!/usr/bin/env python3
"""
Signal analysis module for phys2cvr.

Attributes
----------
LGR :
    Logger
"""

import logging
from pathlib import Path, PosixPath

import numpy as np
import scipy.interpolate as spint
import scipy.stats as sct
from scipy.signal import butter, sosfiltfilt

LGR = logging.getLogger(__name__)
LGR.setLevel(logging.INFO)


def spc(ts):
    """
    Compute signal percentage change over time series (ts).

    Timeseries are divided by the mean.
    Timeseries that have a mean of 0 are divided by 1 instead.

    Parameters
    ----------
    ts : numpy.ndarray
        A timeseries or set of timeseries; last dimension is assumed to be time.

    Returns
    -------
    numpy.ndarray
        Signal percentage change version of the input ts.
    """
    m = ts.mean(axis=-1)[..., np.newaxis]
    md = np.where(m == 0, 1.0, m)
    ts = (ts - m) / md
    np.nan_to_num(ts, copy=False, nan=0.0)
    return ts


def create_hrf(freq=40):
    """
    Create a canonical haemodynamic response function sampled at the given freq.

    Parameters
    ----------
    freq : float
        Sampling frequency used to resample the haemodynamic response function.

    Returns
    -------
    hrf : np.ndarray
        Haemodynamic response function.
    """
    # Create HRF
    RT = 1 / freq
    fMRI_T = 16
    p = [6, 16, 1, 1, 6, 0, 32]

    # Modelled hemodynamic response function - {mixture of Gammas}
    dt = RT / fMRI_T
    u = np.arange(0, p[6] / dt + 1, 1) - p[5] / dt

    a1 = p[0] / p[2]
    b1 = 1 / p[2]
    a2 = p[1] / p[3]
    b2 = 1 / p[3]

    hrf = (
        sct.gamma.pdf(u * dt, a1, scale=b1) - sct.gamma.pdf(u * dt, a2, scale=b2) / p[4]
    ) / dt

    time_axis = np.arange(0, int(p[6] / RT + 1), 1) * fMRI_T
    hrf = hrf[time_axis]

    min_hrf = 1e-9 * min(hrf[hrf > 10 * np.finfo(float).eps])
    if min_hrf < 10 * np.finfo(float).eps:
        min_hrf = 10 * np.finfo(float).eps

    hrf[hrf == 0] = min_hrf
    hrf = hrf / max(hrf)
    return hrf


def filter_signal(data, tr, lowcut=0.02, highcut=0.04, order=9, axis=-1):
    """
    Create a bandpass filter within lower and upper threshold, then filter.

    Parameters
    ----------
    data : np.ndarray
        Data to filter.
    tr : float
        TR of functional files.
    lowcut : float
        Low frequency threshold.
    highcut : float
        High frequency threshold.
    order : int
        Butterworth filter order.
    axis : int
        The axis along which the filter is applied.

    Returns
    -------
    filt_data : np.ndarray
        Bandpass-filtered data.
    """
    nyq = 0.5 / tr
    low = lowcut / nyq
    high = highcut / nyq
    a, b = butter(int(order), [low, high], btype='band', output='sos')
    return sosfiltfilt(a, b, data, axis=axis)


def endtidal_interpolation(signal, peaks, axis=-1):
    """
    Compute end-tidal interpolation of signal.

    Parameters
    ----------
    signal : np.ndarray-like or list
        The signal that needs to be interpolated
    peaks : 1d array-like or list
        The index of the peaks (or point of interest) to interpolate at.
    axis : int, optional
        The axis of signal to interpolate through. Default is last axis.

    Returns
    -------
    IntETSignal
        The end tidal interpolation of the signal

    """
    peaks = np.unique(peaks)
    nx = np.arange(signal.shape[axis])
    f = spint.interp1d(peaks, signal[peaks], fill_value='extrapolate', axis=axis)
    return f(nx)


def convolve_signal(signal, freq, response_function='hrf', mode='full'):
    """
    Convolve provided signal with provided response function.

    Parameters
    ----------
    signal : 1D np.ndarray or list
        The signal to convolve
    freq : int or float
        The sampling frequency of `signal`
    response_function : {`hrf`, `rrf`, `crf`, `icrf`}, None, str, path, or 1D array-like , optional
        Name of the response function to be used in the convolution of the regressor of
        interest or path to a 1D file containing one. Default is `hrf`.
        For `rrf` and `crf`, `phys2denoise` must be installed (see extra installs).
            - `None` (or a string version of none) will skip the convolution
            - `hrf` will use an internally generated canonical HRF
            - `rrf` will use `phys2denoise`'s RRF
            - `crf` will use `phys2denoise`'s CRF
            - `icrf` will use `phys2denoise`'s iCRF
            - Alternatively, specify a valid path file containing a custom one.
            - Alternatively, specify a custom one directly as np.ndarray or list.
    mode : {'full', 'valid', 'same'} str, optional
        Convolution mode, see numpy.convolve.

    Returns
    -------
    colvolved_signal : 1D np.ndarray
        The convolved version of `signal`.

    Raises
    ------
    ImportError
        If phys2denoise is not installed and RRF, CRF, or iCRF are called.
    NotImplementedError
        If signal has more than 1 dimension.
        If response function is not of a supported type or was not found in system.
    OSError
        If response function is a path but it does not exist.
    ValueError
        If the response function is a list or ndarray but is not numeric.
    """
    signal = signal.squeeze()
    if signal.ndim > 2:
        raise NotImplementedError(
            'Convolution on data that has more than 1 dimension is not supported yet.'
        )
    # Check if RF is a string representing a path or a Nonetype
    if isinstance(response_function, (str, Path)):
        rf_path = Path(response_function)
        if rf_path.exists():
            response_function = rf_path
            LGR.debug(f'{response_function} is a valid file')
        elif isinstance(response_function, str):
            rf_str = response_function.lower()
            response_function = None if rf_str == 'none' else rf_str

    if type(response_function) in [list, np.ndarray]:
        # Check if RF is or should be a 1darray-like
        convolving_function = (
            np.array(response_function)
            if type(response_function) is list
            else response_function
        )
        if not np.issubdtype(convolving_function.dtype, np.number):
            raise ValueError(
                'Provided function is not a numeric ndarray-like variable.'
            )
    elif response_function is None:
        LGR.info('Skipping convolution with response function')
        return signal
    elif response_function == 'hrf':
        convolving_function = create_hrf(freq)
    elif response_function == 'rrf':
        try:
            from phys2denoise.metrics.responses import rrf
        except ImportError:
            raise ImportError(
                'phys2denoise is required for the use of RRF response functions. '
                'Please see install instructions.'
            )
        convolving_function = rrf(1 / freq)
    elif response_function == 'crf':
        try:
            from phys2denoise.metrics.responses import crf
        except ImportError:
            raise ImportError(
                'phys2denoise is required for the use of CRF response functions. '
                'Please see install instructions.'
            )
        convolving_function = crf(1 / freq)
    elif response_function == 'icrf':
        try:
            from phys2denoise.metrics.responses import icrf
        except ImportError:
            raise ImportError(
                'phys2denoise is required for the use of iCRF response functions. '
                'Please see install instructions.'
            )
        convolving_function = icrf(1 / freq)
    elif type(response_function) is PosixPath:
        if not Path(response_function).exists():
            raise OSError(f'{response_function} not found in system.')
        # Local import to avoid circularity
        from .io import load_array  # noqa: ABS101

        convolving_function = load_array(response_function)
    else:
        raise NotImplementedError(
            f'Response function "{response_function}" is not supported yet or file was '
            'not found.'
        )

    convolved_signal = np.convolve(signal, convolving_function, mode=mode)
    convolved_signal = np.interp(
        convolved_signal,
        (convolved_signal.min(), convolved_signal.max()),
        (signal.min(), signal.max()),
    )

    return convolved_signal


def resample_signal_samples(ts, samples, axis=-1):
    """
    Resample timeseries based on desired number of samples.

    Parameters
    ----------
    ts : np.ndarray
        The timeseries to be resampled
    samples : int
        Desired number of samples.
    axis : int
        The axis (dimension) over which the interpolation should be applied - by default it's
        -1, i.e., the last dimension.

    Returns
    -------
    numpy.ndarray
        Resampled timeseries.
    """
    len_tp = ts.shape[axis]
    regr_t = np.linspace(0, len_tp - 1, samples)
    time_t = np.linspace(0, len_tp - 1, len_tp)
    f = spint.interp1d(time_t, ts, fill_value='extrapolate', axis=axis)
    return f(regr_t)


def resample_signal_freqs(ts, freq1, freq2, axis=-1):
    """
    Resample timeseries based on current and desired frequency.

    Parameters
    ----------
    ts : np.ndarray
        The timeseries to be resampled.
    freq1 : float
        Current frequency.
    freq2 : float
        Desired frequency.
    axis : int
        The axis (dimension) over which the interpolation should happen - by default it's
        -1, i.e. the last dimension.

    Returns
    -------
    numpy.ndarray
        Resampled timeseries.
    """
    len_tp = ts.shape[axis]
    len_s = (len_tp - 1) / freq1
    regr_t = np.linspace(0, len_s, int(len_s * freq2) + 1)
    time_t = np.linspace(0, len_s, len_tp)
    f = spint.interp1d(time_t, ts, fill_value='extrapolate', axis=axis)
    return f(regr_t)


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
