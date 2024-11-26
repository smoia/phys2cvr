#!/usr/bin/env python3
"""Tests for io."""

import os

import numpy as np

from phys2cvr import signal

# ## Unit tests


def test_spc():
    # Test case 1: test with a single time series
    ts = np.arange(5)
    expected_result = np.array([-1, -0.5, 0, 0.5, 1])
    assert (signal.spc(ts) == expected_result).all()

    # Test case 2: test with multiple time series
    ts2 = np.stack((ts, ts))
    expected_result = np.array([[-1, -0.5, 0, 0.5, 1], [-1, -0.5, 0, 0.5, 1]])
    assert (signal.spc(ts2) == expected_result).all()

    # Test case 3: test with input ts having mean 0
    ts = np.array([[-1, 0, 1], [-1, 0, 1]])
    expected_result = np.array([[-1, 0, 1], [-1, 0, 1]])
    assert (signal.spc(ts) == expected_result).all()

    # Test case 4: test with input ts having nan values
    ts = np.array([[1, 2, np.nan], [-1, np.nan, -3]])
    expected_result = np.array([[-0.33333333, 0.33333333, 0], [-0.5, 0, 0.5]])
    assert np.allclose(signal.spc(ts), expected_result)


def test_create_hrf():
    # Test with default frequency
    expected_length = 1281
    assert len(signal.create_hrf()) == expected_length

    # Test with custom frequency
    freq = 100
    hrf = signal.create_hrf(freq)
    expected_length = 3201
    assert len(hrf) == expected_length
    expected_peak_time = 500
    actual_peak_time = np.argmax(hrf)
    assert actual_peak_time == expected_peak_time

    # Test maximum value of HRF is 1
    assert hrf.max() == 1


def test_filter_signal():
    # Generate random data
    np.random.seed(0)
    data = np.random.rand(10, 20, 30)

    # Test with default filter parameters
    filtered_data = signal.filter_signal(data, tr=2)
    assert filtered_data.shape == (10, 20, 30)

    # Test with custom filter parameters
    filtered_data = signal.filter_signal(data, tr=2, lowcut=0.1, highcut=0.2)
    assert filtered_data.shape == (10, 20, 30)

    # Test that the filtered data is different from the original data
    assert not np.allclose(data, filtered_data)

    # Test that the filtered data has a lower variance than the original data
    assert np.var(filtered_data) < np.var(data)


# #!# Tolerance is too high, need to clean the files!
def test_compute_petco2hrf(testdir):
    # Create test data
    co2 = np.zeros(100)
    co2[20] = 1
    pidx = np.array([18, 19, 20, 21, 22])
    freq = 10
    outname = str(os.path.join(testdir, "test"))

    # Call function
    co2_conv = signal.compute_petco2hrf(co2, pidx, freq, outname)

    # Check output
    assert co2_conv.shape == (420,)
    assert np.isclose(co2_conv.max(), [0.788080])
    assert co2_conv.argmax() == 66

    # Check files were created
    assert os.path.isfile(os.path.join(testdir, "test_petco2.png"))
    assert os.path.isfile(os.path.join(testdir, "test_petco2.1D"))
    assert os.path.isfile(os.path.join(testdir, "test_petco2hrf.png"))
    assert os.path.isfile(os.path.join(testdir, "test_petco2hrf.1D"))


def test_resample_signal_freqs():
    # Create a sinusoidal signal at 10 Hz
    t = np.arange(0, 1, 0.001)
    signal1 = np.sin(2 * np.pi * 10 * t)

    # Resample the signal to 5 Hz
    signal2 = signal.resample_signal_freqs(signal1, freq1=10, freq2=5)

    # Calculate expected signal at 5 Hz
    t2 = np.arange(0, 1, 0.002)
    expected_signal2 = np.sin(2 * np.pi * 10 * t2)

    # Check if the resampled signal matches the expected signal
    assert np.allclose(signal2, expected_signal2)
