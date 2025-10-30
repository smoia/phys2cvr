#!/usr/bin/env python3
"""Tests for io."""

import os

import numpy as np

from phys2cvr import signal

# ## Unit tests


def test_spc():
    # Test case 1: test with a single time series
    ts = np.array([1, 2, 3, 4, 5])
    expected_result = np.array([-0.5, -0.25, 0, 0.25, 0.5])
    assert np.allclose(signal.spc(ts), expected_result)

    # Test case 2: test with multiple time series
    ts = np.array([[1, 2, 3, 4, 5], [5, 4, 3, 2, 1]])
    expected_result = np.array(
        [[-0.5, -0.25, 0, 0.25, 0.5], [0.5, 0.25, 0, -0.25, -0.5]]
    )
    assert np.allclose(signal.spc(ts), expected_result)

    # Test case 3: test with input ts having mean 0
    ts = np.array([[1, 2, 3], [-1, -2, -3]])
    expected_result = np.array([[1, 2, 3], [-1, -2, -3]])
    assert np.allclose(signal.spc(ts), expected_result)

    # Test case 4: test with input ts having nan values
    ts = np.array([[1, 2, np.nan], [-1, np.nan, -3]])
    expected_result = np.array([[1 / 3, 2 / 3, 0], [-1 / 2, 0, -1]])
    assert np.allclose(signal.spc(ts), expected_result)


def test_create_hrf():
    # Test case 1: test with default frequency
    expected_length = 528
    assert len(signal.create_hrf()) == expected_length

    # Test case 2: test with custom frequency
    freq = 100
    expected_length = 1057
    assert len(signal.create_hrf(freq)) == expected_length

    # Test case 3: test that the HRF is non-negative
    assert np.all(signal.create_hrf() >= 0)

    # Test case 4: test that the maximum value of HRF is 1
    assert signal.create_hrf().max() == 1

    # Test case 5: test that the sum of HRF is equal to 1
    assert np.isclose(signal.create_hrf().sum(), 1)

    # Test with custom frequency of 5
    freq = 5
    hrf = signal.create_hrf(freq)

    # Test that the length of the HRF is correct
    expected_length = 132
    assert len(hrf) == expected_length

    # Test that the HRF is non-negative
    assert np.all(hrf >= 0)

    # Test that the first peak of the HRF is at the correct time point
    expected_peak_time = 5
    actual_peak_time = np.argmax(hrf)
    assert actual_peak_time == expected_peak_time


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
def test_convolve_petco2(testdir):
    # Create test data
    co2 = np.zeros(100)
    co2[20] = 1
    co2[40] = 1
    pidx = np.array([20, 40])
    freq = 10
    outname = str(os.path.join(testdir, 'test'))

    # Call function
    co2_conv = signal.convolve_petco2(co2, pidx, freq, outname)

    # Check output
    assert co2_conv.shape == co2.shape
    assert np.allclose(co2_conv[:19], 0)
    assert np.allclose(co2_conv[20], 1, rtol=1e-6)
    assert np.allclose(co2_conv[21:39], 0, rtol=1e-6)
    assert np.allclose(co2_conv[40], 1, rtol=1e-6)
    assert np.allclose(co2_conv[41:], 0)

    # Check files were created
    assert (os.path.join(testdir, 'test_petco2.png')).is_file()
    assert (os.path.join(testdir, 'test_petco2.1D')).is_file()
    assert (os.path.join(testdir, 'test_petco2hrf.png')).is_file()
    assert (os.path.join(testdir, 'test_petco2hrf.1D')).is_file()


def test_resample_signal():
    # Create a sinusoidal signal at 10 Hz
    t = np.linspace(0, 1, 1000)
    signal1 = np.sin(2 * np.pi * 10 * t)

    # Resample the signal to 5 Hz
    signal2 = signal.resample_signal(signal1, freq1=10, freq2=5)

    # Calculate expected signal at 5 Hz
    t2 = np.linspace(0, 1, 500)
    expected_signal2 = np.sin(2 * np.pi * 10 * t2)

    # Check if the resampled signal matches the expected signal
    assert np.allclose(signal2, expected_signal2, decimal=6)
