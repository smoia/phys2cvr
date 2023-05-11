#!/usr/bin/env python3
"""Tests for io."""

import os

import numpy as np

from phys2cvr import stats

# ## Unit tests


def test_x_corr():
    # Define two signals, cosine and sine
    t = np.linspace(0, 1, 1000)
    func = np.cos(2 * np.pi * 5 * t)
    co2 = np.sin(2 * np.pi * 5 * t) + 0.1 * np.random.randn(len(t))

    # Calculate cross correlation
    corr, index, xcorr = stats.x_corr(func, co2, lastrep=100, offset=50)

    # Check that the highest correlation is close to 1
    assert corr > 0.95

    # Check that the index of the highest correlation is close to the true offset
    assert np.abs(index - 50) < 5

    # Check that the length of the cross correlation is correct
    assert len(xcorr) == 50

    # Check that the cross correlation is symmetric
    assert np.allclose(xcorr, xcorr[::-1])

    # Generate signals
    t = np.linspace(0, 10, 1000)
    func = np.sin(t)
    co2 = np.tile(np.cos(t), 2)
    co2 = np.tile(co2, 2)

    # Compute x_corr
    max_corr, max_index, xcorr = stats.x_corr(func, co2, abs_xcorr=True)

    # Compute expected result
    expected_max_corr = np.abs(np.corrcoef(func, co2[: len(func)].T)[1, 0])
    expected_max_index = np.argmax(np.abs(xcorr))

    # Check results
    assert np.isclose(max_corr, expected_max_corr, rtol=1e-6)
    assert max_index == expected_max_index
