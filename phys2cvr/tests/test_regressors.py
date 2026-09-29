#!/usr/bin/env python3
"""Tests for io."""

import os
from unittest.mock import patch

import numpy as np
import pytest

from phys2cvr import regressors

# ## Unit tests


def test_create_legendre():
    L = regressors.create_legendre(3, 10)
    assert L.shape == (10, 4)

    L = regressors.create_legendre(1, 5)
    assert np.allclose(L[:, 1], np.linspace(-1, 1, 5))

    L = regressors.create_legendre(0, 7)
    assert L.shape == (7, 1)
    assert np.all(L[:, 0] == 1)


def test_fourier_basis_shape_and_dtype():
    """Test output shape and data type of generated Fourier basis."""
    trial_len = 60.0
    func_len = 300.0
    lag_range = 10.0
    order = 3
    sample_interval = 0.01

    fourier_mat = regressors.fourier_basis(
        trial_len=trial_len,
        func_len=func_len,
        order=order,
        sample_interval=sample_interval,
    )

    # Expected timepoints: func_len / sample_interval = 300 / 0.01 = 30000
    expected_samples = int(np.round(func_len / sample_interval))
    expected_columns = 2 * order  # sin + cos for each harmonic order

    assert fourier_mat.shape == (expected_samples, expected_columns)
    assert fourier_mat.dtype == np.float32

    fourier_mat = regressors.fourier_basis(
        trial_len=trial_len,
        func_len=func_len,
        sample_interval=sample_interval,
        lag_range=lag_range,
    )

    expected_samples = int(np.round((func_len + lag_range) / sample_interval))
    assert fourier_mat.shape[0] == expected_samples

    total_samples = int(np.round(func_len / sample_interval))
    t = np.linspace(0, func_len, total_samples, endpoint=False, dtype=np.float32)
    w0 = 2 * np.pi / trial_len

    expected_sin = np.sin(w0 * t)
    expected_cos = np.cos(w0 * t)

    np.testing.assert_allclose(fourier_mat[:, 0], expected_sin, atol=1e-5)
    np.testing.assert_allclose(fourier_mat[:, 1], expected_cos, atol=1e-5)


@patch(
    'phys2cvr.regressors.endtidal_interpolation', return_value=np.array([1, 2, 3, 4])
)
@patch('phys2cvr.regressors.plot_timeseries')
@patch('phys2cvr.regressors.convolve_signal', return_value=np.array([9, 8, 7, 6]))
def test_compute_petco2hrf(mock_conv, mock_plot, mock_endtidal, testdir):
    out = os.path.join(testdir, 'Lune')
    co2 = np.array([1, 2, 3, 4])
    pidx = np.array([0, 2])
    r = regressors.compute_petco2hrf(co2, pidx, 1.0, out)
    assert np.all(r == np.array([9, 8, 7, 6]))
    os.remove(f'{out}_petco2.1D')


def test_compute_petco2hrf_skip(testdir):
    out = os.path.join(testdir, 'Catherine')
    r = regressors.compute_petco2hrf(
        np.array([1, 2, 3]),
        np.array([0]),
        1.0,
        out,
        comp_endtidal=False,
        response_function=None,
    )
    assert np.allclose(r, np.array([-1, 0, 1]))  # returns petco2 demeaned


@patch('phys2cvr.regressors.x_corr', return_value=(None, 3, np.array([1, 2, 3])))
@patch('phys2cvr.regressors.plot_xcorr')
def test_compute_bulk_shift(mock_plotx, mock_xcorr, testdir):
    out = os.path.join(testdir, 'Monoco')
    func = np.arange(10.0)
    pet = np.arange(20.0)
    s = regressors.compute_bulk_shift(func, pet, 1.0, out)
    assert s == 3
    assert os.path.exists(f'{out}_optshift.1D')
    os.remove(f'{out}_optshift.1D')


@patch('phys2cvr.regressors.export_regressor', return_value=np.zeros((5, 5)))
def test_create_fine_shift_regressors(mock_export, testdir):
    out = os.path.join(testdir, 'Esquie')
    pet = np.arange(50.0)
    r = regressors.create_fine_shift_regressors(pet, 3, 2, 1.0, 10, 20, out)
    assert r.shape == (5, 5)


@patch('phys2cvr.regressors.export_regressor', return_value=np.zeros((3, 3)))
def test_create_fine_shift_regressors_padding(mock_export, testdir):
    out = os.path.join(testdir, 'Maelle')
    pet = np.arange(10.0)
    r = regressors.create_fine_shift_regressors(pet, 9, 4, 1.0, 5, 8, out)
    assert r.shape == (3, 3)


@patch('phys2cvr.regressors.resample_signal_freqs', return_value=np.arange(20.0))
@patch('phys2cvr.regressors.plot_timeseries')
@patch('phys2cvr.regressors.export_regressor', return_value=np.arange(10.0))
def test_create_physio_regressor(mock_exp, mock_plot, mock_bs, testdir):
    out = os.path.join(testdir, 'Sciel')
    func = np.arange(10.0)
    pet = np.arange(30.0)
    d, lags = regressors.create_physio_regressor(
        func, pet, 1.0, 1.0, out, lag_max=2, skip_xcorr=True
    )
    assert d.shape == (10,)
    assert lags is not None


@patch('phys2cvr.regressors.resample_signal_freqs', return_value=np.arange(20.0))
@patch('phys2cvr.regressors.export_regressor', return_value=np.arange(10.0))
def test_create_physio_regressor_no_lag_max(mock_exp, mock_bs, testdir):
    out = os.path.join(testdir, 'Francois')
    func = np.arange(10.0)
    pet = np.arange(30.0)
    d, lags = regressors.create_physio_regressor(
        func, pet, 1.0, 1.0, out, lag_max=None, skip_xcorr=True
    )
    assert d.shape == (10,)
    assert lags is None
    os.remove(f'{out}_petco2hrf_vs_avgroi.png')


# ## Break tests


def test_break_compute_petco2hrf():
    with pytest.raises(NotImplementedError) as errorinfo:
        regressors.compute_petco2hrf(np.zeros((2, 2)), np.array([0]), 1.0, 'x')
    assert 'Arrays with more' in str(errorinfo.value)


@pytest.mark.parametrize(
    'order,trial_len,sample_interval',
    [
        (0, 60.0, 0.01),
        (-1, 60.0, 0.01),
        (3, 0.0, 0.01),
        (3, -10.0, 0.01),
        (3, 60.0, 0.0),
        (3, 60.0, -0.01),
    ],
)
def test_break_fourier_basis(order, trial_len, sample_interval):
    """Test that ValueError is raised for non-positive input parameters."""
    with pytest.raises(ValueError):
        regressors.fourier_basis(
            trial_len=trial_len,
            func_len=100.0,
            order=order,
            sample_interval=sample_interval,
        )
