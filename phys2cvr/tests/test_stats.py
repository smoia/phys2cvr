#!/usr/bin/env python3
"""Tests for io."""

import numpy as np
import pytest

from phys2cvr import stats


# ## Unit tests
def test_x_corr():
    func = np.array([1.0, -2.0, 3.0])
    co2 = np.array([0.0, 1.0, -2.0, 3.0, 4.0])

    # Manually compute expected correlations
    # sco2 windows: [0,1,-2], [1,-2,3], [-2,3,4]
    # zscores:
    # func z: [-1.224745, 0, 1.224745]
    # window z, then dot product / N
    val, idx, arr = stats.x_corr(func, co2)

    sco2 = np.array([[0.0, 1.0, -2.0], [1.0, -2.0, 3.0], [-2.0, 3.0, 4.0]])
    # Pre-computed expected values
    expected = np.corrcoef(func, sco2)[1:, 0]

    assert isinstance(val, float)
    assert isinstance(idx, (int | np.int64 | np.int32))
    assert isinstance(arr, np.ndarray)
    assert arr.shape[0] == co2.size - func.size + 1
    assert np.allclose(arr, expected, atol=1e-6)
    assert val == pytest.approx(expected.max())
    assert idx == expected.argmax()

    func = np.array([1.0, -2.0, 1.0])
    co2 = np.array([1.0, -1.0, 2.0, -1.0, 1.0])
    sco2 = np.array([[1, -1, 2], [-1, 2, -1], [2, -1, 1]])

    val, idx, arr = stats.x_corr(func, co2, abs_xcorr=True)
    expected_corr = np.corrcoef(func, sco2)[1:, 0]

    assert val >= 0
    assert arr.size == 3
    assert np.allclose(arr, expected_corr)
    assert val == pytest.approx(expected.max())
    assert idx == 1


def test_compute_fstat():
    """Test F-stat calculation against a known linear model setup."""
    np.random.seed(42)
    n_samples = 100
    n_voxels = 5

    # Target regressors (u_cols = 2)
    x_target = np.random.randn(n_samples, 2)
    # Denoising regressor (1 column)
    x_denoise = np.random.randn(n_samples, 1)

    # Full design matrix X
    x_full = np.hstack([x_denoise, x_target])

    # Construct synthetic data Y where target regressors have real signal
    true_betas = np.array([[0.5], [2.0], [-1.5]])  # 1 denoise + 2 target
    noise = np.random.randn(n_samples, n_voxels) * 0.5
    y = x_full @ true_betas + noise

    f_stat = stats.compute_fstat(y=y, x_full=x_full, u_cols=2)

    assert f_stat.shape == (n_voxels,)
    assert np.all(f_stat > 0)

    # Y is pure noise independent of x_target
    y = np.random.randn(n_samples, n_voxels)

    f_stat = stats.compute_fstat(y=y, x_full=x_full, u_cols=2)

    # Average F-statistic under null hypothesis should be around 1.0
    assert np.isclose(np.mean(f_stat), 1.0, atol=0.1)

    y = x_full @ np.array([[1.0], [2.0], [3.0]])
    y = np.repeat(y, n_voxels, axis=1)

    f_stat = stats.compute_fstat(y=y, x_full=x_full, u_cols=2)

    # Should safely return 0.0 (via nan_to_num) rather than NaN or Inf
    assert not np.any(np.isnan(f_stat))
    assert not np.any(np.isinf(f_stat))


def test_ols():
    Y = np.array([1, 2, 3, 4], float)
    X = np.column_stack([np.ones(4), np.arange(4)])
    betas, t, r = stats.ols(Y, X)
    assert betas.shape == (2, 1)
    assert t.shape == (2, 1)
    assert r.shape == (1,)
    assert np.isfinite(r).all()

    # Fit Y = b0*1 + b1*x  → exact solution
    X = np.column_stack([np.ones(5), np.arange(5)])
    Y = np.array([0.0, 1.0, 2.0, 3.0, 4.0])

    betas, tstats, r = stats.ols(Y, X)

    assert np.allclose(betas.flatten(), [0.0, 1.0])
    assert r == pytest.approx(1.0)

    # t-stats for perfect fit → huge (but finite after clipping)
    assert tstats[1, 0] > 1000

    X = np.column_stack([np.ones(4), np.arange(4)])
    Y = np.array([1.0, 2.0, 1.0, 2.0])  # variance around intercept only

    _, _, r = stats.ols(Y, X, r2model='intercept')

    # r^2 intercept = R^2 of X vs intercept-only model
    # Should be >0 because x explains alternating pattern
    assert r[0] > 0.2


def test_ols_residuals():
    X = np.column_stack([np.ones(3), np.arange(3)])
    Y = np.array([1.0, 3.0, 5.0])
    # Perfect Y = 1 + 2x
    res = stats.ols(Y, X, residuals=True)

    assert res.shape == (3, 1)
    assert np.allclose(res.flatten(), 0.0)


def test_regression():
    data = np.random.randn(4, 4, 4, 10)
    regr = np.random.randn(10, 1)
    b, t, r = stats.regression(data, regr)
    assert b.shape == data.shape[:-1]
    assert t.shape == data.shape[:-1]
    assert r.shape == data.shape[:-1]

    mask = np.zeros((4, 4, 4))
    mask[1, 1, 1] = 1
    b, t, r = stats.regression(data, regr, mask=mask)
    assert b[1, 1, 1] != 0
    assert b.sum() != 0
    assert r.shape == data.shape[:-1]

    data = np.zeros((1, 1, 1, 5))
    data[0, 0, 0] = np.array([1.0, 2.0, 3.0, 4.0, 5.0])

    regr = np.arange(5)[:, None]  # perfect predictor

    b, t, r = stats.regression(data, regr)

    assert b[0, 0, 0] == pytest.approx(1.0)
    assert r[0, 0, 0] == pytest.approx(1.0)
    assert t[0, 0, 0] > 1000

    data = np.zeros((2, 2, 2, 4))
    data[..., :] = np.array([1.0, 2.0, 3.0, 4.0])

    regr = np.arange(4)[:, None]

    mask = np.zeros((2, 2, 2))
    mask[0, 0, 0] = 1

    b, t, r = stats.regression(data, regr, mask=mask)

    assert b[0, 0, 0] == pytest.approx(1.0)
    assert r[0, 0, 0] == pytest.approx(1.0)
    assert b.sum() == b[0, 0, 0]
    assert r.sum() == r[0, 0, 0]


# ## Break tests


def test_break_x_corr():
    func = np.ones(5)
    co2 = np.ones(6)
    with pytest.raises(ValueError):
        stats.x_corr(func, co2, offset=5)

    with pytest.raises(NotImplementedError):
        stats.x_corr(func, co2, offset=-1)


def test_break_ols():
    Y = np.arange(5)
    X = np.column_stack([np.ones(5), np.arange(5)])
    with pytest.raises(ValueError):
        stats.ols(Y, X, r2model='invalid')

    Y = np.zeros((2, 2, 2))
    X = np.ones((2, 3))
    with pytest.raises(NotImplementedError):
        stats.ols(Y, X)


def test_break_regression():
    data = np.random.randn(4, 4, 4, 10)
    regr = np.random.randn(10, 1)
    denoise = np.random.randn(9, 2)
    with pytest.raises(ValueError):
        stats.regression(data, regr, denoise_mat=denoise)
