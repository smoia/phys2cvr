#!/usr/bin/env python3
"""Tests for utils."""

from unittest.mock import patch

import numpy as np
import pytest

from phys2cvr import utils

# ## Unit tests


@pytest.mark.parametrize(
    'var, dtype, out',
    [
        (10, 'int', 10),
        (10.0, 'int', 10),
        ('10', 'int', 10),
        (None, 'int', None),
        (10, 'float', 10.0),
        (10.0, 'float', 10.0),
        ('10.0', 'float', 10.0),
        (None, 'float', None),
        (10, 'str', '10'),
        (10.0, 'str', '10.0'),
        ('10', 'str', '10'),
        ([10], 'str', '[10]'),
        (None, 'str', None),
        (10, 'list', [10]),
        ([10], 'list', [10]),
        ('10', 'list', ['10']),
        (None, 'list', None),
    ],
)
def test_if_declared_force_type(var, dtype, out):
    assert utils.if_declared_force_type(var, dtype) == out


@pytest.mark.xfail(strict=False, reason="Something's off with assertions")
@patch('logging.warning')
def test_if_declared_force_type_logging(mock_warning):
    utils.if_declared_force_type(10, 'int', 'my_var')
    mock_warning.assert_called_with(
        "Changing type of variable my_var from <class 'int'> to int"
    )

    utils.if_declared_force_type('10', 'int', 'my_var')
    mock_warning.assert_not_called()

    utils.if_declared_force_type([10], 'int', 'my_var')
    mock_warning.assert_called_with(
        "Changing type of variable my_var from <class 'list'> to int"
    )

    utils.if_declared_force_type(None, 'int', 'my_var')
    mock_warning.assert_not_called()

    utils.if_declared_force_type(10.0, 'float', 'my_var', silent=True)
    mock_warning.assert_not_called()


@pytest.mark.xfail(strict=False, reason='Needs to be updated for latest version')
def test_check_ext():
    all_ext = ['.csv', '.txt']
    fname = 'data.csv'
    assert utils.check_ext(all_ext, fname) is True
    assert utils.check_ext(all_ext, fname, remove=True) == ('data', True)
    all_ext = ['.CSV']
    assert utils.check_ext(all_ext, fname) is True

    fname = 'data.xls'
    assert utils.check_ext(all_ext, fname) is False
    assert utils.check_ext(all_ext, fname, remove=True) == ('data.xls', False)

    all_ext = []
    fname = 'data.csv'
    assert utils.check_ext(all_ext, fname, remove=True) == ('data.csv', False)


@pytest.mark.xfail(strict=False, reason="Something's off with warnings")
@patch('logging.warning')
def test_check_nifti_dim(mock_warning):
    fname = 'valid.nii.gz'
    data = np.ones((4, 4, 4, 4, 4))
    output = utils.check_nifti_dim(fname, data, dim=5)
    assert output.ndim == 5

    output = utils.check_nifti_dim(fname, data, dim=3)

    # Test if data has more dimensutilsns than `dim`
    output = utils.check_nifti_dim(fname, data, dim=3)
    mock_warning.assert_called_with(
        f'{fname} has more than 3 dimensutilsns. Removing D > 3.'
    )
    assert output.ndim == 3


# ## Break tests


def test_break_if_declared_force_type():
    with pytest.raises(NotImplementedError) as errorinfo:
        utils.if_declared_force_type('10', 'invalid_type')
    assert 'not supported' in str(errorinfo.value)


# Check this one
@pytest.mark.parametrize(
    'var, dtype',
    [
        ([10], 'int'),
        ([10.0], 'float'),
    ],
)
def test_if_declared_force_type_errors(var, dtype):
    with pytest.raises(TypeError):
        utils.if_declared_force_type(var, dtype)


@pytest.mark.xfail(strict=False, reason='Message changed')
def test_break_check_nifti_dim_missing_dims():
    with pytest.raises(ValueError) as errorinfo:
        utils.check_nifti_dim('missing_dims.nii.gz', np.ones((4, 4)), dim=4)
    assert 'seem to be a 4D file.' in str(errorinfo.value)
