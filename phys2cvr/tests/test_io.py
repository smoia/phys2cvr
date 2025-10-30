#!/usr/bin/env python3
"""Tests for io."""

import os
from unittest.mock import patch

import nibabel as nib
import numpy as np
import pytest

from phys2cvr import io

# ## Unit tests


@pytest.mark.parametrize(
    'var, dtype, out',
    [
        (10, 'int', 10),
        (10.0, 'int', 10),
        ('10', 'int', 10),
        ([10], 'int', [10]),
        (None, 'int', None),
        (10, 'float', 10.0),
        (10.0, 'float', 10.0),
        ('10.0', 'float', 10.0),
        ([10.0], 'float', [10.0]),
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
    assert io.if_declared_force_type(var, dtype) == out


@patch('logging.warning')
@pytest.mark.parametrize(
    'var, dtype, varname, expected_warning',
    [
        (
            10,
            'int',
            'my_var',
            "Changing type of variable my_var from <class 'int'> to int",
        ),
        ('10', 'int', 'my_var', ''),
        (
            [10],
            'int',
            'my_var',
            "Changing type of variable my_var from <class 'list'> to int",
        ),
        (None, 'int', 'my_var', ''),
        (10.0, 'float', 'my_var', ''),
    ],
)
def test_if_declared_force_type_logging(mock_warning):
    io.if_declared_force_type(10, 'int', 'my_var')
    mock_warning.assert_called_with(
        "Changing type of variable my_var from <class 'int'> to int"
    )

    io.if_declared_force_type('10', 'int', 'my_var')
    mock_warning.assert_not_called()

    io.if_declared_force_type([10], 'int', 'my_var')
    mock_warning.assert_called_with(
        "Changing type of variable my_var from <class 'list'> to int"
    )

    io.if_declared_force_type(None, 'int', 'my_var')
    mock_warning.assert_not_called()

    io.if_declared_force_type(10.0, 'float', 'my_var', silent=True)
    mock_warning.assert_not_called()


def test_check_ext():
    all_ext = ['.csv', '.txt']
    fname = 'data.csv'
    assert io.check_ext(all_ext, fname) is True
    assert io.check_ext(all_ext, fname, remove=True) == ('data', True)
    all_ext = ['.CSV']
    assert io.check_ext(all_ext, fname) is True

    fname = 'data.xls'
    assert io.check_ext(all_ext, fname) is False
    assert io.check_ext(all_ext, fname, remove=True) == ('data.xls', False)

    all_ext = []
    fname = 'data.csv'
    assert io.check_ext(all_ext, fname, remove=True) == ('data.csv', False)


@patch('logging.warning')
def test_check_nifti_dim(mock_warning):
    fname = 'valid.nii.gz'
    data = np.ones((4, 4, 4, 4, 4))
    output = io.check_nifti_dim(fname, data, dim=5)
    assert output.ndim == 5

    output = io.check_nifti_dim(fname, data, dim=3)

    # Test if data has more dimensions than `dim`
    output = io.check_nifti_dim(fname, data, dim=3)
    mock_warning.assert_called_with(
        f'{fname} has more than 3 dimensions. Removing D > 3.'
    )
    assert output.ndim == 3


def test_load_nifti_get_mask(nifti_data, nifti_mask):
    # Test loading a nifti file with default arguments
    fname, expected_data = nifti_data
    data, mask, img = io.load_nifti_get_mask(fname)
    assert np.array_equal(data, expected_data)
    assert np.array_equal(mask, expected_data.any(axis=-1).squeeze())
    assert isinstance(img, nib.nifti1.Nifti1Image)

    # Test loading a nifti mask file with is_mask=True
    fname, expected_data = nifti_mask
    data, mask, img = io.load_nifti_get_mask(fname, is_mask=True)
    assert np.array_equal(data, expected_data)
    assert np.array_equal(mask, expected_data != 0)
    assert isinstance(img, nib.nifti1.Nifti1Image)


def test_export_regressor(testdir):
    petco2hrf_shift = np.random.rand(10)
    freq = 1
    tr = 2
    outname = os.path.join(testdir, 'test_regressor')
    suffix = 'petco2hrf'
    ext = '.1D'

    petco2hrf_demean = io.export_regressor(
        petco2hrf_shift, freq, tr, outname, suffix=suffix, ext=ext
    )

    pco2 = (petco2hrf_shift - petco2hrf_shift.mean())[::2]
    assert np.allclose(petco2hrf_demean, pco2)

    # Check if file was saved and has the correct content
    saved_file = np.loadtxt(f'{outname}_{suffix}{ext}')
    assert np.allclose(saved_file, pco2)


def test_export_nifti(testdir):
    # create some test data
    data = np.random.rand(10, 10, 10)
    affine = np.eye(4)
    header = nib.Nifti1Header()
    img = nib.Nifti1Image(data, affine, header)

    # create a temporary directory to store the output file
    fname = os.path.join(testdir, 'test_output.nii.gz')
    io.export_nifti(data, img, fname)

    # check if the output file exists
    assert os.path.exists(fname)

    # load the output file and check if the data matches the input
    out_img = nib.load(fname)
    assert out_img.get_fdata() == data
    assert out_img.affine == affine
    assert out_img.header == header


# ## Break tests


def test_break_if_declared_force_type():
    with pytest.raises(NotImplementedError) as errorinfo:
        io.if_declared_force_type('10', 'invalid_type')
    assert 'not supported' in str(errorinfo.value)


def test_break_check_nifti_dim_missing_dims():
    with pytest.raises(ValueError) as errorinfo:
        io.check_nifti_dim('missing_dims.nii.gz', np.ones((4, 4)), dim=4)
    assert 'seem to be a 4D file.' in str(errorinfo.value)
