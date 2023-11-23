#!/usr/bin/env python3
"""Tests for io."""

import logging
import os

import nibabel as nib
import numpy as np
import pytest

from phys2cvr import io

LGR = logging.getLogger(__name__)


# ## Unit tests
@pytest.mark.parametrize(
    "var, dtype, out",
    [
        (10, "int", 10),
        (10.0, "int", 10),
        ("10", "int", 10),
        (None, "int", None),
        (10, "float", 10.0),
        (10.0, "float", 10.0),
        ("10.0", "float", 10.0),
        (None, "float", None),
        (10, "str", "10"),
        (10.0, "str", "10.0"),
        ("10", "str", "10"),
        ([10], "str", "[10]"),
        (None, "str", None),
        (10, "list", [10]),
        ([10], "list", [10]),
        ("10", "list", ["10"]),
        (None, "list", None),
    ],
)
def test_if_declared_force_type(var, dtype, out):
    assert io.if_declared_force_type(var, dtype) == out


def test_check_ext():
    all_ext = [".csv", ".txt"]
    fname = "data.csv"
    assert io.check_ext(all_ext, fname) is True
    assert io.check_ext(all_ext, fname, remove=True) == ("data", True)
    all_ext = [".CSV"]
    assert io.check_ext(all_ext, fname) is True

    fname = "data.xls"
    assert io.check_ext(all_ext, fname) is False
    assert io.check_ext(all_ext, fname, remove=True) == ("data.xls", False)

    all_ext = []
    fname = "data.csv"
    assert io.check_ext(all_ext, fname, remove=True) == ("data.csv", False)


def test_check_nifti_dim():
    fname = "valid.nii.gz"
    data = np.ones((4, 4, 4, 4, 4))
    output = io.check_nifti_dim(fname, data, dim=5)
    assert output.ndim == 5

    output = io.check_nifti_dim(fname, data, dim=3)

    # Test if data has more dimensions than `dim`
    output = io.check_nifti_dim(fname, data, dim=3)
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
    petco2hrf_shift = np.arange(10)
    ntp = 5
    outname = os.path.join(testdir, "test_regressor")
    suffix = "petco2hrf"
    ext = ".1D"

    petco2hrf_demean = io.export_regressor(
        petco2hrf_shift, ntp, outname, suffix=suffix, ext=ext
    )

    pco2 = np.asarray([-4.5, -2.25, 0, 2.25, 4.5])
    assert np.allclose(petco2hrf_demean, pco2)

    # Check if file was saved and has the correct content
    saved_file = np.loadtxt(f"{outname}_{suffix}{ext}")
    assert np.allclose(saved_file, pco2)


def test_export_nifti(testdir):
    # create some test data
    data = np.random.rand(10, 10, 10)
    affine = np.eye(4)
    header = nib.Nifti1Header()
    img = nib.Nifti1Image(data, affine, header)

    # create a temporary directory to store the output file
    fname = os.path.join(testdir, "test_output.nii.gz")
    io.export_nifti(data, img, fname)

    # check if the output file exists
    assert os.path.exists(fname)

    # load the output file and check if the data matches the input
    out_img = nib.load(fname)
    assert np.allclose(out_img.get_fdata(), data)
    assert (out_img.affine == affine).all()

def test_array_is_2d():
    # Test case: Valid 2D array
    input_array = np.random.rand(3, 4)
    output_array = array_is_2d(input_array)
    assert output_array.shape == (3, 4)

    # Test case: 1D array
    input_array = np.random.rand(5)
    output_array = array_is_2d(input_array)
    assert output_array.shape == (5, 1)

def test_load_regressor_matrices(tmp_path): 
    # Test for checking if the files upload is correct
    file_path = tmp_path / "regressors.txt"
    np.savetxt(file_path, np.random.rand(10, 3))

    # Test loading a single file
    result = load_regressor_matrices(file_path)
    assert result.shape == (10, 3)
    
    # Test for concatenation of additional matrices
    file_path = tmp_path / "regressors.txt"
    np.savetxt(file_path, np.random.rand(10, 3))
    additional_matrix = np.ones((10, 1))

    result = load_regressor_matrices(file_path, additional_matrix=additional_matrix)
    assert result.shape == (10, 4)

    # Test for manipulation of temporal dimension (ntp)
    file_path = tmp_path / "regressors.txt"
    np.savetxt(file_path, np.random.rand(10, 3))

    result = load_regressor_matrices(file_path, ntp=10)
    assert result.shape == (10, 3)

def test_array_is_2d():
    # Test case: Valid 2D array
    input_array = np.random.rand(3, 4)
    output_array = array_is_2d(input_array)
    assert output_array.shape == (3, 4)

    # Test case: 1D array
    input_array = np.random.rand(5)
    output_array = array_is_2d(input_array)
    assert output_array.shape == (5, 1)


def test_load_regressor_matrices(tmp_path):
    # Test for checking if the files upload is correct
    file_path = tmp_path / "regressors.txt"
    np.savetxt(file_path, np.random.rand(10, 3))

    # Test loading a single file
    result = load_regressor_matrices(file_path)
    assert result.shape == (10, 3)

    # Test for concatenation of additional matrices
    file_path = tmp_path / "regressors.txt"
    np.savetxt(file_path, np.random.rand(10, 3))
    additional_matrix = np.ones((10, 1))

    result = load_regressor_matrices(file_path, additional_matrix=additional_matrix)
    assert result.shape == (10, 4)

    # Test for manipulation of temporal dimension (ntp)
    file_path = tmp_path / "regressors.txt"
    np.savetxt(file_path, np.random.rand(10, 3))

    result = load_regressor_matrices(file_path, ntp=10)
    assert result.shape == (10, 3)


# ## Break tests
def test_break_if_declared_force_type():
    with pytest.raises(NotImplementedError) as errorinfo:
        io.if_declared_force_type("10", "invalid_type")
    assert "not supported" in str(errorinfo.value)


def test_break_check_nifti_dim_missing_dims():
    with pytest.raises(ValueError) as errorinfo:
        io.check_nifti_dim("missing_dims.nii.gz", np.ones((4, 4)), dim=4)
    assert "seem to be a 4D file." in str(errorinfo.value)


def test_array_is_2d():
    # Test case: 0D array (scalar)
    input_array = np.random.rand()
    with pytest.raises(
        ValueError, match="Files with 0 dimensions are not supported yet"
    ):
        array_is_2d(input_array)
    # Test case: 3D array
    input_array = np.random.rand(2, 3, 4)
    with pytest.raises(
        ValueError, match="Files with 3 dimensions are not supported yet"
    ):
        array_is_2d(input_array)
    # Test case: Empty array
    input_array = np.array([])
    with pytest.raises(
        ValueError, match="Files with 0 dimensions are not supported yet"
    ):
        array_is_2d(input_array)


def test_load_regressor_matrices_dimension_error(tmp_path):
    file_path = tmp_path / "regressors.txt"
    np.savetxt(file_path, np.random.rand(10, 3))
    additional_matrix = np.ones((5, 2))
    with pytest.raises(
        ValueError,
        match="Loaded matrix has shape .* but additional matrix has shape .*",
    ):
        load_regressor_matrices(file_path, additional_matrix=additional_matrix)
