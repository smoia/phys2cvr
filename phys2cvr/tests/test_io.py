#!/usr/bin/env python3
"""Tests for io."""

import os
import sys

import nibabel as nib
import numpy as np
import peakdet as pk
import pymatreader
import pytest
from scipy.io import savemat

from phys2cvr import io

# ## Unit tests


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


@pytest.mark.parametrize(
    'extension, delimiter',
    [
        ['.csv', ','],
        ['.csv.gz', ','],
        ['.tsv', '\t'],
        ['.tsv.gz', '\t'],
        ['.txt', ' '],
        ['.1d', ' '],
        ['.par', ' '],
        ['', ' '],
    ],
)
def test_load_txt(extension, delimiter, testdir):
    """Test load_txt."""
    a = np.random.randn(5)
    n = os.path.join(testdir, f'zoe{extension}')
    np.savetxt(n, a, delimiter=delimiter)
    b = io.load_txt(n)

    assert (a == b).all()
    os.remove(n)


def test_load_mat():
    """Test load_mat."""
    a = np.random.randn(5)
    b = 'likealeaf'
    n = 'inthewind'

    savemat(n, {'data': a, 'other': b})

    c = io.load_mat(n)

    assert (a == c).all()
    os.remove(n)


def test_load_array():
    """Test load_mat."""
    a = np.random.randn(5)
    b = 'likealeaf'
    n = 'inthewind.mat'
    m = 'zoe.txt'
    np.savetxt(m, a)

    savemat(n, {'data': a, 'other': b})

    d = io.load_array(n)
    e = io.load_array(m)

    assert (a == d).all()
    assert (a == e).all()
    os.remove(n)
    os.remove(m)


def test_load_physio(tmp_path):
    data = np.random.randn(40)
    peaks = np.array([5, 15, 25, 35])
    d = pk.Physio(
        data,
        fs=10,
        metadata={'peaks': peaks, 'troughs': np.array([]), 'reject': np.array([])},
    )
    n = str(tmp_path / 'kaylee.phys')
    pk.save_physio(n, d)
    p = io.load_physio(n)

    assert (p[0] == d.data).all()
    assert (p[1] == d.peaks).all()
    assert p[2] == 10
    os.remove(n)


def test_export_regressor(testdir):
    petco2hrf_shift = np.random.randn(10)
    ntp = 10
    outname = os.path.join(testdir, 'test_regressor')
    suffix = 'petco2hrf'
    ext = '.1D'

    petco2hrf_demean = io.export_regressor(
        petco2hrf_shift, ntp, outname, suffix=suffix, ext=ext
    )

    pco2 = petco2hrf_shift - petco2hrf_shift.mean()
    assert np.allclose(petco2hrf_demean, pco2, atol=1e-6, rtol=0)

    # Check if file was saved and has the correct content
    saved_file = np.loadtxt(f'{outname}_{suffix}{ext}')
    assert np.allclose(saved_file, pco2, atol=1e-6, rtol=0)
    os.remove(f'{outname}_{suffix}{ext}')


def test_export_nifti(testdir):
    # create some test data
    data = np.random.randn(10, 10, 10)
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
    assert np.allclose(out_img.get_fdata(), data, atol=1e-6, rtol=0)
    assert np.allclose(out_img.affine, affine, atol=1e-6, rtol=0)
    assert out_img.header.__dict__.keys() == header.__dict__.keys()
    os.remove(fname)
    # Eventually check that headers have the same content (although they don't now
    # for good reasons!)


# ## Break tests
def test_break_load_mat():
    """Break load_mat."""
    sys.modules['pymatreader'] = None
    with pytest.raises(ImportError) as errorinfo:
        io.load_mat('simon')
    assert 'is required' in str(errorinfo.value)
    sys.modules['pymatreader'] = pymatreader

    a = 'heart'
    n = 'ofgold'
    savemat(n, {'data': a})

    with pytest.raises(EOFError) as errorinfo:
        io.load_mat(n)
    assert f'{n} does not seem' in str(errorinfo.value)
    os.remove(n)


def test_break_load_xls():
    """Break load_xls."""
    with pytest.raises(NotImplementedError) as errorinfo:
        io.load_xls('firefly')
    assert 'loading is not' in str(errorinfo.value)


def test_break_load_array():
    """Break load_xls."""
    with pytest.raises(NotImplementedError) as errorinfo:
        io.load_array('firefly.leaf')
    assert 'file extension' in str(errorinfo.value)
