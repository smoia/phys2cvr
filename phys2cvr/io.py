#!/usr/bin/env python3
"""
I/O and related utils.

Attributes
----------
LGR :
    logger
"""

import logging

import nibabel as nib
import numpy as np
import scipy.interpolate as spint


LGR = logging.getLogger(__name__)
LGR.setLevel(logging.INFO)


def check_ext(all_ext, fname):
    """
    Check which extension a file has.

    Parameters
    ----------
    all_ext: list
        All possibel extensions to check within
    fname: str
        The filename to check

    Returns
    -------
    has_ext: boolean
        True if the extension is found, false otherwise.
    """
    has_ext = False
    for ext in all_ext:
        if fname.endswith(ext):
            has_ext = True
            break

    return has_ext


def check_nifti_dim(fname, data, dim=4):
    """
    Remove extra dimensions.
    """
    if len(data.shape) < dim:
        raise Exception(f'{fname} does not seem to be a {dim}D file. '
                        f'Plase provide a {dim}D nifti file.')
    if len(data.shape) > dim:
        LGR.warning(f'{fname} has more than {dim} dimensions. Removing D > {dim}.')
        for ax in range(dim, len(data.shape)):
            data = np.delete(data, np.s_[1:], axis=ax)

    return np.squeeze(data)


def load_nifti_get_mask(fname, is_mask=False):
    """
    Load a nifti file and returns its data, its image, and a 3d mask.

    Parameters
    ----------
    fname: str
        The filename to read in
    dim: int
        The number of dimensions expected

    Returns
    -------
    data: np.ndarray
        Data from nifti file.
    mask: np.ndarray

    """
    img = nib.load(fname)
    data = img.get_fdata()

    if is_mask:
        data = check_nifti_dim(fname, data, dim=3)
        mask = (data < 0) + (data > 0)
    else:
        data = check_nifti_dim(fname, data)
        mask = np.squeeze(np.any(data, axis=-1))

    return data, mask, img


def export_regressor(regr_x, petco2hrf_shift, func_x, outname, suffix='petco2hrf', ext='.1D'):
    f = spint.interp1d(regr_x, petco2hrf_shift, fill_value='extrapolate')
    petco2hrf_demean = f(func_x)
    petco2hrf_demean = petco2hrf_demean - petco2hrf_demean.mean()
    np.savetxt(f'{outname}_{suffix}{ext}', petco2hrf_demean, fmt='%.6f')

    return petco2hrf_demean


def export_nifti(data, img, fname):
    """
    Export a nifti file.
    """
    if fname.endswith('.nii.gz'):
        fname = fname[:-7]

    out_img = nib.Nifti1Image(data, img.affine, img.header)
    out_img.to_filename(f'{fname}.nii.gz')
