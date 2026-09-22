#!/usr/bin/env python3
"""
I/O and related utils.

Attributes
----------
EXT_1D : list
    List of supported TXT/1D file extensions, in lower case.
EXT_MAT : list
    List of supported matlab file extensions, in lower case.
EXT_ARRAY : list
    List of supported 1D and 2D array-like file extensions, in lower case.
EXT_PHYS : list
    List of supported physiopy file extensions, in lower case.
EXT_NIFTI : list
    List of supported nifti file extensions, in lower case.
EXT_GIFTI : list
    List of supported gifti file extensions, in lower case.
EXT_NIMG : list
    List of supported neuroimaging file extensions, in lower case.
EXT_ALL : list
    List of ALL supported file extensions, in lower case.
LGR :
    Logger
"""

import logging

import nibabel as nib
import numpy as np
from peakdet.io import load_physio as load_pk_physio

from phys2cvr import utils

EXT_1D = ['.txt', '.csv', '.tsv', '.1d', '.par', '.tsv.gz', '.csv.gz', '.mcdat']
EXT_MAT = ['.mat']
# EXT_XLS = [".xls"]
EXT_ARRAY = EXT_1D + EXT_MAT  # + EXT_XLS
EXT_PHYS = ['.phys']
EXT_NIFTI = ['.nii', '.nii.gz']
EXT_GIFTI = ['.gii', '.gii.gz']
EXT_NIMG = EXT_NIFTI + EXT_GIFTI
EXT_ALL = EXT_ARRAY + EXT_PHYS + EXT_NIMG


LGR = logging.getLogger(__name__)
LGR.setLevel(logging.INFO)


def load_nifti_get_mask(fname, is_mask=False, dim=3):
    """
    Load a nifti-like file and returns its data, its image, and a 3d mask.

    Support all nibabel supported filetypes

    Parameters
    ----------
    fname : str
        The filename to read in
    is_mask : bool, optional
        If the file contains a mask.
        Default: False
    dim : int
        The number of dimensions expected in fname

    Returns
    -------
    data : np.ndarray
        Data from nifti file.
    mask : np.ndarray
        If `is_mask` is False, np.ndarray of one dimension less than data,
        in which any element that has at least a value different from zero
        in the last dimension of `data` is True.
        If `is_mask` is True, mask is a boolean representation of data.
    img : nib.img
        Image object from nibabel.
    """
    img = nib.load(fname)
    LGR.info(f'Loading {fname}')

    if utils.check_ext(EXT_GIFTI, fname)[0]:
        data = img.agg_data().transpose()
    else:
        data = img.get_fdata()

    data = utils.check_nifti_dim(fname, data, dim=dim)
    mask = (data != 0) if is_mask else data.any(axis=-1).squeeze()

    return data, mask, img


def load_txt(fname, shape=None):
    """Read files in textual format.

    Parameters
    ----------
    fname : str | os.PathLike
        Path to the txt file.
    shape : None | ``'square'`` | ``'rectangle'``
        Shape of matrix, if empty, skip check.

    Returns
    -------
    mtx : numpy.ndarray
        Data matrix.

    See Also
    --------
    utils.check_array_dim
    """
    LGR.info(f'Loading {fname}.')
    _, _, ext = utils.check_ext(EXT_1D, fname, scan=True, remove=True)

    delimiter_map = {
        '.csv': ',',
        '.csv.gz': ',',
        '.tsv': '\t',
        '.tsv.gz': '\t',
        '.txt': ' ',
        '.1d': ' ',
        '.par': ' ',
        '.mcdat': None,
        '': ' ',
    }

    mtx = np.genfromtxt(fname, delimiter=delimiter_map.get(ext))

    mtx = mtx[:, 1:7] if ext == '.mcdat' else mtx

    return utils.check_array_dim(fname, mtx, shape)


def load_mat(fname, shape=None):
    """Read files in MATLAB format.

    Assumes the existence of a matrix/vector in the mat file, rendered as
    a numpy.ndarray. If there is more than a matrix, the one with the largest
    size will be selected.

    Parameters
    ----------
    fname : str | os.PathLike
        Path to the ``.mat`` file.
    shape : None | ``'square'`` | ``'rectangle'``
        Shape of matrix, if empty, skip check.

    Returns
    -------
    mtx : numpy.ndarray
        Data matrix.

    Notes
    -----
    Requires module ``pymatreader`` to work.

    See Also
    --------
    utils.check_array_dim

    Raises
    ------
    EOFError
        If the mat file does not contain matrix or vectors.
    ImportError
        If pymatreader is not installed or can't be read.
    """
    try:
        from pymatreader import read_mat
    except ImportError:
        raise ImportError(
            'pymatreader is required to import mat files. Please see install instructions.'
        )

    LGR.info(f'Loading {fname}.')
    data = read_mat(fname)

    data_keys = []
    for k in data.keys():
        # Check data key only if it's not hidden
        # (skip '__header__', '__version__', '__global__')
        if '__' not in k:
            LGR.info(
                f'Checking {fname} key {str(k)} content for data '
                '(float array/matrices in MATLAB).'
            )
            if type(data[k]) is np.ndarray:
                data_keys.append(k)

    if not data_keys:
        raise EOFError(f'{fname} does not seem to contain a numeric matrix.')
    if len(data_keys) > 1:
        LGR.warning('Found multiple possible arrays to load. Selecting the biggest.')

    key = max(data_keys, key=lambda k: data[k].size)
    LGR.info(f'Selected data from MATLAB variable {key}')

    return utils.check_array_dim(fname, data[key], shape)


def load_xls(fname, shape=''):
    """Read files in xls format.

    Parameters
    ----------
    fname : str | os.PathLike
        Path to the xls file.
    shape : None | ``'square'`` | ``'rectangle'``
        Shape of matrix, if empty, skip check.

    See Also
    --------
    utils.check_array_dim

    Raises
    ------
    NotImplementedError
        Spreadheet loading is not implemented yet.
    """
    raise NotImplementedError('Spreadsheet loading is not implemented yet')


def load_array(fname, shape=''):
    """Read files in text-like format.

    Parameters
    ----------
    fname : str | os.PathLike
        Path to txt-like file.
    shape : None | ``'square'`` | ``'rectangle'``
        Shape of matrix, if empty, skip check.

    See Also
    --------
    utils.check_array_dim

    Raises
    ------
    NotImplementedError
        Spreadheet loading is not implemented yet.
    """
    _, _, ext = utils.check_ext(EXT_ARRAY, fname, scan=True, remove=True)

    if ext:
        if ext.lower() in EXT_1D:
            return load_txt(fname, shape=shape)
        if ext.lower() in EXT_MAT:
            return load_mat(fname, shape=shape)

    raise NotImplementedError(
        f'{fname} file extension {ext} was not found or is not supported yet'
    )


def load_physio(fname):
    """Read peakdet and physiopy objects.

    Parameters
    ----------
    fname : str | os.PathLike
        Path to the xls file.

    Returns
    -------
    np.ndarray
        The physiological data
    np.ndarray
        The indexes of peaks
    np.float
        The sampling frequency
    """
    phys = load_pk_physio(fname, allow_pickle=True)

    return phys.data, phys.peaks, phys.fs


def export_regressor(
    regressors_matrix, ntp, outprefix, suffix='petco2hrf', ext='.1D', axis=-1
):
    """
    Export generated regressors for fMRI analysis.

    Parameters
    ----------
    regressors_matrix : np.ndarray
        The regressors that needs to be exported, in its original sample
    ntp : int
        The number of fMRI timepoints
    outprefix : str or path
        Prefix of the output file - can contain a path.
    suffix : str, optional
        The suffix of the output file.
    ext : str, optional
        The extension of the output file.
    axis : int, optional
        The axis along which to perform the operation. Default is -1.

    Returns
    -------
    regressors_demeaned : np.ndarray
        Interpolated and demeaned version of `regressors_matrix` in the sampling of the
        fMRI data.
    """
    # Local import to avoid circularity
    from .signal import resample_signal_samples  # noqa: ABS101

    regressors_matrix = resample_signal_samples(regressors_matrix, ntp, axis=axis)
    regressors_demeaned = regressors_matrix - regressors_matrix.mean(
        axis=axis, keepdims=True
    )
    np.savetxt(f'{outprefix}_{suffix}{ext}', regressors_demeaned, fmt='%.6f')
    return regressors_demeaned


def export_nifti(data, img, fname):
    """
    Export a nifti file.

    Parameters
    ----------
    data : np.ndarray
        Data to be exported
    img : nib.img
        Nibabel image object
    fname : str or path
        Name of the output file
    """
    klass = img.__class__
    out_img = klass(data, img.affine, img.header)
    out_img.to_filename(fname)


"""
Copyright 2021-2026, Stefano Moia & phys2cvr contributors.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""
