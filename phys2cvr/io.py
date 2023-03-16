#!/usr/bin/env python3
"""
I/O and related utils.

Attributes
----------
EXT_1D : list
    List of supported TXT/1D file extensions, in lower case.
EXT_NIFTI : list
    List of supported nifti file extensions, in lower case.
FIGSIZE : tuple
    Figure size
SET_DPI : int
    DPI of the figure
LGR :
    Logger
"""

import logging
from os.path import exists
from pathlib import Path

import nibabel as nib
import numpy as np

from phys2cvr import signal

SET_DPI = 100
FIGSIZE = (18, 10)
EXT_1D = ['.txt', '.csv', '.tsv', '.1d', '.par', '.tsv.gz']
EXT_MAT = ['.mat']
# EXT_XLS = [".xls"]
EXT_ARRAY = EXT_1D + EXT_MAT  # + EXT_XLS
EXT_NIFTI = ['.nii', '.nii.gz']
EXT_GIFTI = ['.gii', '.gii.gz']
EXT_ALL = EXT_ARRAY + EXT_NIFTI + EXT_GIFTI


LGR = logging.getLogger(__name__)
LGR.setLevel(logging.INFO)


def if_declared_force_type(var, dtype, varname='an input variable', silent=False):
    """
    Make sure `var` is of type `dtype`.

    Parameters
    ----------
    var : str, int, or float
        Variable to change type of
    dtype : str
        Type to change `var` to
    varname : str, optional
        The name of the variable
    silent : bool, optional
        If True, don't return any message

    Returns
    -------
    int, float, str, list, or var
        The given `var` in the given `dtype`, or `var` if '' or None

    Raises
    ------
    NotImplementedError
        If dtype is not 'int', 'float', 'str', or 'list'
    """
    if var:
        if dtype == 'int':
            tmpvar = int(var)
        elif dtype == 'float':
            tmpvar = float(var)
        elif dtype == 'str':
            tmpvar = str(var)
        elif dtype == 'list':
            if type(var) is list:
                tmpvar = var
            else:
                tmpvar = [var]
        else:
            raise NotImplementedError(f'Type {dtype} not supported')

        if not silent:
            if type(tmpvar) is not type(var):
                if varname != 'an input variable':
                    varname = f'variable {varname}'

                LGR.warning(f'Changing type of {varname} from {type(var)} to {dtype}')

        return tmpvar

    else:
        return var


def check_ext(all_ext, fname, scan=False, remove=False):
    """Check which extension a file has, and possibly remove it.

    Parameters
    ----------
    all_ext : list
        All possible extensions to check within.
    fname : str or os.PathLike
        The filename to check.
    scan : bool, optional
        Scan the given path to see if there is a file with that extension
        If True and no path declared, check if fname has a path, if not scan '.'
        If False, don't scan any folder.
    remove : bool, optional
        Remove the extension from fname if it has one.

    Returns
    -------
    obj_return : Uses a list to return variable amount of options.
        has_ext : boolean
            True if the extension is found, false otherwise.
        fname : str or os.PathLike
            If ``remove`` is True, return (extensionless) fname.
        ext : str
            If both ``remove`` and ``has_ext`` are True, returns also found extension.
    """
    all_ext = if_declared_force_type(all_ext, 'list', silent=True)

    ext = ''.join(Path(fname).suffixes)

    LGR.debug(f'{fname} ends with extension {ext}')

    has_ext = True if ext in all_ext else False

    if not has_ext and scan:
        for ext in all_ext:
            if exists(f'{fname}{ext}'):
                fname = f'{fname}{ext}'
                LGR.warning(f'Found {fname}{ext}, using it as input henceforth')
                has_ext = True
                break

    obj_return = [has_ext]

    if remove:
        obj_return += [
            fname[: -len(ext)],
            None if ext == '' else ext,
        ]  # case insensitive solution
    else:
        obj_return += [fname]

    return obj_return[:]


def check_nifti_dim(fname, data, dim=4):
    """
    Remove extra dimensions.

    Parameters
    ----------
    fname : str
        The name of the file representing `data`
    data : np.ndarray
        The data which dimensionality needs to be checked
    dim : int, optional
        The amount of dimensions expected/desired in the data.

    Returns
    -------
    np.ndarray
        If `len(data.shape)` = `dim`, returns data.
        If `len(data.shape)` > `dim`, returns a version of data without the
        dimensions above `dim`.

    Raises
    ------
    ValueError
        If `data` has less dimensions than `dim`
    """
    if len(data.shape) < dim:
        raise ValueError(
            f'A {dim}D nifti file is required, but {fname} has {data.ndim}D. Please '
            'check the input file.'
        )
    if len(data.shape) > dim:
        LGR.warning(f'{fname} has more than {dim} dimensions. Removing D > {dim}.')
        for ax in range(dim, len(data.shape)):
            data = np.delete(data, np.s_[1:], axis=ax)

    return np.squeeze(data)


def check_array_dim(fname, data, shape=None):
    """Check dimensions of a matrix.

    For future 3D implementation, check MIPLabCH/nigsp's check_array_dim.

    Parameters
    ----------
    fname : str
        The name of the file representing ``data``.
    data : np.ndarray
        The data which dimensionality needs to be checked.
    shape : None | ``'square'`` | ``'rectangle'``
        Shape of matrix, if empty, skip shape check.

    Returns
    -------
    np.ndarray
        If ``data.ndim = 2``, returns data.
        If ``data.ndim = 1`` and ``shape == 'rectangle'``, returns data with added empty
        axis.

    Raises
    ------
    NotImplementedError
        If ``data`` has more than 2 dimensions.
    ValueError
        If ``data`` is empty
        If ``shape == 'square'`` and ``data`` dimensions have different lengths.
    """
    data = data.squeeze()
    LGR.info('Checking data shape.')

    if data.shape[0] == 0:
        raise ValueError(f'{fname} is empty!')
    if data.ndim > 2:
        raise NotImplementedError(
            f'Only matrices up to 2D are supported, but given matrix is {data.ndim}D.'
        )
    if shape is not None:
        if data.ndim == 1 and shape == 'rectangle':
            data = data[..., np.newaxis]
            LGR.warning(
                f'Rectangular matrix required, but {fname} is a vector. '
                'Adding empty dimension.'
            )
        if shape == 'square' and data.shape[0] != data.shape[1]:
            raise ValueError(
                f'Square matrix required, but {fname} matrix has shape {data.shape}.'
            )

    return data


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

    if check_ext(EXT_GIFTI, fname)[0]:
        data = img.agg_data().transpose()
    else:
        data = img.get_fdata()

    data = check_nifti_dim(fname, data, dim=dim)

    if is_mask:
        mask = data != 0
    else:
        mask = data.any(axis=-1).squeeze()

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
    check_array_dim
    """
    LGR.info(f'Loading {fname}.')

    _, _, ext = check_ext(EXT_1D, fname, scan=True, remove=True)

    if ext in ['.csv', '.csv.gz']:
        delimiter = ','
    elif ext in ['.tsv', '.tsv.gz']:
        delimiter = '\t'
    elif ext in ['.txt', '.1d', '.par']:
        delimiter = ' '
    else:
        delimiter = None

    mtx = np.genfromtxt(fname, delimiter=delimiter)

    mtx = check_array_dim(fname, mtx, shape)

    return mtx


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
    check_array_dim

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
            'pymatreader is required to import mat files. '
            'Please see install instructions.'
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

    if len(data_keys) < 1:
        raise EOFError(f'{fname} does not seem to contain a numeric matrix.')
    elif len(data_keys) > 1:
        LGR.warning(
            'Found multiple possible arrays to load. '
            'Selecting the biggest (highest pythonic size).'
        )

    key = data_keys[0]
    for k in data_keys[1:]:
        if data[k].size > data[key].size:
            key = k

    LGR.info(f'Selected data from MATLAB variable {key}')
    mtx = data[key]
    mtx = check_array_dim(fname, mtx, shape)

    return mtx


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
    check_array_dim

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
    check_array_dim

    Raises
    ------
    NotImplementedError
        Spreadheet loading is not implemented yet.
    """
    _, _, ext = check_ext(EXT_ARRAY, fname, scan=True, remove=True)

    if ext in EXT_1D:
        mtx = load_txt(fname, shape=shape)
    elif ext in EXT_MAT:
        mtx = load_mat(fname, shape=shape)
    # elif ext in EXT_XLS:
    #     mtx = load_xls(fname, shape=shape)
    else:
        raise NotImplementedError(
            f'{fname} file extension {ext} was not found or is not supported yet'
        )

    return mtx


def export_regressor(petco2hrf_shift, freq, tr, outname, suffix='petco2hrf', ext='.1D'):
    """
    Export generated regressors for fMRI analysis.

    Parameters
    ----------
    petco2hrf_shift : np.ndarray
        The regressors that needs to be exported, in its original sample
    ntp : int
        The number of fMRI timepoints
    outname : str or path
        Prefix of the output file - can contain a path.
    suffix : str, optional
        The suffix of the output file.
    ext : str, optional
        The extension of the output file.

    Returns
    -------
    petco2hrf_demean : np.ndarray
        Interpolated version of `petco2hrf_shift` in the sampling of the fMRI data.
    """
    petco2hrf_shift = signal.resample_signal_samples(petco2hrf_shift, ntp)
    petco2hrf_demean = petco2hrf_shift - petco2hrf_shift.mean(axis=-1)
    if petco2hrf_demean.ndim > 1:
        for i in range(petco2hrf_demean.shape[-1]):
            np.savetxt(
                f"{outname}_{suffix}{i:04g}{ext}", petco2hrf_demean[:, i], fmt="%.6f"
            )
    else:
        np.savetxt(f"{outname}_{suffix}{ext}", petco2hrf_demean, fmt="%.6f")

    return petco2hrf_demean


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
Copyright 2021, Stefano Moia.

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
