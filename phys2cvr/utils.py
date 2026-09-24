#!/usr/bin/env python3
"""
Utils functions for phys2cvr.

Attributes
----------
LGR :
    Logger
"""

import datetime
import logging
import os
import sys
from os.path import exists
from pathlib import Path

import numpy as np

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
    if var is None or var == '':
        return var

    converters = {
        'int': int,
        'float': float,
        'str': str,
        'list': lambda v: v if isinstance(v, list) else [v],
    }

    if dtype not in converters:
        raise NotImplementedError(f'Type {dtype} not supported')

    tmpvar = converters[dtype](var)

    if not silent and type(tmpvar) is not type(var):
        name = f'variable {varname}' if varname != 'an input variable' else varname
        LGR.warning(f'Changing type of {name} from {type(var)} to {dtype}')

    return tmpvar


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
    all_ext = set(if_declared_force_type(all_ext, 'list', silent=True))
    ext = ''.join(Path(fname).suffixes)
    LGR.debug(f'{fname} ends with extension {ext}')

    has_ext = ext.lower() in all_ext

    if not has_ext and scan:
        for ext_candidate in all_ext:
            candidate_path = Path(f'{fname}{ext_candidate}')
            if candidate_path.exists():
                fname = str(candidate_path)
                ext = ext_candidate
                LGR.warning(f'Found {fname}, using it as input henceforth')
                has_ext = True
                break

    ext = ext if has_ext else ''

    obj_return = [has_ext]

    if remove:
        obj_return += [
            fname[: -len(ext)],
            None if ext == '' else ext,
        ]
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
    if data.ndim < dim:
        raise ValueError(
            f'A {dim}D nifti file is required, but {fname} has {data.ndim}D. '
            'Please check the input file.'
        )

    if data.ndim > dim:
        LGR.warning(f'{fname} has more than {dim} dimensions. Removing D > {dim}.')
        slc = [slice(None)] * dim + [0] * (data.ndim - dim)
        data = data[tuple(slc)]

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

    if data.size == 0:
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


def load_and_stack_matrices(matrix_files, base_matrix=None, label=None):
    """Load and stack (nuisance) matrices.

    Parameters
    ----------
    matrix_files : list of str or path.Path-like objects
        Path(s) to files containing the matrix to load.
    base_matrix : np.ndarray or None
        Matrix to stack the loaded matrix with. Optional.
    label : None or str
        Type of data tha is being loaded.

    Returns
    -------
    np.ndarray
        Stacked loaded matrices

    """
    # Local import to avoid circularity
    from .io import load_array  # noqa: ABS101

    if not matrix_files:
        return base_matrix

    matrix_files = if_declared_force_type(matrix_files, 'list', silent=True)
    loaded = [base_matrix] if base_matrix is not None else []

    label = 'data' if label is None else label
    for mat_path in matrix_files:
        LGR.info(f'Read {label} from {mat_path}')
        loaded.append(load_array(mat_path))

    return np.hstack(loaded) if len(loaded) > 1 else loaded[0]


def save_bash_call(fname, outdir):
    """
    Save the bash call into file `p2d_call.sh`.

    Parameters
    ----------
    fname : str or path
        Name of or path to functional file
    outdir : str or path
        Output directory
    """
    arg_str = ' '.join(sys.argv[1:])
    call_str = f'phys2cvr {arg_str}'

    outdir = Path(outdir).resolve() if outdir else Path(fname).parent / 'phys2cvr'

    log_path = outdir / 'logs'
    log_path.mkdir(parents=True, exist_ok=True)

    isotime = datetime.datetime.now().strftime('%Y-%m-%dT%H%M%S')
    _, fname, _ = check_ext('.nii.gz', os.path.basename(fname), remove=True)

    script_file = log_path / f'p2c_call_{fname}_{isotime}.sh'

    with open(script_file, 'a', encoding='utf-8') as f:
        f.write(f'#!bin/bash \n{call_str}')


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
