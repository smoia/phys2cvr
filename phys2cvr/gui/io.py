"""Shared GUI data loading and result-folder discovery."""

import glob
import logging
import os
import shlex
from pathlib import Path

import nibabel as nib
import numpy as np

from phys2cvr.io import load_nifti_get_mask

LGR = logging.getLogger(__name__)


def read_spatial_image(path, *, dimensions):
    """Read a nonempty image with one of the requested dimensionalities."""
    image = nib.load(path)
    data = np.asarray(image.get_fdata())
    if data.ndim not in dimensions or any(size == 0 for size in data.shape):
        sizes = ' or '.join(f'{dim}D' for dim in dimensions)
        raise ValueError(f'{path} must be a non-empty {sizes} image.')
    return image, data


def read_log_and_images(path):
    """Read a log and find PNGs under its result folder."""
    path = Path(path)
    text = path.read_text(encoding='utf-8')
    folder = path.parent.parent if path.parent.name == 'logs' else path.parent
    images = sorted(p for p in folder.rglob('*.png') if p.is_file())
    return text, images


def latest_file(folder, pattern):
    """Find the most recently modified matching file, or None."""
    files = (path for path in Path(folder).glob(pattern) if path.is_file())
    return max(
        files, key=lambda path: (path.stat().st_mtime_ns, str(path)), default=None
    )


def read_latest_calls(folder):
    """Read the latest shell and Python calls from a log directory."""
    calls = []
    for pattern in ('p2c_call_*.sh', 'p2c_call_*.py'):
        path = latest_file(folder, pattern)
        calls.append(path.read_text(encoding='utf-8') if path else '')
    return tuple(calls)


def load_timeseries(path):
    """Load a NIfTI file and return data, spatial voxel sizes and TR (seconds)."""
    data, _, img = load_nifti_get_mask(path, dim=4)
    zooms = img.header.get_zooms()
    pixdim = tuple(zooms[:3]) if len(zooms) >= 3 else (1.0, 1.0, 1.0)
    tr = float(zooms[3]) if len(zooms) > 3 else 1.0
    if tr <= 0 or np.isnan(tr):
        tr = 1.0
    return data, pixdim, tr


def load_result_folder(
    folder_path, viewer, comparison_version=None, *, diffviewer=None, logviewer=None
):
    """Scan and load available files and parameters from a phys2cvr result folder."""
    if not os.path.isdir(folder_path):
        raise NotADirectoryError(f'Not a valid directory: {folder_path}')

    LGR.warning('Scanning phys2cvr result folder: %s', folder_path)

    # 1. Parse logs for phys2cvr command arguments if available
    logs_dir = os.path.join(folder_path, 'logs')
    args_p2c = None
    if os.path.isdir(logs_dir):
        latest_sh = latest_file(logs_dir, 'p2c_call_*.sh')
        if latest_sh:
            LGR.warning('Found latest call log script: %s', latest_sh)
            try:
                with open(latest_sh) as f:
                    lines = f.readlines()
                    if lines:
                        last_line = lines[-1].strip()
                        LGR.warning('Parsed command line: %s', last_line)
                        parts = shlex.split(last_line)
                        from phys2cvr.cli.run import _get_parser as get_p2c_parser

                        p2c_parser = get_p2c_parser()

                        def parser_error(message):
                            raise ValueError(message)

                        p2c_parser.error = parser_error
                        args_p2c = p2c_parser.parse_args(parts[1:])
            except Exception as err:
                LGR.warning('Could not parse call script %s: %s', latest_sh, err)

    # Determine base directory for resolving relative paths from args_p2c
    # If outdir is specified and relative, paths are relative to the caller's working directory
    # (which corresponds to folder_path if outdir points to it).
    if args_p2c is not None and getattr(args_p2c, 'outdir', None):
        outdir = args_p2c.outdir
        if not os.path.isabs(outdir):
            # If outdir is relative, we can treat folder_path as the resolved root or relative base
            LGR.warning('Parsed outdir from command args: %s (relative)', outdir)
        else:
            LGR.warning('Parsed outdir from command args: %s (absolute)', outdir)

    def resolve_path(p):
        if not p:
            return None
        if os.path.isabs(p):
            return p

        # Check if it exists directly relative to folder_path
        cand = os.path.join(folder_path, p)
        if os.path.exists(cand):
            return cand
        # Check if it exists directly relative to parent of folder_path (if folder_path is outdir)
        cand_parent = os.path.join(os.path.dirname(folder_path), p)
        if os.path.exists(cand_parent):
            return cand_parent
        return p

    def set_map(kind, path):
        setter = getattr(viewer, f'set_{kind}')
        if comparison_version is None:
            setter(path)
        else:
            setter(path, comparison_version)

    # Apply parameters extracted from logs
    if args_p2c is not None:
        func_path = resolve_path(getattr(args_p2c, 'fname_func', None))
        if (
            comparison_version is None
            and func_path
            and str(func_path).lower().endswith(('.nii', '.nii.gz'))
            and os.path.exists(func_path)
        ):
            LGR.warning('Loading timeseries from log args: %s', func_path)
            viewer.set_timeseries(func_path)

        mask_path = resolve_path(getattr(args_p2c, 'fname_mask', None))
        if mask_path and os.path.exists(mask_path):
            LGR.warning('Loading mask from log args: %s', mask_path)
            set_map('mask', mask_path)

        if comparison_version is None and getattr(args_p2c, 'freq', None) is not None:
            viewer.fs = float(args_p2c.freq)
            LGR.warning('Loaded sampling frequency (fs) from log args: %s', viewer.fs)

        lag_max = getattr(args_p2c, 'lag_max', None)
        lag_min = getattr(args_p2c, 'lag_min', None)
        if comparison_version is None and lag_max is not None:
            l_max = float(lag_max)
            l_min = float(lag_min) if lag_min is not None else -l_max
            viewer.lag_range = (l_min, l_max)
            LGR.warning('Loaded lag range from log args: %s', viewer.lag_range)
        if comparison_version is None and getattr(args_p2c, 'ls', None) is not None:
            viewer.lag_step = float(args_p2c.ls)
            LGR.warning('Loaded lag step from log args: %s', viewer.lag_step)

    # 2. Find CVR map (_cvr.* or fallback to _cvr_simple.*)
    cvr_files = sorted(glob.glob(os.path.join(folder_path, '*_cvr.*')))
    cvr_simple_files = sorted(glob.glob(os.path.join(folder_path, '*_cvr_simple.*')))
    cvr_target = None
    for f in cvr_files:
        if not f.endswith('_cvr_simple.nii.gz'):
            cvr_target = f
            break
    if not cvr_target and cvr_simple_files:
        cvr_target = cvr_simple_files[0]

    if cvr_target and os.path.exists(cvr_target):
        LGR.warning('Loading CVR map: %s', cvr_target)
        set_map('cvr', cvr_target)
    else:
        LGR.warning('No CVR map found in result folder.')

    # 3. Find Lag map
    lag_files = sorted(glob.glob(os.path.join(folder_path, '*_lag.*')))
    log_lag_map = resolve_path(getattr(args_p2c, 'lag_map', None)) if args_p2c else None

    if lag_files and os.path.exists(lag_files[0]):
        LGR.warning('Loading Lag map: %s', lag_files[0])
        set_map('lag', lag_files[0])
    elif log_lag_map and os.path.exists(log_lag_map):
        LGR.warning('Loading Lag map from log args: %s', log_lag_map)
        set_map('lag', log_lag_map)
    else:
        LGR.warning('No Lag map found in result folder.')

    # 4. Find R² map (_r_square_all.* or fallback to _r_square_simple.*)
    r2_files = sorted(glob.glob(os.path.join(folder_path, '*_r_square_all.*')))
    r2_simple_files = sorted(
        glob.glob(os.path.join(folder_path, '*_r_square_simple.*'))
    )
    r2_target = None
    for f in r2_files:
        if not f.endswith('_r_square_simple.nii.gz'):
            r2_target = f
            break
    if not r2_target and r2_simple_files:
        r2_target = r2_simple_files[0]

    if r2_target and os.path.exists(r2_target):
        LGR.warning('Loading R² map: %s', r2_target)
        set_map('r2', r2_target)
    else:
        LGR.warning('No R² map found in result folder.')

    # 5. Find regressor in "regr" subfolder (or via args_p2c.regr_dir)
    if comparison_version is None:
        regr_dir = (
            resolve_path(getattr(args_p2c, 'regr_dir', None)) if args_p2c else None
        )
        if not regr_dir or not os.path.isdir(regr_dir):
            regr_dir = os.path.join(folder_path, 'regr')

        if os.path.isdir(regr_dir):
            all_regr = sorted(glob.glob(os.path.join(regr_dir, '*')))
            if all_regr:
                LGR.warning('Loading regressor file: %s', all_regr[0])
                viewer.set_regressors(all_regr[0])
            else:
                LGR.warning('Regressor directory "regr" (%s) is empty.', regr_dir)
        else:
            LGR.warning('No regr directory found in result folder.')

    if diffviewer is not None and comparison_version is None:
        load_result_folder(folder_path, diffviewer, comparison_version='base')
        diffviewer.refresh()
    if logviewer is not None:
        logs = glob.glob(os.path.join(logs_dir, '*.tsv'))
        if logs:
            latest = max(logs, key=lambda path: (os.path.getmtime(path), path))
            logviewer.load(latest)
        else:
            logviewer.clear()
            logviewer.load_calls(logs_dir)
