#!/usr/bin/env python3
"""Application shell for the phys2cvr GUI."""

import argparse
import glob
import logging
import os
import sys

import bootstack as bs
from bootstack.widgets.sidebar_toggle import SidebarToggle

from phys2cvr import __version__
from phys2cvr.gui.pages.diffviewer import DiffViewer
from phys2cvr.gui.pages.guiargset import build_workflow_setter
from phys2cvr.gui.pages.logviewer import build_log_viewer
from phys2cvr.gui.pages.voxelviewer import VoxelViewer

LGR = logging.getLogger(__name__)


def _get_parser():
    """Create the viewer argument parser."""
    parser = argparse.ArgumentParser(description='phys2cvr viewer.')
    parser.add_argument(
        '-ts',
        '--timeseries',
        default=None,
        help='Path to fMRI file (i.e. the input to phys2cvr)',
    )
    parser.add_argument(
        '-r',
        '--regressor',
        default=None,
        help='Path to file containing all the lagged versions of the regressor of interest',
    )
    parser.add_argument(
        '-r2',
        '--r2',
        default=None,
        help='Path to 4D R² map volume across computed lag GLMs',
    )
    parser.add_argument('-cvr', '--cvr', default=None, help='Path to 3D CVR map')
    parser.add_argument('-lag', '--lag', default=None, help='Path to 3D Lag map')
    parser.add_argument('-m', '--mask', default=None, help='Path to 3D brain mask file')
    parser.add_argument(
        '-fs',
        '--fs',
        type=float,
        default=None,
        help='Regressor sampling frequency in Hz (inferred from the lag range if omitted)',
    )
    parser.add_argument(
        '-lr',
        '--lag_range',
        type=float,
        nargs=2,
        default=None,
        metavar=('START', 'END'),
        help='Set lag range to explore (will be overridden when loading a lag map)',
    )
    parser.add_argument(
        '-rfld',
        '--result_folder',
        default='/home/nemo/Scrivania/Workbench/physiopy/Tests/phys2cvr/gui/p2c_stnd',
        # default=None,
        help='Load a phys2cvr output folder and let the GUI load everything for you!',
    )
    return parser


def load_result_folder(folder_path, viewer):
    """Scan and load available files and parameters from a phys2cvr result folder."""
    if not os.path.isdir(folder_path):
        raise NotADirectoryError(f'Not a valid directory: {folder_path}')

    LGR.warning('Scanning phys2cvr result folder: %s', folder_path)

    # 1. Parse logs for phys2cvr command arguments if available
    logs_dir = os.path.join(folder_path, 'logs')
    args_p2c = None
    if os.path.isdir(logs_dir):
        sh_files = sorted(glob.glob(os.path.join(logs_dir, 'p2c_call_*.sh')))
        if sh_files:
            latest_sh = sh_files[-1]
            LGR.warning('Found latest call log script: %s', latest_sh)
            try:
                with open(latest_sh) as f:
                    lines = f.readlines()
                    if lines:
                        last_line = lines[-1].strip()
                        LGR.warning('Parsed command line: %s', last_line)
                        parts = last_line.split()
                        from phys2cvr.cli.run import _get_parser as get_p2c_parser

                        p2c_parser = get_p2c_parser()
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

    # Apply parameters extracted from logs
    if args_p2c is not None:
        func_path = resolve_path(getattr(args_p2c, 'fname_func', None))
        if func_path and os.path.exists(func_path):
            LGR.warning('Loading timeseries from log args: %s', func_path)
            viewer.set_timeseries(func_path)

        mask_path = resolve_path(getattr(args_p2c, 'fname_mask', None))
        if mask_path and os.path.exists(mask_path):
            LGR.warning('Loading mask from log args: %s', mask_path)
            viewer.set_mask(mask_path)

        if getattr(args_p2c, 'freq', None) is not None:
            viewer.fs = float(args_p2c.freq)
            LGR.warning('Loaded sampling frequency (fs) from log args: %s', viewer.fs)

        lag_max = getattr(args_p2c, 'lag_max', None)
        lag_min = getattr(args_p2c, 'lag_min', None)
        if lag_max is not None:
            l_max = float(lag_max)
            l_min = float(lag_min) if lag_min is not None else -l_max
            viewer.lag_range = (l_min, l_max)
            LGR.warning('Loaded lag range from log args: %s', viewer.lag_range)
        if getattr(args_p2c, 'ls', None) is not None:
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
        viewer.set_cvr(cvr_target)
    else:
        LGR.warning('No CVR map found in result folder.')

    # 3. Find Lag map
    lag_files = sorted(glob.glob(os.path.join(folder_path, '*_lag.*')))
    log_lag_map = resolve_path(getattr(args_p2c, 'lag_map', None)) if args_p2c else None

    if lag_files and os.path.exists(lag_files[0]):
        LGR.warning('Loading Lag map: %s', lag_files[0])
        viewer.set_lag(lag_files[0])
    elif log_lag_map and os.path.exists(log_lag_map):
        LGR.warning('Loading Lag map from log args: %s', log_lag_map)
        viewer.set_lag(log_lag_map)
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
        viewer.set_r2(r2_target)
    else:
        LGR.warning('No R² map found in result folder.')

    # 5. Find regressor in "regr" subfolder (or via args_p2c.regr_dir)
    regr_dir = resolve_path(getattr(args_p2c, 'regr_dir', None)) if args_p2c else None
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


def build_app(voxelviewer):
    """Build the application shell and its viewer page."""
    with bs.AppShell(
        title=f'fMRI Voxel & Regressor Viewer, phys2cvr v{__version__}',
        theme='catppuccin-dark',
        dark_theme='catppuccin-dark',
        light_theme='catppuccin-light',
        remember_window_state=True,
        size=(1200, 850),
        min_size=(1200, 850),
        sidebar_mode='compact',
    ) as shell:
        # --- Status bar -----------------------------------------------------
        ssb = shell.statusbar
        SidebarToggle(
            shell=shell,
            collapse='compact',
            icon='layout-sidebar',
            expand_icon='layout-sidebar',
            collapse_icon='layout-sidebar-inset',
            parent=ssb,
        )
        ssb.add_text(f'phys2cvr v{__version__}')
        ssb.add_spacer()
        bs.ThemeToggle(parent=ssb)

        # --- Pages ----------------------------------------------------------
        with shell.page_nav() as nav:
            with nav.add_page(
                'voxelviewer', text='Voxel viewer', icon='clipboard2-pulse'
            ):
                voxelviewer.build()
            with nav.add_page(
                'workflowsetter', text='Run phys2cvr', icon='terminal-fill'
            ):
                build_workflow_setter()
            with nav.add_page('diffviewer', text='Comparison', icon='subtract'):
                diffviewer = DiffViewer()
                diffviewer.build()
            with nav.add_page(
                'logviewer', text='Log viewer', icon='file-earmark-richtext-fill'
            ):
                build_log_viewer()
        shell.navigate('voxelviewer')
    return shell


def _main(argv=None):
    args = _get_parser().parse_args(argv)
    if args.lag_range is not None:
        if args.lag_range[0] == args.lag_range[1]:
            raise ValueError(
                'Specified lag range extremes are the same. Provide two different '
                'values or leave blank for automatic estimation.'
            )
        if args.lag_range[0] > args.lag_range[1]:
            raise ValueError('Specified minimum lag range is higher than maximum.')

    voxelviewer = VoxelViewer(fs=args.fs, lag_range=args.lag_range)
    shell = build_app(voxelviewer)

    try:
        if args.result_folder:
            load_result_folder(args.result_folder, voxelviewer)

        if args.timeseries:
            voxelviewer.set_timeseries(args.timeseries)
        for path, setter in (
            (args.mask, voxelviewer.set_mask),
            (args.cvr, voxelviewer.set_cvr),
            (args.lag, voxelviewer.set_lag),
            (args.regressor, voxelviewer.set_regressors),
            (args.r2, voxelviewer.set_r2),
        ):
            if path:
                setter(path)
        voxelviewer.refresh()
    except Exception as err:  # noqa: BLE001
        LGR.error('Could not load input: %s', err)
        sys.exit(1)

    shell.run()


if __name__ == '__main__':
    _main(sys.argv[1:])
