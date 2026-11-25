#!/usr/bin/env python3
"""Application shell for the phys2cvr GUI."""

import argparse
import logging
import sys

import bootstack as bs
from bootstack.widgets.sidebar_toggle import SidebarToggle

from phys2cvr import __version__
from phys2cvr.gui.io import load_result_folder
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
        default=None,
        help='Load a phys2cvr output folder and let the GUI load everything for you!',
    )
    return parser


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
                logviewer = build_log_viewer()
        voxelviewer.linked_diffviewer = diffviewer
        voxelviewer.linked_logviewer = logviewer
        diffviewer.linked_logviewer = logviewer
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
            load_result_folder(
                args.result_folder,
                voxelviewer,
                diffviewer=voxelviewer.linked_diffviewer,
                logviewer=voxelviewer.linked_logviewer,
            )

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
