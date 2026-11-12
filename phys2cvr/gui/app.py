#!/usr/bin/env python3
"""Application shell for the phys2cvr GUI."""

import argparse
import logging
import sys

import bootstack as bs

from phys2cvr import __version__
from phys2cvr.gui.pages.voxelviewer import VoxelViewer

logger = logging.getLogger(__name__)


def _get_parser():
    """Create the viewer argument parser."""
    parser = argparse.ArgumentParser(description='phys2cvr viewer.')
    parser.add_argument('nii_file', nargs='?', help='Path to 4D fMRI NIfTI file')
    parser.add_argument('--matrix_file', help='Path to regressor matrix text/mat file')
    parser.add_argument('--cvr', help='Optional path to 3D CVR NIfTI map')
    parser.add_argument('--lag', help='Optional path to 3D Lag NIfTI map')
    parser.add_argument('--mask', help='Optional path to 3D brain mask NIfTI file')
    parser.add_argument(
        '--fs',
        type=float,
        default=None,
        help='Regressor sampling frequency in Hz (default: 1.0)',
    )
    parser.add_argument(
        '--lag_range',
        type=float,
        nargs=2,
        default=None,
        metavar=('START', 'END'),
        help='Optional explicit lag range tuple (start, end) in seconds',
    )
    return parser


def build_app(viewer):
    """Build the application shell and its viewer page."""
    with bs.AppShell(
        title=f'fMRI Voxel & Regressor Viewer, phys2cvr v{__version__}',
        theme='catppuccin-dark',
        dark_theme='catppuccin-dark',
        light_theme='catppuccin-light',
        remember_window_state=True,
        size=(1200, 850),
        min_size=(1200, 850),
    ) as shell:
        with shell.add_toolbar() as bar:
            with bar.add_menu('File') as menu:
                menu.add_action('Exit', shortcut='Mod+Q', on_click=shell.close)
            bar.add_spacer()
            bar.add_theme_toggle()
        with shell.page_nav() as nav:
            with nav.add_page('viewer', text='Voxel viewer', icon='activity'):
                viewer.build()
        shell.navigate('viewer')
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

    viewer = VoxelViewer(fs=args.fs, lag_range=args.lag_range)
    shell = build_app(viewer)

    try:
        if args.nii_file:
            viewer.set_timeseries(args.nii_file)
        for path, setter in (
            (args.mask, viewer.set_mask),
            (args.cvr, viewer.set_cvr),
            (args.lag, viewer.set_lag),
            (args.matrix_file, viewer.set_regressors),
        ):
            if path:
                setter(path)
        viewer.refresh()
    except Exception as err:  # noqa: BLE001
        logger.error('Could not load input: %s', err)
        sys.exit(1)

    shell.run()


if __name__ == '__main__':
    _main(sys.argv[1:])
