#!/usr/bin/env python3
"""Voxel viewer for phys2cvr, built with bootstack."""

import argparse
import logging
import sys

import bootstack as bs
import matplotlib as mpl
import numpy as np
from matplotlib.figure import Figure
from scipy.stats import zscore

from phys2cvr import __version__
from phys2cvr.io import EXT_ARRAY, EXT_NIMG, load_array, load_nifti_get_mask
from phys2cvr.utils import estimate_lag_range

logger = logging.getLogger(__name__)

NIFTI_TYPES = [
    ('NIfTI files', ' '.join(f'*{ext}' for ext in EXT_NIMG)),
    ('All files', '*.*'),
]

ARRAY_TYPES = [
    ('Array files', ' '.join(f'*{ext}' for ext in EXT_ARRAY)),
    ('All files', '*.*'),
]


def _get_parser():
    """
    Parse command line inputs for this function.

    Returns
    -------
    parser.parse_args() : argparse dict

    """
    parser = argparse.ArgumentParser(description='phys2cvr viewer.')
    parser.add_argument(
        'nii_file', nargs='?', default=None, help='Path to 4D fMRI NIfTI file'
    )
    parser.add_argument(
        '--matrix_file',
        default=None,
        help='Path to regressor matrix text/mat file',
    )
    parser.add_argument('--cvr', default=None, help='Optional path to 3D CVR NIfTI map')
    parser.add_argument('--lag', default=None, help='Optional path to 3D Lag NIfTI map')
    parser.add_argument(
        '--mask', default=None, help='Optional path to 3D brain mask NIfTI file'
    )
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


def load_timeseries(path):
    """Load a NIfTI file and return data, spatial voxel sizes and TR (seconds)."""
    data, _, img = load_nifti_get_mask(path)
    zooms = img.header.get_zooms()
    pixdim = tuple(zooms[:3]) if len(zooms) >= 3 else (1.0, 1.0, 1.0)
    tr = float(zooms[3]) if len(zooms) > 3 else 1.0
    if tr <= 0 or np.isnan(tr):
        tr = 1.0
    return data, pixdim, tr


def load_matrix(path):
    """Load a regressor matrix (lags x timepoints); a single regressor is 1 row."""
    return np.atleast_2d(np.loadtxt(path))


class VoxelViewer:
    """Interactive viewer for fMRI voxels, CVR/lag maps and shifted regressors."""

    def __init__(self, fs=None, lag_range=None):
        self.fs = fs if fs else 1.0
        self.lag_range = lag_range

        self.timeseries = None
        self.data_types = [
            'meanvol',
            'mask',
            'cvr',
            'lag',
            'raw_cvr',
            'raw_lag',
        ]
        self.data = {k: None for k in self.data_types}
        self.regressor_matrix = None
        self.pixdim = (1.0, 1.0, 1.0)
        self.tr = 1.0
        self.time = np.arange(1.0)
        self.nx = self.ny = self.nz = self.nt = 1

        self.curr = [0, 0, 0]
        self.curr_shift_sec = 0.0
        self.min_delay_sec = -1.0
        self.max_delay_sec = 1.0
        self.lag_units = 's'
        self._syncing = False
        self.axes = {}

    # ----- data handling -----

    def set_timeseries(self, path):
        self.timeseries, self.pixdim, self.tr = load_timeseries(path)
        self.nx, self.ny, self.nz, self.nt = self.timeserie.shape
        self.time = np.arange(self.nt) * self.tr
        self.data['meanvol'] = np.mean(self.timeserie, axis=-1)
        self.curr = [self.nx // 2, self.ny // 2, self.nz // 2]
        self._apply_mask()

    def set_regressors(self, path):
        self.regressor_matrix = load_array(path, shape='prerectangle')

    def set_cvr(self, path):
        self.data['raw_cvr'], _, _ = load_nifti_get_mask(path)
        self._apply_mask(['cvr'])

    def set_lag(self, path):
        self.data['raw_lag'], _, _ = load_nifti_get_mask(path)
        self._apply_mask(['lag'])
        _, lag_max, lag_min = estimate_lag_range(self.data['lag'], path)
        self.lag_range = [lag_max, lag_min]
        self._update_delay_range()

    def set_mask(self, path):
        _, self.data['mask'], _ = load_nifti_get_mask(path, is_mask=True)
        self._apply_mask()

    def _apply_mask(self, data_types=('cvr', 'lag')):
        """Mask copies of the raw maps, so that masks never compound."""
        for k in data_types:
            self.data[k] = (
                self.data[f'raw_{k}']
                if None in (self.data[f'raw_{k}'], self.data['mask'])
                else self._masked(self.data[f'raw_{k}'])
            )

    def _masked(self, raw):
        data = raw.copy()
        data[~self.data['mask']] = np.nan
        return data

    def _update_delay_range(self):
        """Slider bounds, in seconds; regressor rows span the whole range."""
        if self.lag_range is not None:
            self.min_delay_sec, self.max_delay_sec = self.lag_range
            self.lag_units = 's'
        elif self.regressor_matrix is not None:
            half_lags = (self.regressor_matrix.shape[0] - 1) / 2.0
            self.min_delay_sec = -half_lags / self.fs
            self.max_delay_sec = half_lags / self.fs
            self.lag_units = 'idx'
        else:
            self.min_delay_sec, self.max_delay_sec = -1.0, 1.0
            self.lag_units = 'idx'

        self.slider.min_value = self.min_delay_sec
        self.slider.max_value = self.max_delay_sec
        self.slider.disabled = self.regressor_matrix is None

    def _regressor_row(self, shift_sec):
        """Return the regressor shifted by `shift_sec` and the actual shift used."""
        num_lags = self.regressor_matrix.shape[0]
        idx = int(round((shift_sec - self.min_delay_sec) * self.fs))
        idx = int(np.clip(idx, 0, num_lags - 1))
        return self.regressor_matrix[idx], idx / self.fs + self.min_delay_sec

    # ----- UI -----

    def build(self):
        def loader(title, file_types, setter):
            return lambda: self._load(title, file_types, setter)

        with bs.App(
            title=f'fMRI Voxel & Regressor Viewer, phys2cvr v{__version__}',
            theme='catppuccin-dark',
            dark_theme='catppuccin-dark',
            light_theme='catppuccin-light',
            remember_window_state=True,
            size=(1200, 850),
            min_size=(1200, 850),
        ) as app:
            with app.add_toolbar() as bar:
                with bar.add_menu('File') as menu:
                    menu.add_action(
                        'Load timeseries',
                        on_click=loader(
                            'Select 4D fMRI NIfTI File',
                            NIFTI_TYPES,
                            self.set_timeseries,
                        ),
                    )
                    menu.add_action(
                        'Load regressors',
                        on_click=loader(
                            'Select Regressor Matrix File',
                            ARRAY_TYPES,
                            self.set_regressors,
                        ),
                    )
                    menu.add_action(
                        'Load CVR map',
                        on_click=loader(
                            'Select 3D CVR Map NIfTI File', NIFTI_TYPES, self.set_cvr
                        ),
                    )
                    menu.add_action(
                        'Load lag map',
                        on_click=loader(
                            'Select 3D Lag Map NIfTI File', NIFTI_TYPES, self.set_lag
                        ),
                    )
                    menu.add_action(
                        'Load mask',
                        on_click=loader(
                            'Select 3D Brain Mask NIfTI File',
                            NIFTI_TYPES,
                            self.set_mask,
                        ),
                    )
                    menu.add_divider()
                    menu.add_action('Exit', shortcut='Mod+Q', on_click=app.close)
                bar.add_spacer()
                bar.add_theme_toggle()

            self.chart = bs.Chart(figure=Figure(), grow=True, horizontal='stretch')

            with bs.Row(gap=16, padding=10, horizontal='stretch', vertical='center'):
                self.x_field = bs.NumberField(0, label='X', min_value=0, width=6)
                self.y_field = bs.NumberField(0, label='Y', min_value=0, width=6)
                self.z_field = bs.NumberField(0, label='Z', min_value=0, width=6)
                with bs.Column(grow=True, horizontal='stretch'):
                    bs.Label(f'CO2 Delay ({self.lag_units})', font='caption')
                    self.slider = bs.Slider(
                        0.0,
                        min_value=self.min_delay_sec,
                        max_value=self.max_delay_sec,
                        horizontal='stretch',
                    )
                self.delay_label = bs.Label(f'+0.00 {self.lag_units}', width=10)

        for field in (self.x_field, self.y_field, self.z_field):
            field.on_change(lambda e: self._on_coord_change())
        self.slider.on_change(lambda e: self._on_slider_move())
        self.refresh()
        return app

    def _load(self, title, file_types, setter):
        path = bs.ask_open_file(title=title, file_types=file_types)
        if not path:
            return
        try:
            setter(path)
        except Exception as err:  # noqa: BLE001
            bs.alert(str(err), title='Could not load file')
            return
        self.refresh()

    def refresh(self):
        """Rebuild the figure and sync controls after the data changed."""
        self._update_delay_range()
        if self.timeseries is None:
            self.axes = {}
            self.chart.figure = Figure()
            return

        self.init_plots()
        self.chart.figure.canvas.mpl_connect('button_press_event', self.on_canvas_click)

        self._syncing = True
        try:
            self.x_field.max_value = self.nx - 1
            self.y_field.max_value = self.ny - 1
            self.z_field.max_value = self.nz - 1
        finally:
            self._syncing = False

        self._move_to(*self.curr)

    def init_plots(self):
        num_cols = 3 + (self.data['cvr'] is not None) + (self.data['lag'] is not None)

        fig = Figure(figsize=(4 * num_cols, 7), layout='constrained')
        gs = fig.add_gridspec(2, num_cols, height_ratios=[1.2, 1])

        dx, dy, dz = self.pixdim
        aspects = {'sag': dz / dy, 'cor': dz / dx, 'axi': dy / dx}
        self.axes = {}
        self.imgs = {}
        self.cross = {}

        def add_view(name, subplot, aspect, cmap, color, xlabel, ylabel, size):
            axis = self.axes[name] = fig.add_subplot(subplot)
            self.imgs[name] = axis.imshow(
                np.zeros(size[::-1]), cmap=cmap, origin='lower', aspect=aspect
            )
            self.cross[name] = (
                axis.axvline(0, color=color, lw=1),
                axis.axhline(0, color=color, lw=1),
            )
            axis.set_xlabel(xlabel)
            axis.set_ylabel(ylabel)
            return axis

        nx, ny, nz = self.nx, self.ny, self.nz
        add_view('sag', gs[0, 0], aspects['sag'], 'gray', 'r', 'Y', 'Z', (ny, nz))
        add_view('cor', gs[0, 1], aspects['cor'], 'gray', 'r', 'X', 'Z', (nx, nz))
        add_view('axi', gs[0, 2], aspects['axi'], 'gray', 'r', 'X', 'Y', (nx, ny))

        col = 3
        for name, data, cmap_name, color in (
            ('cvr', self.data['cvr'], 'inferno', 'cyan'),
            ('lag', self.data['lag'], 'viridis_r', 'red'),
        ):
            if data is None:
                continue
            cmap = mpl.colormaps[cmap_name].copy()
            cmap.set_bad(color='black')
            axis = add_view(
                name, gs[0, col], aspects['axi'], cmap, color, 'X', 'Y', (nx, ny)
            )
            axis.set_facecolor('black')
            fig.colorbar(self.imgs[name], ax=axis, fraction=0.046, pad=0.04)
            col += 1

        ts_ax = self.axes['ts'] = fig.add_subplot(gs[1, :])
        (self.line_voxel,) = ts_ax.plot([], [], label='Voxel TS', color='#1f77b4')
        (self.line_co2,) = ts_ax.plot(
            [], [], label='CO2 Regressor', color='#ff7f0e', linestyle='--'
        )
        ts_ax.set_xlabel('Time (s)')
        ts_ax.set_ylabel('Z-Score')
        ts_ax.grid(True, linestyle=':', alpha=0.5)

        self.chart.figure = fig

    def update_all(self, x, y, z, shift_sec):
        self.curr = [x, y, z]
        self.curr_shift_sec = shift_sec
        ax = self.axes

        views = {
            'sag': (self.data['meanvol'][x], y, z, f'Sagittal (X={x})'),
            'cor': (self.data['meanvol'][:, y], x, z, f'Coronal (Y={y})'),
            'axi': (self.data['meanvol'][:, :, z], x, y, f'Axial (Z={z})'),
        }
        for name, data, label, fmt in (
            ('cvr', self.data['cvr'], 'CVR', '{:.2f}'),
            ('lag', self.data['lag'], 'Lag', '{:.2f}s'),
        ):
            if data is not None:
                val = data[x, y, z]
                title = (
                    f'{label}: {fmt.format(val)}'
                    if not np.isnan(val)
                    else f'{label}: Masked'
                )
                views[name] = (data[:, :, z], x, y, title)

        for name, (slice_data, cx, cy, title) in views.items():
            img = self.imgs[name]
            img.set_data(slice_data.T)
            img.autoscale()
            self.cross[name][0].set_xdata([cx, cx])
            self.cross[name][1].set_ydata([cy, cy])
            ax[name].set_title(title)

        ts_ax = ax['ts']
        self.line_voxel.set_data(self.time, zscore(self.timeseries[x, y, z, :]))
        self.line_voxel.set_label(f'Voxel ({x}, {y}, {z})')

        actual_sec = shift_sec
        if self.regressor_matrix is not None:
            raw_co2, actual_sec = self._regressor_row(shift_sec)
            if len(raw_co2) != self.nt:
                reg_time = np.linspace(0, self.time[-1], len(raw_co2))
                raw_co2 = np.interp(self.time, reg_time, raw_co2)
            self.line_co2.set_data(self.time, zscore(raw_co2))
            self.line_co2.set_label(f'CO2 Shift ({actual_sec:+.2f}{self.lag_units})')
        else:
            self.line_co2.set_data([], [])
            self.line_co2.set_label('_nolegend_')

        ts_ax.set_title(f'Timecourse Comparison — Voxel ({x}, {y}, {z})')
        ts_ax.legend(loc='upper right')
        ts_ax.set_xlim(0, self.time[-1] if self.nt > 1 else 1)
        ts_ax.relim()
        ts_ax.autoscale_view(scalex=False, scaley=True)

        self._syncing = True
        try:
            self.x_field.value = x
            self.y_field.value = y
            self.z_field.value = z
            self.slider.value = float(
                np.clip(actual_sec, self.min_delay_sec, self.max_delay_sec)
            )
            self.delay_label.text = f'{actual_sec:+.2f} {self.lag_units}'
        finally:
            self._syncing = False

        self.chart.draw()

    # ----- callbacks -----

    def _move_to(self, x, y, z):
        """Move to a voxel, following its lag value when there is one."""
        x = int(np.clip(x, 0, self.nx - 1))
        y = int(np.clip(y, 0, self.ny - 1))
        z = int(np.clip(z, 0, self.nz - 1))
        shift = self.curr_shift_sec
        if self.data['lag'] is not None and not np.isnan(self.data['lag'][x, y, z]):
            shift = float(self.data['lag'][x, y, z])
        self.update_all(x, y, z, shift)

    def on_canvas_click(self, event):
        names = {axis: name for name, axis in self.axes.items()}
        name = names.get(event.inaxes)
        if name is None or name == 'ts' or event.xdata is None or event.ydata is None:
            return

        x, y, z = self.curr
        u, v = int(round(event.xdata)), int(round(event.ydata))
        if name == 'sag':
            y, z = u, v
        elif name == 'cor':
            x, z = u, v
        else:
            x, y = u, v
        self._move_to(x, y, z)

    def _on_coord_change(self):
        if self._syncing or self.timeseries is None:
            return
        values = [f.value for f in (self.x_field, self.y_field, self.z_field)]
        if any(v is None for v in values):
            return
        self._move_to(*values)

    def _on_slider_move(self):
        if self._syncing or self.timeseries is None:
            return
        self.update_all(*self.curr, float(self.slider.value))


def _main(argv=None):
    args = _get_parser().parse_args(argv)

    if args.lag_range is not None:
        if args.lag_range[0] == args.lag_range[1]:
            raise ValueError(
                'Specified lag range extremes are the same. Provide two different '
                'values or leave blank for automatic estimation.'
            )
        elif args.lag_range[0] > args.lag_range[1]:
            raise ValueError('Specified minimum lag range is higher than maximum.')

    viewer = VoxelViewer(fs=args.fs, lag_range=args.lag_range)
    app = viewer.build()

    # The mask goes first, so that it is applied to the maps
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
    except Exception as err:  # noqa: BLE001
        logger.error('Could not load input: %s', err)
        sys.exit(1)
    viewer.refresh()

    app.run()


if __name__ == '__main__':
    _main(sys.argv[1:])
