#!/usr/bin/env python3
"""Voxel viewer for phys2cvr, built with bootstack."""

import bootstack as bs
import matplotlib as mpl
import numpy as np
from matplotlib.figure import Figure
from scipy.stats import zscore

from phys2cvr.io import EXT_ARRAY, EXT_NIMG, load_array, load_nifti_get_mask
from phys2cvr.utils import estimate_lag_range

NIFTI_TYPES = [
    ('NIfTI files', ' '.join(f'*{ext}' for ext in EXT_NIMG)),
    ('All files', '*.*'),
]

ARRAY_TYPES = [
    ('Array files', '*.1D ' + ' '.join(f'*{ext}' for ext in EXT_ARRAY)),
    ('All files', '*.*'),
]


def load_timeseries(path):
    """Load a NIfTI file and return data, spatial voxel sizes and TR (seconds)."""
    data, _, img = load_nifti_get_mask(path, dim=4)
    zooms = img.header.get_zooms()
    pixdim = tuple(zooms[:3]) if len(zooms) >= 3 else (1.0, 1.0, 1.0)
    tr = float(zooms[3]) if len(zooms) > 3 else 1.0
    if tr <= 0 or np.isnan(tr):
        tr = 1.0
    return data, pixdim, tr


class VoxelViewer:
    """Interactive viewer for fMRI voxels, CVR/lag maps and shifted regressors."""

    def __init__(self, fs=None, lag_range=None):
        self.fs = fs if fs else 1.0
        self.lag_range = lag_range
        self.lag_step = None

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

    def _masked(self, raw):
        data = raw.copy()
        data[~self.data['mask']] = np.nan
        return data

    def _apply_mask(self, data_types=('cvr', 'lag')):
        """Mask copies of the raw maps, so that masks never compound."""
        for k in data_types:
            self.data[k] = (
                self.data[f'raw_{k}']
                if self.data[f'raw_{k}'] is None or self.data['mask'] is None
                else self._masked(self.data[f'raw_{k}'])
            )

    def _update_delay_range(self):
        """Slider bounds, in seconds; regressor rows span the whole range."""
        if self.lag_range is not None:
            self.min_delay_sec, self.max_delay_sec = self.lag_range
            self.lag_units = 's'
            if self.lag_step is not None:
                self.slider.step = self.lag_step
        elif self.regressor_matrix is not None:
            half_lags = (self.regressor_matrix.shape[0] - 1) / 2.0
            self.min_delay_sec = -half_lags / self.fs
            self.max_delay_sec = half_lags / self.fs
            self.lag_units = 'idx'
            self.slider.step = 1
        else:
            self.min_delay_sec, self.max_delay_sec = -1.0, 1.0
            self.lag_units = 'idx'
            self.slider.step = 1

        self.lagmin_field.value = self.min_delay_sec
        self.lagmax_field.value = self.max_delay_sec
        self.slider.min_value = self.min_delay_sec
        self.slider.max_value = self.max_delay_sec
        self.slider.disabled = self.regressor_matrix is None

    def set_timeseries(self, path):
        """Load the functional timeseries and initialize the selected voxel."""
        self.timeseries, self.pixdim, self.tr = load_timeseries(path)
        self.nx, self.ny, self.nz, self.nt = self.timeseries.shape
        self.time = np.arange(self.nt) * self.tr
        self.data['meanvol'] = np.mean(self.timeseries, axis=-1)
        self.curr = [self.nx // 2, self.ny // 2, self.nz // 2]
        self._apply_mask()

    def set_regressors(self, path):
        """Load the shifted regressor matrix."""
        self.regressor_matrix = load_array(path, shape='prerectangle')

    def set_cvr(self, path):
        """Load and mask the CVR map."""
        self.data['raw_cvr'], _, _ = load_nifti_get_mask(path)
        self._apply_mask(['cvr'])

    def set_lag(self, path):
        """Load the lag map and estimate the available lag range."""
        self.data['raw_lag'], _, _ = load_nifti_get_mask(path)
        self._apply_mask(['lag'])
        finite_lags = self.data['raw_lag'][np.isfinite(self.data['raw_lag'])]
        self.lag_step, lag_max, lag_min = estimate_lag_range(finite_lags, path)
        self.lag_range = [lag_min, lag_max]
        self._update_delay_range()

    def set_mask(self, path):
        """Load a mask and apply it to the maps."""
        _, self.data['mask'], _ = load_nifti_get_mask(path, is_mask=True)
        self._apply_mask()

    def _regressor_row(self, shift_sec):
        """Return the regressor shifted by `shift_sec` and the actual shift used."""
        num_lags = self.regressor_matrix.shape[0]
        idx = int(round((shift_sec - self.min_delay_sec) * self.fs))
        idx = int(np.clip(idx, 0, num_lags - 1))
        return self.regressor_matrix[idx], idx / self.fs + self.min_delay_sec

    # ----- UI -----

    def build(self):
        """Build the viewer inside the current page."""

        def loader(title, file_types, setter):
            return lambda: self._load(title, file_types, setter)

        with bs.Toolbar(horizontal='stretch', surface='chrome') as bar:
            self.toolbar = bar
            with bar.add_menu('Load') as menu:
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
        self.chart = bs.Chart(figure=Figure(), grow=True, horizontal='stretch')

        with bs.Row(gap=16, padding=10, horizontal='stretch', vertical='center'):
            self.x_field = bs.NumberField(0, label='X', min_value=0, width=4)
            self.y_field = bs.NumberField(0, label='Y', min_value=0, width=4)
            self.z_field = bs.NumberField(0, label='Z', min_value=0, width=4)
            self.lagmax_field = bs.NumberField(0, label='Lag max', min_value=0, width=4)
            self.lagmin_field = bs.NumberField(0, label='Lag min', min_value=0, width=4)
            with bs.Column(grow=True, horizontal='stretch'):
                self.CO2_label = bs.Label(
                    f'CO2 Delay ({self.lag_units})', font='caption'
                )
                self.slider = bs.Slider(
                    0.0,
                    min_value=self.min_delay_sec,
                    max_value=self.max_delay_sec,
                    horizontal='stretch',
                )
            self.delay_label = bs.Label(f'+0.00 {self.lag_units}', width=10)

        for field in (self.x_field, self.y_field, self.z_field):
            field.on_change(lambda e: self._on_coord_change())
        for field in (self.lagmax_field, self.lagmin_field):
            field.on_change(lambda e: self._on_lagrange_change())
        self.slider.on_change(lambda e: self._on_slider_move())
        self.refresh()

    def _load(self, title, file_types, setter):
        path = bs.ask_open_file(title=title, file_types=file_types)
        if not path:
            return
        try:
            setter(path)
            self.refresh()
        except Exception as err:  # noqa: BLE001
            bs.alert(str(err), title='Could not load file')

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
        """Create the anatomical, map and timecourse plots."""
        has_cvr = self.data['cvr'] is not None
        has_lag = self.data['lag'] is not None

        num_cols = 5

        fig = Figure(figsize=(4 * num_cols, 9), layout='constrained')
        gs = fig.add_gridspec(3, num_cols, height_ratios=[1.2, 1.2, 1])
        gs = fig.add_gridspec(
            3, num_cols, height_ratios=[1.2, 1.2, 1], wspace=0.05, hspace=0.2
        )

        dx, dy, dz = self.pixdim
        aspects = {
            'sag': dz / dy,
            'cor': dz / dx,
            'axi': dy / dx,
        }
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
            axis.set_xticks([])
            axis.set_yticks([])
            axis.set_xticklabels([])
            axis.set_yticklabels([])
            return axis

        nx, ny, nz = self.nx, self.ny, self.nz

        # --- Row 0: CVR Maps (Sag, Cor, Axi) + Average Volume (cols 3 & 4) ---
        if has_cvr:
            cmap_cvr = mpl.colormaps['inferno'].copy()
            cmap_cvr.set_bad(color='black')
            add_view(
                'cvr_sag',
                gs[0, 0],
                aspects['sag'],
                cmap_cvr,
                'cyan',
                'Y',
                'Z',
                (ny, nz),
            )
            add_view(
                'cvr_cor',
                gs[0, 1],
                aspects['cor'],
                cmap_cvr,
                'cyan',
                'X',
                'Z',
                (nx, nz),
            )
            add_view(
                'cvr_axi',
                gs[0, 2],
                aspects['axi'],
                cmap_cvr,
                'cyan',
                'X',
                'Y',
                (nx, ny),
            )
            for name in ('cvr_sag', 'cvr_cor', 'cvr_axi'):
                self.axes[name].set_facecolor('black')
            fig.colorbar(
                self.imgs['cvr_axi'], ax=self.axes['cvr_axi'], fraction=0.046, pad=0.04
            )

        add_view('mean_cor', gs[0, 3], aspects['cor'], 'gray', 'r', 'X', 'Z', (nx, nz))
        add_view('mean_axi', gs[0, 4], aspects['axi'], 'gray', 'r', 'X', 'Y', (nx, ny))

        # --- Row 1: Lag Maps (Sag, Cor, Axi) + Average Volume (col 3) ---
        if has_lag:
            cmap_lag = mpl.colormaps['viridis_r'].copy()
            cmap_lag.set_bad(color='black')
            add_view(
                'lag_sag', gs[1, 0], aspects['sag'], cmap_lag, 'red', 'Y', 'Z', (ny, nz)
            )
            add_view(
                'lag_cor', gs[1, 1], aspects['cor'], cmap_lag, 'red', 'X', 'Z', (nx, nz)
            )
            add_view(
                'lag_axi', gs[1, 2], aspects['axi'], cmap_lag, 'red', 'X', 'Y', (nx, ny)
            )
            for name in ('lag_sag', 'lag_cor', 'lag_axi'):
                self.axes[name].set_facecolor('black')
            fig.colorbar(
                self.imgs['lag_axi'], ax=self.axes['lag_axi'], fraction=0.046, pad=0.04
            )

        add_view('mean_sag', gs[1, 3], aspects['sag'], 'gray', 'r', 'Y', 'Z', (ny, nz))

        # --- Row 2: Timeseries plot spanning all columns ---
        ts_ax = self.axes['ts'] = fig.add_subplot(gs[2, :])
        (self.line_voxel,) = ts_ax.plot([], [], label='Voxel TS', color='#1f77b4')
        (self.line_co2,) = ts_ax.plot(
            [], [], label='CO2 Regressor', color='#ff7f0e', linestyle='--'
        )
        ts_ax.set_xlabel('Time (s)')
        ts_ax.set_ylabel('Z-Score')
        ts_ax.grid(True, linestyle=':', alpha=0.5)

        self.chart.figure = fig

    def update_all(self, x, y, z, shift_sec):
        """Update plots and controls for a voxel and regressor delay."""
        self.curr = [x, y, z]
        self.curr_shift_sec = shift_sec
        ax = self.axes

        mean_vol = self.data['meanvol']
        cvr_data = self.data['cvr']
        lag_data = self.data['lag']

        # --- Define views for all axes across the new grid ---
        views = {
            # Average volume slices
            'mean_cor': (mean_vol[:, y, :], x, z, f'Coronal (Y={y})'),
            'mean_axi': (mean_vol[:, :, z], x, y, f'Axial (Z={z})'),
            'mean_sag': (mean_vol[x, :, :], y, z, f'Sagittal (X={x})'),
        }

        # CVR 3-axis views (Row 0)
        if cvr_data is not None:
            val_cvr = cvr_data[x, y, z]
            cvr_title = (
                f'CVR: {val_cvr:.2f}' if not np.isnan(val_cvr) else 'CVR: Masked'
            )
            views.update(
                {
                    'cvr_sag': (cvr_data[x, :, :], y, z, f'Sagittal (X={x})'),
                    'cvr_cor': (cvr_data[:, y, :], x, z, f'Coronal (Y={y})'),
                    'cvr_axi': (cvr_data[:, :, z], x, y, cvr_title),
                }
            )

        # Lag 3-axis views (Row 1)
        if lag_data is not None:
            val_lag = lag_data[x, y, z]
            lag_title = (
                f'Lag: {val_lag:.2f}s' if not np.isnan(val_lag) else 'Lag: Masked'
            )
            views.update(
                {
                    'lag_sag': (lag_data[x, :, :], y, z, f'Sagittal (X={x})'),
                    'lag_cor': (lag_data[:, y, :], x, z, f'Coronal (Y={y})'),
                    'lag_axi': (lag_data[:, :, z], x, y, lag_title),
                }
            )

        # Update each image slice, crosshairs, and titles
        for name, (slice_data, cx, cy, title) in views.items():
            if name not in self.imgs:
                continue
            img = self.imgs[name]
            img.set_data(slice_data.T)
            img.autoscale()
            self.cross[name][0].set_xdata([cx, cx])
            self.cross[name][1].set_ydata([cy, cy])
            ax[name].set_title(title)

        # --- Timeseries Plot (Row 2) ---
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

        # --- Sync Controls ---
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
        """Select a voxel from any of the displayed anatomical or map slices."""
        names = {axis: name for name, axis in self.axes.items()}
        name = names.get(event.inaxes)
        if name is None or name == 'ts' or event.xdata is None or event.ydata is None:
            return

        x, y, z = self.curr
        u, v = int(round(event.xdata)), int(round(event.ydata))

        # Determine coordinates based on which view was clicked
        if name in ('cvr_sag', 'lag_sag', 'mean_sag'):
            y, z = u, v
        elif name in ('cvr_cor', 'lag_cor', 'mean_cor'):
            x, z = u, v
        elif name in ('cvr_axi', 'lag_axi', 'mean_axi'):
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
