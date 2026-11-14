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

    def __init__(self, fs=None, lag_range=None, cvr_range=None):
        self.fs = fs if fs else 1.0
        self.cvr_range = cvr_range
        self.lag_range = lag_range
        self.lag_step = None

        self.timeseries = None
        self.data_types = [
            'meanvol',
            'mask',
            'cvr',
            'lag',
            'r2',
            'raw_cvr',
            'raw_lag',
            'raw_r2',
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

        self._syncing = True
        try:
            self.lagmin_field.value = self.min_delay_sec
            self.lagmax_field.value = self.max_delay_sec
            self.slider.min_value = self.min_delay_sec
            self.slider.max_value = self.max_delay_sec
            self.slider.disabled = self.regressor_matrix is None
        finally:
            self._syncing = False

    def set_timeseries(self, path):
        """Load the functional timeseries and initialize the selected voxel."""
        self.timeseries, self.pixdim, self.tr = load_timeseries(path)
        self.nx, self.ny, self.nz, self.nt = self.timeseries.shape
        self.time = np.arange(self.nt) * self.tr
        self.data['meanvol'] = np.mean(self.timeseries, axis=-1)
        self.curr = [self.nx // 2, self.ny // 2, self.nz // 2]

        # If a mask was not provided yet, generate it from non-zero mean volume voxels
        if self.data['mask'] is None:
            self.data['mask'] = self.data['meanvol'] != 0

        self._apply_mask()

    def set_regressors(self, path):
        """Load the shifted regressor matrix."""
        self.regressor_matrix = load_array(path, shape='prerectangle')

    def set_cvr(self, path):
        """Load and mask the CVR map."""
        self.data['raw_cvr'], _, _ = load_nifti_get_mask(path)
        self._apply_mask(['cvr'])
        cvr_vals = self.data['cvr'][np.isfinite(self.data['cvr'])]
        self.cvr_range = (
            np.percentile(cvr_vals, [1, 99]) if cvr_vals.size > 0 else (None, None)
        )

    def set_lag(self, path):
        """Load the lag map and estimate the available lag range."""
        self.data['raw_lag'], _, _ = load_nifti_get_mask(path)
        self._apply_mask(['lag'])
        finite_lags = self.data['raw_lag'][np.isfinite(self.data['raw_lag'])]
        self.lag_step, lag_max, lag_min = estimate_lag_range(finite_lags, path)
        self.lag_range = (lag_min, lag_max)
        self._update_delay_range()

    def set_r2(self, path):
        """Load the 4D R² map volume."""
        self.data['raw_r2'], _, _ = load_nifti_get_mask(path, dim=4)
        if self.data['mask'] is not None:
            self.data['r2'] = self.data['raw_r2'].copy()
            # Mask spatial dimensions while keeping 4th dimension intact
            self.data['r2'][~self.data['mask']] = np.nan
        else:
            self.data['r2'] = self.data['raw_r2']

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
                    'Load result folder',
                    on_click=lambda: self._load_folder(),
                )
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
                    'Load R2 map',
                    on_click=loader(
                        'Select 4D R2 Map NIfTI File', NIFTI_TYPES, self.set_r2
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

    def _load_folder(self):
        folder_path = bs.ask_directory(title='Select phys2cvr Result Folder')
        if not folder_path:
            return
        try:
            from phys2cvr.gui.app import load_result_folder

            load_result_folder(folder_path, self)
            self.refresh()
        except Exception as err:  # noqa: BLE001
            bs.alert(str(err), title='Could not load result folder')

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
        """Create anatomical, CVR, Lag, and timecourse plots with mean volume stacked in the 4th column."""
        fig = Figure(figsize=(16, 9), layout='constrained')
        gs = fig.add_gridspec(3, 4, height_ratios=[1.2, 1.2, 1])

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
            axis.set_xticks([])
            axis.set_yticks([])
            axis.set_xticklabels([])
            axis.set_yticklabels([])
            return axis

        nx, ny, nz = self.nx, self.ny, self.nz

        # --- Row 0: CVR Maps (Sag, Cor, Axi) ---
        cmap_cvr = mpl.colormaps['inferno'].copy()
        cmap_cvr.set_bad(color='black')

        if self.cvr_range is not None:
            cvr_bound = max(abs(self.cvr_range[0]), abs(self.cvr_range[1]))
        elif self.data['cvr'] is not None:
            cvr_vals = self.data['cvr'][np.isfinite(self.data['cvr'])]
            cvr_bound = np.max(np.abs(cvr_vals)) if cvr_vals.size > 0 else 1.0
        else:
            cvr_bound = 1.0

        add_view(
            'cvr_sag', gs[0, 0], aspects['sag'], cmap_cvr, 'cyan', 'Y', 'Z', (ny, nz)
        )
        add_view(
            'cvr_cor', gs[0, 1], aspects['cor'], cmap_cvr, 'cyan', 'X', 'Z', (nx, nz)
        )
        add_view(
            'cvr_axi', gs[0, 2], aspects['axi'], cmap_cvr, 'cyan', 'X', 'Y', (nx, ny)
        )

        for name in ('cvr_sag', 'cvr_cor', 'cvr_axi'):
            self.axes[name].set_facecolor('black')
            self.imgs[name].set_clim(vmin=-cvr_bound, vmax=cvr_bound)

        fig.colorbar(
            self.imgs['cvr_axi'], ax=self.axes['cvr_axi'], fraction=0.046, pad=0.04
        )

        # --- Row 1: Lag Maps (Sag, Cor, Axi) ---
        cmap_lag = mpl.colormaps['viridis_r'].copy()
        cmap_lag.set_bad(color='black')

        if self.lag_range is not None:
            lag_bound = max(abs(self.lag_range[0]), abs(self.lag_range[1]))
        elif self.data['lag'] is not None:
            lag_vals = self.data['lag'][np.isfinite(self.data['lag'])]
            lag_bound = np.max(np.abs(lag_vals)) if lag_vals.size > 0 else 1.0
        else:
            lag_bound = 1.0

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
            self.imgs[name].set_clim(vmin=-lag_bound, vmax=lag_bound)

        fig.colorbar(
            self.imgs['lag_axi'], ax=self.axes['lag_axi'], fraction=0.046, pad=0.04
        )

        # --- Mean Volume: 3 axes stacked vertically in column 3, spanning rows 0 and 1, resized (2/3 scale) ---
        gs_mean = gs[:2, 3].subgridspec(3, 1, hspace=0.15)
        add_view(
            'mean_sag', gs_mean[0], aspects['sag'], 'gray', 'r', 'Y', 'Z', (ny, nz)
        )
        add_view(
            'mean_cor', gs_mean[1], aspects['cor'], 'gray', 'r', 'X', 'Z', (nx, nz)
        )
        add_view(
            'mean_axi', gs_mean[2], aspects['axi'], 'gray', 'r', 'X', 'Y', (nx, ny)
        )

        # --- Row 2: Timeseries plot (left, 4/5 width) and R2 curve plot (right, 1/5 width) ---
        gs_row2 = gs[2, :].subgridspec(1, 2, width_ratios=[4, 1], wspace=0.1)

        ts_ax = self.axes['ts'] = fig.add_subplot(gs_row2[0])
        (self.line_voxel,) = ts_ax.plot([], [], label='Voxel TS', color='#1f77b4')
        (self.line_co2,) = ts_ax.plot(
            [], [], label='CO2 Regressor', color='#ff7f0e', linestyle='--'
        )
        ts_ax.set_title('Timecourse Comparison — Voxel (0, 0, 0)')
        ts_ax.set_xlabel('Time (s)')
        ts_ax.set_ylabel('Z-Score')
        ts_ax.set_xlim(0, 1)
        ts_ax.legend(loc='upper right')
        ts_ax.grid(True, linestyle=':', alpha=0.5)

        r2_ax = self.axes['r2'] = fig.add_subplot(gs_row2[1])
        (self.line_r2,) = r2_ax.plot([], [], label='R²', color='#2ca02c')
        self.r2_bar = r2_ax.axvline(0, color='red', linestyle=':', lw=1.5)
        r2_ax.set_title('R²: N/A')
        r2_ax.set_xlabel(
            'Lag'
            if self.data['lag'] is not None or self.lag_range is not None
            else 'Index'
        )
        r2_ax.set_ylabel('R²')
        r2_ax.grid(True, linestyle=':', alpha=0.5)

        # Format y-axis ticks to always show exactly 2 decimals
        r2_ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.2f'))

        self.chart.figure = fig
        self.chart.figure.canvas.mpl_connect('button_press_event', self.on_canvas_click)

        self.chart.figure = fig

    def update_all(self, x, y, z, shift_sec):
        """Update plots and controls for a voxel and regressor delay."""
        if self._syncing or not hasattr(self, 'imgs') or not self.imgs:
            return
        self.curr = [x, y, z]
        self.curr_shift_sec = shift_sec
        ax = self.axes

        views = {
            'mean_cor': (
                self.data['meanvol'][:, y, :]
                if self.data['meanvol'] is not None
                else np.zeros((self.nx, self.nz)),
                x,
                z,
                f'Coronal (Y={y})',
            ),
            'mean_axi': (
                self.data['meanvol'][:, :, z]
                if self.data['meanvol'] is not None
                else np.zeros((self.nx, self.ny)),
                x,
                y,
                f'Axial (Z={z})',
            ),
            'mean_sag': (
                self.data['meanvol'][x, :, :]
                if self.data['meanvol'] is not None
                else np.zeros((self.ny, self.nz)),
                y,
                z,
                f'Sagittal (X={x})',
            ),
        }

        if self.data['cvr'] is not None:
            val_cvr = self.data['cvr'][x, y, z]
            cvr_title = (
                f'CVR: {val_cvr:.2f}' if not np.isnan(val_cvr) else 'CVR: Masked'
            )
            views.update(
                {
                    'cvr_sag': (self.data['cvr'][x, :, :], y, z, f'Sagittal (X={x})'),
                    'cvr_cor': (self.data['cvr'][:, y, :], x, z, f'Coronal (Y={y})'),
                    'cvr_axi': (self.data['cvr'][:, :, z], x, y, cvr_title),
                }
            )
        else:
            views.update(
                {
                    'cvr_sag': (np.zeros((self.ny, self.nz)), y, z, 'CVR: Sagittal'),
                    'cvr_cor': (np.zeros((self.nx, self.nz)), x, z, 'CVR: Coronal'),
                    'cvr_axi': (np.zeros((self.nx, self.ny)), x, y, 'CVR: Axial'),
                }
            )

        if self.data['lag'] is not None:
            val_lag = self.data['lag'][x, y, z]
            lag_title = (
                f'Lag: {val_lag:.2f}s' if not np.isnan(val_lag) else 'Lag: Masked'
            )
            views.update(
                {
                    'lag_sag': (self.data['lag'][x, :, :], y, z, f'Sagittal (X={x})'),
                    'lag_cor': (self.data['lag'][:, y, :], x, z, f'Coronal (Y={y})'),
                    'lag_axi': (self.data['lag'][:, :, z], x, y, lag_title),
                }
            )
        else:
            views.update(
                {
                    'lag_sag': (np.zeros((self.ny, self.nz)), y, z, 'Lag: Sagittal'),
                    'lag_cor': (np.zeros((self.nx, self.nz)), x, z, 'Lag: Coronal'),
                    'lag_axi': (np.zeros((self.nx, self.ny)), x, y, 'Lag: Axial'),
                }
            )

        for name, (slice_data, cx, cy, title) in views.items():
            if name not in self.imgs:
                continue
            img = self.imgs[name]
            img.set_data(slice_data.T)
            img.autoscale()
            self.cross[name][0].set_xdata([cx, cx])
            self.cross[name][1].set_ydata([cy, cy])
            ax[name].set_title(title)

        ts_ax = ax['ts']
        if self.timeseries is not None:
            self.line_voxel.set_data(self.time, zscore(self.timeseries[x, y, z, :]))
            self.line_voxel.set_label(f'Voxel ({x}, {y}, {z})')
        else:
            self.line_voxel.set_data([], [])

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

        # --- R2 Profile Plot ---
        r2_ax = ax['r2']
        current_r2_val = 0.0
        if self.data['r2'] is not None:
            r2_voxel = self.data['r2'][x, y, z, :]
            num_r2 = len(r2_voxel)

            # Find closest R² value corresponding to current shift/slider
            if self.lag_range is not None:
                r2_x = np.linspace(self.lag_range[0], self.lag_range[1], num_r2)
                r2_xlabel = 'Lag (s)'
            elif self.regressor_matrix is not None:
                half_lags = (self.regressor_matrix.shape[0] - 1) / 2.0
                r2_x = np.linspace(-half_lags / self.fs, half_lags / self.fs, num_r2)
                r2_xlabel = 'Lag (s)'
            else:
                r2_x = np.arange(num_r2)
                r2_xlabel = 'Index'

            self.line_r2.set_data(r2_x, r2_voxel)
            r2_ax.set_xlabel(r2_xlabel)
            r2_ax.set_xlim(
                r2_x[0] if len(r2_x) > 0 else 0, r2_x[-1] if len(r2_x) > 0 else 1
            )
            r2_ax.relim()
            r2_ax.autoscale_view(scalex=False, scaley=True)

            # Get exact R² value at the slider/shift index closest to shift_sec
            idx_closest = int(np.argmin(np.abs(r2_x - shift_sec)))
            if 0 <= idx_closest < len(r2_voxel):
                current_r2_val = r2_voxel[idx_closest]
        else:
            self.line_r2.set_data([], [])

        self.r2_bar.set_xdata([shift_sec, shift_sec])
        r2_ax.set_title(
            f'R²: {current_r2_val:.4f}' if not np.isnan(current_r2_val) else 'R²: NaN'
        )

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
        """Move to a voxel, following its lag value or max R² index when applicable."""
        x = int(np.clip(x, 0, self.nx - 1))
        y = int(np.clip(y, 0, self.ny - 1))
        z = int(np.clip(z, 0, self.nz - 1))
        shift = self.curr_shift_sec

        # If a lag map is loaded and has a valid value at this voxel, use it
        if self.data['lag'] is not None and not np.isnan(self.data['lag'][x, y, z]):
            shift = float(self.data['lag'][x, y, z])
        # Otherwise, if lag map is not loaded but R² data is available, set shift to max R² index
        elif self.data['lag'] is None and self.data['r2'] is not None:
            r2_voxel = self.data['r2'][x, y, z, :]
            if r2_voxel.size > 0 and not np.all(np.isnan(r2_voxel)):
                max_idx = int(np.nanargmax(r2_voxel))
                num_r2 = len(r2_voxel)

                # Map the max R² index back to seconds (or index units) based on current settings
                if self.lag_range is not None:
                    r2_x = np.linspace(self.lag_range[0], self.lag_range[1], num_r2)
                    shift = float(r2_x[max_idx])
                elif self.regressor_matrix is not None:
                    half_lags = (self.regressor_matrix.shape[0] - 1) / 2.0
                    r2_x = np.linspace(
                        -half_lags / self.fs, half_lags / self.fs, num_r2
                    )
                    shift = float(r2_x[max_idx])
                else:
                    shift = float(max_idx)

        self.update_all(x, y, z, shift)

    def on_canvas_click(self, event):
        """Select a voxel from slices or update slider/shift from the R² plot."""
        names = {axis: name for name, axis in self.axes.items()}
        name = names.get(event.inaxes)
        if name is None or event.xdata is None or event.ydata is None:
            return

        if name == 'r2':
            # Clicked on R² plot: update the shift/slider value based on x-coordinate
            clicked_shift = float(event.xdata)
            self.update_all(*self.curr, clicked_shift)
            return

        if name == 'ts':
            return

        x, y, z = self.curr
        u, v = int(round(event.xdata)), int(round(event.ydata))

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
