#!/usr/bin/env python3
"""Voxel viewer for phys2cvr, built with bootstack."""

import logging

import bootstack as bs
import matplotlib as mpl
import numpy as np
from matplotlib.figure import Figure
from scipy.stats import zscore

from phys2cvr.gui.io import load_timeseries
from phys2cvr.gui.utils import (
    ARRAY_TYPES,
    NIFTI_TYPES,
    MapViewerControls,
    build_coordinate_fields,
    build_map_fields,
    clip_position,
    make_orthogonal_figure,
    slice_position,
    update_map_view,
)
from phys2cvr.io import load_array, load_nifti_get_mask
from phys2cvr.utils import estimate_lag_range

LGR = logging.getLogger(__name__)



class VoxelViewer(MapViewerControls):
    """Interactive viewer for fMRI voxels, CVR/lag maps and shifted regressors."""

    def __init__(self, fs=None, lag_range=None, cvr_range=None):
        if fs is not None and (not np.isfinite(fs) or fs <= 0):
            raise ValueError(
                'Regressor sampling frequency must be positive and finite.'
            )
        self.fs = fs
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
        self.tr = None
        self.time = None
        self.nx = self.ny = self.nz = self.nt = 1

        self.curr = [0, 0, 0]
        self.curr_shift_sec = 0.0
        self.min_delay_sec = -1.0
        self.max_delay_sec = 1.0
        self.lag_units = 's'
        self._syncing = False
        self.axes = {}
        self.map_settings = {
            'cvr': {
                'limits': cvr_range
                if cvr_range and None not in cvr_range
                else (0.0, 1.0),
                'cmap': 'inferno',
            },
            'lag': {'limits': (-1.0, 1.0), 'cmap': 'viridis_r'},
            'mean': {'limits': (0.0, 1.0), 'cmap': 'gray'},
        }

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
        elif self.regressor_matrix is not None:
            half_lags = (self.regressor_matrix.shape[0] - 1) / 2.0
            frequency = self.fs if self.fs is not None else 1.0
            self.min_delay_sec = -half_lags / frequency
            self.max_delay_sec = half_lags / frequency
            self.lag_units = 's' if self.fs is not None else 'idx'
        elif self.data['r2'] is not None:
            self.min_delay_sec = 0.0
            self.max_delay_sec = max(1, self.data['r2'].shape[-1] - 1)
            self.lag_units = 'idx'
            self.slider.step = 1
        else:
            self.min_delay_sec, self.max_delay_sec = -1.0, 1.0
            self.lag_units = 'idx'

        if self.regressor_matrix is not None:
            step = 1.0 / self._regressor_fs()
        elif self.lag_range is not None and self.data['r2'] is not None:
            count = self.data['r2'].shape[-1]
            step = (self.max_delay_sec - self.min_delay_sec) / max(count - 1, 1)
        else:
            step = 1.0

        self._syncing = True
        try:
            self.lagmin_field.value = self.min_delay_sec
            self.lagmax_field.value = self.max_delay_sec
            self.slider.min_value = self.min_delay_sec
            self.slider.max_value = self.max_delay_sec
            self.slider.step = step
            self.slider.disabled = (
                self.regressor_matrix is None and self.data['r2'] is None
            )
            self.CO2_label.text = f'CO2 Delay ({self.lag_units})'
        finally:
            self._syncing = False

    def _update_dimensions_from_loaded_data(self):
        """Infer spatial dimensions and pixdim from any available map if timeseries is not loaded."""
        if self.timeseries is not None:
            return

        # Check available 3D/4D maps for shape and zooms
        for k in ('cvr', 'lag', 'r2', 'mask', 'meanvol'):
            img_data = self.data.get(k)
            if img_data is None:
                img_data = self.data.get(f'raw_{k}')

            if img_data is not None:
                self.nx, self.ny, self.nz = img_data.shape[:3]
                if self.curr == [0, 0, 0] or self.curr[0] >= self.nx:
                    self.curr = [self.nx // 2, self.ny // 2, self.nz // 2]
                break

    def _update_zooms_from_img(self, img):
        if img is not None and hasattr(img, 'header'):
            zooms = img.header.get_zooms()
            self.pixdim = tuple(zooms[:3]) if len(zooms) >= 3 else (1.0, 1.0, 1.0)

    def set_timeseries(self, path):
        """Load the functional timeseries and initialize the selected voxel."""
        self.timeseries, self.pixdim, self.tr = load_timeseries(path)
        self.nx, self.ny, self.nz, self.nt = self.timeseries.shape
        self.time = np.arange(self.nt) * self.tr
        self.data['meanvol'] = np.mean(self.timeseries, axis=-1)
        finite = self.data['meanvol'][np.isfinite(self.data['meanvol'])]
        if finite.size:
            self._set_map_limits('mean', float(finite.min()), float(finite.max()))
        self.curr = [self.nx // 2, self.ny // 2, self.nz // 2]

        # If a mask was not provided yet, generate it from non-zero mean volume voxels
        if self.data['mask'] is None:
            self.data['mask'] = self.data['meanvol'] != 0

        self._apply_mask()

    def set_regressors(self, path):
        """Load the shifted regressor matrix."""
        self.regressor_matrix = load_array(path, shape='prerectangle')
        self._update_delay_range()

    def set_cvr(self, path):
        """Load and mask the CVR map."""
        self.data['raw_cvr'], _, _ = load_nifti_get_mask(path)
        self._apply_mask(['cvr'])
        cvr_vals = self.data['cvr'][np.isfinite(self.data['cvr'])]
        if cvr_vals.size > 0:
            self.cvr_vmin, self.cvr_vmax = np.percentile(cvr_vals, [1, 99])
        else:
            self.cvr_vmin, self.cvr_vmax = None, None
        self.cvr_range = (self.cvr_vmin, self.cvr_vmax)
        if self.cvr_vmin is not None:
            self._set_map_limits('cvr', self.cvr_vmin, self.cvr_vmax)
        self._update_dimensions_from_loaded_data()

    def set_lag(self, path):
        """Load the lag map and estimate the available lag range."""
        self.data['raw_lag'], _, _ = load_nifti_get_mask(path)
        self._apply_mask(['lag'])
        finite_lags = self.data['raw_lag'][np.isfinite(self.data['raw_lag'])]
        self.lag_step, lag_max, lag_min = estimate_lag_range(finite_lags, path)
        self.lag_range = (lag_min, lag_max)
        lag_bound = max(abs(lag_min), abs(lag_max))
        self._set_map_limits('lag', -lag_bound, lag_bound)
        self._update_delay_range()
        self._update_dimensions_from_loaded_data()

    def set_r2(self, path):
        """Load the 4D R² map volume."""
        self.data['raw_r2'], _, _ = load_nifti_get_mask(path, dim=4)
        if self.data['mask'] is not None:
            self.data['r2'] = self.data['raw_r2'].copy()
            # Mask spatial dimensions while keeping 4th dimension intact
            self.data['r2'][~self.data['mask']] = np.nan
        else:
            self.data['r2'] = self.data['raw_r2']
        self._update_delay_range()
        self._update_dimensions_from_loaded_data()

    def set_mask(self, path):
        """Load a mask and apply it to the maps."""
        _, self.data['mask'], _ = load_nifti_get_mask(path, is_mask=True)
        self._apply_mask()
        self._update_dimensions_from_loaded_data()

    def _regressor_fs(self):
        """Return the supplied frequency, or infer it from lag range and rows."""
        if self.fs is not None:
            return self.fs
        if self.lag_range is not None and self.regressor_matrix.shape[0] > 1:
            return (self.regressor_matrix.shape[0] - 1) / (
                self.lag_range[1] - self.lag_range[0]
            )
        return 1.0

    def _regressor_row(self, shift_sec):
        """Return the regressor shifted by `shift_sec` and the actual shift used."""
        num_lags = self.regressor_matrix.shape[0]
        frequency = self._regressor_fs()
        idx = int(round((shift_sec - self.min_delay_sec) * frequency))
        idx = int(np.clip(idx, 0, num_lags - 1))
        return self.regressor_matrix[idx], idx / frequency + self.min_delay_sec

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
        self.map_fields = {}
        with bs.Row(grow=True, horizontal='stretch', vertical_items='stretch'):
            with bs.Column(
                width=190,
                surface='card_raised',
                vertical='stretch',
                horizontal_items='stretch',
            ):
                with bs.ScrollView(grow=True, horizontal='stretch'):
                    with bs.Column(padding=8, gap=3, horizontal_items='stretch'):
                        for kind, title in (
                            ('cvr', 'CVR map'),
                            ('lag', 'Lag map'),
                            ('mean', 'Mean volume'),
                        ):
                            self.map_fields[kind] = build_map_fields(
                                self.map_settings[kind], title
                            )
                            bs.Divider(margin_y=(10, 5))

                        self.x_field, self.y_field, self.z_field = (
                            build_coordinate_fields()
                        )
                        bs.Divider(margin_y=(10, 5))
                        bs.Label('Timeseries', font='caption', margin_y=(5, 8))
                        self.zscore_switch = bs.Switch('Z-score', value=True)
                        self.lagmin_field = bs.NumberField(
                            0.0, label='Lag min', step=0.1, width=4
                        )
                        self.lagmax_field = bs.NumberField(
                            0.0, label='Lag max', step=0.1, width=4
                        )
            with bs.Column(grow=True, horizontal='stretch', vertical='stretch'):
                self._build_charts()
                with bs.Row(
                    gap=16, padding=10, horizontal='stretch', vertical='center'
                ):
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
        for kind, fields in self.map_fields.items():
            for field in fields:
                field.on_change(lambda e, kind=kind: self._on_map_settings_change(kind))
        self.zscore_switch.on_change(lambda e: self._update_current_plot())
        self.slider.on_change(lambda e: self._on_slider_move())
        self.refresh()

    def _build_charts(self):
        """Create the five chart containers in the viewer layout."""
        self.charts = {}
        with bs.Column(
            grow=True, horizontal='stretch', weights=[7, 3], horizontal_items='stretch'
        ):
            with bs.Row(weights=[5, 1], vertical_items='stretch'):
                with bs.Column(weights=[1, 1], horizontal_items='stretch'):
                    self.charts['cvr'] = bs.Chart(figure=Figure(figsize=(1, 1)))
                    self.charts['lag'] = bs.Chart(figure=Figure(figsize=(1, 1)))
                self.charts['mean'] = bs.Chart(figure=Figure(figsize=(1, 1)))
            with bs.Row(weights=[4, 1], vertical_items='stretch'):
                self.charts['ts'] = bs.Chart(figure=Figure(figsize=(1, 1)))
                self.charts['r2'] = bs.Chart(
                    figure=Figure(figsize=(1, 1)), margin_x=(20, 50)
                )

    def _update_current_plot(self):
        if not self._syncing and self.axes:
            self.update_all(*self.curr, self.curr_shift_sec)

    def refresh(self):
        """Rebuild the figure and sync controls after the data changed."""
        self._update_delay_range()
        self._update_dimensions_from_loaded_data()

        has_data = (
            self.timeseries is not None
            or self.data['cvr'] is not None
            or self.data['lag'] is not None
            or self.data['r2'] is not None
            or self.regressor_matrix is not None
        )

        if not has_data:
            self.axes = {}
            for chart in self.charts.values():
                chart.figure = Figure(figsize=(1, 1))
            return

        self.init_plots()

        self._set_coordinate_limits()
        self._move_to(*self.curr)

    def init_plots(self):
        """Create separate figures for maps, timecourses and R²."""
        self.axes = {}
        self.imgs = {}
        self.cross = {}

        figures = {}
        for kind, color, show_colorbar in (
            ('cvr', 'cyan', True),
            ('lag', 'red', True),
            ('mean', 'red', False),
        ):
            figures[kind] = make_orthogonal_figure(
                self, kind, color, vertical=kind == 'mean', colorbar=show_colorbar
            )

        ts_fig = figures['ts'] = Figure(figsize=(1, 1), layout='constrained')
        ts_ax = self.axes['ts'] = ts_fig.add_subplot()
        (self.line_voxel,) = ts_ax.plot([], [], label='Voxel TS', color='#1f77b4')
        (self.line_co2,) = ts_ax.plot(
            [], [], label='Regressor', color='#ff7f0e', linestyle='--'
        )
        ts_ax.set_title('Timecourse Comparison — Voxel (0, 0, 0)')
        ts_ax.set_xlabel('Time (s)' if self.timeseries is not None else 'TR')
        ts_ax.set_ylabel('Z-Score')
        ts_ax.set_xlim(0, 1)
        ts_ax.legend(loc='upper right')

        text_color = ts_ax.yaxis.label.get_color()
        legend = ts_ax.legend(loc='upper right', frameon=True)
        frame = legend.get_frame()
        frame.set_facecolor(ts_fig.get_facecolor())
        for text in legend.get_texts():
            text.set_color(text_color)

        r2_fig = figures['r2'] = Figure(figsize=(1, 1), layout='constrained')
        r2_ax = self.axes['r2'] = r2_fig.add_subplot()
        (self.line_r2,) = r2_ax.plot([], [], label='R²', color='#2ca02c')
        self.r2_bar = r2_ax.axvline(0, color='red', linestyle=':', lw=1.5)
        (self.r2_point,) = r2_ax.plot([], [], 'o', color='red')
        r2_ax.set_title('R²: N/A')
        r2_ax.set_xlabel(
            'Lag'
            if self.data['lag'] is not None or self.lag_range is not None
            else 'Index'
        )
        r2_ax.set_ylabel('R²')

        # Format y-axis ticks to always show exactly 2 decimals
        r2_ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.2f'))

        for name, fig in figures.items():
            self.charts[name].figure = fig
            fig.canvas.mpl_connect('button_press_event', self.on_canvas_click)

    def update_all(self, x, y, z, shift_sec):
        """Update plots and controls for a voxel and regressor delay."""
        if self._syncing or not hasattr(self, 'imgs') or not self.imgs:
            return
        shift_sec = float(np.clip(shift_sec, self.min_delay_sec, self.max_delay_sec))
        if self.regressor_matrix is not None:
            _, shift_sec = self._regressor_row(shift_sec)
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
                'Coronal',
            ),
            'mean_axi': (
                self.data['meanvol'][:, :, z]
                if self.data['meanvol'] is not None
                else np.zeros((self.nx, self.ny)),
                x,
                y,
                'Axial',
            ),
            'mean_sag': (
                self.data['meanvol'][x, :, :]
                if self.data['meanvol'] is not None
                else np.zeros((self.ny, self.nz)),
                y,
                z,
                'Sagittal',
            ),
        }

        if self.data['cvr'] is not None:
            val_cvr = self.data['cvr'][x, y, z]
            cvr_title = (
                f'CVR: {val_cvr:.2f}' if not np.isnan(val_cvr) else 'CVR: Masked'
            )
            views.update(
                {
                    'cvr_sag': (self.data['cvr'][x, :, :], y, z, 'Sagittal'),
                    'cvr_cor': (self.data['cvr'][:, y, :], x, z, 'Coronal'),
                    'cvr_axi': (self.data['cvr'][:, :, z], x, y, cvr_title),
                }
            )
        else:
            views.update(
                {
                    'cvr_sag': (
                        np.full((self.ny, self.nz), np.nan),
                        y,
                        z,
                        'Sagittal',
                    ),
                    'cvr_cor': (
                        np.full((self.nx, self.nz), np.nan),
                        x,
                        z,
                        'Coronal',
                    ),
                    'cvr_axi': (
                        np.full((self.nx, self.ny), np.nan),
                        x,
                        y,
                        'Axial',
                    ),
                }
            )

        if self.data['lag'] is not None:
            val_lag = self.data['lag'][x, y, z]
            lag_title = (
                f'Lag: {val_lag:.2f}s' if not np.isnan(val_lag) else 'Lag: Masked'
            )
            views.update(
                {
                    'lag_sag': (self.data['lag'][x, :, :], y, z, 'Sagittal'),
                    'lag_cor': (self.data['lag'][:, y, :], x, z, 'Coronal'),
                    'lag_axi': (self.data['lag'][:, :, z], x, y, lag_title),
                }
            )
        else:
            views.update(
                {
                    'lag_sag': (
                        np.full((self.ny, self.nz), np.nan),
                        y,
                        z,
                        'Sagittal',
                    ),
                    'lag_cor': (
                        np.full((self.nx, self.nz), np.nan),
                        x,
                        z,
                        'Coronal',
                    ),
                    'lag_axi': (
                        np.full((self.nx, self.ny), np.nan),
                        x,
                        y,
                        'Axial',
                    ),
                }
            )

        for name, (slice_data, cx, cy, title) in views.items():
            if name not in self.imgs:
                continue
            kind = name.split('_', 1)[0]
            update_map_view(
                self, name, slice_data, cx, cy, title, self.map_settings[kind]
            )

        # --- Timecourse Comparison Plot ---
        ts_ax = ax['ts']
        ts_ax.set_title(f'Timecourse Comparison — Voxel ({x}, {y}, {z})')

        has_ts = self.timeseries is not None
        has_regr = self.regressor_matrix is not None

        # 1. Update Voxel Timeseries line
        if has_ts:
            voxel_ts = self._normalize(self.timeseries[x, y, z, :])
            self.line_voxel.set_data(self.time, voxel_ts)
            ts_ax.set_xlim(0, self.time[-1] if self.nt > 1 else 1)
            self.line_voxel.set_visible(True)
            ts_ax.set_xlabel('Time (s)')
        else:
            self.line_voxel.set_visible(False)

        # 2. Update Regressor line
        if has_regr:
            row, _ = self._regressor_row(shift_sec)
            regr_trace = self._normalize(row)

            if self.time is not None:
                if regr_trace.size != self.time.size:
                    regr_trace = np.interp(
                        self.time,
                        np.linspace(self.time[0], self.time[-1], regr_trace.size),
                        regr_trace,
                    )
                regr_time = self.time
                ts_ax.set_xlabel('Time (s)')
            else:
                regr_time = np.arange(regr_trace.size) + 1
                ts_ax.set_xlabel('TR')

            self.line_co2.set_data(regr_time, regr_trace)
            self.line_co2.set_label(f'Regressor (shift: {shift_sec:.2f})')

            if not has_ts:
                ts_ax.set_xlim(0, regr_time[-1] if len(regr_time) > 1 else 1)
            self.line_co2.set_visible(True)
        else:
            self.line_co2.set_visible(False)

        ts_ax.relim()
        ts_ax.set_ylabel('Z-Score' if self.zscore_switch.value else 'Demeaned')
        ts_ax.autoscale_view(scalex=False, scaley=True)

        # --- R2 Profile Plot ---
        r2_ax = ax['r2']
        current_r2_val = None
        if self.data['r2'] is not None:
            r2_voxel = self.data['r2'][x, y, z, :]
            num_r2 = len(r2_voxel)

            # Find closest R² value corresponding to current shift/slider
            if self.lag_range is not None:
                r2_x = np.linspace(self.lag_range[0], self.lag_range[1], num_r2)
                r2_xlabel = 'Lag (s)'
            elif self.regressor_matrix is not None:
                r2_x = np.linspace(self.min_delay_sec, self.max_delay_sec, num_r2)
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
                self.r2_point.set_data([r2_x[idx_closest]], [current_r2_val])
        else:
            self.line_r2.set_data([], [])
            self.r2_point.set_data([], [])

        self.r2_bar.set_xdata([shift_sec, shift_sec])
        r2_ax.set_title(
            'R²: N/A' if current_r2_val is None else f'R²: {current_r2_val:.4f}'
        )

        # --- Updates a few things ---
        actual_sec = shift_sec
        self._syncing = True
        try:
            self._sync_coordinates((x, y, z))
            self.slider.value = float(
                np.clip(actual_sec, self.min_delay_sec, self.max_delay_sec)
            )
            self.delay_label.text = f'{actual_sec:+.2f} {self.lag_units}'
        finally:
            self._syncing = False

        for chart in self.charts.values():
            chart.draw()

    # ----- callbacks -----

    def _normalize(self, values):
        if self.zscore_switch.value:
            return zscore(values)
        return values - np.mean(values)

    def _move_to(self, x, y, z):
        """Move to a voxel, following its lag value or max R² index when applicable."""
        x, y, z = clip_position((x, y, z), (self.nx, self.ny, self.nz))
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
                    r2_x = np.linspace(self.min_delay_sec, self.max_delay_sec, num_r2)
                    shift = float(r2_x[max_idx])
                else:
                    shift = float(max_idx)

        self.update_all(x, y, z, shift)

    def on_canvas_click(self, event):
        """Select a voxel from slices or update slider/shift from the R² plot."""
        names = {axis: name for name, axis in self.axes.items()}
        name = names.get(event.inaxes)
        if (
            name is None
            or event.xdata is None
            or event.ydata is None
            or not np.isfinite(event.xdata)
            or not np.isfinite(event.ydata)
        ):
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

        orientation = ('sag', 'cor', 'axi').index(name.rsplit('_', 1)[-1])
        x, y, z = slice_position((x, y, z), orientation, u, v)

        self._move_to(x, y, z)

    def _validate_lag_range(self, lower, upper):
        if not np.isfinite(lower) or not np.isfinite(upper) or lower >= upper:
            raise ValueError('Lag minimum must be finite and less than lag maximum.')
        if self.fs is not None and self.regressor_matrix is not None:
            capacity = (self.regressor_matrix.shape[0] - 1) / self.fs
            if upper - lower > capacity and not np.isclose(upper - lower, capacity):
                raise ValueError(
                    'The lag range is too big for the regressor matrix at the '
                    f'specified sampling frequency (maximum span: {capacity:g} s).'
                )

    def _on_lagrange_change(self):
        if self._syncing:
            return
        lower, upper = self.lagmin_field.value, self.lagmax_field.value
        if lower is None or upper is None:
            return
        try:
            self._validate_lag_range(lower, upper)
        except ValueError as err:
            LGR.warning('Invalid lag range: %s', err)
            self._update_delay_range()
            bs.alert(str(err), title='Invalid lag range', severity='warning')
            return
        self.lag_range = (float(lower), float(upper))
        self._update_delay_range()
        if self.axes:
            self.update_all(*self.curr, self.curr_shift_sec)

    def _on_slider_move(self):
        if self._syncing or not self.axes:
            return
        self.update_all(*self.curr, float(self.slider.value))
