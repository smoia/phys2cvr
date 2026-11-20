#!/usr/bin/env python3
"""Voxel viewer for phys2cvr, built with bootstack."""

import logging

import bootstack as bs
import matplotlib as mpl
import numpy as np
from matplotlib.figure import Figure

from phys2cvr.io import EXT_NIMG, load_nifti_get_mask
from phys2cvr.utils import estimate_lag_range

LGR = logging.getLogger(__name__)


NIFTI_TYPES = [
    ('NIfTI files', ' '.join(f'*{ext}' for ext in EXT_NIMG)),
    ('All files', '*.*'),
]


class DiffViewer:
    """Interactive viewer for fMRI voxels, CVR/lag maps and shifted regressors."""

    def __init__(self):
        self.cvr_range_base = None
        self.lag_range_base = None
        self.cvr_range_comp = None
        self.lag_range_comp = None
        self.cvr_range_diff = None
        self.lag_range_diff = None
        self.lag_step = None

        self.timeseries = None
        self.data_types = [
            'mask_base',
            'raw_cvr_base',
            'raw_lag_base',
            'raw_r2_base',
            'cvr_base',
            'lag_base',
            'r2_base',
            'mask_comp',
            'raw_cvr_comp',
            'raw_lag_comp',
            'raw_r2_comp',
            'cvr_comp',
            'lag_comp',
            'r2_comp',
            'cvr_diff',
            'lag_diff',
            'r2_diff',
            'cvr_reld',
            'lag_reld',
            'r2_reld',
        ]
        self.data = {k: None for k in self.data_types}
        self.pixdim = (1.0, 1.0, 1.0)
        self.nx = self.ny = self.nz = 1

        self.curr = [0, 0, 0]
        self._syncing = False
        self.axes = {}
        self.map_settings = {
            'cvr_base': {'limits': (-1.0, 1.0), 'cmap': 'inferno'},
            'lag_base': {'limits': (-1.0, 1.0), 'cmap': 'viridis_r'},
            'r2_base': {'limits': (-1.0, 1.0), 'cmap': 'inferno'},
            'cvr_comp': {'limits': (-1.0, 1.0), 'cmap': 'inferno'},
            'lag_comp': {'limits': (-1.0, 1.0), 'cmap': 'viridis_r'},
            'r2_comp': {'limits': (-1.0, 1.0), 'cmap': 'inferno'},
            'cvr_diff': {'limits': (-1.0, 1.0), 'cmap': 'inferno'},
            'lag_diff': {'limits': (-1.0, 1.0), 'cmap': 'viridis_r'},
            'r2_diff': {'limits': (-1.0, 1.0), 'cmap': 'inferno'},
            'cvr_reldiff': {'limits': (-1.0, 1.0), 'cmap': 'inferno'},
            'lag_reldiff': {'limits': (-1.0, 1.0), 'cmap': 'viridis_r'},
            'r2_reldiff': {'limits': (-1.0, 1.0), 'cmap': 'inferno'},
        }
        self.diffswitch = {}

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

    def _update_dimensions_from_loaded_data(self):
        """Infer spatial dimensions and pixdim from any available map if timeseries is not loaded."""
        # Check available 3D/4D maps for shape and zooms
        for k in self.data_types:
            img_data = self.data.get(k)
            if img_data is None:
                img_data = self.data.get(f'raw_{k}')

            if img_data is not None:
                self.nx, self.ny, self.nz = img_data.shape[:3]
                if self.curr == [0, 0, 0] or self.curr[0] >= self.nx:
                    self.curr = [self.nx // 2, self.ny // 2, self.nz // 2]
                break

    def compute_difference(self, kind):
        """Compute difference and relative difference between maps."""
        if None not in (self.data[f'{kind}_base'], self.data[f'{kind}_comp']):
            self.data[f'{kind}_diff'] = (
                self.data[f'{kind}_base'] - self.data[f'{kind}_comp']
            )
            self.data[f'{kind}_reld'] = np.zeros_like(self.data[f'{kind}_diff'])
            self.data[f'{kind}_reld'] = np.divide(
                self.data[f'{kind}_diff'] * 100,
                self.data[f'{kind}_base'],
                where=self.data[f'{kind}_base'] != 0,
            )
        else:
            self.data[f'{kind}_diff'] = None
            self.data[f'{kind}_reld'] = None

    def set_cvr(self, path, ver):
        """Load and mask the CVR map."""
        self.data[f'raw_cvr_{ver}'], _, _ = load_nifti_get_mask(path)
        self._apply_mask([f'cvr_{ver}'])
        cvr_vals = self.data[f'cvr_{ver}'][np.isfinite(self.data[f'cvr_{ver}'])]
        if cvr_vals.size > 0:
            self.cvr_vmin, self.cvr_vmax = np.percentile(cvr_vals, [1, 99])
        else:
            self.cvr_vmin, self.cvr_vmax = None, None
        self.cvr_range = (self.cvr_vmin, self.cvr_vmax)
        if self.cvr_vmin is not None:
            self._set_map_limits(f'cvr_{ver}', self.cvr_vmin, self.cvr_vmax)
        self._update_dimensions_from_loaded_data()
        self.compute_difference('cvr')

    def set_lag(self, path, ver):
        """Load the lag map and estimate the available lag range."""
        self.data[f'raw_lag_{ver}'], _, _ = load_nifti_get_mask(path)
        self._apply_mask([f'lag_{ver}'])
        finite_lags = self.data[f'raw_lag_{ver}'][
            np.isfinite(self.data[f'raw_lag_{ver}'])
        ]
        self.lag_step, lag_max, lag_min = estimate_lag_range(finite_lags, path)
        self.lag_range = (lag_min, lag_max)
        lag_bound = max(abs(lag_min), abs(lag_max))
        self._set_map_limits(f'lag_{ver}', -lag_bound, lag_bound)
        self._update_dimensions_from_loaded_data()
        self.compute_difference('lag')

    def set_r2(self, path, ver):
        """Load the 4D R² map volume."""
        self.data[f'raw_r2_{ver}'], _, _ = load_nifti_get_mask(path, dim=4)
        self._apply_mask([f'r2_{ver}'])
        self._update_dimensions_from_loaded_data()
        self.compute_difference('r2')

    def set_mask(self, path, ver):
        """Load a mask and apply it to the maps."""
        _, self.data[f'mask_{ver}'], _ = load_nifti_get_mask(path, is_mask=True)
        self._apply_mask()
        self._update_dimensions_from_loaded_data()

    # ----- UI -----

    def build(self):
        """Build the viewer inside the current page."""

        def loader(title, file_types, setter, ver):
            return lambda: self._load(title, file_types, setter, ver)

        # --- Toolbar ------------------------------------------------------
        with bs.Toolbar(horizontal='stretch', surface='chrome') as bar:
            self.toolbar = bar
            with bar.add_menu('Load base') as menu:
                menu.add_action(
                    'Load result folder',
                    on_click=lambda: self._load_folder(),
                )
                menu.add_action(
                    'Load CVR map',
                    on_click=loader(
                        'Select 3D CVR Map NIfTI File',
                        NIFTI_TYPES,
                        self.set_cvr,
                        ver='base',
                    ),
                )
                menu.add_action(
                    'Load lag map',
                    on_click=loader(
                        'Select 3D Lag Map NIfTI File',
                        NIFTI_TYPES,
                        self.set_lag,
                        ver='base',
                    ),
                )
                menu.add_action(
                    'Load R2 map',
                    on_click=loader(
                        'Select 4D R2 Map NIfTI File',
                        NIFTI_TYPES,
                        self.set_r2,
                        ver='base',
                    ),
                )
                menu.add_action(
                    'Load mask',
                    on_click=loader(
                        'Select 3D Mask NIfTI File',
                        NIFTI_TYPES,
                        self.set_mask,
                        ver='base',
                    ),
                )
            with bar.add_menu('Load comparison') as menu:
                menu.add_action(
                    'Load result folder',
                    on_click=lambda: self._load_folder(),
                )
                menu.add_action(
                    'Load CVR map',
                    on_click=loader(
                        'Select 3D CVR Map NIfTI File',
                        NIFTI_TYPES,
                        self.set_cvr,
                        ver='comp',
                    ),
                )
                menu.add_action(
                    'Load lag map',
                    on_click=loader(
                        'Select 3D Lag Map NIfTI File',
                        NIFTI_TYPES,
                        self.set_lag,
                        ver='comp',
                    ),
                )
                menu.add_action(
                    'Load R2 map',
                    on_click=loader(
                        'Select 4D R2 Map NIfTI File',
                        NIFTI_TYPES,
                        self.set_r2,
                        ver='comp',
                    ),
                )
                menu.add_action(
                    'Load mask',
                    on_click=loader(
                        'Select 3D Mask NIfTI File',
                        NIFTI_TYPES,
                        self.set_mask,
                        ver='comp',
                    ),
                )
        self.map_fields = {}
        with bs.Grid(
            columns=['auto', 1, 1, 1],
            rows=3,
            grow=True,
            horizontal_items='stretch',
            vertical_items='stretch',
        ):
            # --- Settings bar ---------------------------------------------------

            with bs.ScrollView(grow=True, horizontal='stretch', column=0, rowspan=3):
                with bs.Column(padding=8, gap=3, horizontal_items='stretch'):
                    for kind, title in (
                        ('cvr_base', 'CVR map: base'),
                        ('lag_base', 'Lag map: base'),
                        ('r2_base', 'R2 map: base'),
                        ('cvr_comp', 'CVR map: comparison'),
                        ('lag_comp', 'Lag map: comparison'),
                        ('r2_comp', 'R2 map: comparison'),
                        ('cvr_diff', 'CVR map: difference'),
                        ('lag_diff', 'Lag map: difference'),
                        ('r2_diff', 'R2 map: difference'),
                    ):
                        bs.Label(title, font='caption', margin_y=(5, 8))
                        lower, upper = self.map_settings[kind]['limits']
                        vmin = bs.NumberField(
                            float(lower), label='Min', step=0.1, width=4
                        )
                        vmax = bs.NumberField(
                            float(upper), label='Max', step=0.1, width=4
                        )
                        cmap = bs.Select(
                            sorted(mpl.colormaps),
                            value=self.map_settings[kind]['cmap'],
                            label='Colorscale',
                            width=10,
                        )
                        self.map_fields[kind] = (vmin, vmax, cmap)

                        if 'diff' in kind:
                            self.diffswitch['rel'] = bs.Switch(
                                'Relative diff', value=False
                            )
                            self.diffswitch['abs'] = bs.Switch(
                                'Absolute diff', value=False
                            )

                        bs.Divider(margin_y=(10, 5))

                    bs.Label('Voxel position', font='caption', margin_y=(5, 8))
                    self.x_field = bs.NumberField(0, label='X', min_value=0, width=4)
                    self.y_field = bs.NumberField(0, label='Y', min_value=0, width=4)
                    self.z_field = bs.NumberField(0, label='Z', min_value=0, width=4)

                    bs.Spacer()

                    # Add R² plot here.

            # --- Views ----------------------------------------------------------
            self.charts['cvr_base'] = bs.Chart(
                figure=Figure(figsize=(1, 1)), row=0, column=1
            )
            self.charts['cvr_comp'] = bs.Chart(
                figure=Figure(figsize=(1, 1)), row=1, column=1
            )
            self.charts['cvr_diff'] = bs.Chart(
                figure=Figure(figsize=(1, 1)), row=2, column=1
            )
            self.charts['lag_base'] = bs.Chart(
                figure=Figure(figsize=(1, 1)), row=0, column=2
            )
            self.charts['lag_comp'] = bs.Chart(
                figure=Figure(figsize=(1, 1)), row=1, column=2
            )
            self.charts['lag_diff'] = bs.Chart(
                figure=Figure(figsize=(1, 1)), row=2, column=2
            )
            self.charts['r2_base'] = bs.Chart(
                figure=Figure(figsize=(1, 1)), row=0, column=3
            )
            self.charts['r2_comp'] = bs.Chart(
                figure=Figure(figsize=(1, 1)), row=1, column=3
            )
            self.charts['r2_diff'] = bs.Chart(
                figure=Figure(figsize=(1, 1)), row=2, column=3
            )

        for field in (self.x_field, self.y_field, self.z_field):
            field.on_change(lambda e: self._on_coord_change())
        for kind, fields in self.map_fields.items():
            for field in fields:
                field.on_change(lambda e, kind=kind: self._on_map_settings_change(kind))
        self.reldiff_switch.on_change(lambda e: self._update_current_plot())
        self.refresh()

    def _set_map_limits(self, kind, lower, upper):
        lower, upper = float(lower), float(upper)
        if lower == upper:
            delta = max(abs(lower) * 0.01, 0.01)
            lower, upper = lower - delta, upper + delta
        self.map_settings[kind]['limits'] = (lower, upper)
        if hasattr(self, 'map_fields') and kind in self.map_fields:
            syncing = self._syncing
            self._syncing = True
            try:
                self.map_fields[kind][0].value = lower
                self.map_fields[kind][1].value = upper
            finally:
                self._syncing = syncing

    def _on_map_settings_change(self, kind):
        if self._syncing:
            return
        vmin, vmax, cmap = self.map_fields[kind]
        lower, upper = vmin.value, vmax.value
        if lower is None or upper is None:
            return
        if not np.isfinite(lower) or not np.isfinite(upper) or lower >= upper:
            previous = self.map_settings[kind]['limits']
            self._set_map_limits(kind, *previous)
            LGR.warning('Invalid %s color range: %s to %s', kind, lower, upper)
            bs.alert(
                'Minimum must be less than maximum.',
                title='Invalid color range',
                severity='warning',
            )
            return
        self.map_settings[kind] = {
            'limits': (float(lower), float(upper)),
            'cmap': cmap.value,
        }
        self._update_current_plot()

    def _update_current_plot(self):
        if not self._syncing and self.axes:
            self.update_all(*self.curr, self.curr_shift_sec)

    def _load(self, title, file_types, setter, ver):
        path = bs.ask_open_file(title=title, file_types=file_types)
        if not path:
            return
        try:
            setter(path)
            self.refresh()
        except Exception as e:  # noqa: BLE001
            LGR.error(
                f'Fatal error encountered while loading {path}: {e}',
                exc_info=True,
            )
            bs.alert(
                f'Fatal error encountered while loading {path}: {e}',
                title='Could not load file',
                severity='danger',
            )

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
        self._update_dimensions_from_loaded_data()

        has_data = any(self.data.values())

        if not has_data:
            self.axes = {}
            for chart in self.charts.values():
                chart.figure = Figure(figsize=(1, 1))
            return

        self.init_plots()

        self._syncing = True
        try:
            self.x_field.max_value = self.nx - 1
            self.y_field.max_value = self.ny - 1
            self.z_field.max_value = self.nz - 1
        finally:
            self._syncing = False

        self._move_to(*self.curr)

    def init_plots(self):
        """Create separate figures for maps, timecourses and R²."""
        dx, dy, dz = self.pixdim
        aspects = {'sag': dz / dy, 'cor': dz / dx, 'axi': dy / dx}
        self.axes = {}
        self.imgs = {}
        self.cross = {}

        def add_view(
            fig, name, subplot, aspect, cmap, norm, color, xlabel, ylabel, size
        ):
            axis = self.axes[name] = fig.add_subplot(subplot)
            self.imgs[name] = axis.imshow(
                np.zeros(size[::-1]),
                cmap=cmap,
                norm=norm,
                origin='lower',
                aspect=aspect,
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

        figures = {}
        for kind, color in (
            ('cvr', 'cyan'),
            ('lag', 'red'),
            ('r2', 'cyan'),
        ):
            fig = figures[kind] = Figure(figsize=(1, 1), layout='compressed')
            fig.get_layout_engine().set(w_pad=0.02, h_pad=0.02, wspace=0, hspace=0)
            grid = fig.add_gridspec(1, 3)
            settings = self.map_settings[kind]
            cmap = mpl.colormaps[settings['cmap']].with_extremes(bad='black')
            norm = mpl.colors.Normalize(*settings['limits'])
            axes = []
            for index, (name, xlabel, ylabel, size) in enumerate(
                (
                    ('sag', 'Y', 'Z', (ny, nz)),
                    ('cor', 'X', 'Z', (nx, nz)),
                    ('axi', 'X', 'Y', (nx, ny)),
                )
            ):
                axis = add_view(
                    fig,
                    f'{kind}_{name}',
                    grid[0, index],
                    aspects[name],
                    cmap,
                    norm,
                    color,
                    xlabel,
                    ylabel,
                    size,
                )
                axis.set_anchor(('E', 'C', 'W')[index])
                axis.set_facecolor('black')
                axes.append(axis)
            fig.colorbar(self.imgs[f'{kind}_axi'], ax=axes, fraction=0.035, pad=0.02)

        for name, fig in figures.items():
            self.charts[name].figure = fig
            fig.canvas.mpl_connect('button_press_event', self.on_canvas_click)

    def update_all(self, x, y, z):
        """Update plots and controls for a voxel and regressor delay."""
        if self._syncing or not hasattr(self, 'imgs') or not self.imgs:
            return
        self.curr = [x, y, z]
        ax = self.axes

        views = {}

        for ver in ['base', 'comp', 'diff', 'reld']:
            if self.data[f'cvr_{ver}'] is not None:
                val_cvr = self.data[f'cvr_{ver}'][x, y, z]
                cvr_title = (
                    f'CVR: {val_cvr:.2f}' if not np.isnan(val_cvr) else 'CVR: Masked'
                )
                views.update(
                    {
                        f'cvr_{ver}_sag': (
                            self.data[f'cvr_{ver}'][x, :, :],
                            y,
                            z,
                            'Sagittal',
                        ),
                        f'cvr_{ver}_cor': (
                            self.data[f'cvr_{ver}'][:, y, :],
                            x,
                            z,
                            'Coronal',
                        ),
                        f'cvr_{ver}_axi': (
                            self.data[f'cvr_{ver}'][:, :, z],
                            x,
                            y,
                            cvr_title,
                        ),
                    }
                )
            else:
                views.update(
                    {
                        f'cvr_{ver}_sag': (
                            np.full((self.ny, self.nz), np.nan),
                            y,
                            z,
                            'Sagittal',
                        ),
                        f'cvr_{ver}_cor': (
                            np.full((self.nx, self.nz), np.nan),
                            x,
                            z,
                            'Coronal',
                        ),
                        f'cvr_{ver}_axi': (
                            np.full((self.nx, self.ny), np.nan),
                            x,
                            y,
                            'Axial',
                        ),
                    }
                )

            if self.data[f'lag_{ver}'] is not None:
                val_lag = self.data[f'lag_{ver}'][x, y, z]
                lag_title = (
                    f'Lag: {val_lag:.2f}s' if not np.isnan(val_lag) else 'Lag: Masked'
                )
                views.update(
                    {
                        f'lag_{ver}_sag': (
                            self.data[f'lag_{ver}'][x, :, :],
                            y,
                            z,
                            'Sagittal',
                        ),
                        f'lag_{ver}_cor': (
                            self.data[f'lag_{ver}'][:, y, :],
                            x,
                            z,
                            'Coronal',
                        ),
                        f'lag_{ver}_axi': (
                            self.data[f'lag_{ver}'][:, :, z],
                            x,
                            y,
                            lag_title,
                        ),
                    }
                )
            else:
                views.update(
                    {
                        f'lag_{ver}_sag': (
                            np.full((self.ny, self.nz), np.nan),
                            y,
                            z,
                            'Sagittal',
                        ),
                        f'lag_{ver}_cor': (
                            np.full((self.nx, self.nz), np.nan),
                            x,
                            z,
                            'Coronal',
                        ),
                        f'lag_{ver}_axi': (
                            np.full((self.nx, self.ny), np.nan),
                            x,
                            y,
                            'Axial',
                        ),
                    }
                )

            if self.data[f'r2_{ver}'] is not None:
                val_lag = self.data[f'r2_{ver}'][x, y, z]
                lag_title = (
                    f'Lag: {val_lag:.2f}s' if not np.isnan(val_lag) else 'Lag: Masked'
                )
                views.update(
                    {
                        f'r2_{ver}_sag': (
                            self.data[f'r2_{ver}'][x, :, :],
                            y,
                            z,
                            'Sagittal',
                        ),
                        f'r2_{ver}_cor': (
                            self.data[f'r2_{ver}'][:, y, :],
                            x,
                            z,
                            'Coronal',
                        ),
                        f'r2_{ver}_axi': (
                            self.data[f'r2_{ver}'][:, :, z],
                            x,
                            y,
                            lag_title,
                        ),
                    }
                )
            else:
                views.update(
                    {
                        f'r2_{ver}_sag': (
                            np.full((self.ny, self.nz), np.nan),
                            y,
                            z,
                            'Sagittal',
                        ),
                        f'r2_{ver}_cor': (
                            np.full((self.nx, self.nz), np.nan),
                            x,
                            z,
                            'Coronal',
                        ),
                        f'r2_{ver}_axi': (
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
            img = self.imgs[name]
            img.set_data(slice_data.T)

            kind = name.split('_', 1)[0]
            settings = self.map_settings[kind]
            if img.get_cmap().name != settings['cmap']:
                img.set_cmap(mpl.colormaps[settings['cmap']].with_extremes(bad='black'))
            if img.get_clim() != tuple(settings['limits']):
                img.set_clim(*settings['limits'])
            self.cross[name][0].set_xdata([cx, cx])
            self.cross[name][1].set_ydata([cy, cy])
            ax[name].set_title(title)

        self._syncing = True
        try:
            self.x_field.value = x
            self.y_field.value = y
            self.z_field.value = z
        finally:
            self._syncing = False

        for chart in self.charts.values():
            chart.draw()

    # ----- callbacks -----

    def _move_to(self, x, y, z):
        """Move to a voxel, following its lag value or max R² index when applicable."""
        x = int(np.clip(x, 0, self.nx - 1))
        y = int(np.clip(y, 0, self.ny - 1))
        z = int(np.clip(z, 0, self.nz - 1))

        self.update_all(x, y, z)

    def on_canvas_click(self, event):
        """Select a voxel from slices or update slider/shift from the R² plot."""
        names = {axis: name for name, axis in self.axes.items()}
        name = names.get(event.inaxes)
        if name is None or event.xdata is None or event.ydata is None:
            return

        x, y, z = self.curr
        u, v = int(round(event.xdata)), int(round(event.ydata))

        if '_sag' in name:
            y, z = u, v
        elif '_cor' in name:
            x, z = u, v
        elif '_axi' in name:
            x, y = u, v

        self._move_to(x, y, z)

    def _on_coord_change(self):
        if self._syncing or not self.axes:
            return
        values = [f.value for f in (self.x_field, self.y_field, self.z_field)]
        if any(v is None for v in values):
            return
        self._move_to(*values)
