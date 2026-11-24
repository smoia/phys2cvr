#!/usr/bin/env python3
"""Comparison viewer for phys2cvr, built with bootstack."""

import logging

import bootstack as bs
import matplotlib as mpl
import nibabel as nib
import numpy as np
from matplotlib.figure import Figure

from phys2cvr.io import EXT_NIMG

LGR = logging.getLogger(__name__)


NIFTI_TYPES = [
    ('NIfTI files', ' '.join(f'*{ext}' for ext in EXT_NIMG)),
    ('All files', '*.*'),
]


class DiffViewer:
    """Compare CVR, lag, and R² maps at a shared voxel position."""

    def __init__(self):
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
        self._geometry = {}
        self._syncing = False
        self.axes = {}
        self.charts = {}
        self.summary_charts = {}
        self.r2_volumes = {'base': None, 'comp': None}
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
        }
        self.diffswitch = {}

    # ----- data handling -----

    def _read_image(self, path, key, *, r2=False, is_mask=False):
        """Read a map and validate its spatial grid before storing its geometry."""
        img = nib.load(path)
        if not hasattr(img, 'affine') or img.affine is None:
            raise ValueError(f'{path} does not have a 3D spatial affine.')
        data = np.asarray(img.get_fdata())
        if data.ndim not in ((3, 4) if r2 else (3,)) or any(
            size == 0 for size in data.shape
        ):
            raise ValueError(
                f'{path} must be a 3D map' + (' or 4D R² map.' if r2 else '.')
            )
        for other, (shape, affine, _) in self._geometry.items():
            if other != key and (
                data.shape[:3] != shape
                or not np.allclose(img.affine, affine, rtol=0, atol=1e-5)
            ):
                raise ValueError(
                    f'{path} does not match the spatial grid of {other}; '
                    'resample the images to the same grid before comparing.'
                )
        zooms = tuple(float(v) for v in img.header.get_zooms()[:3])
        self._geometry[key] = (data.shape[:3], img.affine.copy(), zooms)
        return np.isfinite(data) & (data != 0) if is_mask else data

    def _apply_mask(self, ver):
        """Reapply one dataset's mask to its original maps."""
        mask = self.data[f'mask_{ver}']
        for kind in ('cvr', 'lag', 'r2'):
            raw = self.data[f'raw_{kind}_{ver}']
            if raw is None:
                self.data[f'{kind}_{ver}'] = None
            elif mask is None:
                self.data[f'{kind}_{ver}'] = raw.copy()
            else:
                self.data[f'{kind}_{ver}'] = np.where(mask, raw, np.nan)
            self.compute_difference(kind)

    def _update_dimensions_from_loaded_data(self):
        """Use the common image grid for all plots."""
        if self._geometry:
            shape, _, self.pixdim = next(iter(self._geometry.values()))
            if (self.nx, self.ny, self.nz) != shape:
                self.nx, self.ny, self.nz = shape
                self.curr = [size // 2 for size in shape]
            else:
                self.curr = [
                    int(np.clip(value, 0, size - 1))
                    for value, size in zip(self.curr, shape, strict=True)
                ]

    def compute_difference(self, kind):
        """Compute difference and relative difference between maps."""
        base, comp = self.data[f'{kind}_base'], self.data[f'{kind}_comp']
        if base is not None and comp is not None:
            diff = self.data[f'{kind}_diff'] = base - comp
            reld = np.full(diff.shape, np.nan, dtype=float)
            np.divide(diff * 100, base, out=reld, where=np.isfinite(base) & (base != 0))
            self.data[f'{kind}_reld'] = reld
        else:
            self.data[f'{kind}_diff'] = None
            self.data[f'{kind}_reld'] = None

    def set_cvr(self, path, ver):
        """Load and mask the CVR map."""
        self.data[f'raw_cvr_{ver}'] = self._read_image(path, f'raw_cvr_{ver}')
        self._apply_mask(ver)
        cvr_vals = self.data[f'cvr_{ver}'][np.isfinite(self.data[f'cvr_{ver}'])]
        if cvr_vals.size:
            self._set_map_limits(f'cvr_{ver}', *np.percentile(cvr_vals, [1, 99]))
        self._update_dimensions_from_loaded_data()

    def set_lag(self, path, ver):
        """Load the lag map."""
        self.data[f'raw_lag_{ver}'] = self._read_image(path, f'raw_lag_{ver}')
        self._apply_mask(ver)
        finite_lags = self.data[f'lag_{ver}'][np.isfinite(self.data[f'lag_{ver}'])]
        if finite_lags.size:
            lag_bound = np.max(np.abs(finite_lags))
            self._set_map_limits(f'lag_{ver}', -lag_bound, lag_bound)
        self._update_dimensions_from_loaded_data()

    def set_r2(self, path, ver):
        """Load a 3D R² map or the per-voxel maximum of a 4D R² volume."""
        data = self._read_image(path, f'raw_r2_{ver}', r2=True)
        self.r2_volumes[ver] = data if data.ndim == 4 else None
        if data.ndim == 4:
            finite = np.isfinite(data)
            maximum = np.max(np.where(finite, data, -np.inf), axis=-1)
            maximum[~np.any(finite, axis=-1)] = np.nan
            data = maximum
        self.data[f'raw_r2_{ver}'] = data
        self._apply_mask(ver)
        self._update_dimensions_from_loaded_data()

    def set_mask(self, path, ver):
        """Load a mask and apply it to the maps."""
        self.data[f'mask_{ver}'] = self._read_image(path, f'mask_{ver}', is_mask=True)
        self._apply_mask(ver)
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
                    on_click=lambda: self._load_folder('base'),
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
                        'Select 3D or 4D R2 Map NIfTI File',
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
                    on_click=lambda: self._load_folder('comp'),
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
                        'Select 3D or 4D R2 Map NIfTI File',
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
            rows=4,
            grow=True,
            horizontal_items='stretch',
            vertical_items='stretch',
        ):
            # --- Settings bar ---------------------------------------------------

            with bs.Column(
                row=0,
                column=0,
                rowspan=4,
                width=190,
                surface='card_raised',
                horizontal='stretch',
                vertical='stretch',
            ):
                with bs.ScrollView(grow=True, horizontal='stretch'):
                    self._build_settings()

            # --- Views ----------------------------------------------------------
            for col, kind in enumerate(('cvr', 'lag', 'r2'), start=1):
                for row, ver in enumerate(('base', 'comp', 'diff')):
                    self.charts[f'{kind}_{ver}'] = bs.Chart(
                        figure=Figure(figsize=(1, 1)), row=row, column=col
                    )
                self.summary_charts[kind] = bs.Chart(
                    figure=Figure(figsize=(1, 1)), row=3, column=col
                )

        for field in (self.x_field, self.y_field, self.z_field):
            field.on_change(lambda e: self._on_coord_change())
        for kind, fields in self.map_fields.items():
            for field in fields:
                field.on_change(lambda e, kind=kind: self._on_map_settings_change(kind))
        for switch in self.diffswitch.values():
            switch.on_change(lambda e: self._update_current_plot())
        self.refresh()

    def _build_settings(self):
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
                vmin = bs.NumberField(float(lower), label='Min', step=0.1, width=4)
                vmax = bs.NumberField(float(upper), label='Max', step=0.1, width=4)
                cmap = bs.Select(
                    sorted(mpl.colormaps),
                    value=self.map_settings[kind]['cmap'],
                    label='Colorscale',
                    width=10,
                )
                self.map_fields[kind] = (vmin, vmax, cmap)

                if 'diff' in kind:
                    self.diffswitch[f'{kind[:-4]}rel'] = bs.Switch(
                        'Relative diff', value=False
                    )
                    self.diffswitch[f'{kind[:-4]}abs'] = bs.Switch(
                        'Absolute diff', value=False
                    )

                bs.Divider(margin_y=(10, 5))

            bs.Label('Voxel position', font='caption', margin_y=(5, 8))
            self.x_field = bs.NumberField(0, label='X', min_value=0, width=4)
            self.y_field = bs.NumberField(0, label='Y', min_value=0, width=4)
            self.z_field = bs.NumberField(0, label='Z', min_value=0, width=4)

            bs.Spacer()

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
            self.update_all(*self.curr)

    def _load(self, title, file_types, setter, ver):
        path = bs.ask_open_file(title=title, file_types=file_types)
        if not path:
            return
        try:
            setter(path, ver)
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

    def _load_folder(self, ver):
        folder_path = bs.ask_directory(title='Select phys2cvr Result Folder')
        if not folder_path:
            return
        try:
            from phys2cvr.gui.app import load_result_folder

            load_result_folder(folder_path, self, comparison_version=ver)
            self.refresh()
        except Exception as err:  # noqa: BLE001
            LGR.error(
                'Could not load result folder %s: %s', folder_path, err, exc_info=True
            )
            self.refresh()
            bs.alert(str(err), title='Could not load result folder')

    def refresh(self):
        """Rebuild the figure and sync controls after the data changed."""
        self._update_dimensions_from_loaded_data()

        has_data = any(value is not None for value in self.data.values())

        if not has_data:
            self.axes = {}
            for chart in self.charts.values():
                chart.figure = Figure(figsize=(1, 1))
            self._draw_summaries()
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
        """Create one figure per map, each with three orthogonal views."""
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
        for map_name in self.charts:
            kind = map_name.split('_')[0]
            color = 'red' if kind == 'lag' else 'cyan'
            fig = figures[map_name] = Figure(figsize=(1, 1), layout='compressed')
            fig.get_layout_engine().set(w_pad=0.02, h_pad=0.02, wspace=0, hspace=0)
            grid = fig.add_gridspec(1, 3)
            settings = self.map_settings[map_name]
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
                    f'{map_name}_{name}',
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
            fig.colorbar(
                self.imgs[f'{map_name}_axi'], ax=axes, fraction=0.035, pad=0.02
            )

        for name, fig in figures.items():
            self.charts[name].figure = fig
            fig.canvas.mpl_connect('button_press_event', self.on_canvas_click)

    def update_all(self, x, y, z):
        """Update every map and coordinate control for the selected voxel."""
        if self._syncing or not getattr(self, 'imgs', None):
            return
        self.curr = [x, y, z]
        for map_name in self.charts:
            kind, ver = map_name.split('_')
            data_key = map_name
            if ver == 'diff' and self.diffswitch[f'{kind}_rel'].value:
                data_key = f'{kind}_reld'
            data = self.data[data_key]
            if (
                data is not None
                and ver == 'diff'
                and self.diffswitch[f'{kind}_abs'].value
            ):
                data = np.abs(data)
            if data is None:
                slices = (
                    np.full((self.ny, self.nz), np.nan),
                    np.full((self.nx, self.nz), np.nan),
                    np.full((self.nx, self.ny), np.nan),
                )
                value_title = 'Axial'
            else:
                slices = (data[x, :, :], data[:, y, :], data[:, :, z])
                value = data[x, y, z]
                label = {'cvr': 'CVR', 'lag': 'Lag', 'r2': 'R²'}[kind]
                suffix = (
                    '%'
                    if data_key.endswith('_reld')
                    else ('s' if kind == 'lag' else '')
                )
                value_title = (
                    f'{label}: {value:.2f}{suffix}'
                    if np.isfinite(value)
                    else f'{label}: Masked'
                )

            settings = self.map_settings[map_name]
            for name, slice_data, cx, cy, title in zip(
                ('sag', 'cor', 'axi'),
                slices,
                (y, x, x),
                (z, z, y),
                ('Sagittal', 'Coronal', value_title),
                strict=True,
            ):
                key = f'{map_name}_{name}'
                img = self.imgs[key]
                img.set_data(slice_data.T)
                if img.get_cmap().name != settings['cmap']:
                    img.set_cmap(
                        mpl.colormaps[settings['cmap']].with_extremes(bad='black')
                    )
                if img.get_clim() != tuple(settings['limits']):
                    img.set_clim(*settings['limits'])
                self.cross[key][0].set_xdata([cx, cx])
                self.cross[key][1].set_ydata([cy, cy])
                self.axes[key].set_title(title)

        self._syncing = True
        try:
            self.x_field.value = x
            self.y_field.value = y
            self.z_field.value = z
        finally:
            self._syncing = False
        for chart in self.charts.values():
            chart.draw()
        self._draw_summaries()

    def _draw_summaries(self):
        """Draw masked distributions and the selected voxel's 4D R² profiles."""
        for kind, chart in self.summary_charts.items():
            fig = Figure(figsize=(1, 1), layout='tight')
            axis = fig.add_subplot()
            if kind in ('cvr', 'lag'):
                values = {}
                for ver in ('base', 'comp'):
                    data = self.data[f'{kind}_{ver}']
                    if data is not None:
                        values[ver] = data[np.isfinite(data)]
                populated = [value for value in values.values() if value.size]
                if populated:
                    edges = np.histogram_bin_edges(np.concatenate(populated), bins=30)
                    for ver, color in (('base', 'tab:blue'), ('comp', 'tab:orange')):
                        if ver in values and values[ver].size:
                            axis.hist(
                                values[ver], bins=edges, histtype='step', color=color,
                                label='Base' if ver == 'base' else 'Comparison',
                            )
                axis.set_xlabel('CVR' if kind == 'cvr' else 'Lag (s)')
                axis.set_ylabel('Voxels')
                axis.set_title('Masked distribution')
            else:
                for ver, color in (('base', 'tab:blue'), ('comp', 'tab:orange')):
                    volume = self.r2_volumes[ver]
                    mask = self.data[f'mask_{ver}']
                    if volume is not None and (mask is None or mask[tuple(self.curr)]):
                        profile = volume[tuple(self.curr)]
                        axis.plot(
                            np.arange(profile.size), profile, color=color,
                            label='Base' if ver == 'base' else 'Comparison',
                        )
                axis.set_xlabel('Lag index')
                axis.set_ylabel('R²')
                axis.set_title(f'R² at voxel {tuple(self.curr)}')
            if axis.get_legend_handles_labels()[0]:
                axis.legend()
            chart.figure = fig

    # ----- callbacks -----

    def _move_to(self, x, y, z):
        """Move to a voxel within the common spatial grid."""
        x = int(np.clip(x, 0, self.nx - 1))
        y = int(np.clip(y, 0, self.ny - 1))
        z = int(np.clip(z, 0, self.nz - 1))

        self.update_all(x, y, z)

    def on_canvas_click(self, event):
        """Select a voxel from any orthogonal slice."""
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
