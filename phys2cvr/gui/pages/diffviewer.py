#!/usr/bin/env python3
"""Comparison viewer for phys2cvr, built with bootstack."""

import bootstack as bs
import nibabel as nib
import numpy as np
from matplotlib.figure import Figure

from phys2cvr.gui.utils import (
    NIFTI_TYPES,
    MapViewerControls,
    build_coordinate_fields,
    build_map_fields,
    clip_position,
    make_orthogonal_figure,
    orthogonal_slices,
    slice_position,
    update_map_view,
)


class DiffViewer(MapViewerControls):
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
                self.map_fields[kind] = build_map_fields(self.map_settings[kind], title)

                if 'diff' in kind:
                    self.diffswitch[f'{kind[:-4]}rel'] = bs.Switch(
                        'Relative diff', value=False
                    )
                    self.diffswitch[f'{kind[:-4]}abs'] = bs.Switch(
                        'Absolute diff', value=False
                    )

                bs.Divider(margin_y=(10, 5))

            self.x_field, self.y_field, self.z_field = build_coordinate_fields()

            bs.Spacer()

    def _update_current_plot(self):
        if not self._syncing and self.axes:
            self.update_all(*self.curr)

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

        self._set_coordinate_limits()
        self._draw_summaries(histograms=True)
        self._move_to(*self.curr)

    def init_plots(self):
        """Create one figure per map, each with three orthogonal views."""
        self.axes = {}
        self.imgs = {}
        self.cross = {}

        figures = {}
        for map_name in self.charts:
            kind = map_name.split('_')[0]
            color = 'red' if kind == 'lag' else 'cyan'
            figures[map_name] = make_orthogonal_figure(self, map_name, color)

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
                slices = orthogonal_slices(data, (x, y, z))
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
                update_map_view(self, key, slice_data, cx, cy, title, settings)

        self._sync_coordinates((x, y, z))
        for chart in self.charts.values():
            chart.draw()
        self._draw_summaries(histograms=False)

    def _draw_summaries(self, *, histograms=True):
        """Draw masked distributions and the selected voxel's 4D R² profiles."""
        for kind, chart in self.summary_charts.items():
            if kind != 'r2' and not histograms:
                continue
            fig = chart.figure
            axis = fig.axes[0] if fig.axes else fig.add_subplot()
            axis.clear()
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
                                values[ver],
                                bins=edges,
                                histtype='step',
                                color=color,
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
                            np.arange(profile.size),
                            profile,
                            color=color,
                            label='Base' if ver == 'base' else 'Comparison',
                        )
                axis.set_xlabel('Lag index')
                axis.set_ylabel('R²')
                axis.set_title(f'R² at voxel {tuple(self.curr)}')
            if axis.get_legend_handles_labels()[0]:
                axis.legend()
            chart.draw()

    # ----- callbacks -----

    def _move_to(self, x, y, z):
        """Move to a voxel within the common spatial grid."""
        x, y, z = clip_position((x, y, z), (self.nx, self.ny, self.nz))

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

        orientation = ('sag', 'cor', 'axi').index(name.rsplit('_', 1)[-1])
        x, y, z = slice_position((x, y, z), orientation, u, v)

        self._move_to(x, y, z)
