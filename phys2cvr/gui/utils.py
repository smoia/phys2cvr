"""Shared controls and plotting helpers for the GUI pages."""

import logging

import bootstack as bs
import matplotlib as mpl
import numpy as np
from matplotlib.figure import Figure

from phys2cvr.io import EXT_ARRAY, EXT_NIMG

NIFTI_TYPES = [
    ('NIfTI files', ' '.join(f'*{ext}' for ext in EXT_NIMG)),
    ('All files', '*.*'),
]
ARRAY_TYPES = [
    ('Array files', '*.1D ' + ' '.join(f'*{ext}' for ext in EXT_ARRAY)),
    ('All files', '*.*'),
]


def orthogonal_slices(data, position):
    """Return sagittal, coronal, and axial slices at a voxel."""
    x, y, z = position
    return data[x, :, :], data[:, y, :], data[:, :, z]


def slice_position(position, orientation, horizontal, vertical):
    """Return the voxel selected in an orthogonal slice."""
    position = list(position)
    first, second = ((1, 2), (0, 2), (0, 1))[orientation]
    position[first], position[second] = horizontal, vertical
    return position


def clip_position(position, shape):
    """Constrain voxel coordinates to the spatial grid."""
    return [
        int(np.clip(value, 0, size - 1))
        for value, size in zip(position, shape, strict=True)
    ]


def build_map_fields(settings, title):
    """Create the common map range and color controls."""
    bs.Label(title, font='caption', margin_y=(5, 8))
    lower, upper = settings['limits']
    return (
        bs.NumberField(float(lower), label='Min', step=0.1, width=4),
        bs.NumberField(float(upper), label='Max', step=0.1, width=4),
        bs.Select(
            sorted(mpl.colormaps),
            value=settings['cmap'],
            label='Colorscale',
            width=10,
        ),
    )


def build_coordinate_fields():
    """Create the common voxel coordinate controls."""
    bs.Label('Voxel position', font='caption', margin_y=(5, 8))
    return tuple(
        bs.NumberField(0, label=label, min_value=0, width=4)
        for label in ('X', 'Y', 'Z')
    )


def make_orthogonal_figure(viewer, name, color, *, vertical=False, colorbar=True):
    """Create one map figure and register its images, axes and crosshairs."""
    nx, ny, nz = viewer.nx, viewer.ny, viewer.nz
    dx, dy, dz = viewer.pixdim
    aspects = (dz / dy, dz / dx, dy / dx)
    fig = Figure(figsize=(1, 1), layout='compressed')
    fig.get_layout_engine().set(w_pad=0.02, h_pad=0.02, wspace=0, hspace=0)
    grid = fig.add_gridspec(3 if vertical else 1, 1 if vertical else 3)
    settings = viewer.map_settings[name]
    cmap = mpl.colormaps[settings['cmap']].with_extremes(bad='black')
    norm = mpl.colors.Normalize(*settings['limits'])
    axes = []
    for index, (orientation, xlabel, ylabel, size) in enumerate(
        (
            ('sag', 'Y', 'Z', (ny, nz)),
            ('cor', 'X', 'Z', (nx, nz)),
            ('axi', 'X', 'Y', (nx, ny)),
        )
    ):
        key = f'{name}_{orientation}'
        axis = viewer.axes[key] = fig.add_subplot(
            grid[index, 0] if vertical else grid[0, index]
        )
        viewer.imgs[key] = axis.imshow(
            np.zeros(size[::-1]),
            cmap=cmap,
            norm=norm,
            origin='lower',
            aspect=aspects[index],
        )
        viewer.cross[key] = (
            axis.axvline(0, color=color, lw=1),
            axis.axhline(0, color=color, lw=1),
        )
        axis.set_xlabel(xlabel)
        axis.set_ylabel(ylabel)
        axis.set_xticks([])
        axis.set_yticks([])
        axis.set_xticklabels([])
        axis.set_yticklabels([])
        axis.set_anchor(('S', 'C', 'N')[index] if vertical else ('E', 'C', 'W')[index])
        axis.set_facecolor('black')
        axes.append(axis)
    if colorbar:
        fig.colorbar(viewer.imgs[f'{name}_axi'], ax=axes, fraction=0.035, pad=0.02)
    return fig


def update_map_view(viewer, name, data, x, y, title, settings):
    """Update an orthogonal image, its crosshairs, and its map settings."""
    image = viewer.imgs[name]
    image.set_data(data.T)
    if image.get_cmap().name != settings['cmap']:
        image.set_cmap(mpl.colormaps[settings['cmap']].with_extremes(bad='black'))
    if image.get_clim() != tuple(settings['limits']):
        image.set_clim(*settings['limits'])
    viewer.cross[name][0].set_xdata([x, x])
    viewer.cross[name][1].set_ydata([y, y])
    viewer.axes[name].set_title(title)


class MapViewerControls:
    """Shared map settings, coordinate controls, and file dialogs."""

    def _set_coordinate_limits(self):
        syncing = self._syncing
        self._syncing = True
        try:
            for field, size in zip(
                (self.x_field, self.y_field, self.z_field),
                (self.nx, self.ny, self.nz),
                strict=True,
            ):
                field.max_value = size - 1
        finally:
            self._syncing = syncing

    def _sync_coordinates(self, position):
        syncing = self._syncing
        self._syncing = True
        try:
            for field, value in zip(
                (self.x_field, self.y_field, self.z_field), position, strict=True
            ):
                field.value = value
        finally:
            self._syncing = syncing

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
            logging.getLogger(type(self).__module__).warning(
                'Invalid %s color range: %s to %s', kind, lower, upper
            )
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

    def _on_coord_change(self):
        if self._syncing or not self.axes:
            return
        values = [f.value for f in (self.x_field, self.y_field, self.z_field)]
        if any(v is None for v in values):
            return
        self._move_to(*values)

    def _load(self, title, file_types, setter, *args):
        path = bs.ask_open_file(title=title, file_types=file_types)
        if not path:
            return
        try:
            setter(path, *args)
            self.refresh()
        except Exception as err:  # noqa: BLE001
            logging.getLogger(type(self).__module__).error(
                'Fatal error encountered while loading %s: %s',
                path,
                err,
                exc_info=True,
            )
            bs.alert(
                f'Fatal error encountered while loading {path}: {err}',
                title='Could not load file',
                severity='danger',
            )

    def _load_folder(self, ver=None):
        folder = bs.ask_directory(title='Select phys2cvr Result Folder')
        if not folder:
            return
        try:
            from phys2cvr.gui.app import load_result_folder

            if ver is None:
                load_result_folder(folder, self)
            else:
                load_result_folder(folder, self, comparison_version=ver)
            self.refresh()
        except Exception as err:  # noqa: BLE001
            logging.getLogger(type(self).__module__).error(
                'Could not load result folder %s: %s', folder, err, exc_info=True
            )
            self.refresh()
            bs.alert(str(err), title='Could not load result folder')
