"""Comparison viewer data and interaction tests."""

import importlib.util
import os
import sys
from types import SimpleNamespace
from unittest.mock import patch

import nibabel as nib
import numpy as np
import pytest
from matplotlib.axes import Axes
from matplotlib.backend_bases import MouseEvent

pytestmark = pytest.mark.skipif(
    importlib.util.find_spec('bootstack') is None,
    reason='Comparison viewer requires bootstack',
)

if importlib.util.find_spec('bootstack') is not None:
    import bootstack as bs

    from phys2cvr.gui.app import load_result_folder
    from phys2cvr.gui.pages.diffviewer import DiffViewer


def save_image(path, data, affine=None):
    nib.save(
        nib.Nifti1Image(
            np.asarray(data, dtype=float), np.eye(4) if affine is None else affine
        ),
        path,
    )
    return path


def test_comparison_data_masks_and_geometry(tmp_path):
    viewer = DiffViewer()
    base = np.arange(24, dtype=float).reshape(2, 3, 4)
    comp = base + 2
    viewer.set_cvr(save_image(tmp_path / 'base.nii.gz', base), 'base')
    viewer.set_cvr(save_image(tmp_path / 'comp.nii.gz', comp), 'comp')
    np.testing.assert_array_equal(viewer.data['cvr_diff'], base - comp)
    assert np.isnan(viewer.data['cvr_reld'][0, 0, 0])
    assert viewer.data['cvr_reld'][1, 0, 0] == pytest.approx(-200 / 12)

    mask = np.ones_like(base)
    mask[1, 0, 0] = 0
    viewer.set_mask(save_image(tmp_path / 'mask.nii.gz', mask), 'comp')
    assert np.isnan(viewer.data['cvr_diff'][1, 0, 0])
    assert np.isnan(viewer.data['cvr_base'][1, 0, 0])
    viewer.set_mask(
        save_image(tmp_path / 'full-mask.nii.gz', np.ones_like(base)), 'comp'
    )
    assert viewer.data['cvr_diff'][1, 0, 0] == -2

    r2 = np.stack((base, base + 3), axis=-1)
    r2[0, 0, 0, :] = np.nan
    r2[0, 0, 1, 0] = np.nan
    viewer.set_r2(save_image(tmp_path / 'r2.nii.gz', r2), 'base')
    np.testing.assert_array_equal(viewer.r2_volumes['base'], r2)
    assert np.isnan(viewer.data['r2_base'][0, 0, 0])
    assert np.isnan(viewer.data['r2_base'][0, 0, 1])
    viewer.set_r2(save_image(tmp_path / 'r2-3d.nii.gz', base), 'comp')
    assert viewer.data['r2_comp'].shape == base.shape
    assert viewer.r2_volumes['comp'] is None

    mismatched = save_image(tmp_path / 'wrong-shape.nii.gz', np.ones((3, 3, 4)))
    with pytest.raises(ValueError, match='spatial grid'):
        viewer.set_cvr(mismatched, 'comp')
    np.testing.assert_array_equal(viewer.data['cvr_comp'], comp)
    affine = np.eye(4)
    affine[0, 3] = 2
    shifted = save_image(tmp_path / 'wrong-affine.nii.gz', comp, affine)
    with pytest.raises(ValueError, match='spatial grid'):
        viewer.set_cvr(shifted, 'comp')


@pytest.mark.parametrize('masked_ver', ['base', 'comp'])
def test_shared_mask_updates_both_datasets(tmp_path, masked_ver):
    viewer = DiffViewer()
    shape = (2, 2, 2)
    paths = {}
    for ver, value in (('base', 3), ('comp', 1)):
        paths[ver] = save_image(tmp_path / f'{ver}.nii.gz', np.full(shape, value))
        viewer.set_cvr(paths[ver], ver)
        viewer.set_lag(paths[ver], ver)
    mask = np.ones(shape)
    mask[0, 0, 0] = 0
    viewer.set_mask(save_image(tmp_path / 'mask.nii.gz', mask), masked_ver)
    for kind in ('cvr', 'lag'):
        for ver in ('base', 'comp'):
            assert np.isnan(viewer.data[f'{kind}_{ver}'][0, 0, 0])
        assert np.isnan(viewer.data[f'{kind}_diff'][0, 0, 0])
        assert np.isnan(viewer.data[f'{kind}_reld'][0, 0, 0])

    mask[0, 0, 0] = 1
    mask[1, 0, 0] = 0
    viewer.set_mask(save_image(tmp_path / 'replacement.nii.gz', mask), masked_ver)
    for kind in ('cvr', 'lag'):
        assert viewer.data[f'{kind}_diff'][0, 0, 0] == 2
        assert np.isnan(viewer.data[f'{kind}_diff'][1, 0, 0])

    other = 'comp' if masked_ver == 'base' else 'base'
    other_mask = np.ones(shape)
    other_mask[0, 1, 0] = 0
    viewer.set_mask(save_image(tmp_path / 'other.nii.gz', other_mask), other)
    assert np.isnan(viewer.data[f'cvr_{masked_ver}'][1, 0, 0])
    assert np.isfinite(viewer.data[f'cvr_{other}'][1, 0, 0])
    assert np.isnan(viewer.data[f'cvr_{other}'][0, 1, 0])
    assert np.isfinite(viewer.data[f'cvr_{masked_ver}'][0, 1, 0])


@pytest.mark.parametrize('ver', ['base', 'comp'])
@pytest.mark.parametrize('ndim', [3, 4])
def test_r2_creates_mask_and_remasks_loaded_maps(tmp_path, ver, ndim):
    viewer = DiffViewer()
    shape = (2, 2, 2)
    base = np.full(shape, 4.0)
    comp = np.full(shape, 2.0)
    viewer.set_cvr(save_image(tmp_path / 'base.nii.gz', base), 'base')
    viewer.set_cvr(save_image(tmp_path / 'comp.nii.gz', comp), 'comp')
    r2 = np.ones(shape)
    r2[0, 0, 0] = 0
    r2[0, 0, 1] = np.nan
    if ndim == 4:
        r2 = np.stack((r2, r2), axis=-1)
        r2[0, 1, 0] = [-1, 1]
    viewer.set_r2(save_image(tmp_path / 'r2.nii.gz', r2), ver)
    mask = viewer.data[f'mask_{ver}']
    assert not mask[0, 0, 0]
    assert not mask[0, 0, 1]
    if ndim == 4:
        assert not mask[0, 1, 0]
    for key in ('cvr_base', 'cvr_comp', 'cvr_diff', 'cvr_reld', f'r2_{ver}'):
        assert np.isnan(viewer.data[key][0, 0, 0])
    assert viewer.data['cvr_diff'][1, 1, 1] == 2

    full_mask = np.ones(shape)
    viewer.set_mask(save_image(tmp_path / 'full.nii.gz', full_mask), ver)
    assert viewer.data['cvr_diff'][0, 0, 0] == 2
    assert viewer.data[f'r2_{ver}'][0, 0, 0] == 0
    viewer.set_r2(save_image(tmp_path / 'r2-again.nii.gz', r2), ver)
    assert viewer.data[f'mask_{ver}'][0, 0, 0]


def test_generated_mask_prevents_incompatible_r2_replacement(tmp_path):
    viewer = DiffViewer()
    initial = np.ones((2, 2, 2, 2))
    viewer.set_r2(save_image(tmp_path / 'r2.nii.gz', initial), 'base')
    with pytest.raises(ValueError, match='spatial grid'):
        viewer.set_r2(
            save_image(tmp_path / 'wrong.nii.gz', np.ones((3, 3, 3, 2))), 'base'
        )
    np.testing.assert_array_equal(viewer.r2_volumes['base'], initial)
    np.testing.assert_array_equal(viewer.data['mask_base'], np.ones((2, 2, 2)))


def test_comparison_folder_loading(tmp_path):
    base_dir = tmp_path / 'base'
    comp_dir = tmp_path / 'comp'
    base_dir.mkdir()
    comp_dir.mkdir()
    for folder, value in ((base_dir, 4), (comp_dir, 6)):
        save_image(folder / 'test_cvr.nii.gz', np.full((2, 3, 4), value))
        save_image(folder / 'test_lag.nii.gz', np.full((2, 3, 4), value))
        save_image(folder / 'test_r_square_simple.nii.gz', np.full((2, 3, 4), value))
        (folder / 'regr').mkdir()
        (folder / 'regr' / 'unused.txt').write_text('0')
        (folder / 'logs').mkdir()
        mask = np.ones((2, 3, 4))
        if folder == comp_dir:
            mask[0, 0, 0] = 0
        save_image(folder / 'mask.nii.gz', mask)
        (folder / 'logs' / 'p2c_call_001.sh').write_text(
            'phys2cvr -i unused.nii.gz -m mask.nii.gz\n'
        )
    viewer = DiffViewer()
    load_result_folder(str(base_dir), viewer, comparison_version='base')
    load_result_folder(str(comp_dir), viewer, comparison_version='comp')
    for kind in ('cvr', 'lag', 'r2'):
        difference = viewer.data[f'{kind}_diff']
        assert np.isnan(difference[0, 0, 0])
        np.testing.assert_array_equal(difference[1], -2)


@pytest.mark.skipif(
    sys.platform == 'linux' and not os.environ.get('DISPLAY'),
    reason='Requires a display',
)
def test_comparison_gui_interactions(tmp_path):
    viewer = DiffViewer()
    menus = []
    add_menu = bs.Toolbar.add_menu

    def capture_menu(toolbar, *args, **kwargs):
        menu = add_menu(toolbar, *args, **kwargs)
        menus.append(menu)
        return menu

    with bs.AppShell(title='Comparison test', size=(1600, 1000)) as shell:
        with shell.page_nav() as nav:
            with nav.add_page('comparison', text='Comparison'):
                with patch.object(bs.Toolbar, 'add_menu', capture_menu):
                    viewer.build()
        shell.navigate('comparison')
    shell.show()
    shell._internal.update()
    try:
        assert len(viewer.charts) == 9
        assert set(viewer.summary_charts) == {'cvr', 'lag', 'r2'}
        for col, kind in enumerate(('cvr', 'lag', 'r2'), start=1):
            chart = viewer.summary_charts[kind]._internal
            assert chart.grid_info()['row'] == 3
            assert chart.grid_info()['column'] == col
        assert set(viewer.diffswitch) == {
            f'{kind}_{mode}' for kind in ('cvr', 'lag', 'r2') for mode in ('rel', 'abs')
        }
        shape = (3, 4, 5)
        base = np.arange(1, 61, dtype=float).reshape(shape)
        comp = base + 3
        base_path = save_image(tmp_path / 'base.nii.gz', base)
        comp_path = save_image(tmp_path / 'comp.nii.gz', comp)
        viewer.set_cvr(base_path, 'base')
        viewer.set_lag(base_path, 'base')
        viewer.set_r2(
            save_image(tmp_path / 'r2-base-4d.nii.gz', np.stack((base, base + 1), -1)),
            'base',
        )
        viewer.refresh()
        assert len(viewer.summary_charts['cvr'].figure.axes[0].patches) == 1
        assert len(viewer.summary_charts['lag'].figure.axes[0].patches) == 1
        r2_axis = viewer.summary_charts['r2'].figure.axes[0]
        np.testing.assert_array_equal(
            r2_axis.lines[0].get_ydata(), [base[1, 2, 2], base[1, 2, 2] + 1]
        )
        shared_mask = np.ones(shape)
        shared_mask[0, 0, 0] = 0
        viewer.set_mask(save_image(tmp_path / 'shared-mask.nii.gz', shared_mask), 'base')
        viewer.set_cvr(comp_path, 'comp')
        viewer.set_lag(comp_path, 'comp')
        viewer.refresh()
        assert viewer.data['mask_comp'] is None
        assert np.isnan(viewer.data['cvr_comp'][0, 0, 0])
        viewer._move_to(0, 0, 0)
        assert np.ma.is_masked(viewer.imgs['cvr_comp_axi'].get_array()[0, 0])
        assert not viewer.summary_charts['r2'].figure.axes[0].lines

        viewer.set_r2(
            save_image(tmp_path / 'r2-comp-4d.nii.gz', np.stack((comp, comp + 1), -1)),
            'comp',
        )
        viewer.set_cvr(comp_path, 'comp')
        viewer.set_lag(comp_path, 'comp')
        mask = np.ones(shape)
        mask[0, 0, 0] = 0
        viewer.set_mask(save_image(tmp_path / 'mask.nii.gz', mask), 'base')
        histogram = Axes.hist
        samples = []

        def capture_hist(axis, values, *args, **kwargs):
            samples.append(np.asarray(values).copy())
            return histogram(axis, values, *args, **kwargs)

        with patch.object(Axes, 'hist', capture_hist):
            viewer.refresh()
        assert [len(values) for values in samples] == [59, 60, 59, 60]
        assert base[0, 0, 0] not in samples[0]
        assert len(viewer.summary_charts['cvr'].figure.axes[0].patches) == 2
        assert len(viewer.summary_charts['lag'].figure.axes[0].patches) == 2
        hist_figures = {
            kind: viewer.summary_charts[kind].figure for kind in ('cvr', 'lag')
        }
        with patch.object(Axes, 'hist', side_effect=AssertionError('recomputed')):
            viewer._move_to(0, 0, 0)
        assert all(
            viewer.summary_charts[kind].figure is figure
            for kind, figure in hist_figures.items()
        )
        r2_axis = viewer.summary_charts['r2'].figure.axes[0]
        assert len(r2_axis.lines) == 1
        np.testing.assert_array_equal(
            r2_axis.lines[0].get_ydata(), [comp[0, 0, 0], comp[0, 0, 0] + 1]
        )
        viewer._move_to(1, 1, 1)
        assert len(viewer.summary_charts['r2'].figure.axes[0].lines) == 2
        viewer.set_r2(base_path, 'base')
        viewer.refresh()
        assert len(viewer.summary_charts['r2'].figure.axes[0].lines) == 1
        for menu, path in zip(menus, (base_path, comp_path), strict=True):
            with patch.object(bs, 'ask_open_file', return_value=str(path)):
                for item in menu.items[1:4]:
                    item.on_click()
        shell._internal.update()
        assert len(viewer.axes) == 27
        assert len({id(chart.figure) for chart in viewer.charts.values()}) == 9
        for col, kind in enumerate(('cvr', 'lag', 'r2'), start=1):
            for row, ver in enumerate(('base', 'comp', 'diff')):
                chart = viewer.charts[f'{kind}_{ver}']._internal
                assert chart.grid_info()['row'] == row
                assert chart.grid_info()['column'] == col
                assert chart.winfo_width() > 50
                assert chart.winfo_height() > 50

        chart = viewer.charts['lag_comp']
        chart.figure.canvas.draw()
        axis = viewer.axes['lag_comp_sag']
        px, py = axis.transData.transform((1, 2))
        click = MouseEvent('button_press_event', chart.figure.canvas, px, py, button=1)
        chart.figure.canvas.callbacks.process('button_press_event', click)
        assert viewer.curr[1:] == [1, 2]

        with patch.object(bs, 'ask_directory', return_value=str(tmp_path)):
            with patch('phys2cvr.gui.app.load_result_folder') as load_folder:
                menus[0].items[0].on_click()
                menus[1].items[0].on_click()
        assert [
            call.kwargs['comparison_version'] for call in load_folder.call_args_list
        ] == [
            'base',
            'comp',
        ]

        # Exercise events from each chart and each orthogonal orientation.
        for name, axis in viewer.axes.items():
            orientation = name.rsplit('_', 1)[-1]
            viewer._move_to(0, 0, 0)
            viewer.on_canvas_click(SimpleNamespace(inaxes=axis, xdata=1, ydata=2))
            expected = {'sag': [0, 1, 2], 'cor': [1, 0, 2], 'axi': [1, 2, 0]}[
                orientation
            ]
            assert viewer.curr == expected
            assert [
                viewer.x_field.value,
                viewer.y_field.value,
                viewer.z_field.value,
            ] == expected
            for other, (vertical, horizontal) in viewer.cross.items():
                view = other.rsplit('_', 1)[-1]
                cx, cy = {
                    'sag': (expected[1], expected[2]),
                    'cor': (expected[0], expected[2]),
                    'axi': (expected[0], expected[1]),
                }[view]
                assert list(vertical.get_xdata()) == [cx, cx]
                assert list(horizontal.get_ydata()) == [cy, cy]

        viewer._move_to(1, 1, 1)
        for kind in ('cvr', 'lag', 'r2'):
            name = f'{kind}_diff_axi'
            x, y, z = viewer.curr
            expected = viewer.data[f'{kind}_diff'][:, :, z].T
            np.testing.assert_array_equal(viewer.imgs[name].get_array(), expected)
            viewer.diffswitch[f'{kind}_rel'].toggle()
            np.testing.assert_allclose(
                viewer.imgs[name].get_array(), viewer.data[f'{kind}_reld'][:, :, z].T
            )
            viewer.diffswitch[f'{kind}_abs'].toggle()
            np.testing.assert_allclose(
                viewer.imgs[name].get_array(),
                np.abs(viewer.data[f'{kind}_reld'][:, :, z].T),
            )
            viewer.diffswitch[f'{kind}_rel'].toggle()
            np.testing.assert_allclose(viewer.imgs[name].get_array(), np.abs(expected))
            viewer.diffswitch[f'{kind}_abs'].toggle()
            np.testing.assert_array_equal(viewer.imgs[name].get_array(), expected)
            assert not any(switch.value for switch in viewer.diffswitch.values())

        viewer.diffswitch['cvr_rel'].toggle()
        viewer.diffswitch['lag_abs'].toggle()
        viewer.diffswitch['r2_rel'].toggle()
        viewer.diffswitch['r2_abs'].toggle()
        np.testing.assert_allclose(
            viewer.imgs['cvr_diff_axi'].get_array(), viewer.data['cvr_reld'][:, :, z].T
        )
        np.testing.assert_allclose(
            viewer.imgs['lag_diff_axi'].get_array(),
            np.abs(viewer.data['lag_diff'][:, :, z].T),
        )
        np.testing.assert_allclose(
            viewer.imgs['r2_diff_axi'].get_array(),
            np.abs(viewer.data['r2_reld'][:, :, z].T),
        )

        viewer.x_field.value = 2
        viewer.y_field.value = 3
        viewer.z_field.value = 4
        viewer.z_field._entry_widget().event_generate('<<Change>>')
        assert viewer.curr == [2, 3, 4]
        viewer.map_fields['cvr_comp'][0].value = -5
        viewer.map_fields['cvr_comp'][0]._entry_widget().event_generate('<<Change>>')
        assert viewer.imgs['cvr_comp_axi'].get_clim()[0] == -5
        assert viewer.imgs['cvr_base_axi'].get_clim()[0] != -5
        viewer.map_fields['lag_diff'][2].value = 'plasma'
        viewer._on_map_settings_change('lag_diff')
        assert viewer.imgs['lag_diff_sag'].get_cmap().name == 'plasma'
        assert viewer.imgs['lag_base_sag'].get_cmap().name != 'plasma'
        old_limits = viewer.map_settings['cvr_comp']['limits']
        viewer.map_fields['cvr_comp'][0].value = old_limits[1] + 1
        with patch.object(bs, 'alert') as alert:
            viewer._on_map_settings_change('cvr_comp')
        alert.assert_called_once()
        assert viewer.map_settings['cvr_comp']['limits'] == old_limits
        assert viewer.map_fields['cvr_comp'][0].value == old_limits[0]
    finally:
        shell.destroy()
