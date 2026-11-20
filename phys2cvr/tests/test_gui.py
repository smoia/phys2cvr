"""Tests for the optional GUI and its page-scoped controls."""

import importlib.util
import os
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import nibabel as nib
import numpy as np
from matplotlib.backend_bases import MouseEvent

GUI_AVAILABLE = importlib.util.find_spec('bootstack') is not None
HAS_DISPLAY = sys.platform != 'linux' or bool(os.environ.get('DISPLAY'))

if GUI_AVAILABLE:
    import bootstack as bs

    from phys2cvr.gui import app
    from phys2cvr.gui.pages.voxelviewer import VoxelViewer


@unittest.skipUnless(
    GUI_AVAILABLE and HAS_DISPLAY, 'GUI requires bootstack and a display'
)
class TestViewer(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.viewer = VoxelViewer()
        cls.shell = app.build_app(cls.viewer)
        with cls.shell.page_nav() as nav:
            with nav.add_page('other', text='Other'):
                bs.Label('Other page')
        cls.shell.show()
        cls.shell._internal.update()

    @classmethod
    def tearDownClass(cls):
        cls.shell.destroy()

    def setUp(self):
        self.viewer.fs = None
        self.viewer.zscore_switch.value = True

    def test_viewer_menu_follows_page(self):
        toolbar = self.viewer.toolbar._internal
        self.assertTrue(toolbar.winfo_ismapped())
        self.shell.navigate('other')
        self.shell._internal.update()
        self.assertFalse(toolbar.winfo_ismapped())
        self.shell.navigate('voxelviewer')
        self.shell._internal.update()
        self.assertTrue(toolbar.winfo_ismapped())

    def test_loading_and_interaction(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            rng = np.random.default_rng(42)
            timeseries = rng.normal(size=(4, 5, 6, 8))
            cvr = rng.normal(size=(4, 5, 6))
            lag = np.tile(np.arange(-1, 2), 40).reshape(4, 5, 6).astype(float)
            mask = np.ones((4, 5, 6))
            mask[0, 0, 0] = 0
            for name, data in (
                ('timeseries', timeseries),
                ('cvr', cvr),
                ('lag', lag),
                ('mask', mask),
            ):
                nib.save(nib.Nifti1Image(data, np.eye(4)), tmp / f'{name}.nii.gz')
            np.savetxt(tmp / 'matrix.txt', rng.normal(size=(3, 8)))
            self.viewer.set_timeseries(tmp / 'timeseries.nii.gz')
            self.viewer.set_mask(tmp / 'mask.nii.gz')
            self.viewer.set_cvr(tmp / 'cvr.nii.gz')
            self.viewer.set_lag(tmp / 'lag.nii.gz')
            self.viewer.set_regressors(tmp / 'matrix.txt')
            self.viewer.refresh()
            self.assertEqual(self.viewer.lag_range, (-1, 1))
            self.assertTrue(np.isnan(self.viewer.data['cvr'][0, 0, 0]))
            self.assertEqual(self.viewer.curr, [2, 2, 3])
            event = SimpleNamespace(
                inaxes=self.viewer.axes['cvr_axi'],
                xdata=1,
                ydata=3,
            )
            self.viewer.on_canvas_click(event)
            self.assertEqual(self.viewer.curr, [1, 3, 3])
            self.viewer.slider.value = 1
            self.viewer._on_slider_move()
            self.assertEqual(self.viewer.delay_label.text, '+1.00 s')

    def test_application_shell(self):
        self.assertEqual(self.shell.current, 'voxelviewer')

    def test_slider_without_timeseries(self):
        viewer = self.viewer
        viewer.timeseries = None
        viewer.time = None
        viewer.data = dict.fromkeys(viewer.data_types)
        viewer.lag_range = None
        viewer.regressor_matrix = np.array(
            [[1.0, 3.0, 2.0], [2.0, 3.0, 1.0], [3.0, 1.0, 2.0]]
        )
        viewer.refresh()
        viewer.slider.value = 1
        viewer._on_slider_move()
        np.testing.assert_allclose(
            viewer.line_co2.get_ydata(), [1.22474487, -1.22474487, 0], atol=1e-6
        )
        self.assertFalse(viewer.slider.disabled)

        viewer.regressor_matrix = None
        viewer.data['r2'] = np.tile([0.1, 0.4, 0.8], (4, 5, 6, 1))
        viewer.curr = [0, 0, 0]
        viewer.refresh()
        viewer.slider.value = 1
        viewer._on_slider_move()
        self.assertEqual(viewer.axes['r2'].get_title(), 'R²: 0.4000')
        self.assertFalse(viewer.slider.disabled)

    def test_regressor_slider_uses_frequency_and_inferred_frequency(self):
        viewer = self.viewer
        viewer.data = dict.fromkeys(viewer.data_types)
        viewer.timeseries = None
        viewer.time = None
        viewer.regressor_matrix = np.arange(27.0).reshape(9, 3)
        viewer.lag_range = (-2.0, 2.0)
        for frequency in (None, 2.0):
            viewer.fs = frequency
            viewer.refresh()
            self.assertAlmostEqual(viewer.slider.step, 0.5)
            self.assertAlmostEqual(viewer._regressor_fs(), 2.0)
            np.testing.assert_array_equal(viewer._regressor_row(-0.5)[0], [9, 10, 11])
            viewer.slider.value = -0.5
            viewer._on_slider_move()
            self.assertAlmostEqual(viewer.curr_shift_sec, -0.5)

    def test_r2_uses_closest_sample_to_regressor_step(self):
        viewer = self.viewer
        viewer.data = dict.fromkeys(viewer.data_types)
        viewer.data['r2'] = np.tile([0.1, 0.4, 0.8], (4, 5, 6, 1))
        viewer.curr = [1, 1, 1]
        viewer.timeseries = None
        viewer.time = None
        viewer.fs = None
        viewer.regressor_matrix = np.arange(27.0).reshape(9, 3)
        viewer.lag_range = (-2.0, 2.0)
        viewer.refresh()
        self.assertAlmostEqual(viewer.slider.step, 0.5)
        viewer.slider.value = 1.5
        viewer._on_slider_move()
        self.assertAlmostEqual(viewer.curr_shift_sec, 1.5)
        self.assertEqual(viewer.axes['r2'].get_title(), 'R²: 0.8000')
        np.testing.assert_array_equal(viewer.r2_point.get_xdata(), [2])
        np.testing.assert_array_equal(viewer.r2_point.get_ydata(), [0.8])
        viewer.update_all(*viewer.curr, 0.7)
        self.assertAlmostEqual(viewer.curr_shift_sec, 0.5)

    def test_edit_lag_range_infers_frequency_and_rejects_excess_span(self):
        viewer = self.viewer
        viewer.data = dict.fromkeys(viewer.data_types)
        viewer.timeseries = None
        viewer.time = None
        viewer.regressor_matrix = np.arange(27.0).reshape(9, 3)
        viewer.fs = None
        viewer.lag_range = (-2.0, 2.0)
        viewer.refresh()

        def change_range(lower, upper):
            viewer._syncing = True
            viewer.lagmin_field.value = lower
            viewer.lagmax_field.value = upper
            viewer._syncing = False
            viewer._on_lagrange_change()

        change_range(-1, 1)
        self.assertEqual(viewer.lag_range, (-1, 1))
        self.assertAlmostEqual(viewer.slider.step, 0.25)
        np.testing.assert_array_equal(viewer._regressor_row(1)[0], [24, 25, 26])

        viewer.fs = 2.0
        change_range(-1, 1)
        self.assertAlmostEqual(viewer.slider.step, 0.5)
        np.testing.assert_array_equal(viewer._regressor_row(1)[0], [12, 13, 14])
        with patch('phys2cvr.gui.pages.voxelviewer.bs.alert') as alert:
            change_range(-2, 3)
        alert.assert_called_once()
        self.assertIn('too big', alert.call_args.args[0])
        self.assertEqual(viewer.lag_range, (-1, 1))
        self.assertEqual((viewer.slider.min_value, viewer.slider.max_value), (-1, 1))
        self.assertEqual(
            (viewer.lagmin_field.value, viewer.lagmax_field.value), (-1, 1)
        )

    def test_sidebar_map_settings_and_normalization(self):
        viewer = self.viewer
        viewer.data = dict.fromkeys(viewer.data_types)
        viewer.timeseries = np.tile([2.0, 4.0, 6.0], (2, 2, 2, 1))
        viewer.nx = viewer.ny = viewer.nz = 2
        viewer.nt = 3
        viewer.time = np.arange(3)
        viewer.data['meanvol'] = np.mean(viewer.timeseries, axis=-1)
        viewer.data['cvr'] = np.arange(8.0).reshape(2, 2, 2)
        viewer.data['lag'] = np.arange(8.0).reshape(2, 2, 2)
        viewer.regressor_matrix = np.array([[10.0, 14.0, 18.0]] * 3)
        viewer.lag_range = (-1, 1)
        viewer.refresh()

        for kind in ('cvr', 'lag', 'mean'):
            lower, upper, cmap = viewer.map_fields[kind]
            viewer._syncing = True
            lower.value, upper.value, cmap.value = -2.25, 5.75, 'plasma'
            viewer._syncing = False
            viewer._on_map_settings_change(kind)
            for orientation in ('sag', 'cor', 'axi'):
                image = viewer.imgs[f'{kind}_{orientation}']
                self.assertEqual(image.get_clim(), (-2.25, 5.75))
                self.assertEqual(image.get_cmap().name, 'plasma')
        viewer.update_all(0, 0, 0, 0)
        np.testing.assert_allclose(
            viewer.line_voxel.get_ydata(), [-1.22474487, 0, 1.22474487], atol=1e-6
        )
        viewer.zscore_switch._internal.invoke()
        self.shell._internal.update()
        self.assertFalse(viewer.zscore_switch.value)
        np.testing.assert_allclose(viewer.line_voxel.get_ydata(), [-2, 0, 2])
        np.testing.assert_allclose(viewer.line_co2.get_ydata(), [-4, 0, 4])
        self.assertEqual(viewer.axes['ts'].get_ylabel(), 'Demeaned')

    def test_separate_chart_canvases_and_shared_map_colorbars(self):
        viewer = self.viewer
        viewer.data = dict.fromkeys(viewer.data_types)
        viewer.data['cvr'] = np.arange(120.0).reshape(4, 5, 6)
        viewer.data['lag'] = np.zeros((4, 5, 6))
        viewer.data['meanvol'] = np.ones((4, 5, 6))
        viewer.data['r2'] = np.tile([0.1, 0.4, 0.8], (4, 5, 6, 1))
        viewer.timeseries = None
        viewer.time = None
        viewer.curr = [2, 2, 3]
        viewer.lag_range = (-1, 1)
        viewer.regressor_matrix = np.arange(9.0).reshape(3, 3)
        viewer.refresh()
        self.shell._internal.update()

        self.assertEqual(set(viewer.charts), {'cvr', 'lag', 'mean', 'ts', 'r2'})
        self.assertEqual(len({id(c.figure) for c in viewer.charts.values()}), 5)
        for kind in ('cvr', 'lag', 'mean'):
            fig = viewer.charts[kind].figure
            images = [
                viewer.imgs[f'{kind}_{orientation}']
                for orientation in ('sag', 'cor', 'axi')
            ]
            self.assertEqual(len(fig.axes), 3 if kind == 'mean' else 4)
            self.assertTrue(all(image.axes.figure is fig for image in images))
            self.assertTrue(all(image.norm is images[0].norm for image in images))
            if kind == 'mean':
                self.assertIsNone(images[-1].colorbar)
            else:
                self.assertIs(images[-1].colorbar.mappable, images[-1])

            axis = viewer.axes[f'{kind}_axi']
            fig.canvas.draw()
            px, py = axis.transData.transform((1, 1))
            event = MouseEvent('button_press_event', fig.canvas, px, py, button=1)
            fig.canvas.callbacks.process('button_press_event', event)
            self.assertEqual(viewer.curr[:2], [1, 1])
            self.assertEqual(list(viewer.cross['lag_axi'][0].get_xdata()), [1, 1])

        fig = viewer.charts['r2'].figure
        fig.canvas.draw()
        px, py = viewer.axes['r2'].transData.transform((1, 0.4))
        fig.canvas.callbacks.process(
            'button_press_event',
            MouseEvent('button_press_event', fig.canvas, px, py, button=1),
        )
        self.assertEqual(viewer.curr_shift_sec, 1)
        self.assertEqual(viewer.axes['r2'].get_title(), 'R²: 0.8000')


if __name__ == '__main__':
    unittest.main()
