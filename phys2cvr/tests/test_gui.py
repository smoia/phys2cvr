"""Tests for the optional GUI and its page-scoped controls."""

import importlib.util
import os
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import nibabel as nib
import numpy as np

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
        cls.shell.close()

    def test_viewer_menu_follows_page(self):
        toolbar = self.viewer.toolbar._internal
        self.assertTrue(toolbar.winfo_ismapped())
        self.shell.navigate('other')
        self.shell._internal.update()
        self.assertFalse(toolbar.winfo_ismapped())
        self.shell.navigate('viewer')
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
        self.assertEqual(self.shell.current, 'viewer')

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


if __name__ == '__main__':
    unittest.main()
