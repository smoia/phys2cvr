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

    import app
    from guiplot import VoxelViewer


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
            self.assertEqual(self.viewer.lag_range, [-1, 1])
            self.assertTrue(np.isnan(self.viewer.data['cvr'][0, 0, 0]))
            self.assertEqual(self.viewer.curr, [2, 2, 3])
            event = SimpleNamespace(
                inaxes=self.viewer.axes['axi'],
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


if __name__ == '__main__':
    unittest.main()
