"""Tests for the log viewer page."""

import importlib.util
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from PIL import Image

GUI_AVAILABLE = importlib.util.find_spec('bootstack') is not None
HAS_DISPLAY = sys.platform != 'linux' or bool(os.environ.get('DISPLAY'))

if GUI_AVAILABLE:
    import bootstack as bs

    from phys2cvr.gui.pages.logviewer import MAX_VISIBLE_LINES, build_log_viewer


@unittest.skipUnless(
    GUI_AVAILABLE and HAS_DISPLAY, 'GUI requires bootstack and a display'
)
class TestLogViewer(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        menus = []
        add_menu = bs.Toolbar.add_menu

        def capture_menu(toolbar, *args, **kwargs):
            menu = add_menu(toolbar, *args, **kwargs)
            menus.append(menu)
            return menu

        with bs.AppShell(title='Log viewer test', size=(1200, 850)) as cls.shell:
            with cls.shell.page_nav() as nav:
                with nav.add_page('log', text='Log'):
                    with patch.object(bs.Toolbar, 'add_menu', capture_menu):
                        cls.viewer = build_log_viewer()
                        cls.log_area = cls.viewer.log_area
            cls.shell.navigate('log')
        cls.load_log = staticmethod(menus[0].items[0].on_click)
        cls.shell.show()
        cls.shell._internal.update()

    @classmethod
    def tearDownClass(cls):
        cls.shell.destroy()

    def test_load_and_resize(self):
        text_widget = self.log_area._text_widget()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'run.log'
            for text, height in (
                ('first line\nsecond line\n', 2),
                ('single line', 1),
                ('', 1),
                ('long line ' * 200 + '\n' * 50, MAX_VISIBLE_LINES),
            ):
                path.write_text(text, encoding='utf-8')
                with patch.object(bs, 'ask_open_file', return_value=str(path)):
                    self.load_log()
                self.shell._internal.update()
                self.assertEqual(self.log_area.value, text)
                self.assertEqual(int(text_widget.cget('height')), height)
                self.assertEqual(text_widget.cget('state'), 'disabled')

            with patch.object(bs, 'ask_open_file', return_value=''):
                self.load_log()
            self.assertEqual(self.log_area.value, text)

            with (
                patch.object(bs, 'ask_open_file', return_value=str(path / 'missing')),
                patch.object(bs, 'alert') as alert,
                self.assertLogs('phys2cvr.gui.pages.logviewer', level='ERROR'),
            ):
                self.load_log()
            alert.assert_called_once()
            self.assertEqual(self.log_area.value, text)

            path.write_bytes(b'\xff')
            with (
                patch.object(bs, 'ask_open_file', return_value=str(path)),
                patch.object(bs, 'alert') as alert,
                self.assertLogs('phys2cvr.gui.pages.logviewer', level='ERROR'),
            ):
                self.load_log()
            alert.assert_called_once()
            self.assertEqual(self.log_area.value, text)

        for width in (1000, 1400):
            self.shell._internal.geometry(f'{width}x850')
            self.shell._internal.update()
            grid = self.log_area._internal.master
            area = self.log_area._internal
            self.assertGreater(area.winfo_width(), 200)
            self.assertLessEqual(area.winfo_height(), grid.winfo_height())
            self.assertTrue(self.log_area._internal._core._vscroll.winfo_ismapped())
            self.assertTrue(self.log_area._internal._core._hscroll.winfo_ismapped())

    def test_carousel_and_font(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            logs = root / 'logs'
            logs.mkdir()
            log = logs / 'run.tsv'
            log.write_text('log content', encoding='utf-8')
            for name in ('one.png', 'two.png'):
                Image.new('RGB', (12, 12), 'red').save(root / name)
            self.viewer.load(log)
            self.shell._internal.update()
            self.assertEqual(self.viewer.carousel.count, 2)
            self.assertEqual(self.viewer.carousel.current['name'], 'one.png')
            self.viewer.carousel.next()
            self.assertEqual(self.viewer.carousel.current['name'], 'two.png')
            self.viewer.font_slider.value = 20
            self.shell._internal.update()
            self.assertEqual(self.viewer.font.cget('size'), 20)
            other = root / 'other'
            other.mkdir()
            other_log = other / 'log.tsv'
            other_log.write_text('new content', encoding='utf-8')
            self.viewer.load(other_log)
            self.assertEqual(self.viewer.carousel.count, 0)
