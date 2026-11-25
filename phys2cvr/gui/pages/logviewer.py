"""Log and workflow image viewer."""

import logging
import tkinter.font as tkfont

import bootstack as bs
from bootstack.data import MemoryDataSource

from phys2cvr.gui.io import read_log_and_images

LGR = logging.getLogger(__name__)
MIN_VISIBLE_LINES = 4
MAX_VISIBLE_LINES = 30


class LogViewer:
    """Display a log and images from its result directory."""

    def build(self):
        """Create the log page."""
        with bs.Toolbar(horizontal='stretch', surface='chrome') as bar:
            with bar.add_menu('Load') as menu:
                menu.add_action('Load logger', on_click=self._choose_log)
        with bs.Grid(
            columns=['auto', 1, 1], rows=1, grow=True,
            horizontal_items='stretch', vertical_items='stretch', gap=8,
        ):
            with bs.Column(column=0, row=0, vertical='center', padding=8, height=180):
                bs.Label('Font size')
                self.font_slider = bs.Slider(
                    12, min_value=6, max_value=30, step=1, orient='vertical',
                    grow=True,
                )
            self.log_area = bs.TextArea(
                row=0, column=1, height=MIN_VISIBLE_LINES, label='Log',
                read_only=True, wrap=False, scrollbars='auto',
            )
            self.images = MemoryDataSource()
            self.carousel = bs.Carousel(
                data_source=self.images, caption_field='name', transition='none',
                row=0, column=2,
            )
        text = self.log_area._text_widget()
        self.font = tkfont.Font(root=text, font=text.cget('font'))
        text.configure(font=self.font)
        self.font_slider.signal.subscribe(
            lambda value: self.font.configure(size=int(value)), immediate=True
        )
        return self

    def _choose_log(self):
        path = bs.ask_open_file(title='Select log file', file_types=[('All files', '*.*')])
        if not path:
            return
        try:
            self.load(path)
        except (OSError, UnicodeError) as err:
            LGR.error('Could not load log file %s: %s', path, err, exc_info=True)
            bs.alert(
                f'Could not load log file {path}: {err}',
                title='Could not load file', severity='danger',
            )

    def load(self, path):
        """Load a log and the PNGs from its result folder."""
        text, images = read_log_and_images(path)
        self.log_area.value = text
        lines = max(1, min(MAX_VISIBLE_LINES, len(text.splitlines())))
        self.log_area._text_widget().configure(height=lines)
        self.images.load([{'image': str(p), 'name': p.name} for p in images])
        self.carousel.reload()

    def clear(self):
        """Clear the log and images when no log exists in the result folder."""
        self.log_area.value = ''
        self.log_area._text_widget().configure(height=MIN_VISIBLE_LINES)
        self.images.load([])
        self.carousel.reload()


def build_log_viewer():
    """Build the log viewer page."""
    return LogViewer().build()
