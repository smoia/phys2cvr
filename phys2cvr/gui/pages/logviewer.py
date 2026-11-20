"""Log Viewer page."""

import logging

import bootstack as bs

LGR = logging.getLogger(__name__)

MIN_VISIBLE_LINES = 4
MAX_VISIBLE_LINES = 30


def build_log_viewer():
    """Build the Log Viewer page."""

    def _load(title):
        def run():
            path = bs.ask_open_file(title=title, file_types=[('All files', '*.*')])
            if not path:
                return
            try:
                with open(path, encoding='utf-8') as f:
                    text = f.read()
            except (OSError, UnicodeError) as err:
                LGR.error('Could not load log file %s: %s', path, err, exc_info=True)
                bs.alert(
                    f'Could not load log file {path}: {err}',
                    title='Could not load file',
                    severity='danger',
                )
                return

            log_area.value = text
            lines = max(1, min(MAX_VISIBLE_LINES, len(text.splitlines())))
            # Bootstack's TextArea exposes height only at construction.
            log_area._text_widget().configure(height=lines)

        return run

    with bs.Toolbar(horizontal='stretch', surface='chrome') as bar:
        with bar.add_menu('Load') as menu:
            menu.add_action(
                'Load logger',
                on_click=_load('Select log file'),
            )

    with bs.Grid(horizontal='stretch', grow=True):
        log_area = bs.TextArea(
            relx=0.25,
            relwidth=0.5,
            height=MIN_VISIBLE_LINES,
            label='Log',
            read_only=True,
            wrap=False,
            scrollbars='auto',
        )
    return log_area
