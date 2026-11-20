"""Log Viewer page."""

import logging

import bootstack as bs

LGR = logging.getLogger(__name__)


def build_log_viewer():
    """Build the Log Viewer page."""
    logtext = bs.Signal('')
    num_lines = bs.Signal(4)

    def _load(title):
        def run():
            path = bs.ask_open_file(title=title, file_types=[('All files', '*.*')])
            if not path:
                return
            try:
                with open(path) as f:
                    num_lines.set(sum(1 for _ in f))
                    logtext.set(f.read())
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

        return run

    with bs.Toolbar(horizontal='stretch', surface='chrome') as bar:
        with bar.add_menu('Load') as menu:
            menu.add_action(
                'Load logger',
                on_click=_load('Select log file'),
            )

    with bs.Grid(columns=[1, 2, 1], horizontal_items='center', gap=8):
        bs.TextArea(
            column=1,
            textsignal=logtext,
            height=num_lines,
            label='Log',
            read_only=True,
            wrap=False,
        )
