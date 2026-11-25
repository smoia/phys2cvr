"""Tests for GUI result and log loading without a display."""

import os
from unittest.mock import Mock

from phys2cvr.gui.io import load_result_folder, read_log_and_images


def test_read_log_and_images_searches_result_tree(tmp_path):
    logs = tmp_path / 'logs'
    logs.mkdir()
    log = logs / 'run.tsv'
    log.write_text('run log', encoding='utf-8')
    (tmp_path / 'z.png').touch()
    (logs / 'a.png').touch()
    (tmp_path / 'ignored.txt').touch()

    text, images = read_log_and_images(log)
    assert text == 'run log'
    assert images == sorted([tmp_path / 'z.png', logs / 'a.png'])


def test_result_folder_loads_newest_log_and_clears_missing(tmp_path):
    logs = tmp_path / 'logs'
    logs.mkdir()
    newer = logs / 'a.tsv'
    older = logs / 'z.tsv'
    newer.touch()
    older.touch()
    os.utime(older, (100, 100))
    os.utime(newer, (200, 200))
    viewer = Mock()
    logviewer = Mock()

    load_result_folder(tmp_path, viewer, logviewer=logviewer)
    logviewer.load.assert_called_once_with(str(newer))
    newer.unlink()
    older.unlink()
    load_result_folder(tmp_path, viewer, logviewer=logviewer)
    logviewer.clear.assert_called_once_with()
