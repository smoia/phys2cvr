"""Tests for GUI result and log loading without a display."""

import os
from unittest.mock import Mock

from phys2cvr.gui.io import load_result_folder, read_latest_calls, read_log_and_images


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
    logviewer.load_calls.assert_called_once_with(str(logs))


def test_latest_calls_use_modification_time_and_clear_absent_files(tmp_path):
    old = tmp_path / 'p2c_call_z.sh'
    new = tmp_path / 'p2c_call_a.sh'
    python = tmp_path / 'p2c_call_a.py'
    old.write_text('old shell', encoding='utf-8')
    new.write_text('new shell', encoding='utf-8')
    python.write_text('python call', encoding='utf-8')
    os.utime(old, (100, 100))
    os.utime(new, (200, 200))
    assert read_latest_calls(tmp_path) == ('new shell', 'python call')
    python.unlink()
    assert read_latest_calls(tmp_path) == ('new shell', '')


def test_result_folder_parses_quoted_input_and_skips_non_nifti(tmp_path):
    import shlex

    logs = tmp_path / 'logs'
    logs.mkdir()
    viewer = Mock()
    source = tmp_path / 'input file.nii.gz'
    source.touch()
    script = logs / 'p2c_call_a.sh'
    script.write_text(
        '#!/bin/bash\n' + shlex.join(['phys2cvr', '--input-func', str(source)]),
        encoding='utf-8',
    )
    load_result_folder(tmp_path, viewer)
    viewer.set_timeseries.assert_called_once_with(str(source))
    viewer.reset_mock()
    source = tmp_path / 'input file.txt'
    source.touch()
    script.write_text(
        '#!/bin/bash\n' + shlex.join(['phys2cvr', '--input-func', str(source)]),
        encoding='utf-8',
    )
    load_result_folder(tmp_path, viewer)
    viewer.set_timeseries.assert_not_called()
    script.write_text('#!/bin/bash\nphys2cvr --bad-option\n', encoding='utf-8')
    load_result_folder(tmp_path, viewer)
