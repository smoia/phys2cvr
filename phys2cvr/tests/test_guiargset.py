"""Tests for configuring and running a workflow from the GUI."""

import argparse
import importlib.util
import os
import sys
import time
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import nibabel as nib
import numpy as np
import pytest

GUI_AVAILABLE = importlib.util.find_spec('bootstack') is not None
HAS_DISPLAY = sys.platform != 'linux' or bool(os.environ.get('DISPLAY'))
pytestmark = pytest.mark.skipif(
    not GUI_AVAILABLE or not HAS_DISPLAY, reason='GUI requires bootstack and a display'
)

if GUI_AVAILABLE:
    import bootstack as bs

    from phys2cvr.cli.run import _get_parser
    from phys2cvr.gui.pages.guiargset import build_workflow_setter


@pytest.fixture(scope='module')
def page():
    with bs.AppShell(title='Workflow test', size=(1100, 900)) as shell:
        with shell.page_nav() as nav:
            with nav.add_page('workflow', text='Workflow'):
                form = build_workflow_setter()
        shell.navigate('workflow')
    shell.show()
    shell._internal.update()
    yield form, shell
    shell.destroy()


@pytest.fixture(autouse=True)
def reset_form(page):
    form, _ = page
    for dest, action in form.actions.items():
        field = form.fields[dest]
        if dest in ('response_function', 'r2model', 'workflow_config'):
            field.value = ''
        elif isinstance(
            action, (argparse._StoreTrueAction, argparse._StoreFalseAction)
        ):
            field.value = False
        elif action.type in (int, float) and action.default is not None:
            field.value = action.default
        else:
            field.clear()
    form.response_path.clear()


def test_parser_fields_and_default_values(page):
    form, _ = page
    parser = _get_parser()
    expected = {
        action.dest
        for action in parser._actions
        if not isinstance(action, (argparse._HelpAction, argparse._VersionAction))
    }
    assert set(form.fields) == expected
    with pytest.raises(ValueError, match='required'):
        form.arguments()
    form.fields['fname_func'].value = '/tmp/functional input.nii.gz'
    assert vars(parser.parse_args(form.arguments())) == vars(
        parser.parse_args(['--input-func', '/tmp/functional input.nii.gz'])
    )
    assert form.fields['tr'].value is None


def test_argument_groups_use_groupboxes(page):
    form, _ = page
    expected = [
        group.title
        for group in form.parser._action_groups
        if any(
            not isinstance(action, (argparse._HelpAction, argparse._VersionAction))
            for action in group._group_actions
        )
    ]
    assert len(form.groups) == len(expected)
    assert all(isinstance(group, bs.GroupBox) for group in form.groups)


def test_input_previews_and_navigation(page, tmp_path):
    form, shell = page
    preview = form.preview
    shell._internal.update()
    left = form.form_scroll._internal.winfo_width()
    right = preview.container._internal.winfo_width()
    assert abs(left - right) <= 10
    for key in ('fname_func', 'fname_mask', 'fname_roi'):
        assert len(preview.axes[key]) == 3
        assert len(preview.charts[key].figure.axes) == 3

    shape = (3, 4, 5)
    volume = np.arange(np.prod(shape) * 2, dtype=float).reshape(*shape, 2)
    mask = np.ones(shape)
    mask[0, 0, 0] = 0
    paths = {}
    for key, values in (
        ('fname_func', volume), ('fname_mask', mask), ('fname_roi', mask * 2),
    ):
        paths[key] = tmp_path / f'{key}.nii.gz'
        nib.save(nib.Nifti1Image(values, np.eye(4)), paths[key])
        form.fields[key].value = str(paths[key])
        form.fields[key]._entry_widget().event_generate('<<Change>>')
    np.testing.assert_array_equal(preview.data['fname_func'], volume.mean(axis=-1))
    np.testing.assert_array_equal(preview.data['fname_mask'], mask.astype(bool))
    np.testing.assert_array_equal(preview.data['fname_roi'], mask.astype(bool))
    event = SimpleNamespace(
        name='button_press_event', inaxes=preview.axes['fname_func'][0],
        button=1, xdata=1, ydata=3,
    )
    preview._navigate('fname_func', event)
    assert preview.positions['fname_func'] == [1, 1, 3]
    np.testing.assert_array_equal(
        preview.axes['fname_func'][0].images[0].get_array(), volume.mean(axis=-1)[1].T
    )
    event.name, event.step = 'scroll_event', -1
    preview._navigate('fname_func', event)
    assert preview.positions['fname_func'] == [0, 1, 3]

    regressor = tmp_path / 'co2.txt'
    np.savetxt(regressor, [1, 3, 2])
    with patch.dict(
        'bootstack.widgets._impl.composites.pathentry._MODE_TO_DIALOG',
        {'open': lambda **kwargs: str(regressor)},
    ):
        form.fields['fname_co2']._internal._show_file_chooser()
    np.testing.assert_array_equal(preview.axes['fname_co2'][0].lines[0].get_ydata(), [1, 3, 2])
    form.fields['freq'].value = 2
    form.fields['freq']._entry_widget().event_generate('<<Change>>')
    np.testing.assert_array_equal(preview.axes['fname_co2'][0].lines[0].get_xdata(), [0, .5, 1])
    regressor2 = tmp_path / 'co2-new.txt'
    np.savetxt(regressor2, [6, 8])
    form.fields['fname_co2'].value = str(regressor2)
    form.fields['fname_co2']._entry_widget().event_generate('<<Change>>')
    np.testing.assert_array_equal(preview.axes['fname_co2'][0].lines[0].get_ydata(), [6, 8])
    form.fields['fname_co2'].clear()
    form.fields['fname_co2']._entry_widget().event_generate('<<Change>>')
    assert not preview.axes['fname_co2'][0].lines

    form.fields['fname_roi'].value = str(tmp_path / 'missing.nii.gz')
    form.fields['fname_roi']._entry_widget().event_generate('<<Change>>')
    assert preview.data['fname_roi'] is None
    assert preview.errors['fname_roi']


def test_all_argument_types_and_validation(page):
    form, _ = page
    form.fields['fname_func'].value = '/tmp/input file.nii.gz'
    form.fields['outdir'].value = '/tmp/output folder'
    form.fields['tr'].value = 1.5
    form.fields['n_jobs'].value = -2
    form.fields['run_regression'].value = True
    form.fields['apply_filter'].value = True
    form.fields['denoise_matrix_file'].value = '/tmp/a b.txt\n/tmp/c.txt'
    form.fields['response_function'].value = '--response-function-path'
    form.response_path.value = '/tmp/response file.txt'
    form.fields['r2model'].value = '--r2adjpartial'
    form.fields['workflow_config'].value = '--brightspin'
    args = _get_parser().parse_args(form.arguments())
    assert args.fname_func == '/tmp/input file.nii.gz'
    assert args.outdir == '/tmp/output folder'
    assert args.tr == 1.5 and args.n_jobs == -2
    assert not args.run_regression and args.apply_filter
    assert args.denoise_matrix_file == ['/tmp/a b.txt', '/tmp/c.txt']
    assert args.response_function == '/tmp/response file.txt'
    assert args.r2model == 'adj_partial'
    assert args.workflow_config == 'brightspin'

    entry = form.fields['tr']._entry_widget()
    entry.delete(0, 'end')
    entry.insert(0, '2.75')
    assert _get_parser().parse_args(form.arguments()).tr == 2.75
    entry.delete(0, 'end')
    entry.insert(0, 'invalid')
    with pytest.raises(ValueError, match='number'):
        form.arguments()
    form.fields['tr'].value = 1.5

    form.fields['n_jobs'].value = 1.5
    with pytest.raises(ValueError, match='integer'):
        form.arguments()
    form.fields['n_jobs'].value = 1
    form.response_path.clear()
    with pytest.raises(ValueError, match='custom response'):
        form.arguments()
    form.fields['response_function'].value = '--no-response-function'
    form.fields['workflow_config'].value = ''
    form.fields['lag_min'].value = -2
    with pytest.raises(ValueError, match='minimum lag'):
        form.arguments()
    form.fields['lag_min'].clear()


def test_multiple_file_picker(page):
    form, _ = page
    field = form.fields['extra_matrix_file']
    field.value = '/tmp/first matrix.txt'
    with patch.object(bs, 'ask_open_files', return_value=['/tmp/second matrix.txt']):
        form._add_files(field)
    assert field.value.splitlines() == [
        '/tmp/first matrix.txt',
        '/tmp/second matrix.txt',
    ]
    with patch.object(bs, 'ask_open_files', return_value=[]):
        form._add_files(field)
    assert len(field.value.splitlines()) == 2


def test_mutually_exclusive_options(page):
    form, _ = page
    form.fields['fname_func'].value = '/tmp/input.nii.gz'
    form.fields['lag_min'].clear()
    parser = _get_parser()
    for dest in ('response_function', 'r2model', 'workflow_config'):
        for action in parser._actions:
            if dest != action.dest or not isinstance(
                action, argparse._StoreConstAction
            ):
                continue
            option = next(opt for opt in action.option_strings if opt.startswith('--'))
            form.fields[dest].value = option
            assert getattr(parser.parse_args(form.arguments()), dest) == action.const
        form.fields[dest].value = ''


def test_invalid_input_and_process_start_failure(page):
    form, _ = page
    form.fields['fname_func'].clear()
    with patch.object(bs, 'alert') as alert, patch('subprocess.Popen') as popen:
        form.run()
    alert.assert_called_once()
    popen.assert_not_called()
    assert form.process is None

    form.fields['fname_func'].value = '/tmp/input.nii.gz'
    with (
        patch.object(bs, 'alert') as alert,
        patch('subprocess.Popen', side_effect=OSError('cannot start')),
    ):
        form.run()
    assert 'cannot start' in alert.call_args.args[0]
    assert form.process is None
    assert not form.run_button.disabled


def wait_for_workflow(form, shell):
    deadline = time.monotonic() + 30
    while form.process is not None and time.monotonic() < deadline:
        shell._internal.update()
        time.sleep(0.05)
    assert form.process is None, 'Workflow did not finish in time'


def test_real_workflow_and_error_reporting(page, tmp_path):
    form, shell = page
    signal = tmp_path / 'signal.txt'
    np.savetxt(signal, np.sin(np.arange(100) / 10) + 10)
    output = tmp_path / 'workflow output'
    form.fields['fname_func'].value = str(signal)
    form.fields['outdir'].value = str(output)
    form.fields['tr'].value = 1
    form.fields['run_regression'].value = True
    form.fields['lagged_regression'].value = True
    form.fields['skip_xcorr'].value = True
    form.fields['comp_endtidal'].value = True
    form.fields['apply_filter'].value = False
    form.fields['response_function'].value = '--no-response-function'
    for dest in (
        'denoise_matrix_file',
        'extra_matrix_file',
        'orthogonalised_matrix_file',
    ):
        form.fields[dest].clear()
    form.fields['r2model'].value = ''
    with patch.object(bs, 'alert') as alert:
        form.run_button._internal.invoke()
        assert form.process is not None
        assert form.run_button.disabled
        wait_for_workflow(form, shell)
    alert.assert_not_called()
    assert form.status.text == 'Workflow completed.'
    assert list((output / 'logs').glob('p2c_call_*.sh'))
    assert not form.run_button.disabled

    form.fields['fname_func'].value = str(tmp_path / 'missing.txt')
    with patch.object(bs, 'alert') as alert:
        form.run_button._internal.invoke()
        wait_for_workflow(form, shell)
    alert.assert_called_once()
    assert 'failed' in form.status.text.lower()
    assert not form.run_button.disabled


def test_window_destruction_stops_active_workflow(page):
    _, shell = page
    with shell.page_nav() as nav:
        with nav.add_page('workflow-stop', text='Stop test'):
            form = build_workflow_setter()
    process = MagicMock()
    output = MagicMock()
    form.process = process
    form.output = output
    form.run_button.destroy()
    process.terminate.assert_called_once()
    process.wait.assert_called_once_with(timeout=5)
    output.close.assert_called_once()
    assert form.process is None
