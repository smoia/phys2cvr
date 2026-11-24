"""Configure and run a phys2cvr workflow."""

import argparse
import logging
import math
import subprocess
import sys
import tempfile

import bootstack as bs

from phys2cvr.cli.run import _check_opt_conf, _get_parser

LGR = logging.getLogger(__name__)

_GROUPED = {'response_function', 'r2model', 'workflow_config'}
_DIRECTORIES = {'outdir', 'regr_dir'}
_MULTI_FILES = {
    'denoise_matrix_file',
    'orthogonalised_matrix_file',
    'extra_matrix_file',
}


def _option(action):
    return next(
        (name for name in action.option_strings if name.startswith('--')),
        action.option_strings[0],
    )


def _label(action):
    return _option(action).lstrip('-').replace('-', ' ').replace('_', ' ').capitalize()


class WorkflowSetter:
    """Build a form from the CLI parser and run the workflow in a child process."""

    def __init__(self):
        self.parser = _get_parser()
        self.fields = {}
        self.actions = {}
        self.groups = []
        self.process = None
        self.output = None
        self._tooltips = []

    def build(self):
        """Create controls for all workflow arguments."""
        with bs.Row(horizontal='stretch', gap=12, padding=8):
            self.run_button = bs.Button('Run phys2cvr', on_click=self.run)
            self.status = bs.Label('Ready')

        with bs.ScrollView(grow=True, horizontal='stretch'):
            with bs.Column(padding=12, gap=8, horizontal_items='stretch'):
                bs.Label(
                    'Configure the workflow. Hover over a field for CLI help.',
                    font='body',
                )
                for group in self.parser._action_groups:
                    actions = [
                        action
                        for action in group._group_actions
                        if not isinstance(
                            action, (argparse._HelpAction, argparse._VersionAction)
                        )
                    ]
                    if not actions:
                        continue
                    with bs.GroupBox(
                        group.title, padding=8, horizontal='stretch'
                    ) as group_box:
                        self.groups.append(group_box)
                        with bs.Grid(
                            columns=2, gap=8, horizontal='stretch',
                            horizontal_items='stretch', vertical_items='top',
                        ):
                            for action in actions:
                                if action.dest in self.actions:
                                    continue
                                self.actions[action.dest] = action
                                with bs.Column(
                                    gap=4, horizontal_items='stretch',
                                    horizontal='stretch',
                                ):
                                    self._build_field(action, actions)
        self.run_button.on_destroy(lambda event: self._stop())
        return self

    def _build_field(self, action, group_actions):
        dest = action.dest
        if dest in _GROUPED:
            choices = [item for item in group_actions if item.dest == dest]
            options = [('Default', '')]
            for item in choices:
                text = (
                    'Custom file'
                    if item.const is None and item.type is str
                    else (
                        'No response function'
                        if item.const is None
                        else str(item.const)
                    )
                )
                options.append((text, _option(item)))
            field = bs.Select(
                options, value='', label=dest.replace('_', ' ').capitalize()
            )
            if dest == 'response_function':
                self.response_path = bs.PathField(label='Custom response function file')
            help_text = '\n\n'.join(f'{_option(item)}: {item.help}' for item in choices)
        elif isinstance(
            action, (argparse._StoreTrueAction, argparse._StoreFalseAction)
        ):
            field = bs.Switch(_label(action), value=False)
            help_text = action.help
        elif action.nargs == '*':
            if dest not in _MULTI_FILES:
                raise ValueError(f'Unsupported multiple-value option: {dest}')
            field = bs.TextArea(label=f'{_label(action)} (one path per line)', height=3)
            bs.Button('Add files', on_click=lambda field=field: self._add_files(field))
            help_text = action.help
        elif action.type in (int, float):
            field = bs.NumberField(
                value=action.default if action.default is not None else 0,
                label=_label(action),
                step=1 if action.type is int else 0.1,
            )
            if action.default is None:
                field.clear()
            help_text = action.help
        elif action.type is str and action.nargs is None:
            field = bs.PathField(
                label=_label(action),
                mode='directory' if dest in _DIRECTORIES else 'open',
                required=action.required,
            )
            help_text = action.help
        else:
            raise ValueError(f'Unsupported parser action: {action}')
        self.fields[dest] = field
        if help_text:
            self._tooltips.append(
                bs.Tooltip(
                    field, help_text.replace('%(prog)s', 'phys2cvr'), wrap_width=450
                )
            )

    @staticmethod
    def _add_files(field):
        paths = bs.ask_open_files(title='Select matrix files')
        if paths:
            existing = field.value.splitlines()
            field.value = '\n'.join(existing + list(paths))

    def arguments(self):
        """Return CLI arguments after validating the form and parser constraints."""
        argv = []
        for dest, action in self.actions.items():
            field = self.fields[dest]
            value = field.value
            if dest in _GROUPED:
                if value:
                    argv.append(value)
                    if value == '--response-function-path':
                        path = self.response_path.text.strip()
                        if not path:
                            raise ValueError('Select a custom response function file.')
                        argv.append(path)
            elif isinstance(
                action, (argparse._StoreTrueAction, argparse._StoreFalseAction)
            ):
                if value:
                    argv.append(_option(action))
            elif dest in _MULTI_FILES:
                paths = [line.strip() for line in value.splitlines() if line.strip()]
                if paths:
                    argv.extend([_option(action), *paths])
            elif action.type in (int, float):
                text = field.text.strip()
                if text:
                    try:
                        number = action.type(text)
                    except ValueError:
                        expected = 'an integer' if action.type is int else 'a number'
                        raise ValueError(
                            f'{_label(action)} must be {expected}.'
                        ) from None
                    if action.type is float and not math.isfinite(number):
                        raise ValueError(f'{_label(action)} must be finite.')
                    argv.extend([_option(action), str(number)])
            elif field.text.strip():
                argv.extend([_option(action), field.text.strip()])
            elif action.required:
                raise ValueError(f'{_label(action)} is required.')

        def parser_error(message):
            raise ValueError(message)

        parser = _get_parser()
        parser.error = parser_error
        _check_opt_conf(parser.parse_args(argv))
        return argv

    def run(self):
        """Validate inputs and start the workflow without blocking the GUI."""
        if self.process is not None:
            return
        try:
            argv = self.arguments()
            output = tempfile.TemporaryFile(mode='w+t', encoding='utf-8')
            try:
                process = subprocess.Popen(
                    [sys.executable, '-m', 'phys2cvr.workflows', *argv],
                    stdout=output,
                    stderr=subprocess.STDOUT,
                )
            except OSError:
                output.close()
                raise
        except (ValueError, OSError) as err:
            LGR.warning('Could not start workflow: %s', err)
            bs.alert(str(err), title='Could not start workflow', severity='warning')
            return
        self.output = output
        self.process = process
        self.run_button.disabled = True
        self.status.text = 'Running phys2cvr…'
        self.run_button.schedule.delay(250, self._poll)

    def _poll(self):
        code = self.process.poll()
        if code is None:
            self.run_button.schedule.delay(250, self._poll)
            return
        self.output.seek(0)
        output = self.output.read()
        self.output.close()
        self.output = None
        self.process = None
        self.run_button.disabled = False
        if code == 0:
            self.status.text = 'Workflow completed.'
        else:
            self.status.text = f'Workflow failed (exit code {code}).'
            LGR.error('Workflow failed (exit code %s):\n%s', code, output)
            bs.alert(
                f'Workflow failed (exit code {code}).\n\n{output[-4000:]}',
                title='Workflow failed',
                severity='danger',
            )

    def _stop(self):
        if self.process is not None:
            self.process.terminate()
            try:
                self.process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                self.process.kill()
                self.process.wait()
            self.process = None
            self.output.close()
            self.output = None


def build_workflow_setter():
    """Build the workflow configuration page."""
    return WorkflowSetter().build()
