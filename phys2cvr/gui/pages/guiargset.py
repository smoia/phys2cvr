"""Configure and run a phys2cvr workflow."""

import argparse
import logging
import math
import subprocess
import sys
import tempfile

import bootstack as bs
import nibabel as nib
import numpy as np
from matplotlib.figure import Figure

from phys2cvr.cli.run import _check_opt_conf, _get_parser
from phys2cvr.gui.utils import clip_position, orthogonal_slices, slice_position
from phys2cvr.io import EXT_ARRAY, load_array, load_physio

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

        with bs.Row(
            grow=True,
            horizontal='stretch',
            grow_items=True,
            vertical_items='stretch',
            gap=8,
        ):
            with bs.ScrollView(
                horizontal='stretch', scroll_direction='both'
            ) as self.form_scroll:
                self._build_form()
            self.preview = WorkflowPreview()
            self.preview.build()
        for dest in ('fname_func', 'fname_mask', 'fname_roi', 'fname_co2'):
            self.fields[dest].on_change(
                lambda event, dest=dest: self.preview.load(
                    dest, self.fields[dest].text.strip()
                )
            )
        self.fields['freq'].on_change(
            lambda event: self.preview.set_frequency(self.fields['freq'].value)
        )
        self.run_button.on_destroy(lambda event: self._stop())
        return self

    def _build_form(self):
        """Build the scrollable argument controls."""
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
                    group.title,
                    padding=8,
                    horizontal='stretch',
                    accent='warning' if any(action.required for action in actions) else None,
                ) as group_box:
                    self.groups.append(group_box)
                    with bs.Grid(
                        columns=2,
                        gap=8,
                        horizontal='stretch',
                        horizontal_items='stretch',
                        vertical_items='top',
                    ):
                        for action in actions:
                            if action.dest in self.actions:
                                continue
                            self.actions[action.dest] = action
                            with bs.Column(
                                gap=4,
                                horizontal_items='stretch',
                                horizontal='stretch',
                            ):
                                self._build_field(action, actions)

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


class WorkflowPreview:
    """Preview workflow inputs without changing the configured arguments."""

    _titles = {
        'fname_func': 'Mean functional volume',
        'fname_mask': 'Mask',
        'fname_roi': 'ROI mask',
        'fname_co2': 'Regressor',
    }

    def __init__(self):
        self.data = dict.fromkeys(self._titles)
        self.positions = {}
        self.pixdims = {}
        self.frequency = None
        self.file_frequency = None
        self.charts = {}
        self.axes = {}
        self.errors = {}

    def build(self):
        """Create empty plots in the order used by the workflow."""
        with bs.Column(
            horizontal='stretch',
            horizontal_items='stretch',
            weights=[3, 2, 3, 3],
        ) as self.container:
            for key in ('fname_func', 'fname_co2', 'fname_mask', 'fname_roi'):
                fig = Figure(figsize=(1, 1), layout='tight')
                self.charts[key] = bs.Chart(figure=fig)
                self.axes[key] = (
                    [fig.add_subplot(111)]
                    if key == 'fname_co2'
                    else [fig.add_subplot(1, 3, i + 1) for i in range(3)]
                )
                if key != 'fname_co2':
                    fig.canvas.mpl_connect(
                        'button_press_event',
                        lambda event, key=key: self._navigate(key, event),
                    )
                    fig.canvas.mpl_connect(
                        'scroll_event',
                        lambda event, key=key: self._navigate(key, event),
                    )
                self._draw(key)

    def load(self, key, path):
        """Load an input preview; report errors without blocking form editing."""
        self.data[key] = None
        self.errors.pop(key, None)
        if key == 'fname_co2':
            self.file_frequency = None
        if path:
            try:
                if key == 'fname_co2':
                    if path.lower().endswith('.phys'):
                        values, _, self.file_frequency = load_physio(path)
                        if self.file_frequency is not None and (
                            not np.isfinite(self.file_frequency)
                            or self.file_frequency <= 0
                        ):
                            self.file_frequency = None
                    else:
                        values = load_array(path)
                    values = np.asarray(values).squeeze()
                    if values.ndim != 1 or not values.size:
                        raise ValueError('Regressor must be a non-empty vector.')
                    self.data[key] = values
                elif key == 'fname_func' and any(
                    path.lower().endswith(ext.lower()) for ext in EXT_ARRAY
                ):
                    self.errors[key] = 'No volume preview for a 1D input.'
                else:
                    image = nib.load(path)
                    values = np.asarray(image.get_fdata())
                    expected = 4 if key == 'fname_func' else 3
                    if values.ndim != expected or any(n == 0 for n in values.shape):
                        raise ValueError(f'Expected a non-empty {expected}D image.')
                    self.data[key] = (
                        np.mean(values, axis=-1)
                        if key == 'fname_func'
                        else np.isfinite(values) & (values != 0)
                    )
                    self.positions[key] = [n // 2 for n in values.shape[:3]]
                    self.pixdims[key] = tuple(
                        float(v) if np.isfinite(v) and v > 0 else 1.0
                        for v in image.header.get_zooms()[:3]
                    )
            except (
                OSError,
                ValueError,
                EOFError,
                ImportError,
                NotImplementedError,
                nib.filebasedimages.ImageFileError,
            ) as err:
                LGR.warning('Could not preview %s: %s', path, err)
                self.errors[key] = str(err)
        self._draw(key)

    def set_frequency(self, frequency):
        """Set the sampling frequency for the regressor time axis."""
        self.frequency = (
            frequency
            if frequency and np.isfinite(frequency) and frequency > 0
            else None
        )
        self._draw('fname_co2')

    def _draw(self, key):
        axes = self.axes[key]
        data = self.data[key]
        for axis in axes:
            axis.clear()
        if key == 'fname_co2':
            axis = axes[0]
            frequency = self.frequency or self.file_frequency
            if data is not None:
                time = np.arange(data.size)
                if frequency:
                    time = time / frequency
                axis.plot(time, data)
            axis.set_xlabel('Time (s)' if frequency else 'Sample')
            axis.set_title(self._titles[key])
        else:
            if data is not None:
                x, y, z = self.positions[key]
                dx, dy, dz = self.pixdims[key]
                for axis, image, aspect, point in zip(
                    axes,
                    orthogonal_slices(data, (x, y, z)),
                    (dz / dy, dz / dx, dy / dx),
                    ((y, z), (x, z), (x, y)),
                    strict=True,
                ):
                    axis.imshow(image.T, origin='lower', cmap='gray', aspect=aspect)
                    axis.axvline(point[0], color='red', lw=0.7)
                    axis.axhline(point[1], color='red', lw=0.7)
            for axis, orientation in zip(
                axes, ('Sagittal', 'Coronal', 'Axial'), strict=True
            ):
                axis.set_title(f'{self._titles[key]}: {orientation}')
                axis.set_xticks([])
                axis.set_yticks([])
        if key in self.errors:
            axes[0].text(
                0.5,
                0.5,
                self.errors[key],
                transform=axes[0].transAxes,
                ha='center',
                va='center',
                wrap=True,
            )
        self.charts[key].draw()

    def _navigate(self, key, event):
        data = self.data[key]
        if data is None or event.inaxes not in self.axes[key]:
            return
        orientation = self.axes[key].index(event.inaxes)
        position = self.positions[key].copy()
        if event.name == 'scroll_event':
            position[orientation] += int(np.sign(event.step))
        elif (
            event.button == 1
            and event.xdata is not None
            and event.ydata is not None
            and np.isfinite(event.xdata)
            and np.isfinite(event.ydata)
        ):
            position = slice_position(
                position, orientation, round(event.xdata), round(event.ydata)
            )
        else:
            return
        self.positions[key] = clip_position(position, data.shape)
        self._draw(key)


def build_workflow_setter():
    """Build the workflow configuration page."""
    return WorkflowSetter().build()
