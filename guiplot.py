#!/usr/bin/env python3

import argparse
import sys

import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
from matplotlib.widgets import Slider, TextBox
from scipy.stats import zscore


def load_and_prep_data(nii_path, matrix_path):
    # Load 4D fMRI NIfTI file
    img = nib.load(nii_path)
    func_data = img.get_fdata()

    # Load matrix file (assumes space-separated or CSV text format)
    # Each row corresponds to a lagged regressor version over time
    regressor_matrix = np.loadtxt(matrix_path)

    return func_data, regressor_matrix


def launch_gui(func_data, regressor_matrix, init_coords=None):
    nx, ny, nz, nt = func_data.shape
    max_lag = regressor_matrix.shape[0] - 1

    # Initial voxel coordinates (middle of volume) and initial shift
    # Set initial voxel coordinates (use provided or fallback to center)
    if init_coords is not None:  # <-- updated block start
        x_init, y_init, z_init = init_coords
    else:
        x_init, y_init, z_init = nx // 2, ny // 2, nz // 2  # <-- updated block end

    shift_init = 0

    # Create figure and main axis
    fig, ax = plt.subplots(figsize=(10, 5))
    plt.subplots_adjust(bottom=0.3)

    # Initial plot setup
    voxel_ts = zscore(func_data[x_init, y_init, z_init, :])
    co2_ts = zscore(regressor_matrix[shift_init, :])

    (line_voxel,) = ax.plot(
        voxel_ts, label=f'Voxel ({x_init}, {y_init}, {z_init})', color='b'
    )
    (line_co2,) = ax.plot(
        co2_ts, label=f'CO2 Shift ({shift_init})', color='r', linestyle='--'
    )

    ax.set_title(f'Timecourse Comparison — Voxel ({x_init}, {y_init}, {z_init})')
    ax.set_xlabel('Timepoint (TR)')
    ax.set_ylabel('Z-Score')
    ax.legend(loc='upper right')
    ax.grid(True, linestyle=':', alpha=0.6)

    # Add Slider for CO2 Shift
    ax_shift = plt.axes([0.2, 0.15, 0.65, 0.03])
    slider_shift = Slider(
        ax_shift,
        'CO2 Shift',
        0,
        max_lag,
        valinit=shift_init,
        valfmt='%d',
        valstep=1,
    )

    # Add Textboxes for Voxel Coordinates (X, Y, Z)
    ax_box_x = plt.axes([0.2, 0.05, 0.1, 0.04])
    ax_box_y = plt.axes([0.45, 0.05, 0.1, 0.04])
    ax_box_z = plt.axes([0.7, 0.05, 0.1, 0.04])

    text_x = TextBox(ax_box_x, 'X: ', initial=str(x_init))
    text_y = TextBox(ax_box_y, 'Y: ', initial=str(y_init))
    text_z = TextBox(ax_box_z, 'Z: ', initial=str(z_init))

    # Update callback function
    def update(_=None):
        try:
            x = int(text_x.text)
            y = int(text_y.text)
            z = int(text_z.text)
            shift = int(slider_shift.val)

            # Validate coordinate bounds
            if not (0 <= x < nx and 0 <= y < ny and 0 <= z < nz):
                return

            # Extract and z-score
            new_voxel_ts = zscore(func_data[x, y, z, :])
            new_co2_ts = zscore(regressor_matrix[shift, :])

            # Update plot lines
            line_voxel.set_ydata(new_voxel_ts)
            line_voxel.set_label(f'Voxel ({x}, {y}, {z})')

            line_co2.set_ydata(new_co2_ts)
            line_co2.set_label(f'CO2 Shift ({shift})')

            ax.set_title(f'Timecourse Comparison — Voxel ({x}, {y}, {z})')
            ax.legend(loc='upper right')

            # Rescale axis if necessary
            ax.relim()
            ax.autoscale_view(scalex=False, scaley=True)

            fig.canvas.draw_idle()
        except ValueError:
            pass  # Handle non-integer input gracefully

    # Connect callbacks
    slider_shift.on_changed(update)
    text_x.on_submit(update)
    text_y.on_submit(update)
    text_z.on_submit(update)

    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Plot voxel timecourses against a lag-shifted regressor matrix.'
    )
    parser.add_argument('nii_file', help='Path to the 4D fMRI NIfTI file')
    parser.add_argument(
        'matrix_file', help='Path to the regressor matrix text/mat file'
    )
    parser.add_argument(
        '--coords',
        '-c',
        nargs=3,
        type=int,
        metavar=('X', 'Y', 'Z'),
        help='Initial voxel coordinates (e.g., -c 32 32 15)',
    )
    args = parser.parse_args()

    func_data, regressor_matrix = load_and_prep_data(args.nii_file, args.matrix_file)
    launch_gui(func_data, regressor_matrix, init_coords=args.coords)
