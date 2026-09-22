#!/usr/bin/env python3

import argparse
import sys

import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
from matplotlib.widgets import Slider, TextBox
from scipy.stats import zscore


def load_and_prep_data(nii_path, matrix_path):
    img = nib.load(nii_path)
    func_data = img.get_fdata()
    regressor_matrix = np.loadtxt(matrix_path)
    return func_data, regressor_matrix


def launch_gui(func_data, regressor_matrix, init_coords=None):
    nx, ny, nz, nt = func_data.shape
    max_lag = regressor_matrix.shape[0] - 1

    # Mean functional volume for anatomical slice viewing
    mean_vol = np.mean(func_data, axis=-1)

    # Coordinates
    if init_coords is not None:
        curr_x, curr_y, curr_z = init_coords
    else:
        curr_x, curr_y, curr_z = nx // 2, ny // 2, nz // 2

    curr_shift = 0

    # Layout Setup: GridSpec with slices on top, timecourse on bottom
    fig = plt.figure(figsize=(12, 8))
    gs = fig.add_gridspec(3, 3, height_ratios=[1.2, 1, 0.45], hspace=0.35)

    ax_sag = fig.add_subplot(gs[0, 0])
    ax_cor = fig.add_subplot(gs[0, 1])
    ax_axi = fig.add_subplot(gs[0, 2])
    ax_ts = fig.add_subplot(gs[1, :])

    plt.subplots_adjust(bottom=0.18, left=0.08, right=0.95)

    # 1. Plot Slices
    img_sag = ax_sag.imshow(
        mean_vol[curr_x, :, :].T,
        cmap='gray',
        origin='lower',
        aspect='auto',
    )
    (cross_sag_v,) = ax_sag.plot([curr_y, curr_y], [0, nz - 1], 'r-', lw=1)
    (cross_sag_h,) = ax_sag.plot([0, ny - 1], [curr_z, curr_z], 'r-', lw=1)
    ax_sag.set_title(f'Sagittal (X={curr_x})')
    ax_sag.set_xlabel('Y')
    ax_sag.set_ylabel('Z')

    img_cor = ax_cor.imshow(
        mean_vol[:, curr_y, :].T,
        cmap='gray',
        origin='lower',
        aspect='auto',
    )
    (cross_cor_v,) = ax_cor.plot([curr_x, curr_x], [0, nz - 1], 'r-', lw=1)
    (cross_cor_h,) = ax_cor.plot([0, nx - 1], [curr_z, curr_z], 'r-', lw=1)
    ax_cor.set_title(f'Coronal (Y={curr_y})')
    ax_cor.set_xlabel('X')
    ax_cor.set_ylabel('Z')

    img_axi = ax_axi.imshow(
        mean_vol[:, :, curr_z].T,
        cmap='gray',
        origin='lower',
        aspect='auto',
    )
    (cross_axi_v,) = ax_axi.plot([curr_x, curr_x], [0, ny - 1], 'r-', lw=1)
    (cross_axi_h,) = ax_axi.plot([0, nx - 1], [curr_y, curr_y], 'r-', lw=1)
    ax_axi.set_title(f'Axial (Z={curr_z})')
    ax_axi.set_xlabel('X')
    ax_axi.set_ylabel('Y')

    # 2. Timecourse Plot Setup
    voxel_ts = zscore(func_data[curr_x, curr_y, curr_z, :])
    co2_ts = zscore(regressor_matrix[curr_shift, :])

    (line_voxel,) = ax_ts.plot(
        voxel_ts, label=f'Voxel ({curr_x}, {curr_y}, {curr_z})', color='b'
    )
    (line_co2,) = ax_ts.plot(
        co2_ts, label=f'CO2 Shift ({curr_shift})', color='r', linestyle='--'
    )

    ax_ts.set_title(f'Timecourse Comparison — Voxel ({curr_x}, {curr_y}, {curr_z})')
    ax_ts.set_xlabel('Timepoint (TR)')
    ax_ts.set_ylabel('Z-Score')
    ax_ts.legend(loc='upper right')
    ax_ts.grid(True, linestyle=':', alpha=0.6)

    # 3. Widgets Setup
    ax_shift = plt.axes([0.2, 0.08, 0.65, 0.03])
    slider_shift = Slider(
        ax_shift,
        'CO2 Shift',
        0,
        max_lag,
        valinit=curr_shift,
        valfmt='%d',
        valstep=1,
    )

    ax_box_x = plt.axes([0.2, 0.02, 0.1, 0.03])
    ax_box_y = plt.axes([0.45, 0.02, 0.1, 0.03])
    ax_box_z = plt.axes([0.7, 0.02, 0.1, 0.03])

    text_x = TextBox(ax_box_x, 'X: ', initial=str(curr_x))
    text_y = TextBox(ax_box_y, 'Y: ', initial=str(curr_y))
    text_z = TextBox(ax_box_z, 'Z: ', initial=str(curr_z))

    # Master Update Function
    def update_all(x, y, z, shift, update_textboxes=True):
        nonlocal curr_x, curr_y, curr_z, curr_shift
        curr_x, curr_y, curr_z, curr_shift = x, y, z, shift

        # Update Images
        img_sag.set_data(mean_vol[x, :, :].T)
        cross_sag_v.set_xdata([y, y])
        cross_sag_h.set_ydata([z, z])
        ax_sag.set_title(f'Sagittal (X={x})')

        img_cor.set_data(mean_vol[:, y, :].T)
        cross_cor_v.set_xdata([x, x])
        cross_cor_h.set_ydata([z, z])
        ax_cor.set_title(f'Coronal (Y={y})')

        img_axi.set_data(mean_vol[:, :, z].T)
        cross_axi_v.set_xdata([x, x])
        cross_axi_h.set_ydata([y, y])
        ax_axi.set_title(f'Axial (Z={z})')

        # Update Timecourses
        new_voxel_ts = zscore(func_data[x, y, z, :])
        new_co2_ts = zscore(regressor_matrix[shift, :])

        line_voxel.set_ydata(new_voxel_ts)
        line_voxel.set_label(f'Voxel ({x}, {y}, {z})')

        line_co2.set_ydata(new_co2_ts)
        line_co2.set_label(f'CO2 Shift ({shift})')

        ax_ts.set_title(f'Timecourse Comparison — Voxel ({x}, {y}, {z})')
        ax_ts.legend(loc='upper right')

        ax_ts.relim()
        ax_ts.autoscale_view(scalex=False, scaley=True)

        # Sync textboxes without triggering recursive callback
        if update_textboxes:
            text_x.eventson = text_y.eventson = text_z.eventson = False
            text_x.set_val(str(x))
            text_y.set_val(str(y))
            text_z.set_val(str(z))
            text_x.eventson = text_y.eventson = text_z.eventson = True

        fig.canvas.draw_idle()

    # Callbacks
    def on_click(event):
        if event.inaxes not in [ax_sag, ax_cor, ax_axi]:
            return

        x, y, z = curr_x, curr_y, curr_z
        if event.inaxes == ax_sag:
            y, z = int(round(event.xdata)), int(round(event.ydata))
        elif event.inaxes == ax_cor:
            x, z = int(round(event.xdata)), int(round(event.ydata))
        elif event.inaxes == ax_axi:
            x, y = int(round(event.xdata)), int(round(event.ydata))

        x, y, z = (
            np.clip(x, 0, nx - 1),
            np.clip(y, 0, ny - 1),
            np.clip(z, 0, nz - 1),
        )
        update_all(x, y, z, curr_shift)

    def on_textbox_submit(_=None):
        try:
            x = np.clip(int(text_x.text), 0, nx - 1)
            y = np.clip(int(text_y.text), 0, ny - 1)
            z = np.clip(int(text_z.text), 0, nz - 1)
            update_all(x, y, z, curr_shift, update_textboxes=False)
        except ValueError:
            pass

    def on_slider_change(val):
        update_all(curr_x, curr_y, curr_z, int(val))

    fig.canvas.mpl_connect('button_press_event', on_click)
    slider_shift.on_changed(on_slider_change)
    text_x.on_submit(on_textbox_submit)
    text_y.on_submit(on_textbox_submit)
    text_z.on_submit(on_textbox_submit)

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
