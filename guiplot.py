#!/usr/bin/env python3

import argparse
import sys

import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
from matplotlib.widgets import Slider, TextBox
from scipy.stats import zscore


def load_and_prep_data(nii_path, matrix_path, cvr_path=None, lag_path=None):
    img = nib.load(nii_path)
    func_data = img.get_fdata()
    regressor_matrix = np.loadtxt(matrix_path)

    cvr_data = nib.load(cvr_path).get_fdata() if cvr_path else None
    lag_data = nib.load(lag_path).get_fdata() if lag_path else None

    return func_data, regressor_matrix, cvr_data, lag_data


def launch_gui(
    func_data, regressor_matrix, cvr_data=None, lag_data=None, init_coords=None, fs=1.0
):
    nx, ny, nz, nt = func_data.shape
    num_lags = regressor_matrix.shape[0]

    # Calculate lag range in seconds centered around 0
    # Half of total lags mapped to negative seconds, half to positive seconds
    half_lags = (num_lags - 1) / 2.0
    min_delay_sec = -half_lags / fs
    max_delay_sec = half_lags / fs
    sec_per_step = 1.0 / fs

    # Mean anatomical volume
    mean_vol = np.mean(func_data, axis=-1)

    # Initial Coordinates
    if init_coords is not None:
        curr_x, curr_y, curr_z = init_coords
    else:
        curr_x, curr_y, curr_z = nx // 2, ny // 2, nz // 2

    # Initial shift set to middle index (0 seconds)
    curr_shift_idx = int(round(half_lags))

    # Dynamic Column Layout
    num_cols = (
        3 + (1 if cvr_data is not None else 0) + (1 if lag_data is not None else 0)
    )

    fig = plt.figure(figsize=(4 * num_cols, 8))
    gs = fig.add_gridspec(
        3, num_cols, height_ratios=[1.2, 1, 0.45], hspace=0.35, wspace=0.25
    )

    # Orthogonal Anatomical Axes
    ax_sag = fig.add_subplot(gs[0, 0])
    ax_cor = fig.add_subplot(gs[0, 1])
    ax_axi = fig.add_subplot(gs[0, 2])

    # Optional Map Axes
    col_idx = 3
    ax_cvr = None
    if cvr_data is not None:
        ax_cvr = fig.add_subplot(gs[0, col_idx])
        col_idx += 1

    ax_lag = None
    if lag_data is not None:
        ax_lag = fig.add_subplot(gs[0, col_idx])

    # Timecourse Plot Axis
    ax_ts = fig.add_subplot(gs[1, :])

    plt.subplots_adjust(bottom=0.18, left=0.06, right=0.96)

    # 1. Plot Anatomical Slices
    img_sag = ax_sag.imshow(
        mean_vol[curr_x, :, :].T, cmap='gray', origin='lower', aspect='auto'
    )
    (cross_sag_v,) = ax_sag.plot([curr_y, curr_y], [0, nz - 1], 'r-', lw=1)
    (cross_sag_h,) = ax_sag.plot([0, ny - 1], [curr_z, curr_z], 'r-', lw=1)
    ax_sag.set_title(f'Sagittal (X={curr_x})')
    ax_sag.set_xlabel('Y')
    ax_sag.set_ylabel('Z')

    img_cor = ax_cor.imshow(
        mean_vol[:, curr_y, :].T, cmap='gray', origin='lower', aspect='auto'
    )
    (cross_cor_v,) = ax_cor.plot([curr_x, curr_x], [0, nz - 1], 'r-', lw=1)
    (cross_cor_h,) = ax_cor.plot([0, nx - 1], [curr_z, curr_z], 'r-', lw=1)
    ax_cor.set_title(f'Coronal (Y={curr_y})')
    ax_cor.set_xlabel('X')
    ax_cor.set_ylabel('Z')

    img_axi = ax_axi.imshow(
        mean_vol[:, :, curr_z].T, cmap='gray', origin='lower', aspect='auto'
    )
    (cross_axi_v,) = ax_axi.plot([curr_x, curr_x], [0, ny - 1], 'r-', lw=1)
    (cross_axi_h,) = ax_axi.plot([0, nx - 1], [curr_y, curr_y], 'r-', lw=1)
    ax_axi.set_title(f'Axial (Z={curr_z})')
    ax_axi.set_xlabel('X')
    ax_axi.set_ylabel('Y')

    # 2. Plot CVR Map (Axial Slice)
    img_cvr, cross_cvr_v, cross_cvr_h = None, None, None
    if cvr_data is not None:
        img_cvr = ax_cvr.imshow(
            cvr_data[:, :, curr_z].T, cmap='inferno', origin='lower', aspect='auto'
        )
        (cross_cvr_v,) = ax_cvr.plot([curr_x, curr_x], [0, ny - 1], 'cyan', lw=1)
        (cross_cvr_h,) = ax_cvr.plot([0, nx - 1], [curr_y, curr_y], 'cyan', lw=1)
        fig.colorbar(img_cvr, ax=ax_cvr, fraction=0.046, pad=0.04)
        ax_cvr.set_title(f'CVR: {cvr_data[curr_x, curr_y, curr_z]:.2f}')
        ax_cvr.set_xlabel('X')
        ax_cvr.set_ylabel('Y')

    # 3. Plot Lag Map (Axial Slice)
    img_lag, cross_lag_v, cross_lag_h = None, None, None
    if lag_data is not None:
        img_lag = ax_lag.imshow(
            lag_data[:, :, curr_z].T, cmap='turbo', origin='lower', aspect='auto'
        )
        (cross_lag_v,) = ax_lag.plot([curr_x, curr_x], [0, ny - 1], 'black', lw=1)
        (cross_lag_h,) = ax_lag.plot([0, nx - 1], [curr_y, curr_y], 'black', lw=1)
        fig.colorbar(img_lag, ax=ax_lag, fraction=0.046, pad=0.04)
        ax_lag.set_title(f'Lag: {lag_data[curr_x, curr_y, curr_z]:.2f}s')
        ax_lag.set_xlabel('X')
        ax_lag.set_ylabel('Y')

    # 4. Timecourse Plot Setup
    init_sec = (curr_shift_idx - half_lags) / fs
    voxel_ts = zscore(func_data[curr_x, curr_y, curr_z, :])
    co2_ts = zscore(regressor_matrix[curr_shift_idx, :])

    (line_voxel,) = ax_ts.plot(
        voxel_ts, label=f'Voxel ({curr_x}, {curr_y}, {curr_z})', color='b'
    )
    (line_co2,) = ax_ts.plot(
        co2_ts, label=f'CO2 Shift ({init_sec:+.2f}s)', color='r', linestyle='--'
    )

    ax_ts.set_title(f'Timecourse Comparison — Voxel ({curr_x}, {curr_y}, {curr_z})')
    ax_ts.set_xlabel('Timepoint (TR)')
    ax_ts.set_ylabel('Z-Score')
    ax_ts.legend(loc='upper right')
    ax_ts.grid(True, linestyle=':', alpha=0.6)

    # 5. Widgets Setup (Slider in Seconds)
    ax_shift = plt.axes([0.2, 0.08, 0.65, 0.03])
    slider_shift = Slider(
        ax_shift,
        'CO2 Delay (s)',
        min_delay_sec,
        max_delay_sec,
        valinit=0.0,
        valstep=sec_per_step,
        valfmt='%+.2f s',
    )

    ax_box_x = plt.axes([0.2, 0.02, 0.1, 0.03])
    ax_box_y = plt.axes([0.45, 0.02, 0.1, 0.03])
    ax_box_z = plt.axes([0.7, 0.02, 0.1, 0.03])

    text_x = TextBox(ax_box_x, 'X: ', initial=str(curr_x))
    text_y = TextBox(ax_box_y, 'Y: ', initial=str(curr_y))
    text_z = TextBox(ax_box_z, 'Z: ', initial=str(curr_z))

    # Master Update Function
    def update_all(x, y, z, shift_sec, update_textboxes=True):
        nonlocal curr_x, curr_y, curr_z, curr_shift_idx
        curr_x, curr_y, curr_z = x, y, z

        # Convert delay in seconds back to row index in regressor_matrix
        shift_idx = int(round(shift_sec * fs + half_lags))
        curr_shift_idx = np.clip(shift_idx, 0, num_lags - 1)

        # Update Slices
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

        # Update CVR Map
        if cvr_data is not None:
            img_cvr.set_data(cvr_data[:, :, z].T)
            cross_cvr_v.set_xdata([x, x])
            cross_cvr_h.set_ydata([y, y])
            ax_cvr.set_title(f'CVR: {cvr_data[x, y, z]:.2f}')

        # Update Lag Map
        if lag_data is not None:
            img_lag.set_data(lag_data[:, :, z].T)
            cross_lag_v.set_xdata([x, x])
            cross_lag_h.set_ydata([y, y])
            ax_lag.set_title(f'Lag: {lag_data[x, y, z]:.2f}s')

        # Update Timecourse
        new_voxel_ts = zscore(func_data[x, y, z, :])
        new_co2_ts = zscore(regressor_matrix[curr_shift_idx, :])

        actual_sec = (curr_shift_idx - half_lags) / fs

        line_voxel.set_ydata(new_voxel_ts)
        line_voxel.set_label(f'Voxel ({x}, {y}, {z})')

        line_co2.set_ydata(new_co2_ts)
        line_co2.set_label(f'CO2 Shift ({actual_sec:+.2f}s)')

        ax_ts.set_title(f'Timecourse Comparison — Voxel ({x}, {y}, {z})')
        ax_ts.legend(loc='upper right')

        ax_ts.relim()
        ax_ts.autoscale_view(scalex=False, scaley=True)

        # Sync Textboxes
        if update_textboxes:
            text_x.eventson = text_y.eventson = text_z.eventson = False
            text_x.set_val(str(x))
            text_y.set_val(str(y))
            text_z.set_val(str(z))
            text_x.eventson = text_y.eventson = text_z.eventson = True

        fig.canvas.draw_idle()

    # Callbacks
    def on_click(event):
        clickable_axes = [ax_sag, ax_cor, ax_axi]
        if ax_cvr:
            clickable_axes.append(ax_cvr)
        if ax_lag:
            clickable_axes.append(ax_lag)

        if event.inaxes not in clickable_axes:
            return

        x, y, z = curr_x, curr_y, curr_z
        if event.inaxes == ax_sag:
            y, z = int(round(event.xdata)), int(round(event.ydata))
        elif event.inaxes == ax_cor:
            x, z = int(round(event.xdata)), int(round(event.ydata))
        elif event.inaxes in [ax_axi, ax_cvr, ax_lag]:
            x, y = int(round(event.xdata)), int(round(event.ydata))

        x, y, z = (
            np.clip(x, 0, nx - 1),
            np.clip(y, 0, ny - 1),
            np.clip(z, 0, nz - 1),
        )
        update_all(x, y, z, slider_shift.val)

    def on_textbox_submit(_=None):
        try:
            x = np.clip(int(text_x.text), 0, nx - 1)
            y = np.clip(int(text_y.text), 0, ny - 1)
            z = np.clip(int(text_z.text), 0, nz - 1)
            update_all(x, y, z, slider_shift.val, update_textboxes=False)
        except ValueError:
            pass

    def on_slider_change(val_sec):
        update_all(curr_x, curr_y, curr_z, val_sec)

    fig.canvas.mpl_connect('button_press_event', on_click)
    slider_shift.on_changed(on_slider_change)
    text_x.on_submit(on_textbox_submit)
    text_y.on_submit(on_textbox_submit)
    text_z.on_submit(on_textbox_submit)

    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Plot voxel timecourses against a lag-shifted regressor matrix alongside CVR and Lag maps.'
    )
    parser.add_argument('nii_file', help='Path to the 4D fMRI NIfTI file')
    parser.add_argument(
        'matrix_file', help='Path to the regressor matrix text/mat file'
    )
    parser.add_argument('--cvr', help='Optional path to 3D CVR NIfTI map')
    parser.add_argument('--lag', help='Optional path to 3D Lag NIfTI map')
    parser.add_argument(
        '--fs',
        type=float,
        default=1.0,
        help='Sampling frequency of the regressor matrix in Hz (default: 1.0)',
    )
    parser.add_argument(
        '--coords',
        '-c',
        nargs=3,
        type=int,
        metavar=('X', 'Y', 'Z'),
        help='Initial voxel coordinates',
    )

    args = parser.parse_args()

    func_data, regressor_matrix, cvr_data, lag_data = load_and_prep_data(
        args.nii_file, args.matrix_file, args.cvr, args.lag
    )

    launch_gui(
        func_data,
        regressor_matrix,
        cvr_data,
        lag_data,
        init_coords=args.coords,
        fs=args.fs,
    )
