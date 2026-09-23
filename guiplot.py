#!/usr/bin/env python3

import argparse
import sys
import tkinter as tk
from tkinter import filedialog, ttk

import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
from darkdetect import theme
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from scipy.stats import zscore
from sv_ttk import set_theme

from phys2cvr import __version__


def emptyfunc():
    # Momentarily here for menu
    pass


def _get_parser():
    """
    Parse command line inputs for this function.

    Returns
    -------
    parser.parse_args() : argparse dict

    """
    parser = argparse.ArgumentParser(description='phys2cvr viewer.')
    parser.add_argument('nii_file', help='Path to 4D fMRI NIfTI file')
    parser.add_argument('matrix_file', help='Path to regressor matrix text/mat file')
    parser.add_argument('--cvr', default=None, help='Optional path to 3D CVR NIfTI map')
    parser.add_argument('--lag', default=None, help='Optional path to 3D Lag NIfTI map')
    parser.add_argument(
        '--mask', default=None, help='Optional path to 3D brain mask NIfTI file'
    )
    parser.add_argument(
        '--fs',
        type=float,
        default=1.0,
        help='Regressor sampling frequency in Hz (default: 1.0)',
    )

    return parser


def load_and_prep_data(
    nii_path, matrix_path, cvr_path=None, lag_path=None, mask_path=None
):
    img = nib.load(nii_path)
    func_data = img.get_fdata()

    # Get spatial pixel dimensions (pixdim)
    zooms = img.header.get_zooms()
    pixdim = zooms[:3] if len(zooms) >= 3 else (1.0, 1.0, 1.0)

    # Get TR in seconds from NIfTI header (pixdim[4]), fallback to 1.0s if invalid
    tr = zooms[3] if len(zooms) > 3 else 1.0
    if tr <= 0 or np.isnan(tr):
        tr = 1.0

    regressor_matrix = np.loadtxt(matrix_path)

    # Load mask if provided
    mask_data = (
        nib.load(mask_path).get_fdata().astype(bool)
        if mask_path and mask_path != ''
        else None
    )

    # Optional CVR and Lag Maps
    raw_cvr_data = (
        nib.load(cvr_path).get_fdata() if cvr_path and cvr_path != '' else None
    )
    raw_lag_data = (
        nib.load(lag_path).get_fdata() if lag_path and lag_path != '' else None
    )

    # Apply mask: set voxels outside the mask to NaN
    cvr_data = raw_cvr_data.copy() if raw_cvr_data is not None else None
    lag_data = raw_lag_data.copy() if raw_lag_data is not None else None

    if mask_data is not None:
        if cvr_data is not None:
            cvr_data[~mask_data] = np.nan
        if lag_data is not None:
            lag_data[~mask_data] = np.nan

    return (
        func_data,
        regressor_matrix,
        cvr_data,
        lag_data,
        tr,
        pixdim,
        raw_cvr_data,
        raw_lag_data,
    )


class VoxelViewerApp(tk.Tk):
    def __init__(
        self,
        func_data,
        regressor_matrix,
        cvr_data=None,
        lag_data=None,
        fs=1.0,
        tr=None,
        pixdim=(1.0, 1.0, 1.0),
        mask_data=None,
        raw_cvr_data=None,
        raw_lag_data=None,
    ):
        super().__init__()

        self.func_data = func_data
        self.regressor_matrix = regressor_matrix
        self.cvr_data = cvr_data
        self.lag_data = lag_data
        self.raw_cvr_data = raw_cvr_data if raw_cvr_data is not None else cvr_data
        self.raw_lag_data = raw_lag_data if raw_lag_data is not None else lag_data
        self.mask_data = mask_data
        self.fs = fs
        self.tr = tr
        self.pixdim = pixdim

        self.nx, self.ny, self.nz, self.nt = func_data.shape
        self.num_lags = regressor_matrix.shape[0]

        # Calculate time vector for full timecourses in seconds
        self.time = (
            np.arange(self.nt) * self.tr if self.tr is not None else np.arange(self.nt)
        )

        # Calculate delay parameters for the slider
        self.half_lags = (self.num_lags - 1) / 2.0
        self.min_delay_sec = -self.half_lags / self.fs
        self.max_delay_sec = self.half_lags / self.fs

        self.mean_vol = np.mean(func_data, axis=-1)

        # Default to center voxel
        self.curr_x = self.nx // 2
        self.curr_y = self.ny // 2
        self.curr_z = self.nz // 2
        self.curr_shift_sec = 0.0

        # Check if initial voxel has a valid lag value
        if self.lag_data is not None:
            init_lag = self.lag_data[self.curr_x, self.curr_y, self.curr_z]
            if not np.isnan(init_lag):
                self.curr_shift_sec = float(init_lag)

        # Configure Window
        self.title(f'fMRI Voxel & Regressor Viewer, phys2cvr v{__version__}')
        self.geometry('1200x850')

        # Apply darkdetect + sv_ttk theme
        current_theme = theme()
        self.current_theme_mode = (
            'dark' if (current_theme and current_theme.lower() == 'dark') else 'light'
        )
        set_theme(self.current_theme_mode)
        if self.current_theme_mode == 'dark':
            plt.style.use('dark_background')
        else:
            plt.style.use('seaborn-v0_8-colorblind')

        # Build UI Components
        self.create_menu()
        self.create_widgets()
        self.init_plots()
        self.update_all(
            self.curr_x,
            self.curr_y,
            self.curr_z,
            self.curr_shift_sec,
            update_controls=True,
        )

    def set_app_theme(self, theme_mode):
        self.current_theme_mode = theme_mode
        set_theme(theme_mode)

        # Toggle Matplotlib style dynamically
        style_name = (
            'dark_background' if theme_mode == 'dark' else 'seaborn-v0_8-colorblind'
        )
        plt.style.use(style_name)

        bg_color = '#000000' if theme_mode == 'dark' else '#ffffff'
        fg_color = '#ffffff' if theme_mode == 'dark' else '#000000'

        self.fig.patch.set_facecolor(bg_color)
        for ax in self.fig.axes:
            if ax not in [self.ax_cvr, self.ax_lag]:
                ax.set_facecolor(bg_color)
            ax.title.set_color(fg_color)
            ax.xaxis.label.set_color(fg_color)
            ax.yaxis.label.set_color(fg_color)
            ax.tick_params(colors=fg_color, which='both')
            for spine in ax.spines.values():
                spine.set_color(fg_color)

        self.canvas.draw_idle()

    def create_menu(self):
        # Disable menu tear-offs globally
        self.option_add('*tearOff', tk.FALSE)

        # Root menu bar
        menu_master = tk.Menu(self)
        self.config(menu=menu_master)

        # File menu
        menu_file = tk.Menu(menu_master, tearoff=0)
        menu_master.add_cascade(label='File', menu=menu_file)
        menu_file.add_command(
            label='Load timeseries', command=self.load_timeseries_file
        )
        menu_file.add_command(
            label='Load regressors', command=self.load_regressors_file
        )
        menu_file.add_command(label='Load CVR map', command=self.load_cvr_file)
        menu_file.add_command(label='Load lag map', command=self.load_lag_file)
        menu_file.add_command(label='Load mask', command=self.load_mask_file)
        menu_file.add_separator()
        menu_file.add_command(label='Exit', command=self.destroy)

        # Options menu
        menu_options = tk.Menu(menu_master, tearoff=0)
        menu_master.add_cascade(label='Options', menu=menu_options)
        menu_options.add_command(label='About', command=lambda: emptyfunc())
        menu_options.add_command(label='Documentation', command=lambda: emptyfunc())
        menu_options.add_separator()
        menu_theme = tk.Menu(menu_options, tearoff=0)
        menu_options.add_cascade(label='Theme', menu=menu_theme)
        menu_theme.add_command(
            label='Light', command=lambda: self.set_app_theme('light')
        )
        menu_theme.add_command(label='Dark', command=lambda: self.set_app_theme('dark'))

    def load_timeseries_file(self):
        file_path = filedialog.askopenfilename(
            title='Select 4D fMRI NIfTI File',
            filetypes=[('NIfTI files', '*.nii *.nii.gz'), ('All files', '*.*')],
        )
        if not file_path:
            return

        img = nib.load(file_path)
        self.func_data = img.get_fdata()

        zooms = img.header.get_zooms()
        self.pixdim = zooms[:3] if len(zooms) >= 3 else (1.0, 1.0, 1.0)

        tr = zooms[3] if len(zooms) > 3 else 1.0
        if tr <= 0 or np.isnan(tr):
            tr = 1.0
        self.tr = tr

        self.nx, self.ny, self.nz, self.nt = self.func_data.shape
        self.time = (
            np.arange(self.nt) * self.tr if self.tr is not None else np.arange(self.nt)
        )
        self.mean_vol = np.mean(self.func_data, axis=-1)

        self.curr_x = np.clip(self.curr_x, 0, self.nx - 1)
        self.curr_y = np.clip(self.curr_y, 0, self.ny - 1)
        self.curr_z = np.clip(self.curr_z, 0, self.nz - 1)

        self.reinit_plots()

    def load_regressors_file(self):
        file_path = filedialog.askopenfilename(
            title='Select Regressor Matrix File',
            filetypes=[('Text/Mat files', '*.txt *.mat *.1D'), ('All files', '*.*')],
        )
        if not file_path:
            return

        self.regressor_matrix = np.loadtxt(file_path)
        self.num_lags = self.regressor_matrix.shape[0]

        self.half_lags = (self.num_lags - 1) / 2.0
        self.min_delay_sec = -self.half_lags / self.fs
        self.max_delay_sec = self.half_lags / self.fs

        self.slider.config(from_=self.min_delay_sec, to=self.max_delay_sec)
        self.update_all(
            self.curr_x,
            self.curr_y,
            self.curr_z,
            self.curr_shift_sec,
            update_controls=True,
        )

    def load_cvr_file(self):
        file_path = filedialog.askopenfilename(
            title='Select 3D CVR Map NIfTI File',
            filetypes=[('NIfTI files', '*.nii *.nii.gz'), ('All files', '*.*')],
        )
        if not file_path:
            return

        self.raw_cvr_data = nib.load(file_path).get_fdata()
        cvr_data = self.raw_cvr_data.copy()
        if self.mask_data is not None:
            cvr_data[~self.mask_data] = np.nan
        self.cvr_data = cvr_data

        self.reinit_plots()

    def load_lag_file(self):
        file_path = filedialog.askopenfilename(
            title='Select 3D Lag Map NIfTI File',
            filetypes=[('NIfTI files', '*.nii *.nii.gz'), ('All files', '*.*')],
        )
        if not file_path:
            return

        self.raw_lag_data = nib.load(file_path).get_fdata()
        lag_data = self.raw_lag_data.copy()
        if self.mask_data is not None:
            lag_data[~self.mask_data] = np.nan
        self.lag_data = lag_data

        self.reinit_plots()

    def load_mask_file(self):
        file_path = filedialog.askopenfilename(
            title='Select 3D Brain Mask NIfTI File',
            filetypes=[('NIfTI files', '*.nii *.nii.gz'), ('All files', '*.*')],
        )
        if not file_path:
            return

        self.mask_data = nib.load(file_path).get_fdata().astype(bool)

        # Apply mask on raw data copies to avoid compounding mask operations
        if self.raw_cvr_data is not None:
            cvr_data = self.raw_cvr_data.copy()
            cvr_data[~self.mask_data] = np.nan
            self.cvr_data = cvr_data

        if self.raw_lag_data is not None:
            lag_data = self.raw_lag_data.copy()
            lag_data[~self.mask_data] = np.nan
            self.lag_data = lag_data

        self.reinit_plots()

    def reinit_plots(self):
        if hasattr(self, 'canvas'):
            self.canvas.get_tk_widget().destroy()
        plt.close(self.fig)

        self.init_plots()
        self.set_app_theme(self.current_theme_mode)
        self.update_all(
            self.curr_x,
            self.curr_y,
            self.curr_z,
            self.curr_shift_sec,
            update_controls=True,
        )

    def create_widgets(self):
        # Main Layout Frames
        self.plot_frame = ttk.Frame(self)
        self.plot_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True, pady=5)

        self.control_frame = ttk.Frame(self, padding=10)
        self.control_frame.pack(side=tk.BOTTOM, fill=tk.X)

        # 1. Coordinate Inputs
        coord_frame = ttk.LabelFrame(
            self.control_frame, text=' Voxel Coordinates ', padding=5
        )
        coord_frame.pack(side=tk.LEFT, padx=10)

        ttk.Label(coord_frame, text='X:').grid(row=0, column=0, padx=2)
        self.entry_x = ttk.Entry(coord_frame, width=5)
        self.entry_x.grid(row=0, column=1, padx=4)

        ttk.Label(coord_frame, text='Y:').grid(row=0, column=2, padx=2)
        self.entry_y = ttk.Entry(coord_frame, width=5)
        self.entry_y.grid(row=0, column=3, padx=4)

        ttk.Label(coord_frame, text='Z:').grid(row=0, column=4, padx=2)
        self.entry_z = ttk.Entry(coord_frame, width=5)
        self.entry_z.grid(row=0, column=5, padx=4)

        btn_update = ttk.Button(
            coord_frame, text='Set', command=self.on_coord_entry_submit
        )
        btn_update.grid(row=0, column=6, padx=6)

        # Bind Enter key to coordinate fields
        self.entry_x.bind('<Return>', lambda e: self.on_coord_entry_submit())
        self.entry_y.bind('<Return>', lambda e: self.on_coord_entry_submit())
        self.entry_z.bind('<Return>', lambda e: self.on_coord_entry_submit())

        # 2. Regressor Shift Slider
        slider_frame = ttk.LabelFrame(
            self.control_frame, text=' CO2 Delay (seconds) ', padding=5
        )
        slider_frame.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=10)

        self.slider_var = tk.DoubleVar(value=0.0)
        self.slider = ttk.Scale(
            slider_frame,
            from_=self.min_delay_sec,
            to=self.max_delay_sec,
            variable=self.slider_var,
            orient=tk.HORIZONTAL,
            command=self.on_slider_move,
        )
        self.slider.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=5)

        self.lbl_delay = ttk.Label(
            slider_frame, text='+0.00 s', width=10, anchor='center'
        )
        self.lbl_delay.pack(side=tk.RIGHT, padx=5)

    def init_plots(self):
        num_cols = (
            3
            + (1 if self.cvr_data is not None else 0)
            + (1 if self.lag_data is not None else 0)
        )

        self.fig = plt.figure(figsize=(4 * num_cols, 7), tight_layout=True)
        gs = self.fig.add_gridspec(
            2, num_cols, height_ratios=[1.2, 1], hspace=0.3, wspace=0.25
        )

        # Calculate voxel aspect ratios from pixdim
        dx, dy, dz = self.pixdim
        aspect_sag = dz / dy
        aspect_cor = dz / dx
        aspect_axi = dy / dx

        # Axes setup
        self.ax_sag = self.fig.add_subplot(gs[0, 0])
        self.ax_cor = self.fig.add_subplot(gs[0, 1])
        self.ax_axi = self.fig.add_subplot(gs[0, 2])

        col_idx = 3
        self.ax_cvr = None
        if self.cvr_data is not None:
            self.ax_cvr = self.fig.add_subplot(gs[0, col_idx])
            col_idx += 1

        self.ax_lag = None
        if self.lag_data is not None:
            self.ax_lag = self.fig.add_subplot(gs[0, col_idx])

        self.ax_ts = self.fig.add_subplot(gs[1, :])

        # Render Base Anatomical Slices with pixdim scaling
        self.img_sag = self.ax_sag.imshow(
            self.mean_vol[self.curr_x, :, :].T,
            cmap='gray',
            origin='lower',
            aspect=aspect_sag,
        )
        (self.cross_sag_v,) = self.ax_sag.plot(
            [self.curr_y, self.curr_y], [0, self.nz - 1], 'r-', lw=1
        )
        (self.cross_sag_h,) = self.ax_sag.plot(
            [0, self.ny - 1], [self.curr_z, self.curr_z], 'r-', lw=1
        )
        self.ax_sag.set_xlabel('Y')
        self.ax_sag.set_ylabel('Z')

        self.img_cor = self.ax_cor.imshow(
            self.mean_vol[:, self.curr_y, :].T,
            cmap='gray',
            origin='lower',
            aspect=aspect_cor,
        )
        (self.cross_cor_v,) = self.ax_cor.plot(
            [self.curr_x, self.curr_x], [0, self.nz - 1], 'r-', lw=1
        )
        (self.cross_cor_h,) = self.ax_cor.plot(
            [0, self.nx - 1], [self.curr_z, self.curr_z], 'r-', lw=1
        )
        self.ax_cor.set_xlabel('X')
        self.ax_cor.set_ylabel('Z')

        self.img_axi = self.ax_axi.imshow(
            self.mean_vol[:, :, self.curr_z].T,
            cmap='gray',
            origin='lower',
            aspect=aspect_axi,
        )
        (self.cross_axi_v,) = self.ax_axi.plot(
            [self.curr_x, self.curr_x], [0, self.ny - 1], 'r-', lw=1
        )
        (self.cross_axi_h,) = self.ax_axi.plot(
            [0, self.nx - 1], [self.curr_y, self.curr_y], 'r-', lw=1
        )
        self.ax_axi.set_xlabel('X')
        self.ax_axi.set_ylabel('Y')

        # CVR Map (Black background for masked-out NaNs)
        if self.cvr_data is not None:
            self.ax_cvr.set_facecolor('black')
            cmap_cvr = plt.cm.inferno.copy()
            cmap_cvr.set_bad(color='black')

            self.img_cvr = self.ax_cvr.imshow(
                self.cvr_data[:, :, self.curr_z].T,
                cmap=cmap_cvr,
                origin='lower',
                aspect=aspect_axi,
            )
            (self.cross_cvr_v,) = self.ax_cvr.plot(
                [self.curr_x, self.curr_x], [0, self.ny - 1], 'cyan', lw=1
            )
            (self.cross_cvr_h,) = self.ax_cvr.plot(
                [0, self.nx - 1], [self.curr_y, self.curr_y], 'cyan', lw=1
            )
            self.fig.colorbar(self.img_cvr, ax=self.ax_cvr, fraction=0.046, pad=0.04)
            self.ax_cvr.set_xlabel('X')
            self.ax_cvr.set_ylabel('Y')

        # Lag Map (Inverted Viridis + Black background for masked-out NaNs)
        if self.lag_data is not None:
            self.ax_lag.set_facecolor('black')
            cmap_lag = plt.cm.viridis_r.copy()  # Inverted viridis
            cmap_lag.set_bad(color='black')

            self.img_lag = self.ax_lag.imshow(
                self.lag_data[:, :, self.curr_z].T,
                cmap=cmap_lag,
                origin='lower',
                aspect=aspect_axi,
            )
            (self.cross_lag_v,) = self.ax_lag.plot(
                [self.curr_x, self.curr_x], [0, self.ny - 1], 'red', lw=1
            )
            (self.cross_lag_h,) = self.ax_lag.plot(
                [0, self.nx - 1], [self.curr_y, self.curr_y], 'red', lw=1
            )
            self.fig.colorbar(self.img_lag, ax=self.ax_lag, fraction=0.046, pad=0.04)
            self.ax_lag.set_xlabel('X')
            self.ax_lag.set_ylabel('Y')

        # Full Timecourse setup (X-axis in seconds)
        (self.line_voxel,) = self.ax_ts.plot([], [], label='Voxel TS', color='#1f77b4')
        (self.line_co2,) = self.ax_ts.plot(
            [], [], label='CO2 Regressor', color='#ff7f0e', linestyle='--'
        )
        self.ax_ts.set_xlabel('Time (s)')
        self.ax_ts.set_ylabel('Z-Score')
        self.ax_ts.set_xlim(0, self.time[-1])
        self.ax_ts.grid(True, linestyle=':', alpha=0.5)

        # Canvas Embedding
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.plot_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # Canvas Mouse Click Binding
        self.canvas.mpl_connect('button_press_event', self.on_canvas_click)

    def update_all(self, x, y, z, shift_sec, update_controls=False):
        self.curr_x, self.curr_y, self.curr_z = x, y, z
        self.curr_shift_sec = shift_sec

        # Compute regressor matrix row index
        shift_idx = int(round(shift_sec * self.fs + self.half_lags))
        shift_idx = np.clip(shift_idx, 0, self.num_lags - 1)
        actual_sec = (shift_idx - self.half_lags) / self.fs

        # Update Images & Crosshairs
        self.img_sag.set_data(self.mean_vol[x, :, :].T)
        self.cross_sag_v.set_xdata([y, y])
        self.cross_sag_h.set_ydata([z, z])
        self.ax_sag.set_title(f'Sagittal (X={x})')

        self.img_cor.set_data(self.mean_vol[:, y, :].T)
        self.cross_cor_v.set_xdata([x, x])
        self.cross_cor_h.set_ydata([z, z])
        self.ax_cor.set_title(f'Coronal (Y={y})')

        self.img_axi.set_data(self.mean_vol[:, :, z].T)
        self.cross_axi_v.set_xdata([x, x])
        self.cross_axi_h.set_ydata([y, y])
        self.ax_axi.set_title(f'Axial (Z={z})')

        if self.cvr_data is not None:
            self.img_cvr.set_data(self.cvr_data[:, :, z].T)
            self.cross_cvr_v.set_xdata([x, x])
            self.cross_cvr_h.set_ydata([y, y])
            val = self.cvr_data[x, y, z]
            self.ax_cvr.set_title(
                f'CVR: {val:.2f}' if not np.isnan(val) else 'CVR: Masked'
            )

        if self.lag_data is not None:
            self.img_lag.set_data(self.lag_data[:, :, z].T)
            self.cross_lag_v.set_xdata([x, x])
            self.cross_lag_h.set_ydata([y, y])
            val = self.lag_data[x, y, z]
            self.ax_lag.set_title(
                f'Lag: {val:.2f}s' if not np.isnan(val) else 'Lag: Masked'
            )

        # Plot full voxel and regressor timecourses against time in seconds
        full_voxel_ts = zscore(self.func_data[x, y, z, :])
        full_co2_ts = zscore(self.regressor_matrix[shift_idx, :])

        self.line_voxel.set_data(self.time, full_voxel_ts)
        self.line_voxel.set_label(f'Voxel ({x}, {y}, {z})')

        self.line_co2.set_data(self.time, full_co2_ts)
        self.line_co2.set_label(f'CO2 Shift ({actual_sec:+.2f}s)')

        self.ax_ts.set_title(f'Timecourse Comparison — Voxel ({x}, {y}, {z})')
        self.ax_ts.legend(loc='upper right')
        self.ax_ts.set_xlim(0, self.time[-1])
        self.ax_ts.relim()
        self.ax_ts.autoscale_view(scalex=False, scaley=True)

        # Update Tkinter Controls
        if update_controls:
            self.entry_x.delete(0, tk.END)
            self.entry_x.insert(0, str(x))
            self.entry_y.delete(0, tk.END)
            self.entry_y.insert(0, str(y))
            self.entry_z.delete(0, tk.END)
            self.entry_z.insert(0, str(z))

            self.slider_var.set(actual_sec)

        self.lbl_delay.config(text=f'{actual_sec:+.2f} s')
        self.canvas.draw_idle()

    # Callbacks
    def on_canvas_click(self, event):
        clickable = [self.ax_sag, self.ax_cor, self.ax_axi]
        if self.ax_cvr:
            clickable.append(self.ax_cvr)
        if self.ax_lag:
            clickable.append(self.ax_lag)

        if event.inaxes not in clickable or event.xdata is None or event.ydata is None:
            return

        x, y, z = self.curr_x, self.curr_y, self.curr_z

        if event.inaxes == self.ax_sag:
            y, z = int(round(event.xdata)), int(round(event.ydata))
        elif event.inaxes == self.ax_cor:
            x, z = int(round(event.xdata)), int(round(event.ydata))
        elif event.inaxes in [self.ax_axi, self.ax_cvr, self.ax_lag]:
            x, y = int(round(event.xdata)), int(round(event.ydata))

        x = np.clip(x, 0, self.nx - 1)
        y = np.clip(y, 0, self.ny - 1)
        z = np.clip(z, 0, self.nz - 1)

        shift_sec = self.curr_shift_sec
        if self.lag_data is not None:
            lag_val = self.lag_data[x, y, z]
            if not np.isnan(lag_val):
                shift_sec = float(lag_val)

        self.update_all(x, y, z, shift_sec, update_controls=True)

    def on_coord_entry_submit(self):
        try:
            x = np.clip(int(self.entry_x.get()), 0, self.nx - 1)
            y = np.clip(int(self.entry_y.get()), 0, self.ny - 1)
            z = np.clip(int(self.entry_z.get()), 0, self.nz - 1)

            shift_sec = self.curr_shift_sec
            if self.lag_data is not None:
                lag_val = self.lag_data[x, y, z]
                if not np.isnan(lag_val):
                    shift_sec = float(lag_val)

            self.update_all(x, y, z, shift_sec, update_controls=True)
        except ValueError:
            pass

    def on_slider_move(self, val):
        self.update_all(
            self.curr_x,
            self.curr_y,
            self.curr_z,
            float(val),
            update_controls=False,
        )


def _main(argv=None):
    args = _get_parser().parse_args(argv)

    (
        func_data,
        regressor_matrix,
        cvr_data,
        lag_data,
        tr,
        pixdim,
        raw_cvr_data,
        raw_lag_data,
    ) = load_and_prep_data(
        args.nii_file, args.matrix_file, args.cvr, args.lag, args.mask
    )

    mask_data = (
        nib.load(args.mask).get_fdata().astype(bool)
        if args.mask and args.mask != ''
        else None
    )

    app = VoxelViewerApp(
        func_data=func_data,
        regressor_matrix=regressor_matrix,
        cvr_data=cvr_data,
        lag_data=lag_data,
        mask_data=mask_data,
        raw_cvr_data=raw_cvr_data,
        raw_lag_data=raw_lag_data,
        fs=args.fs,
        tr=tr,
        pixdim=pixdim,
    )
    app.mainloop()


if __name__ == '__main__':
    _main(sys.argv[1:])
