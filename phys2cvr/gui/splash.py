"""Startup splash screen."""

import tkinter as tk
from pathlib import Path

from PIL import Image, ImageTk

CITATION = 'Please cite https://doi.org/10.5281/zenodo.5559756'


class SplashScreen:
    """Display the project logo while the application is being initialized."""

    def __init__(self, parent):
        self.window = tk.Toplevel(parent)
        self.window.overrideredirect(True)
        self.window.configure(background='white')
        with Image.open(Path(__file__).parent / 'assets' / 'logo.png') as logo:
            logo.thumbnail((440, 180))
            self.image = ImageTk.PhotoImage(logo, master=self.window)
        tk.Label(self.window, image=self.image, background='white').pack(
            padx=20, pady=8
        )
        tk.Label(self.window, text='Loading phys2cvr...', background='white').pack()
        tk.Label(self.window, text=CITATION, background='white').pack(padx=20, pady=12)
        self.window.update_idletasks()
        width, height = self.window.winfo_reqwidth(), self.window.winfo_reqheight()
        x = (self.window.winfo_screenwidth() - width) // 2
        y = (self.window.winfo_screenheight() - height) // 2
        self.window.geometry(f'{width}x{height}+{x}+{y}')
        self.window.update()

    def close(self):
        """Close the splash."""
        self.window.destroy()
