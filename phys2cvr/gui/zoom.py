"""Live application font and chart zoom."""

from tkinter import font as tkfont
from weakref import WeakKeyDictionary

from bootstack.style.typography import _DERIVED_FONTS, TK_FONT_OVERRIDES, Typography
from matplotlib.text import Text


def ensure_fonts(root):
    """Initialize Bootstack's named fonts when reopening an app in one process."""
    if 'body' not in root.tk.splitlist(root.tk.call('font', 'names')):
        Typography.initialize(root)


class Zoom:
    """Scale GUI fonts and plot text relative to their original sizes."""

    def __init__(self, shell, charts, logviewer):
        self.shell = shell
        self.charts = tuple(charts)
        self.logviewer = logviewer
        self.level = 0
        self.base_sizes = {
            name: Typography.get_token(name).size for name in Typography.token_names()
        }
        self.text_sizes = WeakKeyDictionary()
        for chart in self.charts:
            original_draw = chart.draw

            def draw(chart=chart, original_draw=original_draw):
                self._scale_figure(chart.figure)
                original_draw()

            chart.draw = draw
        original_destroy = shell._internal.destroy

        def destroy():
            self._apply_font_sizes(1.0)
            original_destroy()

        shell._internal.destroy = destroy

    @property
    def factor(self):
        """Current zoom as a multiplier."""
        return 1.1**self.level

    def change(self, steps):
        """Change zoom within a useful range."""
        self.level = max(-4, min(6, self.level + steps))
        self._apply_font_sizes(self.factor)
        self.logviewer.set_zoom(self.factor)
        for chart in self.charts:
            chart.draw()

    def _apply_font_sizes(self, factor):
        sizes = {
            name: max(1, round(size * factor)) for name, size in self.base_sizes.items()
        }
        # Preserve Bootstack's owned named-font objects while updating its tokens.
        Typography._fonts = Typography._fonts._replace(
            **{
                name.replace('-', '_'): Typography.get_token(name)._replace(size=size)
                for name, size in sizes.items()
            }
        )
        for name, size in sizes.items():
            Typography.get_named_token_font(name).configure(size=size)
        for tk_name, token in TK_FONT_OVERRIDES.items():
            tkfont.nametofont(tk_name, root=self.shell._internal).configure(
                size=sizes[token]
            )
        # Bootstack's modifier fonts are cached separately from their base tokens.
        for spec, font in _DERIVED_FONTS.items():
            token, size, delta, *_ = spec
            base = sizes[token]
            font.configure(size=max(1, round(size * factor) if size else base + delta))

    def _scale_figure(self, figure):
        for text in figure.findobj(match=Text):
            if text not in self.text_sizes:
                self.text_sizes[text] = text.get_fontsize()
            text.set_fontsize(self.text_sizes[text] * self.factor)
