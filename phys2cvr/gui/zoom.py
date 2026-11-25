"""Live application font and chart zoom."""

from weakref import WeakKeyDictionary

from bootstack.style import update_font_token
from bootstack.style.typography import _DERIVED_FONTS, Typography
from matplotlib.text import Text


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
        shell._internal.bind('<Destroy>', self._on_destroy, add='+')

    @property
    def factor(self):
        """Current zoom as a multiplier."""
        return 1.1**self.level

    def change(self, steps):
        """Change zoom within a useful range."""
        self.level = max(-4, min(6, self.level + steps))
        for name, size in self.base_sizes.items():
            update_font_token(name, size=max(1, round(size * self.factor)))
        # Bootstack's modifier fonts are cached separately from their base tokens.
        for spec, font in _DERIVED_FONTS.items():
            token, size, delta, *_ = spec
            base = Typography.get_token(token).size
            font.configure(
                size=max(1, round(size * self.factor) if size else base + delta)
            )
        self.logviewer.set_zoom(self.factor)
        for chart in self.charts:
            chart.draw()

    def _scale_figure(self, figure):
        for text in figure.findobj(match=Text):
            if text not in self.text_sizes:
                self.text_sizes[text] = text.get_fontsize()
            text.set_fontsize(self.text_sizes[text] * self.factor)

    def _on_destroy(self, event):
        if event.widget is self.shell._internal:
            for name, size in self.base_sizes.items():
                update_font_token(name, size=size)
