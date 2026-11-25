Use/run phys2cvr
================

You can run the `phys2cvr` workflow in a shell session (or in your code) - just follow the help:
```shell
$ phys2cvr --help
```

Alternatively, you can use phys2cvr as a module in a python session (or within your python script):
```python
import phys2cvr

phys2cvr.__version__
```

## Graphical voxel viewer

The optional GUI requires Python 3.12 or newer and Tk. Install it with
`pip install 'phys2cvr[gui]'`, then start the viewer with:

```shell
python -m phys2cvr.gui.app
```

Use the viewer's **Load** menu to open a 4D fMRI timeseries, regressor matrix,
CVR map, lag map, R² volume, mask, or result folder. These inputs can be viewed
independently. The left settings panel controls voxel coordinates, map color
limits and colormaps, lag bounds, and whether timecourses are z-scored or
demeaned. The lag slider moves by one regressor row; when the sampling frequency
is unspecified, it is inferred from the lag range and matrix size.

The **Comparison** page loads base and comparison result folders or individual
CVR, lag, and R² maps. Inputs must share the same spatial grid and affine;
resample them before loading if necessary. A 4D R² volume is reduced to its
maximum across lags at each voxel. Each column's difference controls select
signed or absolute values and either base minus comparison or the relative
difference `100 * (base - comparison) / base`. Relative differences are masked
where the base value is zero. Clicking any slice updates the shared voxel
position across all maps.
The bottom row compares masked CVR and lag distributions and, for 4D R²
inputs, plots the selected voxel's R² values against lag index. Each dataset
uses its own mask; if only one mask exists, it is shared by both datasets.
Loading an R² image without a corresponding mask creates one from nonzero,
finite values of the R² map (or its mean across lags for 4D input). A later
loaded mask replaces the generated mask. 3D R² maps have no lag profile.

You can also pass paths and lag settings at startup;
run `python -m phys2cvr.gui.app --help` for the available options. From a
repository checkout, `bootstack run` also launches the GUI using
`bootstack.toml`.

The **Run phys2cvr** page provides the workflow's command-line options as
grouped controls. Hover over a control for its CLI help. Matrix inputs accept
one file path per line or multiple files selected with the browse button.
Optional numeric fields may be left blank to use the CLI default. Select a
custom response-function file using the response-function selector and its
file field. Selecting a workflow preset applies its command-line overrides to
the form so the resulting settings are visible. The preview column shows the mean functional volume, regressor,
mask, and ROI mask as their paths are selected. Click a slice to move its
crosshairs or scroll over it to change slices; a 1D functional input has no
volume preview. The regressor preview marks peaks from a `.phys` file or an
explicit peaks file. Loading an image or `.phys` file fills the TR or frequency
from its metadata. Controls unavailable under the selected options are disabled.
The run button starts a separate process, keeping the GUI responsive;
the page reports completion or an error. Closing the application stops an
active workflow.

The **Log viewer** loads a log and displays PNG images found in its result
folder in a navigable carousel. Use the slider beside the log to change its
font size.

To run the GUI tests with pytest, install the test and GUI extras and use a
graphical display. On headless Linux, use Xvfb:

```shell
pip install -e '.[test,gui]'
xvfb-run -a python -m pytest phys2cvr/tests/test_gui.py \
  phys2cvr/tests/test_diffviewer.py phys2cvr/tests/test_logviewer.py \
  phys2cvr/tests/test_guiargset.py
```

Display-dependent tests are skipped when no display is available.

## In-depth tutorial coming soon...
