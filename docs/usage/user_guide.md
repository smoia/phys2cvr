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
CVR map, lag map, and mask. You can also pass paths and lag settings at startup;
run `python -m phys2cvr.gui.app --help` for the available options. From a
repository checkout, `bootstack run` also launches the GUI using
`bootstack.toml`.

## In-depth tutorial coming soon...
