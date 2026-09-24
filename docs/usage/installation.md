Installation
============

Install on any `*nix` system using python and pip, or clone this repository and install locally (run setup.py or pip).
`phys2cvr` supports python versions 3.6+. However, please note that no tests are currently run.

## Install with `pip` (recommended)

:exclamation::exclamation::exclamation: Please note that some systems might require to use `pip3` instead of `pip`.

## Basic installation:
For basic installation, simply run:
```shell
$ pip install phys2cvr
```

### Richer installation
To install the dependencies to enable more features, you can append labels to `phys2cvr`, e.g:
```shell
$ pip install phys2cvr[all]
```

The possible features are:

-  `[responses]`: to use physiological response functions from physiopy's `phys2denoise`.
-  `[matlab]`: to load and export MATLAB (`.mat`) files.
-  `[parallel]`: to run L-GLMs in parallel and drop execution process time (at the cost of extra resources).
-  `[all]`: to install all of the above.

## Clone from Github / install without `pip`

:exclamation::exclamation::exclamation: Please note that `phys2cvr` is continuously deployed, i.e. the latest feature available are immediately released on PyPI.
To install `phys2cvr` from Github, clone the repository first, then move to the cloned folder and run:
```shell
$ python setup.py install
```

Note that in this case, requirements for extra features need to be installed manually.

Alternatively, `pip` can be used too:
```shell
$ pip install .
```

## Run/use `phys2cvr`

You can run the `phys2cvr` workflow in a shell session (or in your code) - just follow the help:
```shell
$ phys2cvr --help
```

Alternatively, you can use `phys2cvr` as a module in a python session (or within your python script):

```python
import phys2cvr as p2c

p2c.__version__
```

## Developer installation

Follow the developer installation [pt. 1](../developers/how_to_contribute.html#linux-mac-and-windows-developer-installation) and [pt. 2](developers/how_to_contribute.html#full-developer-installation) in this documentation, [including installing pre-commit](../developers/how_to_contribute.html#install-pre-commit-hooks)
