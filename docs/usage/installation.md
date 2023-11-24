Installation
============

Install on any `*nix` system using python and pip, or clone this repository and install locally (run setup.py or pip).
`phys2cvr` supports python versions 3.6+. However, please note that no tests are currently run.

## Install with `pip` (recommended)

:exclamation::exclamation::exclamation: Please note that some systems might require to use `pip3` instead of `pip`.

### Basic installation:
For basic installation, simply run:
```shell
$ pip install phys2cvr
```

## Clone from Github / install without `pip`

:exclamation::exclamation::exclamation: Please note that `phys2cvr` is continuously deployed, i.e. the latest feature available are immediately released on PyPI.
To install `phys2cvr` from Github, clone the repository first, then move to the cloned folder and run:
```shell
$ python setup.py install
```

Alternatively, `pip` can be used too:
```shell
$ pip install .
```

## Developer installation
To be sure you have everything installed to develop (and test) `phys2cvr``, fork smoia/phys2cvr to your repository, then clone it locally and move inside the cloned folder. Finally, run the following commands from within the repository main folder:
```shell
# Add upstream remote
$ git remote add upstream git@github.com:smoia/phys2cvr.git

# Fetch everything, tags included
$ git fetch --all --tags

# Checkout master (the main development branch) and make it track upstream
$ git checkout master
$ git branch --set-upstream-to=upstream/master

# !!! VERY IMPORTANT !!!
# Set the default push to origin, in order NOT to push by mistake to upstream.
$ git config remote.pushDefault origin

# Install package with pip using the developer mode and the `[dev]` label
# You might need to use pip3 depending on how you set up your system
$ pip install -e .[dev]
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
