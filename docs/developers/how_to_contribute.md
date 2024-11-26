Getting Started!
================

First of all: thank you!

Contributions can be made in different ways, not only code! As we follow
the
[all-contributors](https://github.com/all-contributors/all-contributors)
specification, any contribution will be recognised accordingly.

Follow these steps to get started:

1.  Have a look at the [contributor guide](contributor_guide.md) page as
    well as the [code of conduct](code_of_conduct.md).
2.  Make sure that you have a GitHub account. You can set up a [free
    GitHub account](https://github.com/); here are some
    [instructions](https://help.github.com/articles/signing-up-for-a-new-github-account).
3.  If you intend to contribute code and/or use the `phys2cvr` packages
    in any way, check that you have `git` and `pip` installed on your
    system. Then install the package as a developer. This will let you
    run the program with the latest modifications, without requiring to
    re-install it every time.

!!! note
    The following instructions are provided assuming that python 3 is
    **not** your default version of python. If it is, you might need to use
    `pip` instead of `pip3`, although some OSs do adopt `pip3` anyway. If
    you want to check your OS behaviour, type `python --version` in a terminal.


## Linux, Mac, and Windows developer installation

Be sure to have `git` and `pip` installed. Fork the `phys2cvr` repository
in GitHub, then open a terminal and run the following code to clone the
forked repository and set it as your *origin*:

```shell
$ git clone https://github.com/{username}/phys2cvr.git
# or in case you prefer to use ssh:
$ git clone git@github.com:{username}/phys2cvr.git
```

We also recommend to set up the [`smoia/phys2cvr`](https://github.com/smoia/phys2cvr) repository as
*upstream*. In this way you can always keep your main branch
up to date with the command *git pull upstream master*:

```shell
$ cd phys2cvr
$ git remote add upstream https://github.com/smoia/phys2cvr.git
$ git pull upstream master
```

## Full developer installation

If it's your first experience as a python developer, or you just want
to be sure that you have everything you need to work on `phys2cvr`, you
can install it with all the other packages that are frequently
used during development in one step!

Go to the `phys2cvr` repository folder and execute the command:

```shell
$ cd phys2cvr
$ pip3 install -e .[dev]
```

This will install:

- `phys2cvr` as an editable package, which means that you can
 modify the program and run it without having to reinstall it every
 time!
- All `phys2cvr` required dependencies.
- All `phys2cvr` optional dependencies:
    + All packages used for **filetypes I/O** (`pymatreader`, `scipy`, `nibabel`)
    + All packages used for **visualisation and plotting** (`matplotlib`, `nilearn`)
- All **documentation** modules (`mkdocs` based), so that you can
 build the docs locally before submitting them.
- All **test** modules (`pytest`, `coverage`), in order for you to test your
 code locally before committing it!
- `pre-commit`, for commit hooks.

## Install pre-commit hooks

`pre-commit` is a tool that allows you to _automagically_ run smaller CI operations locally when committing code with git. For instance, it will check for merge conflicts and black the code, improving the quality of your PRs.

Go to the `phys2cvr` repository folder and execute the command:

```shell
$ cd phys2cvr
$ pre-commit install
```

To check that pre-commits are correctly installed, run:

```shell
$ pre-commit run --all-files
```

The output should look like:

```
trim trailing whitespace.................................................Passed
fix end of files.........................................................Passed
check yaml...............................................................Passed
check for added large files..............................................Passed
check for case conflicts.................................................Passed
check for merge conflicts................................................Passed
black....................................................................Passed
```

## Check your installation!

Type the commands:

```shell
$ cd phys2cvr
$ pytest
```

This will execute the tests locally and check that your phys2bids
installation works properly - it should look like this:

```
==================================== test session starts ====================================
platform linux -- Python 3.7.13, pytest-7.1.1, pluggy-1.0.0
rootdir: /home/nemo/Scrivania/gitlab/phys2cvr, configfile: setup.cfg
plugins: cov-3.0.0
collected 56 items

phys2cvr/tests/test_graph.py ..                                                           [  3%]
phys2cvr/tests/test_integration.py .                                                      [  5%]
phys2cvr/tests/test_io.py ..........................                                      [ 51%]
phys2cvr/tests/test_nifti.py ...........                                                  [ 71%]
phys2cvr/tests/test_utils.py ................                                             [100%]

============================== 56 passed in 96.57s (0:01:36) ================================

```

Do **not** worry if there is a xfail error in the log. This happens when
we know that a test will fail for known reasons, and we are probably
working to fix it (see
[here](https://docs.pytest.org/en/latest/skipping.html#xfail-mark-test-functions-as-expected-to-fail)).
However, if you do encounter any other error, check that you are connected to internet, you have all the extra dependencies installed, and their version meets `phys2cvr`
requirements.
