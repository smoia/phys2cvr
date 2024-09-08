<a name="readme"></a>
<!-- <img alt="Phys2BIDS" src="https://github.com/physiopy/phys2bids/blob/master/docs/_static/phys2bids_logo1280×640.png" height="150"> -->

phys2cvr
========

[![Latest version](https://img.shields.io/pypi/v/phys2cvr?style=flat&logo=pypi)](https://pypi.org/project/phys2cvr/)
[![Release date](https://img.shields.io/github/release-date/MIPLabCH/nigsp?style=flat&logo=github)](https://github.com/MIPLabCH/nigsp/releases)
[![Auto Release](https://img.shields.io/badge/release-auto.svg?style=flat&colorA=888888&colorB=9B065A&label=auto&logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABQAAAAUCAYAAACNiR0NAAACzElEQVR4AYXBW2iVBQAA4O+/nLlLO9NM7JSXasko2ASZMaKyhRKEDH2ohxHVWy6EiIiiLOgiZG9CtdgG0VNQoJEXRogVgZYylI1skiKVITPTTtnv3M7+v8UvnG3M+r7APLIRxStn69qzqeBBrMYyBDiL4SD0VeFmRwtrkrI5IjP0F7rjzrSjvbTqwubiLZffySrhRrSghBJa8EBYY0NyLJt8bDBOtzbEY72TldQ1kRm6otana8JK3/kzN/3V/NBPU6HsNnNlZAz/ukOalb0RBJKeQnykd7LiX5Fp/YXuQlfUuhXbg8Di5GL9jbXFq/tLa86PpxPhAPrwCYaiorS8L/uuPJh1hZFbcR8mewrx0d7JShr3F7pNW4vX0GRakKWVk7taDq7uPvFWw8YkMcPVb+vfvfRZ1i7zqFwjtmFouL72y6C/0L0Ie3GvaQXRyYVB3YZNE32/+A/D9bVLcRB3yw3hkRCdaDUtFl6Ykr20aaLvKoqIXUdbMj6GFzAmdxfWx9iIRrkDr1f27cFONGMUo/gRI/jNbIMYxJOoR1cY0OGaVPb5z9mlKbyJP/EsdmIXvsFmM7Ql42nEblX3xI1BbYbTkXCqRnxUbgzPo4T7sQBNeBG7zbAiDI8nWfZDhQWYCG4PFr+HMBQ6l5VPJybeRyJXwsdYJ/cRnlJV0yB4ZlUYtFQIkMZnst8fRrPcKezHCblz2IInMIkPzbbyb9mW42nWInc2xmE0y61AJ06oGsXL5rcOK1UdCbEXiVwNXsEy/6+EbaiVG8eeEAfxvaoSBnCH61uOD7BS1Ul8ESHBKWxCrdyd6EYNKihgEVrwOAbQruoytuBYIFfAc3gVN6iawhjKyNCEpYhVJXgbOzARyaU4hCtYizq5EI1YgiUoIlT1B7ZjByqmRWYbwtdYjoWoN7+LOIQefIqKawLzK6ID69GGpQgwhhEcwGGUzfEPAiPqsCXadFsAAAAASUVORK5CYII=)](https://github.com/intuit/auto)

<!-- [![See the documentation at: https://nigsp.readthedocs.io](https://img.shields.io/badge/docs-read%20latest-informational?style=flat&logo=readthedocs)](https://nigsp.readthedocs.io/en/latest/?badge=latest) -->
[![Latest DOI](https://zenodo.org/badge/357980417.svg)](https://zenodo.org/badge/latestdoi/357980417)
[![Licensed Apache 2.0](https://img.shields.io/github/license/smoia/phys2cvr?style=flat)](https://github.com/smoia/phys2cvr/blob/master/LICENSE)

<!-- [![Codecov](https://img.shields.io/codecov/c/gh/MIPlabCH/nigsp?style=flat&label=codecov&logo=codecov)](https://codecov.io/gh/MIPLabCH/nigsp)
[![Build Status](https://img.shields.io/circleci/build/github/MIPLabCH/nigsp?style=flat&label=circleci&logo=circleci)](https://circleci.com/gh/MIPLabCH/nigsp)
[![Documentation Status](https://img.shields.io/readthedocs/nigsp?style=flat&label=readthedocs&logo=readthedocs)](https://nigsp.readthedocs.io/en/latest/?badge=latest) -->

[![Latest version](https://img.shields.io/pypi/v/phys2cvr?style=flat&logo=pypi&logoColor=white)](https://pypi.org/project/phys2cvr/)
[![Supports python version](https://img.shields.io/pypi/pyversions/phys2cvr?style=shield&logo=python)](https://pypi.org/project/phys2cvr/)

<!-- ALL-CONTRIBUTORS-BADGE:START - Do not remove or modify this section -->
[![All Contributors](https://img.shields.io/badge/all_contributors-3-orange.svg?style=flat)](#contributors)
<!-- ALL-CONTRIBUTORS-BADGE:END -->

A python-based tool to generate regressor for and/or estimate CVR maps and their lag.

**The project is currently under development stage alpha**.
Any suggestion/bug report is welcome! Feel free to open an issue.

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification. Contributions of any kind welcome!

Documentation
=============

Full documentation coming soon!

Cite
----

If you use `phys2cvr` in your work, please cite either the all-time Zenodo DOI [![general Zenodo DOI](https://zenodo.org/badge/5559756.svg)](https://doi.org/10.5281/zenodo.5559756) or the Zenodo DOI related to the version you are using.
Please cite the following paper(s) too:
>Moia, S., Stickland, R. C., Ayyagari, A., Termenon, M., Caballero-Gaudes, C., & Bright, M. G. (2020). Voxelwise optimization of hemodynamic lags to improve regional CVR estimates in breath-hold fMRI. In 2020 42nd Annual International Conference of the IEEE Engineering in Medicine & Biology Society (EMBC) (pp. 1489–1492). Montreal, QC, Canada: IEEE. https://doi.org/10.1109/EMBC44109.2020.9176225

If you are using the `--brightspin` configuration option:
>Moia, S., Termenon, M., Uruñuela, E., Chen, G., Stickland, R. C., Bright, M. G., & Caballero-Gaudes, C. (2021). ICA-based denoising strategies in breath-hold induced cerebrovascular reactivity mapping with multi echo BOLD fMRI. NeuroImage, 233, 117914. https://doi.org/10.1016/j.neuroimage.2021.117914

If you are using the `--brightspin-clinical` configuration option:
>Stickland, R. C., Zvolanek, K. M., Moia, S., Ayyagari, A., & Bright, M. G. (2021). A practical modification to a resting state fMRI protocol for improved characterization of cerebrovascular function. Supplementary Material. Neuroimage.

If you are using the `--baltimore-lag` configuration option:
>Liu, P., Li, Y., Pinho, M., Park, D. C., Welch, B. G., & Lu, H. (2017). Cerebrovascular reactivity mapping without gas challenges. NeuroImage, 146(November 2016), 320–326. https://doi.org/10.1016/j.neuroimage.2016.11.054

If you are using the `--baltimore` configuration option, please cite only the Zenodo DOI and the last listed paper.

Installation
------------

Install on any `*nix` system using python and pip, or clone this repository and install locally (run setup.py or pip).
`phys2cvr` supports python versions 3.6+. However, please note that no tests are currently run.

### Install with `pip` (recommended)

:exclamation::exclamation::exclamation: Please note that some systems might require to use `pip3` instead of `pip`.

#### Basic installation:
For basic installation, simply run:
```bash
pip install phys2cvr
```

### Clone from Github / install without `pip`

:exclamation::exclamation::exclamation: Please note that `phys2cvr` is continuously deployed, i.e. the latest feature available are immediately released on PyPI.
To install `phys2cvr` from Github, clone the repository first, then move to the cloned folder and run:
```bash
python setup.py install
```

Alternatively, `pip` can be used too:
```bash
pip install .
```

### Developer installation

To be sure you have everything installed to develop (and test) `phys2cvr`, **fork** `smoia/phys2cvr` to your repository, then clone it locally and move inside the cloned folder. Finally, run the following commands from within the repository main folder:
```bash
# Add upstream remote
git remote add upstream git@github.com:smoia/phys2cvr.git

# Fetch everything, tags included
git fetch --all --tags

# Checkout master (the main development branch) and make it track upstream
git checkout master
git branch --set-upstream-to=upstream/master

# !!! VERY IMPORTANT !!!
# Set the default push to origin, in order NOT to push by mistake to upstream.
git config remote.pushDefault origin

# Install package with pip using the developer mode and the `[dev]` label
# You might need to use pip3 depending on how you set up your system
pip install -e .[dev]
```
If you make changes that you consider fundamental/interesting for the whole community, feel free to open a PR!

Run/use `phys2cvr`
---------------

You can run the `phys2cvr` workflow in a shell session (or in your code) - just follow the help:
```shell
phys2cvr --help
```

Alternatively, you can use phys2cvr as a module in a python session (or within your python script):
```python
import phys2cvr as p2c

p2c.__version__
```

Full API coming soon.


## Contributors ✨

Thanks goes to these wonderful people ([emoji key](https://allcontributors.org/docs/en/emoji-key)):

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tr>
    <td align="center"><a href="https://github.com/smoia"><img src="https://avatars3.githubusercontent.com/u/35300580?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Stefano Moia</b></sub></a><br /><a href="https://github.com/smoia/phys2cvr/commits?author=smoia" title="Code">💻</a> <a href="#ideas-smoia" title="Ideas, Planning, & Feedback">🤔</a> <a href="#infra-smoia" title="Infrastructure (Hosting, Build-Tools, etc)">🚇</a> <a href="#projectManagement-smoia" title="Project Management">📆</a></td>
    <td align="center"><a href="https://github.com/kristinazvolanek"><img src="https://avatars3.githubusercontent.com/u/54590158?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Kristina Zvolanek</b></sub></a><br /><a href="https://github.com/smoia/phys2cvr/commits?author=kristinazvolanek" title="Code">💻</a> <a href="https://github.com/smoia/phys2cvr/issues?q=author%3Akristinazvolanek" title="Bug reports">🐛</a></td>
    <td align="center"><a href="https://github.com/avigotsky"><img src="https://avatars.githubusercontent.com/u/904218?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Andrew Vigotsky</b></sub></a><br /><a href="https://github.com/smoia/phys2cvr/commits?author=avigotsky" title="Code">💻</a></td>
  </tr>
</table>

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->


License
-------

Copyright 2021, Stefano Moia.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
