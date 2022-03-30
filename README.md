<a name="readme"></a>
<!-- <img alt="Phys2BIDS" src="https://github.com/physiopy/phys2bids/blob/master/docs/_static/phys2bids_logo1280×640.png" height="150"> -->

phys2cvr
========

[![Latest version](https://img.shields.io/pypi/v/phys2cvr?style=flat&logo=pypi)](https://pypi.org/project/phys2cvr/)
[![Latest DOI](https://zenodo.org/badge/357980417.svg)](https://zenodo.org/badge/latestdoi/357980417)
[![Licensed Apache 2.0](https://img.shields.io/github/license/smoia/phys2cvr?style=flat)](https://github.com/smoia/phys2cvr/blob/master/LICENSE)

[![Auto Release](https://img.shields.io/badge/release-auto.svg?colorA=888888&colorB=9B065A&label=auto)](https://github.com/intuit/auto)
[![Supports python version](https://img.shields.io/pypi/pyversions/phys2cvr?style=shield&logo=python)](https://pypi.org/project/phys2cvr/)

<!-- ALL-CONTRIBUTORS-BADGE:START - Do not remove or modify this section -->
[![All Contributors](https://img.shields.io/badge/all_contributors-1-orange.svg?style=flat)](#contributors)
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

If you use `phys2cvr` in your work, please cite either the all-time Zenodo DOI [![general Zenodo DOI](https://zenodo.org/badge/110845855.svg)](https://zenodo.org/badge/latestdoi/110845855) or the Zenodo DOI related to the version you are using.
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

To be sure you have everything installed to develop (and test) `phys2cvr`, **fork** `smoia/phys2cvr` to your repository, then clone it locally and move inside the cloned folder. Finally, install with `pip` using the developer mode and the `[all]` label:
```bash
pip install -e .[all]
```


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


<!-- ## Contributors ✨

Thanks goes to these wonderful people ([emoji key](https://allcontributors.org/docs/en/emoji-key)): -->

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->


<!-- markdownlint-enable -->
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
