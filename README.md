<a name="readme"></a>
<!-- <img alt="Phys2BIDS" src="https://github.com/physiopy/phys2bids/blob/master/docs/_static/phys2bids_logo1280×640.png" height="150"> -->

phys2cvr
========

[![Latest version](https://img.shields.io/pypi/v/phys2cvr?style=flat&logo=pypi)](https://pypi.org/project/phys2cvr/)
[![Release date](https://img.shields.io/github/release-date/smoia/phys2cvr?style=flat&logo=github)](https://github.com/smoia/phys2cvr/releases)
[![Auto Release](https://img.shields.io/badge/release-auto.svg?style=flat&colorA=888888&colorB=9B065A&label=auto&logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABQAAAAUCAYAAACNiR0NAAACzElEQVR4AYXBW2iVBQAA4O+/nLlLO9NM7JSXasko2ASZMaKyhRKEDH2ohxHVWy6EiIiiLOgiZG9CtdgG0VNQoJEXRogVgZYylI1skiKVITPTTtnv3M7+v8UvnG3M+r7APLIRxStn69qzqeBBrMYyBDiL4SD0VeFmRwtrkrI5IjP0F7rjzrSjvbTqwubiLZffySrhRrSghBJa8EBYY0NyLJt8bDBOtzbEY72TldQ1kRm6otana8JK3/kzN/3V/NBPU6HsNnNlZAz/ukOalb0RBJKeQnykd7LiX5Fp/YXuQlfUuhXbg8Di5GL9jbXFq/tLa86PpxPhAPrwCYaiorS8L/uuPJh1hZFbcR8mewrx0d7JShr3F7pNW4vX0GRakKWVk7taDq7uPvFWw8YkMcPVb+vfvfRZ1i7zqFwjtmFouL72y6C/0L0Ie3GvaQXRyYVB3YZNE32/+A/D9bVLcRB3yw3hkRCdaDUtFl6Ykr20aaLvKoqIXUdbMj6GFzAmdxfWx9iIRrkDr1f27cFONGMUo/gRI/jNbIMYxJOoR1cY0OGaVPb5z9mlKbyJP/EsdmIXvsFmM7Ql42nEblX3xI1BbYbTkXCqRnxUbgzPo4T7sQBNeBG7zbAiDI8nWfZDhQWYCG4PFr+HMBQ6l5VPJybeRyJXwsdYJ/cRnlJV0yB4ZlUYtFQIkMZnst8fRrPcKezHCblz2IInMIkPzbbyb9mW42nWInc2xmE0y61AJ06oGsXL5rcOK1UdCbEXiVwNXsEy/6+EbaiVG8eeEAfxvaoSBnCH61uOD7BS1Ul8ESHBKWxCrdyd6EYNKihgEVrwOAbQruoytuBYIFfAc3gVN6iawhjKyNCEpYhVJXgbOzARyaU4hCtYizq5EI1YgiUoIlT1B7ZjByqmRWYbwtdYjoWoN7+LOIQefIqKawLzK6ID69GGpQgwhhEcwGGUzfEPAiPqsCXadFsAAAAASUVORK5CYII=)](https://github.com/intuit/auto)

<!-- [![See the documentation at: https://phys2cvr.readthedocs.io](https://img.shields.io/badge/docs-read%20latest-informational?style=flat&logo=readthedocs)](https://phys2cvr.readthedocs.io/en/latest/?badge=latest) -->
[![Latest DOI](https://zenodo.org/badge/357980417.svg)](https://doi.org/10.5281/zenodo.5559756)
[![Licensed Apache 2.0](https://img.shields.io/github/license/smoia/phys2cvr?style=flat)](https://github.com/smoia/phys2cvr/blob/master/LICENSE)

<!-- [![Codecov](https://img.shields.io/codecov/c/gh/smoia/phys2cvr?style=flat&label=codecov&logo=codecov)](https://codecov.io/gh/smoia/phys2cvr)
[![Build Status](https://img.shields.io/circleci/build/github/smoia/phys2cvr?style=flat&label=circleci&logo=circleci)](https://circleci.com/gh/smoia/phys2cvr) -->
[![Documentation Status](https://img.shields.io/readthedocs/phys2cvr?style=flat&label=readthedocs&logo=readthedocs)](https://phys2cvr.readthedocs.io/en/latest/?badge=latest)

[![Latest version](https://img.shields.io/pypi/v/phys2cvr?style=flat&logo=pypi&logoColor=white)](https://pypi.org/project/phys2cvr/)
[![Supports python version](https://img.shields.io/pypi/pyversions/phys2cvr?style=shield&logo=python)](https://pypi.org/project/phys2cvr/)

<!-- ALL-CONTRIBUTORS-BADGE:START - Do not remove or modify this section -->
[![All Contributors](https://img.shields.io/badge/all_contributors-7-orange.svg?style=flat)](#contributors)
<!-- ALL-CONTRIBUTORS-BADGE:END -->

A python-based tool to generate regressor for and/or estimate CVR maps and their lag.

**The project is currently under development stage alpha**.
Any suggestion/bug report is welcome! Feel free to open an issue.

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification. Contributions of any kind welcome!

Documentation
=============

Full documentation [here](https://phys2cvr.readthedocs.io/en/latest/)

Cite
----

If you use `phys2cvr` in your work, please cite either the all-time Zenodo DOI [![general Zenodo DOI](https://zenodo.org/badge/5559756.svg)](https://doi.org/10.5281/zenodo.5559756) or the Zenodo DOI related to the version you are using.
Please cite the following paper(s) too:
>Moia, S., Stickland, R. C., Ayyagari, A., Termenon, M., Caballero-Gaudes, C., & Bright, M. G. (2020). Voxelwise optimization of hemodynamic lags to improve regional CVR estimates in breath-hold fMRI. In 2020 42nd Annual International Conference of the IEEE Engineering in Medicine & Biology Society (EMBC) (pp. 1489–1492). Montreal, QC, Canada: IEEE. https://doi.org/10.1109/EMBC44109.2020.9176225

If you are using the `--brightspin` configuration option:
>Moia, S., Termenon, M., Uruñuela, E., Chen, G., Stickland, R. C., Bright, M. G., & Caballero-Gaudes, C. (2021). ICA-based denoising strategies in breath-hold induced cerebrovascular reactivity mapping with multi echo BOLD fMRI. NeuroImage, 233, 117914. https://doi.org/10.1016/j.neuroimage.2021.117914

If you are using the ``--brightspin-clinical`` configuration option:
>Stickland, R. C., Zvolanek, K. M., Moia, S., Ayyagari, A., & Bright, M. G. (2021). A practical modification to a resting state fMRI protocol for improved characterization of cerebrovascular function. Neuroimage, 239, 118306. https://doi.org/10.1016/j.neuroimage.2021.118306

If you are using the `--baltimore-lag` configuration option:
>Liu, P., Li, Y., Pinho, M., Park, D. C., Welch, B. G., & Lu, H. (2017). Cerebrovascular reactivity mapping without gas challenges. NeuroImage, 146, 320–326. https://doi.org/10.1016/j.neuroimage.2016.11.054

If you are using the `--baltimore` configuration option, please cite only the Zenodo DOI and the last listed paper.

If you are using the ``--fourier-set-lag`` configuration option:
>Pinto, J., Jorge, J., Sousa, I., Vilela, P., & Figueiredo, P. (2016). Fourier modeling of the BOLD response to a breath-hold task: Optimization and reproducibility.    NeuroImage, 135, 223-231. https://doi.org/10.1016/j.neuroimage.2016.02.037

If you are using the ``--fourier-set`` configuration option, please cite only the Zenodo DOI and this last paper.


Installation
------------

Instructions [here](https://phys2cvr.readthedocs.io/en/latest/usage/installation.html#installation)

### Developer installation

(Potential) Contributors, instead see [here](https://phys2cvr.readthedocs.io/en/latest/developers/how_to_contribute.html#linux-mac-and-windows-developer-installation)!

Run/use `phys2cvr`
---------------

You can run the `phys2cvr` workflow in a shell session (or in your code) - just follow the help or see [here](https://phys2cvr.readthedocs.io/en/latest/usage/cli.html):
```shell
phys2cvr --help
```

Alternatively, you can use phys2cvr as a module in a python session (or within your python script):
```python
import phys2cvr as p2c

p2c.__version__
```

Full API [here](https://phys2cvr.readthedocs.io/en/latest/api.html)


## Contributors ✨

Thanks goes to these wonderful people ([emoji key](https://allcontributors.org/docs/en/emoji-key)):

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tbody>
    <tr>
      <td align="center" valign="top" width="20%"><a href="https://github.com/smoia"><img src="https://avatars3.githubusercontent.com/u/35300580?v=4?s=100" width="100px;" alt="Stefano Moia"/><br /><sub><b>Stefano Moia</b></sub></a><br /><a href="https://github.com/smoia/phys2cvr/commits?author=smoia" title="Code">💻</a> <a href="#ideas-smoia" title="Ideas, Planning, & Feedback">🤔</a> <a href="#infra-smoia" title="Infrastructure (Hosting, Build-Tools, etc)">🚇</a> <a href="#projectManagement-smoia" title="Project Management">📆</a></td>
      <td align="center" valign="top" width="20%"><a href="https://github.com/kristinazvolanek"><img src="https://avatars3.githubusercontent.com/u/54590158?v=4?s=100" width="100px;" alt="Kristina Zvolanek"/><br /><sub><b>Kristina Zvolanek</b></sub></a><br /><a href="https://github.com/smoia/phys2cvr/commits?author=kristinazvolanek" title="Code">💻</a> <a href="https://github.com/smoia/phys2cvr/issues?q=author%3Akristinazvolanek" title="Bug reports">🐛</a> <a href="#infra-kristinazvolanek" title="Infrastructure (Hosting, Build-Tools, etc)">🚇</a></td>
      <td align="center" valign="top" width="20%"><a href="https://github.com/beccaclements99"><img src="https://avatars.githubusercontent.com/u/95257798?v=4?s=100" width="100px;" alt="Becca Clements"/><br /><sub><b>Becca Clements</b></sub></a><br /><a href="https://github.com/smoia/phys2cvr/issues?q=author%3Abeccaclements99" title="Bug reports">🐛</a> <a href="https://github.com/smoia/phys2cvr/commits?author=beccaclements99" title="Code">💻</a> <a href="#userTesting-beccaclements99" title="User Testing">📓</a></td>
      <td align="center" valign="top" width="20%"><a href="https://github.com/avigotsky"><img src="https://avatars.githubusercontent.com/u/904218?v=4?s=100" width="100px;" alt="Andrew Vigotsky"/><br /><sub><b>Andrew Vigotsky</b></sub></a><br /><a href="https://github.com/smoia/phys2cvr/commits?author=avigotsky" title="Code">💻</a></td>
      <td align="center" valign="top" width="20%"><a href="https://github.com/ccomellalue"><img src="https://avatars.githubusercontent.com/u/209748720?v=4?s=100" width="100px;" alt="ccomellalue"/><br /><sub><b>ccomellalue</b></sub></a><br /><a href="https://github.com/smoia/phys2cvr/commits?author=ccomellalue" title="Code">💻</a> <a href="#ideas-ccomellalue" title="Ideas, Planning, & Feedback">🤔</a></td>
    </tr>
    <tr>
      <td align="center" valign="top" width="20%"><a href="https://github.com/merelvdthiel"><img src="https://avatars.githubusercontent.com/u/72999546?v=4?s=100" width="100px;" alt="merelvdthiel"/><br /><sub><b>merelvdthiel</b></sub></a><br /><a href="https://github.com/smoia/phys2cvr/commits?author=merelvdthiel" title="Documentation">📖</a></td>
      <td align="center" valign="top" width="20%"><a href="https://github.com/RazkinMalen"><img src="https://avatars.githubusercontent.com/u/147508933?v=4?s=100" width="100px;" alt="razkin"/><br /><sub><b>razkin</b></sub></a><br /><a href="#design-RazkinMalen" title="Design">🎨</a> <a href="https://github.com/smoia/phys2cvr/commits?author=RazkinMalen" title="Documentation">📖</a> <a href="https://github.com/smoia/phys2cvr/commits?author=RazkinMalen" title="Tests">⚠️</a></td>
    </tr>
  </tbody>
</table>

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->


License
-------

Copyright 2021-2026, Stefano Moia & phys2cvr contributors.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
