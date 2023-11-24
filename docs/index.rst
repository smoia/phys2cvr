.. phys2cvr documentation master file, created by
   sphinx-quickstart on Tue Jul 12 11:04:40 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root ``toctree`` directive.

:hide-toc:

phys2cvr
========

|Latest version| |Release date| |Auto Release|

|See the documentation at: https://phys2cvr.readthedocs.io| |Latest DOI|
|Licensed Apache 2.0|

|Codecov| |Build Status| |Documentation Status|

|image1| |Supports python version|

|All Contributors|

A python-based tool to generate regressor for and/or estimate CVR maps and their lag.

**The project is currently under development stage alpha**. Any
suggestion/bug report is welcome! Feel free to open an issue.

This project follows the
`all-contributors <https://github.com/all-contributors/all-contributors>`__
specification. Contributions of any kind welcome!

Cite
----

If you use ``phys2cvr`` in your work, please cite either the all-time
Zenodo DOI [![general Zenodo DOI](https://zenodo.org/badge/5559756.svg)](https://doi.org/10.5281/zenodo.5559756) or the Zenodo DOI related to the version
you are using. Please cite the following paper(s) too:

   Moia, S., Stickland, R. C., Ayyagari, A., Termenon, M., Caballero-Gaudes, C.,
   & Bright, M. G. (2020). *Voxelwise optimization of hemodynamic lags to improve regional CVR estimates in breath-hold fMRI.*
   In 2020 42nd Annual International Conference of the IEEE Engineering in Medicine & Biology Society (EMBC) (pp. 1489–1492).
   Montreal, QC, Canada: IEEE. `https://doi.org/10.1109/EMBC44109.2020.9176225 <https://doi.org/10.1109/EMBC44109.2020.9176225>`__.

If you are using the ``--brightspin`` configuration option:

   Moia, S., Termenon, M., Uruñuela, E., Chen, G., Stickland, R. C., Bright, M. G., & Caballero-Gaudes, C. (2021).
   *ICA-based denoising strategies in breath-hold induced cerebrovascular reactivity mapping with multi echo BOLD fMRI.*
   NeuroImage, 233, 117914.
   `https://doi.org/10.1016/j.neuroimage.2021.117914 <https://doi.org/10.1016/j.neuroimage.2021.117914>`__.

If you are using the ``--brightspin-clinical`` configuration option:
   Stickland, R. C., Zvolanek, K. M., Moia, S., Ayyagari, A., & Bright, M. G. (2021).
   *A practical modification to a resting state fMRI protocol for improved characterization of cerebrovascular function.*
   Supplementary Material. Neuroimage.

If you are using the ``--baltimore-lag`` configuration option:
   Liu, P., Li, Y., Pinho, M., Park, D. C., Welch, B. G., & Lu, H. (2017). *Cerebrovascular reactivity mapping without gas challenges.*
   NeuroImage, 146(November 2016), 320–326. `https://doi.org/10.1016/j.neuroimage.2016.11.054 <https://doi.org/10.1016/j.neuroimage.2016.11.054>`__.

If you are using the ``--baltimore`` configuration option, please cite only the Zenodo DOI and the last listed paper.

.. |Latest version| image:: https://img.shields.io/github/v/release/smoia/phys2cvr?style=flat&logo=github&sort=semver
   :target: https://github.com/smoia/phys2cvr/releases
.. |Release date| image:: https://img.shields.io/github/release-date/smoia/phys2cvr?style=flat&logo=github
   :target: https://github.com/smoia/phys2cvr/releases
.. |Auto Release| image:: https://img.shields.io/badge/release-auto.svg?style=flat&colorA=888888&colorB=9B065A&label=auto&logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABQAAAAUCAYAAACNiR0NAAACzElEQVR4AYXBW2iVBQAA4O+/nLlLO9NM7JSXasko2ASZMaKyhRKEDH2ohxHVWy6EiIiiLOgiZG9CtdgG0VNQoJEXRogVgZYylI1skiKVITPTTtnv3M7+v8UvnG3M+r7APLIRxStn69qzqeBBrMYyBDiL4SD0VeFmRwtrkrI5IjP0F7rjzrSjvbTqwubiLZffySrhRrSghBJa8EBYY0NyLJt8bDBOtzbEY72TldQ1kRm6otana8JK3/kzN/3V/NBPU6HsNnNlZAz/ukOalb0RBJKeQnykd7LiX5Fp/YXuQlfUuhXbg8Di5GL9jbXFq/tLa86PpxPhAPrwCYaiorS8L/uuPJh1hZFbcR8mewrx0d7JShr3F7pNW4vX0GRakKWVk7taDq7uPvFWw8YkMcPVb+vfvfRZ1i7zqFwjtmFouL72y6C/0L0Ie3GvaQXRyYVB3YZNE32/+A/D9bVLcRB3yw3hkRCdaDUtFl6Ykr20aaLvKoqIXUdbMj6GFzAmdxfWx9iIRrkDr1f27cFONGMUo/gRI/jNbIMYxJOoR1cY0OGaVPb5z9mlKbyJP/EsdmIXvsFmM7Ql42nEblX3xI1BbYbTkXCqRnxUbgzPo4T7sQBNeBG7zbAiDI8nWfZDhQWYCG4PFr+HMBQ6l5VPJybeRyJXwsdYJ/cRnlJV0yB4ZlUYtFQIkMZnst8fRrPcKezHCblz2IInMIkPzbbyb9mW42nWInc2xmE0y61AJ06oGsXL5rcOK1UdCbEXiVwNXsEy/6+EbaiVG8eeEAfxvaoSBnCH61uOD7BS1Ul8ESHBKWxCrdyd6EYNKihgEVrwOAbQruoytuBYIFfAc3gVN6iawhjKyNCEpYhVJXgbOzARyaU4hCtYizq5EI1YgiUoIlT1B7ZjByqmRWYbwtdYjoWoN7+LOIQefIqKawLzK6ID69GGpQgwhhEcwGGUzfEPAiPqsCXadFsAAAAASUVORK5CYII=
   :target: https://github.com/intuit/auto
.. |See the documentation at: https://phys2cvr.readthedocs.io| image:: https://img.shields.io/badge/docs-read%20latest-informational?style=flat&logo=readthedocs
   :target: https://phys2cvr.readthedocs.io/en/latest/?badge=latest
.. |Latest DOI| image:: https://zenodo.org/badge/446805866.svg?style=flat&logo=zenodo
   :target: https://zenodo.org/badge/latestdoi/446805866
.. |Licensed Apache 2.0| image:: https://img.shields.io/github/license/smoia/phys2cvr?style=flat&logo=apache
   :target: https://github.com/smoia/phys2cvr/blob/master/LICENSE
.. |Codecov| image:: https://img.shields.io/codecov/c/gh/smoia/phys2cvr?style=flat&label=codecov&logo=codecov
   :target: https://codecov.io/gh/smoia/phys2cvr
.. |Build Status| image:: https://img.shields.io/circleci/build/github/smoia/phys2cvr?style=flat&label=circleci&logo=circleci
   :target: https://circleci.com/gh/smoia/phys2cvr
.. |Documentation Status| image:: https://img.shields.io/readthedocs/phys2cvr?style=flat&label=readthedocs&logo=readthedocs
   :target: https://phys2cvr.readthedocs.io/en/latest/?badge=latest
.. |image1| image:: https://img.shields.io/pypi/v/phys2cvr?style=flat&logo=pypi&logoColor=white
   :target: https://pypi.org/project/phys2cvr/
.. |Supports python version| image:: https://img.shields.io/pypi/pyversions/phys2cvr?style=flat&logo=python&logoColor=white
   :target: https://pypi.org/project/phys2cvr/
.. |All Contributors| image:: https://img.shields.io/badge/all_contributors-1-orange.svg?style=flat
   :target: #contributors
.. |general Zenodo DOI| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.6373436.svg
   :target: https://zenodo.org/badge/latestdoi/446805866

.. toctree::
   :caption: Usage
   :hidden:
   :maxdepth: 2

   Installation <usage/installation>
   User Guide <usage/user_guide>
   Command Line Interface (CLI) <usage/cli>
   Output <usage/output>
   Licence <usage/licence>
   Changelog <https://github.com/smoia/phys2cvr/releases>

.. toctree::
   :caption: API
   :hidden:
   :maxdepth: 1
   :glob:

   API <api>

.. toctree::
   :caption: Cerebrovascular reactivity mapping
   :hidden:
   :maxdepth: 1

   About phys2cvr <about_phys2cvr>

.. toctree::
   :caption: Developers
   :hidden:
   :maxdepth: 1

   How to Contribute <developers/how_to_contribute>
   Contributor Guide <developers/contributor_guide>
   Code of Conduct <developers/code_of_conduct>
