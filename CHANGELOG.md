# 0.30.1 (Thu Sep 24 2026)

#### 🐛 Bug Fix

- Fix OLS function's poor handling of rank deficient matrices that caused R² to hit 1 (thus setting lags at their minimum) and refactor stats [#168](https://github.com/smoia/phys2cvr/pull/168) ([@smoia](https://github.com/smoia))

#### Authors: 1

- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.30.0 (Tue Sep 22 2026)

#### 🚀 Enhancement

- Support I/O of Freesurfer's motion parameters and fix debug volumes header's pixdim [#166](https://github.com/smoia/phys2cvr/pull/166) ([@smoia](https://github.com/smoia))

#### 🐛 Bug Fix

- Fix runner [#167](https://github.com/smoia/phys2cvr/pull/167) ([@smoia](https://github.com/smoia))

#### Authors: 1

- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.29.1 (Tue Sep 22 2026)

#### 🐛 Bug Fix

- Export lag index if debug flag is True, and change header of exported volumes to reflect the lag step timescale rather than the neuroimaging TR [#165](https://github.com/smoia/phys2cvr/pull/165) ([@smoia](https://github.com/smoia))

#### 🏠 Internal

- Bump actions/labeler from 6 to 7 [#164](https://github.com/smoia/phys2cvr/pull/164) ([@dependabot[bot]](https://github.com/dependabot[bot]))
- Bump actions/labeler from 6 to 7 [#163](https://github.com/smoia/phys2cvr/pull/163) ([@dependabot[bot]](https://github.com/dependabot[bot]))
- Bump actions/setup-python from 6 to 7 [#162](https://github.com/smoia/phys2cvr/pull/162) ([@dependabot[bot]](https://github.com/dependabot[bot]))
- [pre-commit.ci] pre-commit autoupdate [#161](https://github.com/smoia/phys2cvr/pull/161) ([@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot]))

#### Authors: 3

- [@dependabot[bot]](https://github.com/dependabot[bot])
- [@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot])
- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.29.0 (Tue Jun 30 2026)

#### 🚀 Enhancement

- Add the option to use an asymmetrical lag range [#155](https://github.com/smoia/phys2cvr/pull/155) ([@beccaclements99](https://github.com/beccaclements99) [@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot]) [@smoia](https://github.com/smoia))

#### ⚠️ Tests

- Bump actions/checkout from 6 to 7 [#159](https://github.com/smoia/phys2cvr/pull/159) ([@dependabot[bot]](https://github.com/dependabot[bot]) [@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot]))

#### 🏠 Internal

- [pre-commit.ci] pre-commit autoupdate [#158](https://github.com/smoia/phys2cvr/pull/158) ([@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot]))

#### Authors: 4

- [@dependabot[bot]](https://github.com/dependabot[bot])
- [@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot])
- Becca Clements ([@beccaclements99](https://github.com/beccaclements99))
- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.28.4 (Fri Mar 13 2026)

#### 🐛 Bug Fix

- fix: Fix unnecessary skip reintroduced in merging PR #154 ([@smoia](https://github.com/smoia))

#### 📝 Documentation

- Update authorship and licencing [#157](https://github.com/smoia/phys2cvr/pull/157) ([@smoia](https://github.com/smoia))

#### Authors: 1

- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.28.3 (Thu Mar 12 2026)

:tada: This release contains work from a new contributor! :tada:

Thank you, Becca Clements ([@beccaclements99](https://github.com/beccaclements99)), for all your work!

#### 🐛 Bug Fix

- Fix usage of previously computed lag maps to run a (new) lagged regression [#154](https://github.com/smoia/phys2cvr/pull/154) ([@beccaclements99](https://github.com/beccaclements99) [@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot]) [@smoia](https://github.com/smoia))
- Give an error message if peaks are not provided AND comp_endtidal is True [#156](https://github.com/smoia/phys2cvr/pull/156) ([@beccaclements99](https://github.com/beccaclements99))

#### ⚠️ Pushed to `master`

- tests: Reorganise tests, remove or mark xfails the failing ones ([@smoia](https://github.com/smoia))
- docs: update licenses ([@smoia](https://github.com/smoia))
- docs: Fix logos ([@smoia](https://github.com/smoia))

#### ⚠️ Tests

- Fix workflow halting on missing peaks when they are not required (i.e. when skipping end-tidal convolution) [#153](https://github.com/smoia/phys2cvr/pull/153) ([@smoia](https://github.com/smoia) [@beccaclements99](https://github.com/beccaclements99) [@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot]))

#### 🏠 Internal

- [pre-commit.ci] pre-commit autoupdate [#152](https://github.com/smoia/phys2cvr/pull/152) ([@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot]))
- Bump actions/checkout from 5 to 6 [#151](https://github.com/smoia/phys2cvr/pull/151) ([@dependabot[bot]](https://github.com/dependabot[bot]))

#### Authors: 4

- [@dependabot[bot]](https://github.com/dependabot[bot])
- [@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot])
- Becca Clements ([@beccaclements99](https://github.com/beccaclements99))
- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.28.2 (Wed Nov 19 2025)

#### 🐛 Bug Fix

- fix: Fix circular imports ([@smoia](https://github.com/smoia))

#### Authors: 1

- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.28.1 (Wed Nov 19 2025)

#### ⚠️ Pushed to `master`

- Refactor compute_petco2hrf to improve modularity and hopefully fix recursive import. ([@smoia](https://github.com/smoia))
- docs: Fix import ([@smoia](https://github.com/smoia))

#### 📝 Documentation

- Try to fix documentation [#147](https://github.com/smoia/phys2cvr/pull/147) ([@smoia](https://github.com/smoia))

#### Authors: 1

- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.28.0 (Wed Nov 19 2025)

:tada: This release contains work from a new contributor! :tada:

Thank you, razkin ([@RazkinMalen](https://github.com/RazkinMalen)), for all your work!
Thank you, Merel ([@merelvdthiel](https://github.com/merelvdthiel)), for all your work!

#### 💥 Breaking Change during development

- Rename internal API (requires re-installation), allow end-tidal interpolation or convolution to be carried out independently, split utility and plotting function in their own modules, small functional refactors, fix extension function [#146](https://github.com/smoia/phys2cvr/pull/146) ([@smoia](https://github.com/smoia))

#### 📝 Documentation

- Documentation and logos for phys2cvr [#145](https://github.com/smoia/phys2cvr/pull/145) ([@RazkinMalen](https://github.com/RazkinMalen) [@merelvdthiel](https://github.com/merelvdthiel) [@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot]) [@smoia](https://github.com/smoia))

#### Authors: 4

- [@merelvdthiel](https://github.com/merelvdthiel)
- [@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot])
- razkin ([@RazkinMalen](https://github.com/RazkinMalen))
- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.27.0 (Tue Nov 11 2025)

#### 💥 Breaking Change during development

- Compute average of ROI signal before eventually filtering it [#144](https://github.com/smoia/phys2cvr/pull/144) ([@smoia](https://github.com/smoia) [@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot]))

#### Authors: 2

- [@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot])
- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.26.1 (Tue Nov 11 2025)

#### 🐛 Bug Fix

- Better handle response functions [#143](https://github.com/smoia/phys2cvr/pull/143) ([@smoia](https://github.com/smoia))

#### Authors: 1

- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.26.0 (Tue Nov 11 2025)

#### 💥 Breaking Change during development

- Split stats models from regressors creation into different files and increase modularity of phys2cvr [#142](https://github.com/smoia/phys2cvr/pull/142) ([@smoia](https://github.com/smoia))

#### Authors: 1

- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.25.1 (Tue Nov 11 2025)

#### 🐛 Bug Fix

- Few minor fixes (including correctly skipping xcorr in brightspin-clinical call) [#141](https://github.com/smoia/phys2cvr/pull/141) ([@smoia](https://github.com/smoia))

#### Authors: 1

- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.25.0 (Tue Nov 11 2025)

#### 💥 Breaking Change during development

- Split resampling between frequencies and samples and fix possible sampling issues in resampling. [#140](https://github.com/smoia/phys2cvr/pull/140) ([@smoia](https://github.com/smoia))

Due to the change of resampling algorithms, this Breaking Change introduces minor differences compared to the previous version in the results, particularly in some configurations' CVR maps and T-stat maps.
CVR maps diff: max=`4.8*10^-7`, mean=`0` (max=`1.2*10^-7`, mean=`0` using `--emat`).
T-stat maps diff: max=`3.81*10^-6`, mean=`10^-8` (no difference using `--emat`).
No difference in lag maps nor in non-lagged maps.
All differences are in absolute value. Having set max tolerance to `10^-5` and mean tolerance to `10^-6`, I consider the differences being within tolerance.

#### Authors: 1

- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.24.0 (Tue Nov 11 2025)

#### 💥 Breaking Change during development

- Use numpy's 1.20+ sliding windows view to fasten up regressors creation [#139](https://github.com/smoia/phys2cvr/pull/139) ([@smoia](https://github.com/smoia))

#### Authors: 1

- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.23.2 (Mon Nov 10 2025)

#### 🐛 Bug Fix

- Allow internal functions to work with 2D matrices and any of their axes. [#138](https://github.com/smoia/phys2cvr/pull/138) ([@smoia](https://github.com/smoia))

#### Authors: 1

- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.23.1 (Wed Nov 05 2025)

#### 🐛 Bug Fix

- Fix petco2hrf call in CLI [#137](https://github.com/smoia/phys2cvr/pull/137) ([@smoia](https://github.com/smoia))

#### Authors: 1

- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.23.0 (Wed Nov 05 2025)

#### 🚀 Enhancement

- Add possibility to select the response function for the PetCO2 signal convolution or to skip it altogether [#135](https://github.com/smoia/phys2cvr/pull/135) ([@smoia](https://github.com/smoia) [@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot]))

#### Authors: 2

- [@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot])
- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.22.0 (Wed Nov 05 2025)

#### 🚀 Enhancement

- Allow specifying a mode for the convolution of petco2, improve convolution plots, and force co2 trace to be 1D [#134](https://github.com/smoia/phys2cvr/pull/134) ([@smoia](https://github.com/smoia) [@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot]))

#### Authors: 2

- [@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot])
- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.21.0 (Tue Nov 04 2025)

#### 💥 Breaking Change during development

- Refactor cross correlation code for simplification. [#133](https://github.com/smoia/phys2cvr/pull/133) ([@smoia](https://github.com/smoia))

#### Authors: 1

- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.20.3 (Tue Nov 04 2025)

#### 🐛 Bug Fix

- Force all new matrices to be created as zeros to avoid background and number formatting issues. [#132](https://github.com/smoia/phys2cvr/pull/132) ([@smoia](https://github.com/smoia) [@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot]))

#### 🏠 Internal

- Bump actions/labeler from 5 to 6 [#131](https://github.com/smoia/phys2cvr/pull/131) ([@dependabot[bot]](https://github.com/dependabot[bot]))

#### Authors: 3

- [@dependabot[bot]](https://github.com/dependabot[bot])
- [@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot])
- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.20.2 (Sat Nov 01 2025)

#### 🐛 Bug Fix

- fix: Fix broken export ([@smoia](https://github.com/smoia))

#### Authors: 1

- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.20.1 (Thu Oct 30 2025)

#### 🐛 Bug Fix

- fix: add extra installation for matlab files and update readme ([@smoia](https://github.com/smoia))

#### Authors: 1

- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.20.0 (Thu Oct 30 2025)

#### 🚀 Enhancement

- Massively increase I/O filetype support [#130](https://github.com/smoia/phys2cvr/pull/130) ([@smoia](https://github.com/smoia))

#### 🏠 Internal

- Add auto labeler workflow [#129](https://github.com/smoia/phys2cvr/pull/129) ([@smoia](https://github.com/smoia))
- Bump actions/setup-python from 5 to 6 [#127](https://github.com/smoia/phys2cvr/pull/127) ([@dependabot[bot]](https://github.com/dependabot[bot]))
- [pre-commit.ci] pre-commit autoupdate [#125](https://github.com/smoia/phys2cvr/pull/125) ([@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot]) [@smoia](https://github.com/smoia))
- Bump actions/checkout from 4 to 5 [#128](https://github.com/smoia/phys2cvr/pull/128) ([@dependabot[bot]](https://github.com/dependabot[bot]))
- Setup Ruff for linting and add bots for automatic updates [#126](https://github.com/smoia/phys2cvr/pull/126) ([@smoia](https://github.com/smoia))

#### Authors: 3

- [@dependabot[bot]](https://github.com/dependabot[bot])
- [@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot])
- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.19.0 (Mon May 05 2025)

#### 🚀 Enhancement

- Update installation to add python 3.12 compatibility and update a few settings [#120](https://github.com/smoia/phys2cvr/pull/120) ([@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot]) [@smoia](https://github.com/smoia))

#### 🏠 Internal

- Updated versioneer to v0.29 [#123](https://github.com/smoia/phys2cvr/pull/123) ([@kristinazvolanek](https://github.com/kristinazvolanek) [@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot]))

#### Authors: 3

- [@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot])
- Kristina Zvolanek ([@kristinazvolanek](https://github.com/kristinazvolanek))
- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.18.6 (Sat Jun 01 2024)

#### 🐛 Bug Fix

- [pre-commit.ci] pre-commit autoupdate [#119](https://github.com/smoia/phys2cvr/pull/119) ([@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot]) [@smoia](https://github.com/smoia))

#### Authors: 2

- [@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot])
- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.18.5 (Mon May 20 2024)

#### 🐛 Bug Fix

- [Snyk] Security upgrade pillow from 9.5.0 to 10.3.0 [#118](https://github.com/smoia/phys2cvr/pull/118) ([@snyk-bot](https://github.com/snyk-bot) [@smoia](https://github.com/smoia))

#### Authors: 2

- Snyk bot ([@snyk-bot](https://github.com/snyk-bot))
- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.18.4 (Mon May 20 2024)

#### 🐛 Bug Fix

- [Snyk] Security upgrade fonttools from 4.38.0 to 4.43.0 [#115](https://github.com/smoia/phys2cvr/pull/115) ([@snyk-bot](https://github.com/snyk-bot) [@smoia](https://github.com/smoia))

#### 🏠 Internal

- [pre-commit.ci] pre-commit autoupdate [#114](https://github.com/smoia/phys2cvr/pull/114) ([@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot]))

#### Authors: 3

- [@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot])
- Snyk bot ([@snyk-bot](https://github.com/snyk-bot))
- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.18.3 (Tue Sep 05 2023)

#### 🐛 Bug Fix

- fix: Make sure butterworth order is given ([@smoia](https://github.com/smoia))

#### Authors: 1

- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.18.2 (Mon Aug 28 2023)

#### 🐛 Bug Fix

- fix: Do not crawl tests folder when importing modules (s.moia@bcbl.eu)

#### ⚠️ Pushed to `master`

- int: Update actions version (s.moia@bcbl.eu)

#### Authors: 1

- Stefano Moia (s.moia@bcbl.eu)

---

# 0.18.1 (Sun Aug 13 2023)

#### 🐛 Bug Fix

- Fix bandpass Butterworth filter order to get a better filter, add option to change order in CLI [#90](https://github.com/smoia/phys2cvr/pull/90) ([@smoia](https://github.com/smoia))
- [pre-commit.ci] pre-commit autoupdate [#88](https://github.com/smoia/phys2cvr/pull/88) ([@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot]) [@smoia](https://github.com/smoia))

#### ⚠️ Pushed to `master`

- int: codespell happy ([@smoia](https://github.com/smoia))
- Start adding a bunch of tests ([@smoia](https://github.com/smoia))
- Wrong library ;) ([@smoia](https://github.com/smoia))
- int: Add codespell ([@smoia](https://github.com/smoia))
- int: update auto-release environment ([@smoia](https://github.com/smoia))
- int: update python-publish environment ([@smoia](https://github.com/smoia))

#### 🏠 Internal

- [pre-commit.ci] pre-commit autoupdate [#87](https://github.com/smoia/phys2cvr/pull/87) ([@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot]))

#### Authors: 2

- [@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot])
- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.18.1 (Thu May 11 2023)

#### ⚠️ Pushed to `master`

- Start adding a bunch of tests ([@smoia](https://github.com/smoia))
- Wrong library ;) ([@smoia](https://github.com/smoia))
- int: Add codespell ([@smoia](https://github.com/smoia))
- int: update auto-release environment ([@smoia](https://github.com/smoia))
- int: update python-publish environment ([@smoia](https://github.com/smoia))

#### 🏠 Internal

- [pre-commit.ci] pre-commit autoupdate [#87](https://github.com/smoia/phys2cvr/pull/87) ([@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot]))

#### Authors: 2

- [@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot])
- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.18.0 (Wed Mar 08 2023)

#### 💥 Breaking Change during development

- Add option to add orthogonalised regressors in phys2cvr [#84](https://github.com/smoia/phys2cvr/pull/84) ([@smoia](https://github.com/smoia) [@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot]))

#### 🐛 Bug Fix

- Improve `auto` settings to include all labels and plugins [#72](https://github.com/smoia/phys2cvr/pull/72) ([@smoia](https://github.com/smoia))

#### ⚠️ Pushed to `master`

- int: Fix issue templates ([@smoia](https://github.com/smoia))

#### Authors: 2

- [@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot])
- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.17.2 (Sat Feb 25 2023)

#### 🐛 Bug Fix

- [pre-commit.ci] pre-commit autoupdate [#64](https://github.com/smoia/phys2cvr/pull/64) ([@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot]))

#### ⚠️ Pushed to `master`

- docs: Update developer installation ([@smoia](https://github.com/smoia))

#### 🏠 Internal

- [pre-commit.ci] pre-commit autoupdate [#63](https://github.com/smoia/phys2cvr/pull/63) ([@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot]))
- [pre-commit.ci] pre-commit autoupdate [#62](https://github.com/smoia/phys2cvr/pull/62) ([@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot]))
- [pre-commit.ci] pre-commit autoupdate [#61](https://github.com/smoia/phys2cvr/pull/61) ([@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot]))
- [pre-commit.ci] pre-commit autoupdate [#59](https://github.com/smoia/phys2cvr/pull/59) ([@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot]))

#### Authors: 2

- [@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot])
- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.17.1 (Mon Dec 12 2022)

#### 🐛 Bug Fix

- [pre-commit.ci] pre-commit autoupdate [#58](https://github.com/smoia/phys2cvr/pull/58) ([@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot]))

#### Authors: 1

- [@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot])

---

# 0.17.0 (Fri Nov 18 2022)

#### 🚀 Enhancement

- Add colormaps to plot lag maps (and sneak in contributors updates) [#56](https://github.com/smoia/phys2cvr/pull/56) ([@smoia](https://github.com/smoia) [@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot]))

#### Authors: 2

- [@pre-commit-ci[bot]](https://github.com/pre-commit-ci[bot])
- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.16.0 (Mon Oct 10 2022)

#### 🚀 Enhancement

- Add pre-commit and format all code [#54](https://github.com/smoia/phys2cvr/pull/54) ([@smoia](https://github.com/smoia))

#### ⚠️ Pushed to `master`

- docs: Fix badges and latest DOI ([@smoia](https://github.com/smoia))

#### Authors: 1

- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.15.3 (Fri Jul 01 2022)

#### 🐛 Bug Fix

- fix: force mask to be an ndarray of booleans in regression ([@smoia](https://github.com/smoia))

#### Authors: 1

- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.15.2 (Fri Jul 01 2022)

#### 🐛 Bug Fix

- fix: Fix plot of xcorr that was called even when skipping xcorr ([@smoia](https://github.com/smoia))

#### Authors: 1

- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.15.1 (Fri Jun 24 2022)

#### 🐛 Bug Fix

- Fix dimensionality check for functional images [#52](https://github.com/smoia/phys2cvr/pull/52) ([@smoia](https://github.com/smoia))

#### Authors: 1

- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.15.0 (Wed Apr 27 2022)

#### 🚀 Enhancement

- Add the possibility to skip the cross correlation step when creating lagged regressors [#50](https://github.com/smoia/phys2cvr/pull/50) ([@smoia](https://github.com/smoia))

#### Authors: 1

- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.14.1 (Tue Apr 19 2022)

#### 🐛 Bug Fix

- fix: Fix nifti import and mask creation ([@smoia](https://github.com/smoia))

#### Authors: 1

- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.14.0 (Tue Apr 19 2022)

#### 🚀 Enhancement

- Allow confounding matrix and mask to be optional in `stats.regression` [#49](https://github.com/smoia/phys2cvr/pull/49) ([@smoia](https://github.com/smoia))

#### Authors: 1

- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.13.0 (Tue Apr 05 2022)

#### 🚀 Enhancement

- Enhance cross correlation to allow both regressor and reference to be shifted and consider max abs corr [#48](https://github.com/smoia/phys2cvr/pull/48) ([@smoia](https://github.com/smoia))

#### Authors: 1

- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.12.1 (Tue Apr 05 2022)

#### 🐛 Bug Fix

- Let signal be resampled over any axis of the input timeseries (default to first dimension) [#47](https://github.com/smoia/phys2cvr/pull/47) ([@smoia](https://github.com/smoia))

#### ⚠️ Pushed to `master`

- docs: fix readme citations ([@smoia](https://github.com/smoia))
- docs: update README with installation instruction and badges ([@smoia](https://github.com/smoia))

#### Authors: 1

- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.12.0 (Tue Mar 22 2022)

#### 💥 Breaking Change during development

- Change behaviour to create petco2hrf fmri proxy due to the SPC of average signal being more robust than average of SPC [#46](https://github.com/smoia/phys2cvr/pull/46) (s.moia@bcbl.eu)

#### Authors: 1

- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.11.2 (Tue Mar 22 2022)

#### 🐛 Bug Fix

- Fix SPC computation by dividing by the transposed mean (not the matrix itself) [#45](https://github.com/smoia/phys2cvr/pull/45) (s.moia@bcbl.eu)

#### Authors: 1

- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.11.1 (Mon Mar 07 2022)

#### 🐛 Bug Fix

- Fix if statement for skipping convolution [#43](https://github.com/smoia/phys2cvr/pull/43) ([@kristinazvolanek](https://github.com/kristinazvolanek))

#### Authors: 1

- Kristina Zvolanek ([@kristinazvolanek](https://github.com/kristinazvolanek))

---

# 0.11.0 (Thu Mar 03 2022)

#### 💥 Breaking Change during development

- Fix workflow to estimate CVR from fMRI only [#42](https://github.com/smoia/phys2cvr/pull/42) (s.moia@bcbl.eu [@smoia](https://github.com/smoia))

#### Authors: 2

- smoia (s.moia@bcbl.eu)
- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.10.0 (Wed Feb 23 2022)

#### 💥 Breaking Change during development

- Force all empty matrices to be float32 [#41](https://github.com/smoia/phys2cvr/pull/41) (s.moia@bcbl.eu [@smoia](https://github.com/smoia))

#### Authors: 2

- smoia (s.moia@bcbl.eu)
- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.9.0 (Wed Feb 23 2022)

#### 💥 Breaking Change during development

- Remove second bulk shift estimation (doesn't make sense). [#40](https://github.com/smoia/phys2cvr/pull/40) (s.moia@bcbl.eu [@smoia](https://github.com/smoia))

#### Authors: 2

- smoia (s.moia@bcbl.eu)
- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.8.0 (Wed Feb 23 2022)

#### 💥 Breaking Change during development

- Include upper boundary of regression ranges (and add a legacy option for pythonic ranges). [#39](https://github.com/smoia/phys2cvr/pull/39) (s.moia@bcbl.eu [@smoia](https://github.com/smoia))

#### Authors: 2

- smoia (s.moia@bcbl.eu)
- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.7.1 (Sun Feb 06 2022)

#### 🐛 Bug Fix

- Add mailmap to provide better "git shortlog -sn" [#37](https://github.com/smoia/phys2cvr/pull/37) (s.moia@bcbl.eu [@smoia](https://github.com/smoia))

#### Authors: 2

- smoia (s.moia@bcbl.eu)
- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.7.0 (Sun Feb 06 2022)

#### 🚀 Enhancement

- Fix lag regression and improve regression debugging options [#36](https://github.com/smoia/phys2cvr/pull/36) (s.moia@bcbl.eu [@smoia](https://github.com/smoia))

#### Authors: 2

- smoia (s.moia@bcbl.eu)
- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.6.2 (Tue Feb 01 2022)

#### ⚠️ Pushed to `master`

- Merge branch 'master' of github.com:smoia/phys2cvr (s.moia@bcbl.eu)
- Add full import (s.moia@bcbl.eu)

#### Authors: 1

- Stefano Moia (s.moia@bcbl.eu)

---

# 0.6.1 (Sat Jan 29 2022)

:tada: This release contains work from a new contributor! :tada:

Thank you, DeepSource Bot ([@deepsourcebot](https://github.com/deepsourcebot)), for all your work!

#### ⚠️ Pushed to `master`

- Add .deepsource.toml ([@deepsourcebot](https://github.com/deepsourcebot))

#### Authors: 1

- DeepSource Bot ([@deepsourcebot](https://github.com/deepsourcebot))

---

# 0.6.0 (Thu Jan 27 2022)

#### 💥 Breaking Change during development

- Updated regression function [#32](https://github.com/smoia/phys2cvr/pull/32) ([@kristinazvolanek](https://github.com/kristinazvolanek) [@smoia](https://github.com/smoia))

#### Authors: 2

- Kristina Zvolanek ([@kristinazvolanek](https://github.com/kristinazvolanek))
- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.5.1 (Thu Jan 13 2022)

#### 🐛 Bug Fix

- Fix astype call on physiological indexes array [#31](https://github.com/smoia/phys2cvr/pull/31) ([@smoia](https://github.com/smoia))

#### Authors: 1

- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.5.0 (Mon Nov 15 2021)

#### 🚀 Enhancement

- Add output of all regression results for debug [#30](https://github.com/smoia/phys2cvr/pull/30) ([@smoia](https://github.com/smoia))

#### Authors: 1

- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.4.2 (Mon Nov 15 2021)

#### 🐛 Bug Fix

- Fix scale factor (divide instead of multiply) for bulk shift output [#29](https://github.com/smoia/phys2cvr/pull/29) ([@smoia](https://github.com/smoia))

#### Authors: 1

- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.4.1 (Thu Nov 04 2021)

#### 🐛 Bug Fix

- Fix peak array loading for index purposes [#28](https://github.com/smoia/phys2cvr/pull/28) ([@kristinazvolanek](https://github.com/kristinazvolanek))

#### Authors: 1

- Kristina Zvolanek ([@kristinazvolanek](https://github.com/kristinazvolanek))

---

# 0.4.0 (Wed Nov 03 2021)

#### 💥 Breaking Change during development

- If scale factor is declared, *divide* the betas by the scale factor [#27](https://github.com/smoia/phys2cvr/pull/27) ([@smoia](https://github.com/smoia))

#### Authors: 1

- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.3.1 (Wed Nov 03 2021)

#### ⚠️ Pushed to `master`

- Extension check are case insensitive ([@smoia](https://github.com/smoia))

#### Authors: 1

- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.3.0 (Sun Oct 31 2021)

#### 💥 Breaking Change during development

- Allow to provide a mask to exclude voxels from the computations, and a separate ROI to use for Xcorr average signal (and lag correction) [#26](https://github.com/smoia/phys2cvr/pull/26) ([@smoia](https://github.com/smoia))

#### Authors: 1

- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.2.4 (Sun Oct 31 2021)

:tada: This release contains work from a new contributor! :tada:

Thank you, Kristina Zvolanek ([@kristinazvolanek](https://github.com/kristinazvolanek)), for all your work!

#### 🐛 Bug Fix

- Fixed output regressor to be the optimally shifted version, rather than the original demeaned regressor [#25](https://github.com/smoia/phys2cvr/pull/25) ([@kristinazvolanek](https://github.com/kristinazvolanek))

#### Authors: 1

- Kristina Zvolanek ([@kristinazvolanek](https://github.com/kristinazvolanek))

---

# 0.2.3 (Sat Oct 30 2021)

#### 🐛 Bug Fix

- Fix denoise_matrix input from CLI that was causing it to be a list of lists [#24](https://github.com/smoia/phys2cvr/pull/24) ([@smoia](https://github.com/smoia))

#### Authors: 1

- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.2.2 (Sat Oct 30 2021)

#### ⚠️ Pushed to `master`

- Fix check_ext call in generating output name ([@smoia](https://github.com/smoia))

#### Authors: 1

- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.2.1 (Sat Oct 30 2021)

#### 🐛 Bug Fix

- Fix outputs paths and names [#23](https://github.com/smoia/phys2cvr/pull/23) ([@smoia](https://github.com/smoia))

#### Authors: 1

- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.2.0 (Fri Oct 15 2021)

#### 💥 Breaking Change during development

- Improve parser options [#22](https://github.com/smoia/phys2cvr/pull/22) ([@smoia](https://github.com/smoia))

#### Authors: 1

- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.1.1 (Sun Oct 10 2021)

#### 🐛 Bug Fix

- Fix zenodo metadata [#21](https://github.com/smoia/phys2cvr/pull/21) ([@smoia](https://github.com/smoia))

#### Authors: 1

- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.1.0 (Sun Oct 10 2021)

#### 💥 Breaking Change during development

- Fix `auto` config file [#20](https://github.com/smoia/phys2cvr/pull/20) ([@smoia](https://github.com/smoia))
- Add configuration files for CD, declare licence in all codefiles, update development status [#19](https://github.com/smoia/phys2cvr/pull/19) ([@smoia](https://github.com/smoia))

#### 🐛 Bug Fix

- Add docstrings for `io.py` [#18](https://github.com/smoia/phys2cvr/pull/18) ([@smoia](https://github.com/smoia))

#### Authors: 1

- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.0.3 (Sun Oct 10 2021)

#### 🐛 Bug Fix

- Improve Exception types [#17](https://github.com/smoia/phys2cvr/pull/17) ([@smoia](https://github.com/smoia))

#### Authors: 1

- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.0.2 (Sat Oct 09 2021)

#### 🐛 Bug Fix

- Invert CLI behaviour for lagged regressors, improve workflow input parameters definition and check [#16](https://github.com/smoia/phys2cvr/pull/16) ([@smoia](https://github.com/smoia))

#### Authors: 1

- Stefano Moia ([@smoia](https://github.com/smoia))

---

# 0.0.1 (Sat Oct 09 2021)

#### 🐛 Bug Fix

- Clean attributes in subfiles [#15](https://github.com/smoia/phys2cvr/pull/15) ([@smoia](https://github.com/smoia))

#### ⚠️ Pushed to `master`

- Add workflows and templates ([@smoia](https://github.com/smoia))
- Parser options all there. ([@smoia](https://github.com/smoia))
- lag_step becomes 0.3 ([@smoia](https://github.com/smoia))
- Rename lag variables ([@smoia](https://github.com/smoia))
- Rename maxlag for coherence ([@smoia](https://github.com/smoia))
- Force run_regression to be true if lagged_regression is! ([@smoia](https://github.com/smoia))
- Remove unused parameter ([@smoia](https://github.com/smoia))
- Start working that parser ([@smoia](https://github.com/smoia))
- Fix relative lag background ([@smoia](https://github.com/smoia))
- Minor text correction and transform comments into issues ([@smoia](https://github.com/smoia))
- Signal percentage change the data outside of regression ([@smoia](https://github.com/smoia))
- When co2 signal is not provided, force skip convolution ([@smoia](https://github.com/smoia))
- Add relative lag export ([@smoia](https://github.com/smoia))
- Reduce butter order to 1 ([@smoia](https://github.com/smoia))
- Improve empty array creation, change lag variable names, change take command due to heavy memory load ([@smoia](https://github.com/smoia))
- Fix f string ([@smoia](https://github.com/smoia))
- Remove breakpoints ([@smoia](https://github.com/smoia))
- Fix beta and tstat assignment ([@smoia](https://github.com/smoia))
- Invert regr shift call and fix assignment of volumes after regression ([@smoia](https://github.com/smoia))
- Close those plot ([@smoia](https://github.com/smoia))
- Skip hrf plot ([@smoia](https://github.com/smoia))
- Better try than have the doubt ([@smoia](https://github.com/smoia))
- Fine shifts exist only if they exist ([@smoia](https://github.com/smoia))
- Fix lagged regressor matrices shape ([@smoia](https://github.com/smoia))
- Fix step ([@smoia](https://github.com/smoia))
- Move output to right folder ([@smoia](https://github.com/smoia))
- stats should be debugged. ([@smoia](https://github.com/smoia))
- Correct regressors hstack ([@smoia](https://github.com/smoia))
- Also say the right matrix! ([@smoia](https://github.com/smoia))
- read denoising matrix, not its name ([@smoia](https://github.com/smoia))
- Invert petco2hrf shifts matrix orientation ([@smoia](https://github.com/smoia))
- Remember to add a plot of the xcorr with both xcorrs if two were used ([@smoia](https://github.com/smoia))
- Shorter tp assignment ([@smoia](https://github.com/smoia))
- join comes from paths... ([@smoia](https://github.com/smoia))
- Fix func import ([@smoia](https://github.com/smoia))
- Fix mask import and average of functional signal ([@smoia](https://github.com/smoia))
- Boolean input variables as booleans rather than empty ([@smoia](https://github.com/smoia))
- Fix general pre-installation stuff ([@smoia](https://github.com/smoia))
- Fix authors ([@smoia](https://github.com/smoia))
- DEBUG ([@smoia](https://github.com/smoia))
- Add multiple denoise matrices ([@smoia](https://github.com/smoia))
- If user specifies lag map, use it! ([@smoia](https://github.com/smoia))
- stats.regression returns all lagged regressors in memory ([@smoia](https://github.com/smoia))
- Update comments ([@smoia](https://github.com/smoia))
- If no outdir is specified, make one in the functional input folder ([@smoia](https://github.com/smoia))
- Log to log folder (skip code folder) ([@smoia](https://github.com/smoia))
- Simplify input ([@smoia](https://github.com/smoia))
- Update inputs accordingly ([@smoia](https://github.com/smoia))
- Check dimension of confounders input ([@smoia](https://github.com/smoia))
- Add R square computation ([@smoia](https://github.com/smoia))
- Add possibility to skip convolution (e.g. if input is already petco2hrf) ([@smoia](https://github.com/smoia))
- Add legendre and OLS regression ([@smoia](https://github.com/smoia))
- Move everything related to statistical models to stats module ([@smoia](https://github.com/smoia))
- Move signal related functions to signal module ([@smoia](https://github.com/smoia))
- Move io function to their own module ([@smoia](https://github.com/smoia))
- Remove extra main ([@smoia](https://github.com/smoia))
- Remove nilearn and fix imports ([@smoia](https://github.com/smoia))
- Start adding lagged regression - but it's totally untested. ([@smoia](https://github.com/smoia))
- Update necessary libraries and import OLSModel ([@smoia](https://github.com/smoia))
- Read phys files from peakdet! ([@smoia](https://github.com/smoia))
- Remove save call from main workflow ([@smoia](https://github.com/smoia))
- Less variables, more understanding ([@smoia](https://github.com/smoia))
- Add support for gzipped tsv ([@smoia](https://github.com/smoia))
- Load phys files ([@smoia](https://github.com/smoia))
- Fix check extension ([@smoia](https://github.com/smoia))
- Continue file import, add TR reading, Add filter application ([@smoia](https://github.com/smoia))
- Start working on get_regressor ([@smoia](https://github.com/smoia))
- Add a filter signal function ([@smoia](https://github.com/smoia))
- Save bash call as a function not in main workflow ([@smoia](https://github.com/smoia))
- Change scipy import ([@smoia](https://github.com/smoia))
- Change author ([@smoia](https://github.com/smoia))
- Add most of file import ([@smoia](https://github.com/smoia))
- Clear up package ([@smoia](https://github.com/smoia))
- Update contributors ([@smoia](https://github.com/smoia))
- Add physiopy logo ([@smoia](https://github.com/smoia))
- Add parser and main workflow ([@smoia](https://github.com/smoia))
- Correct development status ([@smoia](https://github.com/smoia))
- Add infra files ([@smoia](https://github.com/smoia))
- Initial commit ([@smoia](https://github.com/smoia))

#### 📝 Documentation

- Add docstrings [#14](https://github.com/smoia/phys2cvr/pull/14) ([@smoia](https://github.com/smoia))

#### Authors: 1

- Stefano Moia ([@smoia](https://github.com/smoia))
