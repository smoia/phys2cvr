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
