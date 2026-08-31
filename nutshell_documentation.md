# Nutshell explanation of the main function of **phys2cvr**.

**phys2cvr** is a python tool that allows to estimate cerebrovascular reactivity (CVR) in the brain by assuming a linear relationship between BOLD fMRI signal and the level of C02 in the body.

## The program phys2cvr take as **inputs**:

- BOLD fMRI signal
- optionally C02 file (physiological measure)
- optionally you can input:
  - mask → domain where the cvr map is computed
  - roi file → used to compute the average signal that will be used as a regressor when C02 is not available.
  - additional regressors that are computed using e.g. ICA or motion parameters
- Specify the lag range that is the expected physiological delay between the C02 signal and the BOLD fmri signal

## **Core computation:**

The core algorithm performs a voxel-wise general linear model (GLM) between the fMRI signal and either:
- The processed CO₂ trace, or
- The ROI-based surrogate regressor (if no CO₂ is provided).

The CO₂ trace is internally preprocessed including:
- Peak detection
- Demeaning
- Band-pass filtering
- Convolution with a hemodynamic response function (HRF)
- Optional resampling to match the fMRI time series

## The **main output** of the code are two maps (in nifti format):
- CVR map as 3D nifti
- Lag map as 3D nifti
