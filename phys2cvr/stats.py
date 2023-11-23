#!/usr/bin/env python3
"""
Statistical module for phys2cvr.

Attributes
----------
LGR :
    Logger
"""

import logging
import os

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view as swv
from scipy.stats import zscore

from phys2cvr import io

R2MODEL = ["full", "partial", "intercept", "adj_full", "adj_partial", "adj_intercept"]

LGR = logging.getLogger(__name__)
LGR.setLevel(logging.INFO)


def x_corr(func, co2, n_shifts=None, offset=0, abs_xcorr=False):
    """
    Calculates the cross correlation between `func` and `co2`.

    Parameters
    ----------
    func : np.ndarray
        Timeseries of functional data, must be SHORTER than `co2` (or of equal length)
    co2 : np.ndarray
        Timeseries of CO2 (or physiological) data, can be LONGER than `func`
    n_shifts : int or None, optional
        Number of maximal timepoints to shift when cross-correlating func and co2.
        When None (default), consider all possible shifts.
        Each shift consists of one shifted timepoint (one step).
    offset : int, optional
        This is optional. The desired timepoint at which `func` offsets, i.e., the number of timepoints of `co2`
        to exclude from the cross correlation.
    abs_xcorr : bool, optional
        If True, x_corr will find the strongest absolute correlation,
        i.e., max(|corr(func, co2)|), which could be the strongest negative correlation rather than the strongest positive correlation.

    Returns
    -------
    float :
        Strongest cross correlation
    int :
        Index of strongest cross correlation
    xcorr : np.ndarray
        Full array containing all computed cross correlations between `func` and `co2`

    Raises
    ------
    ValueError
        If `offset` is higher than the difference between the length of `co2` and `func`.
    NotImplementedError
        If `offset` < 0
        If `co2` length is smaller than `func` length.
    """
    if offset < 0:
        raise NotImplementedError("Negative offsets are not supported yet.")

    if func.shape[0] + offset > co2.shape[0]:
        if offset > 0:
            raise ValueError(
                f"The specified offset of {offset} is too high."
                f"func of length {func.shape[0]} can not be compared with co2 of "
                f"length {co2.shape[0]}"
            )
        else:
            raise NotImplementedError(
                f"The timeseries has length of {func.shape[0]}, which is longer than the"
                f"length of the given co2 regressor ({co2.shape[0]}). This "
                "is not supported."
            )

    if n_shifts is None:
        n_shifts = co2.shape[0] - (func.shape[0] + offset) + 1
        LGR.info(
            f"Considering all possible shifts of regressor for Xcorr, i.e. {n_shifts}"
        )
    else:
        if n_shifts + offset + func.shape[0] > co2.shape[0]:
            LGR.warning(
                f"The specified number of shifts ({n_shifts}) is too high for the "
                f"length of the regressor ({co2.shape[0]})."
            )
            n_shifts = co2.shape[0] - (func.shape[0] + offset) + 1
            LGR.warning(f"Considering {n_shifts} shifts instead.")

    sco2 = swv(co2, func.shape[0], axis=-1)[offset : n_shifts + offset]

    xcorr = np.dot(zscore(sco2, axis=-1), zscore(func)) / func.shape[0]

    if abs_xcorr:
        return np.abs(xcorr).max(), np.abs(xcorr).argmax() + offset, xcorr
    else:
        return xcorr.max(), xcorr.argmax() + offset, xcorr


def ols(Ymat, Xmat, r2model="full", residuals=False, demean=False):
    """
    Implement Ordinary Least Squares linear regression.

    For both Ymat and Xmat, axis 0 must reflect the time axis.
    This is the barebone OLS implementation. For the full regression step,
    see `stats.regression`.

    Parameters
    ----------
    Ymat : np.ndarray
        Dependent variable, 1D or 2D.
        The columns must represent the variable over time, i.e., time axis must be axis 0.
    Xmat : np.ndarray
        Independent variables, 1D or 2D. The regressor of interest MUST be the last dimension.
        The columns must represent the regressors over time, i.e., time axis must be axis 0.
    r2model : {'full', 'partial', 'intercept', 'adj_full', 'adj_partial', 'adj_intercept'}, optional
        R^2 to report.
        Possible models are:
            - 'full' (default)
                Use all regressors in the model, i.e., compare all regressors versus a baseline of 0
            - 'partial'
                Consider only the first vector of Xmat as regressor of interest in the model, i.e. compare the first vector with
                a baseline which is composed of all other vectors of Xmat beside the first.
            - 'intercept'
                Use all regressors in the model, except for the intercept, i.e., compare all regressors versus a baseline which is the
                intercept (Legendre polynomial order 0, a.k.a.,
                average signal)
            - 'adj_*'
                Adjusted R^2 version of 'full'
        Under normal conditions, although the R^2 values will vary between the different options,
        a lagged regression based on any R^2 model will provide the same results independent of the chosen option.
        This WILL NOT be the case if orthogonalisations are applied between the first vector of Xmat
        and the others vectors. In this case, a lagged regression based on `partial` might hold
        different results from the others.
        Default: 'full'
    residuals : bool, optional
        If True, outputs only the residuals of the model - to be mainly used for orthogonalisation
        If False, outputs betas, tstats, and R^2 (default).
    demean : bool, optional
        If True, demean Xmat before running OLS.
        Default is False.

    Returns
    -------
    betas : np.ndarray
        Beta values
    t_stats : np.ndarray
        T-statistic values
    r_square : np.ndarray
        R^2 values

    Raises
    ------
    NotImplementedError
        If Ymat has more than 2 dimensions
        If Xmat has more than 2 dimensions
        If "poly" R^2 is declared (as it's not implemented yet).
    ValueError
        If a non-valid R^2 value is declared
    """
    if Ymat.ndim > 2:
        raise NotImplementedError(
            "OLS on data with more than 2 dimensions is not implemented yet."
        )
    elif Ymat.ndim < 2:
        Ymat = Ymat[..., np.newaxis]
    if Xmat.ndim > 2:
        raise NotImplementedError(
            "OLS on regressors with more than 2 dimensions is not implemented yet."
        )
    elif Xmat.ndim < 2:
        Xmat = Xmat[..., np.newaxis]

    if demean:
        LGR.info("Demean regressors")
        Xmat = Xmat - Xmat.mean(axis=0)

    try:
        betas, RSS, _, _ = np.linalg.lstsq(Xmat, Ymat, rcond=None)
        if not RSS.any():
            RSS = np.zeros(betas.shape[1])

    except np.linalg.LinAlgError:
        raise ValueError(
            "The given matrices might not be oriented correctly. Try to first transpose the "
            "regressor matrix."
        )

    if residuals:
        return Ymat - (Xmat @ betas)
    else:
        # compute t-values of betas (estimates)
        # first compute number of degrees of freedom
        df = Xmat.shape[0] - Xmat.shape[1]

        # compute sigma:
        # RSS = sum{[mdata - (X * betas)]^2}
        # sigma = RSS / Degrees_of_Freedom
        # RSS = np.sum(np.power(Ymat.T - np.dot(Xmat, betas), 2), axis=0)
        sigma = RSS / df
        sigma = sigma[..., np.newaxis]

        # Copmute std of betas:
        # C = (mat^T * mat)_ii^(-1)
        # std(betas) = sqrt(sigma * C)
        C = np.diag(np.linalg.pinv(np.dot(Xmat.T, Xmat)))
        C = C[..., np.newaxis]
        std_betas = np.sqrt(np.outer(C, sigma))
        tstats = betas / std_betas
        tstats[np.isneginf(tstats)] = -9999
        tstats[np.isposinf(tstats)] = 9999

        # Compute R^2 coefficient of multiple determination!
        r2model_support = False
        for model in R2MODEL:
            if r2model == model:
                r2model_support = True
        if not r2model_support:
            raise ValueError(
                f"{r2model} R^2 not supported. Supported R^2 models are {R2MODEL}"
            )

        r2msg = ""
        # R^2 = 1 - RSS/TSS, (TSS = total sum of square)
        if "full" in r2model:
            # Baseline model is 0 - this is what we're using ATM
            TSS = np.sum(np.power(Ymat, 2), axis=0)
        elif "intercept" in r2model:
            # Baseline model is intercept: TSS is variance*samples
            TSS = Ymat.var(axis=0) * Ymat.shape[0]
        elif "poly" in r2model:
            # Baseline model is legendre polynomials - or others: TSS is RSS of partial matrix
            # Could improve efficiency by moving this fitting step outside the regression loop
            # polynomials = Xmat[:, 0:4]
            # TSS = np.linalg.lstsq(polynomials, Ymat.T, rcond=None)[1]
            raise NotImplementedError("'poly' R^2 not implemented yet")
            r2msg = "polynomial"
        elif "partial" in r2model:
            pass

        if "partial" in r2model:
            # We could also compute PARTIAL R square of regr instead (or on top)
            # See for computation: https://sscc.nimh.nih.gov/sscc/gangc/tr.html
            tstats_square = np.power(tstats, 2)
            r_square = (tstats_square / (tstats_square + df))[-1, :]
        else:
            r_square = np.ones(Ymat.shape[1], dtype="float32") - (RSS / TSS)

        if r2msg == "":
            r2msg = r2model
        if "adj_" in r2model:
            # We could compute ADJUSTED R^2 instead
            r_square = 1 - ((1 - r_square) * (Xmat.shape[0] - 1) / (df - 1))
            r2msg = f"adjusted {r2msg}"

        LGR.info(f"Adopting {r2msg} baseline to compute R^2.")

        return betas, tstats, r_square


def regression(
    data,
    regr,
    denoise_mat=None,
    ortho_mat=None,
    extra_mat=None,
    mask=None,
    r2model="full",
    debug=False,
    x1D="mat.1D",
):
    """
    Estimate regression parameters.

    Parameters
    ----------
    data : np.ndarray
        Dependent variable of the model (i.e., Y), as a 4D volume.
    regr : np.ndarray
        Regressor of interest
    denoise_mat : np.ndarray or None, optional
        Confounding effects (regressors)
    ortho_mat : np.ndarray or None, optional
        Confounding effects (regressors) which will be orthogonalised with respect to `regr`, `denoise_mat`, and
        `extra_mat`
    extra_mat : np.ndarray or None, optional
        Extra factors (regressors) that will be used to orthogonalise `ortho_mat`
    mask : np.ndarray or None, optional
        A 3D mask to reduce the number of voxels to run the regression on.
    r2model : {'full', 'partial', 'intercept', 'adj_full', 'adj_partial', 'adj_intercept'}, optional
        R^2 to report.
        Possible models are:
            - 'full' (default)
                Use all regressors in the model, i.e., compare all regressors versus a baseline of 0
            - 'partial'
                Consider only `regr` in the model, i.e., compare `regr` with a baseline which is
                composed of all other regressors (`denoise_mat`)
            - 'intercept'
                Use all regressors in the model, except for the intercept, i.e., compare all regressors versus a baseline which is the
                intercept (Legendre polynomial order 0, a.k.a.
                average signal)
            - 'adj_*'
                Adjusted R^2 version of 'full'
        Under normal conditions, although the R^2 value will vary between the different options,
        a lagged regression based on any R^2 model will provide the same results independent of the chosen option.
        This WILL NOT be the case if orthogonalisations are applied  between `regr` and `denoise_mat`. In that case,
        a lagged regression based on `partial` might hold different results from the others.
        Default: 'full'

    Returns
    -------
    bout : np.ndarray
        Beta map
    tout : np.ndarray
        T-statics map
    rout : np.ndarray
        R^2 map

    Raises
    ------
    ValueError
        If `denoise_mat` and `regr` do not have at least one common dimension.
    """
    # Obtain data in 2d
    if mask is None:
        mask = np.ones(data.shape[:-1])

    mask = mask.astype("bool")

    Ymat = data[mask]
    # Check that regr has "two" dimensions
    regr = io.array_is_2d(regr)

    if denoise_mat is not None:
        if regr.shape[0] != denoise_mat.shape[0]:
            if regr.shape[0] != denoise_mat.shape[1]:
                raise ValueError(
                    "The provided confounding matrix does not match "
                    "the dimensionality of the PetCO2hrf regressor!"
                )
            else:
                denoise_mat = denoise_mat.T
        # Stack mat
        # Note: Xmat is not currently demeaned within this function, so inputs
        # should already be demeaned
        Xmat = np.hstack([denoise_mat, regr])
        if ortho_mat is not None:
            if ortho_mat.shape[0] != denoise_mat.shape[0]:
                ortho_mat = ortho_mat.T
                if ortho_mat.shape[0] != denoise_mat.shape[0]:
                    raise ValueError(
                        "The provided orthogonalised matrix does not match "
                        "the dimensionality of the PetCO2hrf regressor!"
                    )

            nuisance_mat = Xmat.copy()
            if extra_mat is not None:
                if extra_mat.shape[0] != denoise_mat.shape[0]:
                    extra_mat = extra_mat.T
                    if extra_mat.shape[0] != denoise_mat.shape[0]:
                        raise ValueError(
                            "The provided extra matrix does not match "
                            "the dimensionality of the PetCO2hrf regressor!"
                        )
                nuisance_mat = np.hstack([nuisance_mat, extra_mat])

            ortho_mat = ols(ortho_mat, nuisance_mat, residuals=True)
            Xmat = np.hstack([ortho_mat, Xmat])
    else:
        Xmat = regr

    if debug:
        os.makedirs(os.path.dirname(x1D), exist_ok=True)
        np.savetxt(x1D, Xmat, fmt="%.6f")

    betas, tstats, r_square = ols(
        Ymat.T, Xmat, r2model="full", residuals=False, demean=False
    )

    if debug:
        # debug should export betas, tstats, r_square
        pass

    # Assign betas, Rsquare and tstats to new volume
    bout = np.zeros(mask.shape)
    tout = np.zeros(mask.shape)
    rout = np.zeros(mask.shape)
    bout[mask] = betas[-1, :]
    tout[mask] = tstats[-1, :]
    rout[mask] = r_square
    return bout, tout, rout


"""
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
"""
