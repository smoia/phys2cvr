#!/usr/bin/env python3
"""Parallelisation runs for phys2cvr."""

import os

import numpy as np

from phys2cvr import stats


def _l_glm_lagmap(
    i,
    step,
    regr_shifts,
    outdir,
    func,
    denoise_matrix,
    orthogonalised_matrix,
    extra_matrix,
    mask,
    r2model,
    debug,
    lag_idx,
):
    """Worker function for running GLM per unique lag index."""
    regr = regr_shifts[(i * step), :, np.newaxis]
    x1D = os.path.join(outdir, 'mat', f'mat_{i:04g}.1D')
    idx_mask = lag_idx == i

    b, t, _ = stats.regression(
        func[idx_mask],
        regr,
        denoise_matrix,
        orthogonalised_matrix,
        extra_matrix,
        mask[idx_mask],
        r2model,
        debug,
        x1D,
    )
    return idx_mask, b, t


def _l_glm_range(
    n,
    i,
    regr_shifts,
    outdir,
    func,
    denoise_matrix,
    orthogonalised_matrix,
    extra_matrix,
    mask,
    r2model,
    debug,
):
    """Worker function for running GLM across full volume per lag range iteration."""
    regr = regr_shifts[i, :, np.newaxis]
    x1D = os.path.join(outdir, 'mat', f'mat_{i:04g}.1D')

    b, t, r2 = stats.regression(
        func,
        regr,
        denoise_matrix,
        orthogonalised_matrix,
        extra_matrix,
        mask,
        r2model,
        debug,
        x1D,
    )
    return n, b, t, r2


"""
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
"""
