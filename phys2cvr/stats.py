#!/usr/bin/env python3

import logging
import os

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as spint
import scipy.stats as sct

from phys2cvr import io


SET_DPI = 100
FIGSIZE = (18, 10)
LGR = logging.getLogger(__name__)
LGR.setLevel(logging.INFO)
EXT_1D = ['.txt', '.csv', '.tsv', '.1d', '.par', '.tsv.gz']
EXT_NIFTI = ['.nii', '.nii.gz']


def x_corr(func, co2, lastrep, firstrep=0, offset=0):
    if len(func) + offset > len(co2):
        raise ValueError(f'The specified offset of {offset} is too high to '
                         f'compare func of length {len(func)} with co2 of '
                         f'length {len(co2)}')
    if firstrep + offset < 0:
        firstrep = -offset
    if lastrep + offset + len(func) > len(co2):
        lastrep = len(co2) - offset - len(func)

    xcorr = np.empty(lastrep+firstrep)
    for i in range(firstrep, lastrep):
        xcorr[i] = np.corrcoef(func, co2[0+i+offset:len(func)+i+offset].T)[1, 0]

    return xcorr.max(), (xcorr.argmax() + firstrep + offset), xcorr


def get_regr(func_avg, petco2hrf, tr, freq, outname, maxlag=9, trial_len='',
             n_trials='', no_pad=False, ext='.1D', lagged_regression=True):
    # Setting up some variables
    first_tp = 0
    last_tp = -1

    if trial_len and n_trials:
        # If both are specified, disregard two extreme _trial from matching.
        LGR.info(f'Specified {n_trials} trials lasting {trial_len} seconds')
        if n_trials > 2:
            LGR.info('Ignoring first trial to improve first bulk shift estimation')
            first_tp = int(trial_len*freq)
        else:
            LGR.info('Using all trials for bulk shift estimation')
        if n_trials > 3:
            LGR.info('Ignoring last trial to improve first bulk shift estimation')
            last_tp = first_tp*(n_trials-1)

    elif trial_len and not n_trials:
        LGR.warning('The length of trial was specified, but the number of '
                    'trials was not. Using all trials for bulk shift estimation')
    elif not trial_len and n_trials:
        LGR.warning('The number of trials was specified, but the length of '
                    'trial was not. Using all trials for bulk shift estimation')
    else:
        LGR.info('Using all trials for bulk shift estimation.')

    # Upsample functional signal
    func_len = len(func_avg)
    regr_x = np.arange(0, ((func_len-1) * tr + 1/freq), 1/freq)
    func_x = np.linspace(0, (func_len - 1) * tr, func_len)
    f = spint.interp1d(func_x, func_avg, fill_value='extrapolate')
    func_upsampled = f(regr_x)
    len_upd = len(func_upsampled)

    # Preparing breathhold and CO2 trace for Xcorr
    func_cut = func_upsampled[first_tp:last_tp]
    petco2hrf_cut = petco2hrf[first_tp:]

    nrep = len(petco2hrf_cut) - len(func_cut)

    _, optshift, xcorr = x_corr(func_cut, petco2hrf, nrep)
    LGR.info(f'First cross correlation estimated bulk shift at {optshift/freq} seconds')

    if trial_len and n_trials and n_trials > 2:
        LGR.info('Running second bulk shift estimation')
        if len(func_upsampled) + nrep > len(petco2hrf):
            pad = len(func_upsampled) + nrep - len(petco2hrf)
            petco2hrf_padded = np.pad(petco2hrf, pad, mode='mean')
        else:
            petco2hrf_padded = petco2hrf

        _, optshift, xcorr = x_corr(func_upsampled, petco2hrf_padded, nrep, -nrep, optshift)
        LGR.info(f'Second cross correlation estimated bulk shift at {optshift/freq} seconds')

    # Export estimated optimal shift in seconds
    with open(f'{outname}_optshift.1D', 'w') as f:
        print(f'{(optshift/freq):.4f}', file=f)

    petco2hrf_shift = petco2hrf[optshift:optshift+len_upd]

    # preparing for and exporting figures of shift
    time_axis = np.arange(0, nrep/freq, 1/freq)

    if nrep < len(time_axis):
        time_axis = time_axis[:nrep]
    elif nrep > len(time_axis):
        time_axis = np.pad(time_axis, (0, int(nrep - len(time_axis))), 'linear_ramp')

    plt.figure(figsize=FIGSIZE, dpi=SET_DPI)
    plt.plot(time_axis, xcorr)
    plt.title('optshift')
    plt.savefig(f'{outname}_optshift.png', dpi=SET_DPI)

    plt.figure(figsize=FIGSIZE, dpi=SET_DPI)
    plt.plot(sct.zscore(petco2hrf_shift), '-', sct.zscore(func_upsampled), '-')
    plt.title('GM and shift')
    plt.savefig(f'{outname}_petco2hrf.png', dpi=SET_DPI)

    petco2hrf_demean = io.export_regressor(regr_x, petco2hrf_shift, func_x, outname, 'petco2hrf', ext)

    if lagged_regression:
        outprefix = os.path.join(os.path.split(outname)[0], 'regr', os.path.split(outname)[1])
        os.makedirs(os.path.join(os.path.split(outname)[0], 'regr'), exist_ok=True)

        # Set num of fine shifts: 9 seconds is a bit more than physiologically feasible
        nrep = int(maxlag * freq)

        # Padding regressor for shift, and padding optshift too
        if (optshift - nrep) < 0:
            lpad = nrep - optshift
        else:
            lpad = 0

        if (optshift + nrep + len_upd) > len(petco2hrf):
            rpad = (optshift + nrep + len_upd) - len(petco2hrf)
        else:
            rpad = 0

        petco2hrf_padded = np.pad(petco2hrf, (int(lpad), int(rpad)), 'mean')

        for i in range(-nrep, nrep):
            petco2hrf_shift = petco2hrf_padded[optshift+lpad-i:optshift+lpad-i+len_upd]
            io.export_regressor(regr_x, petco2hrf_shift, func_x, outprefix, f'_{(i + nrep):04g}', ext)

    return petco2hrf_demean
