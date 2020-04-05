#!/usr/bin/env python3

import os

import numpy as np
import peakutils as pk

import matplotlib.pyplot as plt
import scipy.signal as sgn
import scipy.interpolate as spint
import scipy.stats as sct

# matplotlib.use("TkAgg")
# matplotlib.interactive(True)
reject_list = np.array([])
SET_DPI = 100
FIGSIZE = (18, 10)


def create_hrf(newfreq):
	"""
	
	"""
    # Create HRF
    RT = 1/newfreq
    fMRI_T = 16
    p = [6, 16, 1, 1, 6, 0, 32]

    # Modelled hemodynamic response function - {mixture of Gammas}
    dt = RT / fMRI_T
    u = np.arange(0, p[6]/dt+1, 1) - p[5]/dt
    a1 = p[0] / p[2]
    b1 = 1 / p[3]
    a2 = p[1] / p[3]
    b2 = 1 / p[3]
    hrf = (sct.gamma.pdf(u*dt, a1, scale=b1) - sct.gamma.pdf(u*dt, a2, scale=b2)/p[4])/dt
    time_axis = np.arange(0, int(p[6]/RT+1), 1) * fMRI_T
    hrf = hrf[time_axis]
    min_hrf = 1e-9*min(hrf[hrf > 10*np.finfo(float).eps])

    if min_hrf < 10*np.finfo(float).eps:
        min_hrf = 10*np.finfo(float).eps

    hrf[hrf == 0] = min_hrf
    hrf = hrf/max(hrf)

    plt.figure()
    plt.plot(time_axis, hrf)

    return hrf


def decimate_data(filename, newfreq=40):
    data = np.genfromtxt(filename + '.tsv.gz')
    idz = (data[:, 0]>=0).argmax()
    data = data[idz:, ]
    f = spint.interp1d(data[:, 0], data[:, ], axis=0, fill_value='extrapolate')
    data_tdec = np.arange(0, data[-1, 0], 1/newfreq)
    data_dec = f(data_tdec)

    del data, idz

    np.savetxt(filename + '_dec.tsv.gz', data_dec)

    return data_dec


def filter_signal(data_dec, channel=4):
    ba = sgn.butter(7, 2/20, 'lowpass')
    resp_filt = sgn.filtfilt(ba[0], ba[1], data_dec[:, channel])
    plt.figure()
    plt.plot(resp_filt)

    return resp_filt


def get_peaks(co):
    # Finding peaks
    pidx = pk.peak.indexes(co, thres=0.5, min_dist=120).tolist()
    plt.figure()
    plt.plot(co)
    plt.plot(co, 'm*', markevery=pidx)
    # pidx = np.delete(pidx,[56])
    return co, pidx


def get_petco2(co, pidx, hrf, filename):
    # Extract PETco2
    co_lenght = len(co)
    nx = np.linspace(0, co_lenght, co_lenght)
    f = spint.interp1d(pidx, co[pidx], fill_value='extrapolate')
    co_peakline = f(nx)

    # Plot PETco2
    plt.figure(figsize=FIGSIZE, dpi=SET_DPI)
    plt.title('CO2 and regressor')
    plt.plot(co, '-', co_peakline, '-')
    plt.savefig(filename + '_co_peakline.png', dpi=SET_DPI)
    plt.close()

    # Demean and export
    co_peakline = co_peakline - co_peakline.mean()
    np.savetxt(filename + '_co_peakline.1D',co_peakline,fmt='%.18f')

    # Convolve, and then rescale to have same amplitude (?)
    co_conv = np.convolve(co_peakline, hrf)
    co_conv = np.interp(co_conv, (co_conv.min(), co_conv.max()), (co_peakline.min(), co_peakline.max()))

    # co_conv = co_conv[ign_tr:co_lenght]
    plt.figure(figsize=FIGSIZE, dpi=SET_DPI)
    plt.title('regressor and convolved regressor')
    plt.plot(co_conv, '-', co_peakline, '-')
    plt.savefig(filename + '_co_regressors.png', dpi=SET_DPI)
    plt.close()

    np.savetxt(filename + '_co_conv.1D', co_conv, fmt='%.18f')

    return co_conv


def export_regressor(regr_x, co_shift, GM_x, GM_name, suffix='_co_regr'):
    f = spint.interp1d(regr_x, co_shift, fill_value='extrapolate')
    co_tr = f(GM_x)
    co_tr = co_tr - co_tr.mean()
    textname = GM_name + suffix + '.1D'
    np.savetxt(textname, co_tr, fmt='%.18f')


def get_regr(GM_name, co_conv, tr=1.5, newfreq=40, BH_len=58, nBH=8, ext='.1D'):
    GM = np.genfromtxt(GM_name + ext)
    sequence_tps = len(GM)

    regr_x = np.arange(0, ((sequence_tps-1) * tr + 1/newfreq), 1/newfreq)
    GM_x = np.linspace(0, (sequence_tps - 1) * tr, sequence_tps)

    regr_len = len(regr_x)
    BH_len_upsampled = int(BH_len*newfreq)
    nBH = int(nBH)

    f = spint.interp1d(GM_x, GM, fill_value='extrapolate')
    GM_upsampled = f(regr_x)

    # Preparing central breathhold and CO2 trace for Xcorr
    # CO2 trace should have the equivalent of
    # ten tr of bad data at the beginning of the file
    if BH_len:
        last_tp = BH_len_upsampled*(nBH-1)
    else:
        last_tp = -1

    GM_cut = GM_upsampled[BH_len_upsampled:last_tp]
    co_conv_cut = co_conv[BH_len_upsampled:]

    # Detrend GM # Molly hinted it might be better not to
    # GM_dt = sgn.detrend(GM_cut, type='linear', bp=0)

    GM_cut_len = len(GM_cut)
    nrep = len(co_conv_cut) - GM_cut_len
    if BH_len and nrep > BH_len_upsampled:
        nrep = BH_len_upsampled

    GM_co_r = np.zeros(nrep)
    for i in range(0, nrep):
        GM_co_r[i] = np.corrcoef(GM_cut, co_conv[0+i:GM_cut_len+i].T)[1, 0]

    optshift = int(GM_co_r.argmax())
    textname = GM_name + '_optshift.1D'
    # #!#
    optshiftout = np.array((optshift/newfreq,0))
    np.savetxt(textname, optshiftout, fmt='%.4f')
    co_shift = co_conv[optshift:optshift+regr_len]

    # preparing for and exporting figures of shift
    time_axis = np.arange(0, nrep/newfreq, 1/newfreq)
    # #!# I should change to following line but not tested yet
    # time_axis = np.linspace(0, nrep*newfreq, nrep)
    if nrep < len(time_axis):
        time_axis = time_axis[:nrep]
    elif nrep > len(time_axis):
        time_axis = np.pad(time_axis, (0, int(nrep - len(time_axis))), 'linear_ramp')

    plt.figure(figsize=FIGSIZE, dpi=SET_DPI)
    plt.plot(time_axis, GM_co_r)
    plt.title('optshift')
    plt.savefig(GM_name + '_optshift.png', dpi=SET_DPI)

    plt.figure(figsize=FIGSIZE, dpi=SET_DPI)
    plt.plot(sct.zscore(co_shift), '-', sct.zscore(GM_upsampled), '-')
    plt.title('GM and shift')
    plt.savefig(GM_name + '_co_regr.png', dpi=SET_DPI)

    export_regressor(regr_x, co_shift, GM_x, GM_name, '_co_regr')

    # Create folder
    GM_dir = GM_name + '_regr_shift'
    if not os.path.exists(GM_dir):
        os.makedirs(GM_dir)

    # Set num of fine shifts: 9 seconds is a bit more than physiologically feasible
    nrep = int(9 * newfreq)

    # Padding regressor for shift, and padding optshift too
    if nrep > optshift:
        left_pad = nrep - optshift
    else:
        left_pad = 0

    if (optshift + nrep + regr_len) > len(co_conv):
        right_pad = (optshift + nrep + regr_len) - len(co_conv)
    else:
        right_pad = 0

    co_padded = np.pad(co_conv, (int(left_pad), int(right_pad)),'mean')
    optshift_padded = optshift + left_pad

    for i in range(-nrep, nrep):
        co_shift = co_padded[optshift_padded-i:optshift_padded-i+regr_len]
        suffix = '/shift_' + '%04d' % (i + nrep)
        export_regressor(regr_x, co_shift, GM_x, GM_dir, suffix)


def onpick_manualedit(event):
    thisline = event.artist
    xdata = thisline.get_xdata()
    ind = event.ind
    xdataind = xdata[ind]
    print('onpick xind:',xdataind)

    global reject_list
    reject_list = np.append(reject_list,xdataind)


def partone(filename, channel=4, tr=1.5, newfreq=40):
    # data_dec = decimate_data(filename, newfreq)
    data_dec = np.genfromtxt(filename + '_dec.tsv.gz')
    resp_filt = filter_signal(data_dec, channel)
    [co, pidx] = get_peaks(resp_filt)
    # export original peaks
    textname = filename + '_autopeaks.1D'
    np.savetxt(textname, pidx)
    textname = filename + '_co_orig.1D'
    np.savetxt(textname, co)
    return co, pidx


def manualchange(filename, pidx, reject_list):
    # masking pidx and saving manual selection
    pidx = np.array(pidx)
    pmask = np.in1d(pidx, reject_list)
    npidx = list(pidx[~pmask])
    textname = filename + '_manualpeaks.1D'
    np.savetxt(textname, npidx)
    return npidx


# def parttwo(filename):
def parttwo(co, pidx, filename, GM_name, tr=1.5, newfreq=40, BH_len=58, nBH=8):
    hrf = create_hrf(newfreq)
    co_conv = get_petco2(co, pidx, hrf, filename)
    if not os.path.exists('regr'):
        os.makedirs('regr')

    # co_conv = np.genfromtxt('regr/' + filename + '_co_conv.1D')
    #!#
    get_regr(GM_name, co_conv, tr, newfreq, BH_len, nBH)
    plt.close('all')


def _main(argv=None):
    options = _get_parser().parse_args(argv)
    # Reading first data
    # newfreq = 40
    # nrep = 2000
    # tr = 1.5
    # tps=340
    filename = options.filename
    GM_name = options.GM_name
    newfreq = options.newfreq
    tr = options.tr
    channel = options.channel
    BH_len = options.BH_len
    nBH = options.nBH

    co, pidx = partone(filename, channel, tr, newfreq)
    parttwo(co, pidx, filename, GM_name, tr, newfreq, BH_len, nBH)


if __name__ == '__main__':
    _main()