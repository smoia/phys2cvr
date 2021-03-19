#!/usr/bin/env python3

import datetime
import logging
import os
import sys
from copy import deepcopy

import numpy as np
from peakdet.io import load_physio

from phys2cvr import io, signal, stats, _version
from phys2cvr.cli.run import _get_parser


SET_DPI = 100
FIGSIZE = (18, 10)
LGR = logging.getLogger(__name__)
LGR.setLevel(logging.INFO)
EXT_1D = ['.txt', '.csv', '.tsv', '.1d', '.par', '.tsv.gz']
EXT_NIFTI = ['.nii', '.nii.gz']


def save_bash_call(outdir):
    """
    Save the bash call into file `p2d_call.sh`.

    Parameters
    ----------
    metric : function
        Metric function to retrieve arguments for
    metric_args : dict
        Dictionary containing all arguments for all functions requested by the
        user
    """
    arg_str = ' '.join(sys.argv[1:])
    call_str = f'phys2cvr {arg_str}'
    outdir = os.path.abspath(outdir)
    log_path = os.path.join(outdir, 'code', 'logs')
    os.makedirs(log_path, exist_ok=True)
    isotime = datetime.datetime.now().strftime('%Y-%m-%dT%H%M%S')
    f = open(os.path.join(log_path, f'p2c_call_{isotime}.sh'), "a")
    f.write(f'#!bin/bash \n{call_str}')
    f.close()


def phys2cvr(fname_func, fname_co2='', fname_pidx='', fname_mask='', outdir='',
             freq='', tr='', trial_len='', n_trials='', highcut='', lowcut='',
             apply_filter='', no_pad='', no_phys=False, do_regression=False,
             lagged_regression=True, maxlag=9, denoise_matrix='',
             quiet=False, debug=False):
    """
    Run main workflow of phys2cvr.
    """
    # Add logger and suff
    outdir = os.path.abspath(outdir)
    os.makedirs(outdir, exist_ok=True)
    petco2log_path = outdir
    # petco2log_path = os.path.join(outdir, 'code', 'petco2log')
    # os.makedirs(petco2log_path, exist_ok=True)

    # Create logfile name
    basename = 'phys2cvr_'
    extension = 'tsv'
    isotime = datetime.datetime.now().strftime('%Y-%m-%dT%H%M%S')
    logname = os.path.join(petco2log_path, (basename + isotime + '.' + extension))

    # Set logging format
    log_formatter = logging.Formatter(
        '%(asctime)s\t%(name)-12s\t%(levelname)-8s\t%(message)s',
        datefmt='%Y-%m-%dT%H:%M:%S')

    # Set up logging file and open it for writing
    log_handler = logging.FileHandler(logname)
    log_handler.setFormatter(log_formatter)
    sh = logging.StreamHandler()

    if quiet:
        logging.basicConfig(level=logging.WARNING,
                            handlers=[log_handler, sh], format='%(levelname)-10s %(message)s')
    elif debug:
        logging.basicConfig(level=logging.DEBUG,
                            handlers=[log_handler, sh], format='%(levelname)-10s %(message)s')
    else:
        logging.basicConfig(level=logging.INFO,
                            handlers=[log_handler, sh], format='%(levelname)-10s %(message)s')

    version_number = _version.get_versions()['version']
    LGR.info(f'Currently running phys2cvr version {version_number}')
    LGR.info(f'Input file is {fname_func}')

    # Check func type and read it
    func_is_1d = io.check_ext(EXT_1D, fname_func)
    func_is_nifti = io.check_ext(EXT_NIFTI, fname_func)

    if func_is_1d:
        if tr:
            func_avg = np.genfromtxt(fname_func)
            if apply_filter:
                LGR.info('Applying butterworth filter to {fname_func}')
                func_avg = signal.filter_signal(func_avg, tr, lowcut, highcut)
        else:
            raise Exception('Provided functional signal, but no TR specified! '
                            'Rerun specifying the TR')
    elif func_is_nifti:
        func, dmask, img = io.load_nifti_get_mask(fname_func)
        if len(func.shape) < 4:
            raise Exception(f'Provided functional file {fname_func} is not a 4D file!')
        # Read TR or declare its overwriting
        if tr:
            LGR.warning(f'Forcing TR to be {tr} seconds')
        else:
            tr = img.header['pixdim'][4]

        # Read mask if provided
        if fname_mask:
            _, mask, _ = io.load_nifti_get_mask(fname_mask)
            if func.shape[:3] != mask.shape:
                raise Exception(f'{fname_mask} and {fname_func} have different sizes!')
            mask = mask * dmask
        else:
            mask = dmask
            LGR.warning(f'No mask specified, using any voxel different from 0 in {fname_func}')

        if apply_filter:
            LGR.info('Obtaining filtered average of {fname_func}')
            func_filt = signal.filter_signal(func, tr, lowcut, highcut)
            func_avg = func_filt[mask].mean(axis=-1)
        else:
            func_avg = func[mask].mean(axis=-1)

    else:
        raise Exception(f'{fname_func} file type is not supported yet, or '
                        'the extension was not specified.')

    if fname_co2 == '' and not no_phys:
        raise Exception('The pipeline for no physiological files was not selected, and '
                        'no file for physiological regressor was found. Please rerun with '
                        'one option or the other.')
    elif no_phys:
        LGR.info(f'Computing "CVR" maps using {fname_func} only')
        if func_is_1d:
            LGR.warning('Using an average signal only, solution might be unoptimal.')

            if not apply_filter:
                LGR.warning('No filter applied to the input average! You know '
                            'what you are doing, right?')

        petco2hrf = func_avg

    else:
        co2_is_phys = io.check_ext('.phys', fname_co2)
        co2_is_1d = io.check_ext(EXT_1D, fname_co2)

        if co2_is_1d:
            if fname_pidx:
                pidx = np.genfromtxt(fname_pidx)
            else:
                raise Exception(f'{fname_co2} file is a text file, but no '
                                'file containing its peaks was provided. '
                                ' Please provide peak file!')

            if not freq:
                raise Exception(f'{fname_co2} file is a text file, but no '
                                'frequency was specified. Please provide peak '
                                ' file!')

            co2 = np.genfromtxt(fname_co2)
        elif co2_is_phys:
            # Read a phys file!
            phys = load_physio(fname_co2)

            co2 = phys.data
            pidx = phys.peaks
            if freq:
                LGR.warning(f'Forcing CO2 frequency to be {freq} Hz')
            else:
                freq = phys.fs
        else:
            raise Exception(f'{fname_co2} file type is not supported yet, or '
                            'the extension was not specified.')

        # Set output file & path - calling splitext twice cause .gz
        basename_co2 = os.path.splitext(os.path.splitext(os.path.basename(fname_co2))[0])[0]
        outname = os.join(outdir, basename_co2)
        petco2hrf = signal.convolve_petco2(co2, pidx, freq, outname)

    regr = stats.get_regr(func_avg, petco2hrf, tr, freq, outname, maxlag, trial_len,
                          n_trials, no_pad, '.1D', lagged_regression)

    # Add internal regression if required!
    if func_is_nifti and do_regression:
        print('Running regression!')
        pass

        outfuncname = os.path.splitext(os.path.splitext(fname_func)[0])[0]
        newdim = deepcopy(img.header['dim'])
        newdim[0], newdim[4] = 3, 1
        oimg = deepcopy(img)
        oimg.header['dim'] = newdim

        # Read in eventual denoising factors
        if denoise_matrix:
            mat_conf = np.genfromtxt(denoise_matrix)
            if not np.any(np.all(mat_conf == 1, axis=1)):
                mat_conf = np.stack(mat_conf, np.ones_like(regr))
        else:
            mat_conf = np.ones_like(regr)

        mat = np.stack([regr, mat_conf])

        # #!# func has to be SPC!

        model = OLSModel(mat).fit(func)

        beta = model.predicted()

        # #!# beta is not cvr!

        io.export_nifti(beta, img, f'{outfuncname}_cvr_simple')

        if lagged_regression:
            nrep = int(maxlag * freq * 2)
            outprefix = os.path.join(os.path.split(outname)[0], 'regr', os.path.split(outname)[1])
            r_square = np.empty((func.shape[0], func.shape[1], func.shape[2], nrep))
            beta_all = np.empty((func.shape[0], func.shape[1], func.shape[2], nrep))

            for i in range(nrep):
                regr = np.genfromtxt(f'{outprefix}_{(i + nrep):04g}')

                mat = np.stack([regr, mat_conf])

                model = OLSModel(mat).fit(func)

                beta_all[:, :, :, i] = model.predicted()
                r_square[:, :, :, i] = model.r_square()

            lag_idx = np.argmax(r_square, axis=-1)
            lag = lag_idx / freq
            beta = beta_all[:, :, :, lag_idx]

            io.export_nifti(beta, oimg, f'{outfuncname}_cvr')
            io.export_nifti(lag, oimg, f'{outfuncname}_lag')

    elif do_regression:
        LGR.warning('The input file is not a nifti volume. At the moment, '
                    'regression is not supported for other formats.')

    LGR.info('phys2cvr finished! Enjoy your outputs!')


def _main(argv=None):
    options = _get_parser().parse_args(argv)

    save_bash_call(options.outdir)

    phys2cvr(**vars(options))


if __name__ == '__main__':
    _main(sys.argv[1:])
