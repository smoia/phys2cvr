# -*- coding: utf-8 -*-
"""
Parser for phys2cvr.
"""

import argparse

from phys2cvr import __version__


def _get_parser():
    """
    Parses command line inputs for this function

    Returns
    -------
    parser.parse_args() : argparse dict

    Notes
    -----
    # Argument parser follow template provided by RalphyZ.
    # https://stackoverflow.com/a/43456577
    """
    parser = argparse.ArgumentParser()
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('Required Argument:')
    required.add_argument('-i', '--input-func',
                          dest='fname_func',
                          type=str,
                          help=('Complete path (absolute or relative) and name '
                                'of the file containing fMRI signal. This file '
                                'can be a nifti file or a 1D txt file.'),
                          required=True)
    optional.add_argument('-co2', '--input-co2',
                          dest='fname_co2',
                          type=str,
                          help=('Complete path (absolute or relative) and name '
                                'of the file containing CO2 signal (or equivalent '
                                'physiological trace to compute the regressor). '
                                'This file can be a 1D txt file or a .phys file '
                                'from peakdet. If nothing is specified, the '
                                'average timeseries of the mask will be used '
                                'as regressor.'),
                          default='')
    optional.add_argument('-pk', '--input-peaks',
                          dest='fname_pidx',
                          type=str,
                          help=('Complete path (absolute or relative) and name '
                                'of the file containing the peak of the '
                                'physiological trace. Required if the physiological '
                                'trace file is not a .phys file. '
                                'Use this option to overwrite the peaks specified '
                                'in the .phys file.'),
                          default='')
    optional.add_argument('-m', '--input-mask',
                          dest='fname_mask',
                          type=str,
                          help=('Complete path (absolute or relative) and name '
                                'of the file containing a brain mask (nifti file). '
                                'This mask will be used to extract the functional '
                                'signal to run the cross correlation with the '
                                'physiological regressor. Use this option to '
                                'specify a GM mask or overwrite a full brain mask.'
                                'If the functional file is specified and this '
                                'option is not used, or the mask cannot be '
                                'loaded, the program will create a mask using '
                                'any voxel of the functional file constantly '
                                'different from zero.'),
                          default='')
    optional.add_argument('-o', '--output-directory',
                          dest='outdir',
                          type=str,
                          help=('Complete path (absolute or relative) and name '
                                'of the desired output directory. If it does not '
                                'exist, it will be created. If it is not '
                                'specified, a folder named "phys2cvr" will be '
                                'created in the folder containing the functional '
                                'file.'),
                          default='')
    optional.add_argument('-fr', '--frequency',
                          dest='freq',
                          type=float,
                          help=('Frequency of the physiological trace. Required '
                                'if the latter is not passed as a .phys file. '
                                'Use this option to overwrite the frequency of a '
                                '.phys file.'),
                          default='')
    optional.add_argument('-tr', '--repetition-time',
                          dest='tr',
                          type=float,
                          help=('TR of functional data. Required '
                                'if the latter is not passed as a nifti file. '
                                'Use this option to overwrite the frequency of a '
                                'nifti file.'),
                          default='')
    optional.add_argument('-tlen', '--trial-length',
                          dest='trial_len',
                          type=float,
                          help=('Total duration of a respiration trial in seconds. '
                                'Specify this with the number of trials to run '
                                'a double cross-correlation between functional '
                                'signal and physiological regressor to improve '
                                'the detection of the bulk shift.'),
                          default='')
    optional.add_argument('-ntrial', '--trial-number',
                          dest='n_trials',
                          type=int,
                          help=('Number or respiration trials in the sequence. '
                                'Specify this with the duration of trials to run '
                                'a double cross-correlation between functional '
                                'signal and physiological regressor to improve '
                                'the detection of the bulk shift.'),
                          default='')
    optional.add_argument('-af', '--apply-filter',
                          dest='apply_filter',
                          action='store_true',
                          type=bool,
                          help=('Apply a filter to the functional data before '
                                'estimating the bulk shift. The filter will not '
                                'be applied on the data before the GLM computation. '
                                'If you want that, consider applying it before '
                                'running phys2cvr (or running it a second time).'),
                          default=False)
    optional.add_argument('-hf', '--highcut-frequency',
                          dest='highcut',
                          type=float,
                          help=('Higher frequency to use in signal filtering. '
                                'The filter will be applied only to the functional '
                                'data to estimate the bulk shift. This option '
                                'is suggested when only using a functional file.'),
                          default='')
    optional.add_argument('-lf', '--lowcut-frequency',
                          dest='lowcut',
                          type=float,
                          help=('Lower frequency to use in signal filtering. '
                                'The filter will be applied only to the functional '
                                'data to estimate the bulk shift. This option '
                                'is suggested when only using a functional file.'),
                          default='')
    optional.add_argument('-reg', '--run-regression',
                          dest='run_regression',
                          action='store_true',
                          type=bool,
                          help=('Actually run physiological regression(s) '
                                'internally. By default phys2cvr only '
                                'estimates and produces physiological regressors, '
                                'assuming that the regression itself will be '
                                'carried out with other software (e.g. AFNI, FSL, ...).'),
                          default=False)
    optional.add_argument('-skip_lagreg', '--skip-lagged-regression',
                          dest='lagged_regression',
                          action='store_false',
                          type=bool,
                          help=('Skip estimating the lagged regressors, '
                                'estimating only the central one. '
                                'Skip running the lagged regression if the regression is run.'),
                          default=True)
    optional.add_argument('-lm', '--lag-max',
                          dest='lag_max',
                          type=float,
                          help=('Maximum lag to consider during lag regression '
                                'in seconds. The same lag will be considered in '
                                'both directions. Default is ±9 seconds. '
                                'Remember that being this python, the upper limit '
                                'is excluded from the computation, i.e. default is'
                                '[-9, +8.7] or [-9, +9).'),
                          default=9)
    optional.add_argument('-ls', '--lag-step',
                          dest='lag_step',
                          type=float,
                          help=('Lag step to consider during lagged regression '
                                'in seconds. Default is 0.3 seconds.'),
                          default=0.3)
    optional.add_argument('-ldeg', '--legendre-degree',
                          dest='l_degree',
                          type=int,
                          help=('Maximum legendre degree to add to the regression '
                                'matrix as nuisance. Default is 0, to account for '
                                'the degree of freedom lost in computing the SPC.'),
                          default=0)
    optional.add_argument('-dmat', '--denoise-matrix',
                          dest='denoise_matrix',
                          action='append',
                          nargs='*',
                          type=str,
                          help=('Complete path (absolute or relative) and name '
                                'of denoising matrices to add to the regression model. '
                                'This option can be specified multiple times to '
                                'add multiple denoising matrices, but multiple '
                                'denoising matrices can be specified one after '
                                'the other, separated by a space.'),
                          default=None)
    optional.add_argument('-scale', '--scale-factor',
                          dest='scale_factor',
                          type=float,
                          help=('Scale factor for the CVR map output. Since BIDS '
                                'currently does not support mmHg as unit, if '
                                'using CO2 traces check their unit of measure '
                                'and their scaling factor to transform Volts into '
                                'mmHg. Use this option for other standardisations too.'),
                          default='')
    optional.add_argument('-lmap', '--lag-map',
                          dest='lag_map',
                          type=str,
                          help=('Complete path (absolute or relative) and name '
                                'of a previously computed lag map to use in '
                                'lagged regression.'),
                          default=None)
    optional.add_argument('-rdir', '--regr-dir',
                          dest='regr_dir',
                          type=str,
                          help=('Complete path (absolute or relative) and name '
                                'of previously computed lagged regressors to '
                                'use in a new lagged regression.'),
                          default=None)
    optional.add_argument('-sc', '--skip-conv',
                          dest='skip_conv',
                          action='store_true',
                          type=bool,
                          help=('Skip convolution of physiological trace. '
                                'By default phys2cvr convolves the physiological '
                                'trace with a standard HRF. Skip it when using '
                                'fMRI signal only.'),
                          default=False)
    optional.add_argument('-debug', '--debug',
                          dest='debug',
                          action='store_true',
                          help='Only print debugging info to log file. Default is False.',
                          default=False)
    optional.add_argument('-quiet', '--quiet',
                          dest='quiet',
                          action='store_true',
                          help='Only print warnings to log file. Default is False.',
                          default=False)
    optional.add_argument('-v', '--version', action='version',
                          version=('%(prog)s ' + __version__))
    return parser


if __name__ == '__main__':
    raise RuntimeError('phys2cvr/cli/run.py should not be run directly;\n'
                       'Please `pip install` phys2cvr and use the '
                       '`phys2cvr` command')
