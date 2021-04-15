# -*- coding: utf-8 -*-
"""
Parser for phys2cvr
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
    required.add_argument('-f', '--input-func',
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
                                'from peakdet.'),
                          default='')
    optional.add_argument('-pk', '--input-peaks',
                          dest='fname_pidx',
                          type=str,
                          help=('Complete path (absolute or relative) and name '
                                'of the file containing the peak of the CO2 '
                                'signal. Required if the CO2 '
                                'signal file is not a .phys file. '
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



    required.add_argument('-gm', '--gmfile',
                          dest='GM_mask',
                          type=str,
                          help=('Fullpath to the file containing '
                                'the Gray Matter mask'),
                          required=True)
    optional.add_argument('-tr', '--tr',
                          dest='tr',
                          type=float,
                          help='TR of the GM data',
                          default=1.5)
    optional.add_argument('-nf', '--newfreq',
                          dest='newfreq',
                          type=float,
                          help=('Desired frequency to work with. As default, '
                                'it uses the frequency of the BIDS formatted '
                                'file'),
                          default=0)
    optional.add_argument('-bl', '--bh_len',
                          dest='BH_len',
                          type=float,
                          help='Duration of a whole Breathhold trial, in sec',
                          default=58)
    optional.add_argument('-nbh', '--nunmber_bh',
                          dest='n_BH',
                          type=int,
                          help='Number of Breathhold trial repetitions',
                          default=8)
    optional.add_argument('-ch', '--channel',
                          dest='channel',
                          type=int,
                          help=('Channel (column) containing CO2 trace.\n'
                                'Remember that numeration starts from 1!'),
                          default=4)

    parser._action_groups.append(optional)

    return parser


if __name__ == '__main__':
    raise RuntimeError('phys2cvr/cli/run.py should not be run directly;\n'
                       'Please `pip install` phys2cvr and use the '
                       '`phys2cvr` command')
