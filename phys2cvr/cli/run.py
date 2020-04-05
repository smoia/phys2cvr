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
    required.add_argument('-in', '--input-file',
                          dest='filename',
                          type=str,
                          help=('The name of the BIDS formatted '
                                'file containing physiological data, with or '
                                'without extension.'),
                          required=True)
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
