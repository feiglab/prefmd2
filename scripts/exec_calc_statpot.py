#!/usr/bin/env python

import os
import sys
import argparse

import libscore


def main(args=None):

    arg = argparse.ArgumentParser(prog='calc_tmscore')
    arg.add_argument('-l', '--list', dest='pdb_list', help='PDB file lists',
                     required=True)
    arg.add_argument('-o', '--out', dest='out_fp', help='output filepath',
                     default="statpot.dat")
    arg.add_argument('-j', '--cpu', dest='n_proc', help='number of proc.',
                     default=1, type=int)
    arg.add_argument('--rwplus', dest='run_rwplus', help='evaluate RWplus',
                     default=False, action='store_true')
    arg.add_argument('--dfire', dest='run_dfire', help='evaluate dfire',
                     default=False, action='store_true')
    arg = arg.parse_args(args)

    libscore.run_scoring_statpot(input_fpl=arg.pdb_list, out_fp=arg.out_fp,
                                 n_jobs=arg.n_proc,
                                 use_rwplus=arg.run_rwplus,
                                 use_ddfire=arg.run_dfire)


if __name__ == '__main__':
    main()
