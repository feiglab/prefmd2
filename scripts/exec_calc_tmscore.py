#!/usr/bin/env python

import os
import sys
import argparse

import libscore


def main(args=None):
    arg = argparse.ArgumentParser(prog='calc_tmscore')
    arg.add_argument('-r', '--ref', dest='ref_fn', help='reference PDB file',
                     required=True)
    arg.add_argument('-l', '--list', dest='pdb_list', help='PDB file lists',
                     required=True)
    arg.add_argument('-o', '--out', dest='out_fp', help='output filepath',
                     default="qual_init")
    arg.add_argument('-j', '--cpu', dest='n_proc', help='number of proc.',
                     default=1, type=int)
    arg = arg.parse_args(args)

    libscore.run_scoring_tmscore(reference_fp=arg.ref_fn,
                                  input_fpl=arg.pdb_list,
                                  out_fp=arg.out_fp, n_jobs=arg.n_proc)

if __name__=='__main__':
    main()
