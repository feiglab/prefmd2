#!/usr/bin/env python

"""
Module to actually run a full initial 3D model hybridization task:
    - First look for templates and build template-based models with MODELLER
      (using the hybrid.exec_prep_hybrid module).
    - Then hybridize with Rosetta (using the hybrid.exec_run_hybrid module).
"""

import os
import sys
import json
import argparse

import warnings
warnings.filterwarnings("ignore")

import path
import libcommon
from hybrid.libpdb import Sequence, PDB

EXTENSION = os.getenv("ROSETTA_EXTENSION", "linuxgccrelease")


def run(arg):

    model_summary = path.Path("model_s")
    original_cwd = os.getcwd()

    # Prepare hybridization: search for templates using HHsuite, select the
    # templates and build 3D models with MODELLER.
    print(("- Attempting to build template-based models with HHsuite and"
           " MODELLER..."))
    model_NA = path.Path("hybrid/DONE")
    model_list = path.Path("hybrid/init_s")
    cmd = [arg.output_prefix, arg.init_pdb.short(), str(arg.n_proc)]
    libcommon.asystem(module="hybrid.exec_prep_hybrid", args=cmd)

    #
    if model_NA.status():   # Not enough template
        with model_summary.open("wt") as fout:
            fout.write("#")
        return

    # Run hybridization.
    sel_s = path.Path.glob("hybrid/iter_*/sel.out")
    if len(sel_s) < 10:
        print("- Starting hybridization using Rosetta...")
        cmd = [model_list.short(), '%d'%arg.n_proc]
        libcommon.asystem(module="hybrid.exec_run_hybrid", args=cmd)

    os.chdir(original_cwd)

    # Parse hybridization output.
    sel_s = path.Path.glob("hybrid/iter_*/sel.out")
    if len(sel_s) == 0:
        raise NotImplementedError("ERROR!")
        with model_summary.open("wt") as fout:
            fout.write("#")
        return

    # Extract hybridization output.
    print("- Extracting hybridization output...")
    final = sorted(sel_s, key=lambda x: int(x.split("/")[-2].split("_")[-1]))[-1]
    model_home = path.Dir("model", build=True)
    model_home.chdir()
    #
    cmd = ['%s/main/source/bin/extract_pdbs.%s'%(os.environ['ROSETTA_HOME'],
                                                 EXTENSION)]
    cmd.extend(['-in:file:silent', final.short()])
    libcommon.system(cmd, stdout=True)
    #
    fn_s = sorted(path.Path.glob("iter*.*.pdb"), key=lambda x: int(x.split(".")[-2]))
    model_s = []
    for i,fn in enumerate(fn_s):
        out_fn = path.Path("model.%d.pdb"%i)
        with out_fn.open("wt") as fout:
            libcommon.asystem(module="hybrid.exec_match_resNo",
                              args=[arg.init_pdb.short(), fn.short()],
                              outfile=fout)
        model_s.append(out_fn)
    #
    with model_summary.open("wt") as fout:
        for model in model_s:
            fout.write("%s\n"%model)

    os.chdir(original_cwd)


def main(args=None):
    arg = argparse.ArgumentParser(prog='hybrid')
    arg.add_argument(dest='output_prefix')
    arg.add_argument(dest='init_pdb')
    arg.add_argument('-j', '--cpu', dest='n_proc', type=int, default=20)
    #
    arg = arg.parse_args(args)
    arg.init_pdb = path.Path(arg.init_pdb)
    #
    run(arg)

if __name__=='__main__':
    main()
