#!/usr/bin/env python

"""
Module to build template-based models used for the hybridization.
"""


import os
import sys
import numpy as np
import subprocess as sp

import path
import libcommon
import hybrid.libhhsuite

N_TEMPL = 100
N_HYBRID = 10
N_HYBRID_MIN = 2
TM_CUTOFF = 0.6 ; TM_OFFSET = 0.2

VERBOSE=False
EXCLUDE = []
if libcommon.TBM_EXCLUDE is not None:
    with open(libcommon.TBM_EXCLUDE) as fp:
        for line in fp:
            EXCLUDE.append(line.strip())

errfile = None  # '/dev/null'


def calc_tmalign(ref_fn, pdb_s, out_fp, n_proc):
    cmd = []
    if os.getenv("PREFMD2_PYTHON_MPROCS") is None: # or 1:
        _n_proc = str(n_proc)
    else:
        _n_proc = os.getenv("PREFMD2_PYTHON_MPROCS")
    cmd.append('--cpu')
    cmd.append('%s'%_n_proc)
    cmd.append('--ref')
    cmd.append(ref_fn.short())
    cmd.append("--list")
    cmd.append(pdb_s)
    cmd.append("--out")
    cmd.append(out_fp)
    libcommon.asystem(module="hybrid.exec_calc_tmalign", args=cmd)
    with open(out_fp, "r") as i_fh:
        scores = [float(line.strip().split()[0]) for line in i_fh.readlines()]
    return scores

def calc_tmscore(ref_fn, pdb_s, out_fp, n_proc):
    cmd = []
    if os.getenv("PREFMD2_PYTHON_MPROCS") is None:
        _n_proc = str(n_proc)
    else:
        _n_proc = os.getenv("PREFMD2_PYTHON_MPROCS")
    cmd.append('--cpu')
    cmd.append('%s'%n_proc)
    cmd.append('--ref')
    cmd.append(ref_fn.short())
    cmd.append("--list")
    cmd.append(pdb_s)
    cmd.append("--out")
    cmd.append(out_fp)
    libcommon.asystem(module="exec_calc_tmscore", args=cmd)
    with open(out_fp, "r") as i_fh:
        scores = [float(line.strip().split()[1]) for line in i_fh.readlines()]
    return scores

def run_hhsuite(homolog_home, id, fa_fn, input_pdb, n_proc):

    # First run hhblits (on a sequence database) and hhsearch (on the PDB70
    # database).
    hh_fn = homolog_home.fn("%s.vit.local"%id)
    cmd = [fa_fn.short()]
    cmd.extend(['-p', 'simple'])
    cmd.extend(['-d', os.getenv("HHSUITE_SEQ_DB")])
    cmd.extend(['-c', str(n_proc)])
    libcommon.asystem(module="hybrid.exec_run_hhpred", args=cmd,
                      errfile=errfile)

    # Parse the results.
    hh_s = hybrid.libhhsuite.parse_hhr(hh_fn.short())

    # Fetch the PDB files of the templates and confront them with the initial
    # 3D model using TMalign.
    tm_fn = homolog_home.fn("tm.dat")
    if not tm_fn.status():

        # Fetch the .pdb file of the template from the remote PDB.
        pdb_fn_s = []
        for hh in hh_s[:N_TEMPL]:
            pdb_id = hh.pdb_id
            if pdb_id[:4] in EXCLUDE: continue
            #
            pdb_fn = homolog_home.fn("%s.pdb"%pdb_id)
            if not pdb_fn.status():
                libcommon.asystem(module="hybrid.exec_pdb_get",
                                  args=['-f', pdb_id],
                                  stdout=True, errfile=errfile)
            if pdb_fn.status():
                pdb_fn_s.append(pdb_fn)
        with homolog_home.fn("pdb_s").open("wt") as fout:
            for pdb_fn in pdb_fn_s:
                fout.write(pdb_fn.short() + '\n')

        # Compare the templates with the initial 3D model.
        tm_s = calc_tmalign(input_pdb, "pdb_s", tm_fn.path(), n_proc)
        selected = []
        with tm_fn.open("wt") as fout:
            for pdb_fn,tm in zip(pdb_fn_s, tm_s):
                fout.write("%s  "%pdb_fn.name())
                fout.write(" %6.4f"%tm)
                fout.write("\n")
                if tm > TM_CUTOFF:
                    if pdb_fn.name() not in selected:
                        selected.append(pdb_fn.name())
                pdb_fn.remove()
        with homolog_home.fn("selected").open("wt") as fout:
            for pdb_id in selected:
                fout.write("%s\n"%pdb_id)

    else:
        selected = []
        with homolog_home.fn("selected").open() as fp:
            for line in fp:
                selected.append(line.strip())

    return selected

def run_modeller(homolog_home, id, fa_fn, input_pdb, selected, n_proc):
    model_fn = homolog_home.fn("model_s")
    if not model_fn.status():

        # Actually builds the 3D template-based models with MODELLER.
        cmd = []
        cmd.append(fa_fn.short())
        cmd.append("-p")
        cmd.append("build")
        cmd.append("--include")
        cmd.append("selected")
        cmd.append("--force")
        cmd.append("--alt")
        cmd.append("3")
        cmd.append("--cpu")
        cmd.append("%s"%n_proc)
        cmd.append('-d')
        cmd.append(os.getenv("HHSUITE_SEQ_DB"))
        cmd.append("--hhdb")
        cmd.append(os.getenv("HHSUITE_PDB_DB"))
        if len(selected) > 0:
            libcommon.asystem(module="hybrid.exec_run_hhpred", args=cmd,
                              stdout=True, errfile=errfile)

        # Modify the numbers assigned to residues of the MODELLER models to
        # match the initial model.
        model_s = []
        for pdb_id in selected:
            model = path.Path.glob("%s-%s*.model.pdb"%(id, pdb_id))
            for m in model:
                name = '.'.join(m.fname().split(".")[:-2])
                out_fn = homolog_home.fn("%s.pdb"%name)
                if not out_fn.status():
                    cmd = []
                    cmd.append(input_pdb.short())
                    cmd.append(m.short())
                    with out_fn.open("wt") as fout:
                        libcommon.asystem(module="hybrid.exec_match_resNo",
                                          args=cmd, outfile=fout,
                                          stdout=True, errfile=errfile)
                if out_fn not in model_s:
                    model_s.append(out_fn)

        with model_fn.open("wt") as fout:
            for model in model_s:
                fout.write('%s\n'%model.short())
    else:
        model_s = []
        with model_fn.open() as fp:
            for line in fp:
                model_s.append(path.Path(line.strip()))

    # Compare the templated-based models with the initial model using the
    # TMscore program.
    tm_fn = homolog_home.fn("model.dat")
    if not tm_fn.status():
        tm_s = calc_tmscore(input_pdb, "model_s", tm_fn.path(), n_proc)
        with tm_fn.open("wt") as fout:
            for pdb_fn,tm in zip(model_s, tm_s):
                fout.write("%6.4f %s\n"%(tm, pdb_fn))
    else:
        tm_s = []
        with tm_fn.open() as fp:
            for line in fp:
                x = line.strip().split()
                tm_s.append(float(x[0]))
    return model_s, tm_s

def select_init(hybrid_home, input_pdb, model_s, tm_s):
    if len(model_s) > 0:
        tm_cutoff = max(max(tm_s)-TM_OFFSET, TM_CUTOFF)
        sorted_index = np.argsort(tm_s)[::-1]
        init_s = [] ; init_s.append(input_pdb)
        for i in sorted_index:
            tm = tm_s[i]
            model = model_s[i]
            if tm > tm_cutoff:
                init_s.append(model)
        init_s = init_s[:N_HYBRID]
    else:
        init_s = [input_pdb]
    #
    hybrid_home.chdir()
    if len(init_s) < N_HYBRID_MIN:
        with hybrid_home.fn("DONE").open("wt") as fout:
            fout.write("Failed to detect templates for hybridization\n")
    else:
        with hybrid_home.fn("init_s").open("wt") as fout:
            n_init = len(init_s)
            for i in range(N_HYBRID):
                fout.write("%s\n"%init_s[i%n_init].short())

def use_selected(hybrid_home, input_pdb, selected_home):
    init_s = [] ; init_s.append(input_pdb)
    init_s.extend(selected_home.glob("*.pdb"))
    #
    hybrid_home.chdir()
    with hybrid_home.fn("init_s").open("wt") as fout:
        n_init = len(init_s)
        for i in range(N_HYBRID):
            fout.write("%s\n"%init_s[i%n_init].short())

def main(args=None):
    if args is not None:
        arg_1 = args[0]
        arg_2 = args[1]
        n_proc = int(args[2])
    else:
        arg_1 = sys.argv[1]
        arg_2 = sys.argv[2]
        n_proc = int(args[3])

    id = arg_1
    input_pdb = path.Path(arg_2)
    #
    selected_home = path.Dir("selected")
    homolog_home = path.Dir("homolog", build=True)
    hybrid_home = path.Dir("hybrid", build=True)
    #
    if not selected_home.status():
        #
        homolog_home.chdir()
        # Prepare the input for HHsuite.
        fa_fn = homolog_home.fn("%s.fa"%id)
        with fa_fn.open("wt") as fout:
            cmd = [input_pdb.short()]
            libcommon.asystem(module="hybrid.exec_pdb_seq", args=cmd,
                              outfile=fout,
                              stdout=True, errfile=errfile)
        # First run HHsuite tools to identify some templates.
        selected = run_hhsuite(homolog_home, id, fa_fn, input_pdb, n_proc)
        # The build template-based models with MODELLER.
        model_s, tm_s = run_modeller(homolog_home, id, fa_fn, input_pdb,
                                     selected, n_proc)
        #
        select_init(hybrid_home, input_pdb, model_s, tm_s)
    else:
        use_selected(hybrid_home, input_pdb, selected_home)

if __name__ == '__main__':
    main()
