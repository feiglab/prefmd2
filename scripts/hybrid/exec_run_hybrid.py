#!/usr/bin/env python

"""
Script to use Rosetta to run the hybridization process of the PREFMD pipeline.
It will search for fragments for the input protein using the 'fragment_picker'
application and then it will perform the actual hybridization using a modified
version of the IterationMaster.py from the Rosetta package.
"""


import os
import sys
import numpy as np
import subprocess as sp

import path
import libcommon

ROSETTA_HOME = os.getenv("ROSETTA_HOME")
EXTENSION = os.getenv("ROSETTA_EXTENSION", "linuxgccrelease")
module_dirpath = os.path.dirname(sys.modules[__name__].__file__)
SCRIPT_HOME = os.path.join(module_dirpath, "iterative_hybridize")
N_ITERATIONS = 10
# PSIPRED_EXEC='%s/apps/psipred/current/run_psipred'%os.getenv("HOME")

def prep_init_models(init_s):
    prep_s = []
    init = path.Path("init.pdb") ; prep_s.append(init)
    if not init.status():
        with init.open("wt") as fout:
            libcommon.system("convpdb.pl -out generic -renumber 1 %s"%init_s[0].short(),
                             outfile=fout)
    #
    disu_s = []
    with init.open() as fp:
        for line in fp:
            if line.startswith("SSBOND"):
                disu_s.append("%s %s\n"%(line[17:21], line[31:35]))
    with open("disulf.def", 'wt') as fout:
        if len(disu_s) > 0:
            fout.writelines(disu_s)
        else:
            fout.write("1 1\n")
    #
    for i,model in enumerate(init_s[1:]):
        prep = path.Path("model.%d.pdb"%(i+1))
        if not prep.status():
            with prep.open("wt") as fout:
                libcommon.asystem(module="hybrid.exec_match_resNo",
                                  args=[init.short(), model.short()],
                                  outfile=fout)
        prep_s.append(prep)
    return prep_s

def build_frag(fa_fn, n_proc):
    """
    Use Rosetta's fragment_picker application to pick fragments for the
    hybridization.
    """

    name = fa_fn.name()
    if not os.path.exists("%s.200.3mers"%name) or not os.path.exists("%s.200.9mers"%name):
        # ss2_fn = path.Path("%s.ss2"%name)
        # if not ss2_fn.status():
        #     cmd = [PSIPRED_EXEC, fa_fn.short()]
        #     system(cmd, errfile='/dev/null')
        # cmd = ['python', "%s/tools/fragment_tools/make_fragments.py"%ROSETTA_HOME]
        # cmd.extend(['-cpus', '%d'%n_proc])
        # cmd.extend(['-psipredfile', ss2_fn.path()])
        # cmd.append(fa_fn.short())

        cmd = ["%s/main/source/bin/fragment_picker.%s" % (ROSETTA_HOME,
                                                          EXTENSION),
               # Input data.
               "-in:file:fasta", "input.fa",
               # Database.
               "-in:file:vall",
               "%s/tools/fragment_tools/vall.jul19.2011.gz" % ROSETTA_HOME,
               # Output.
               "-out:file:frag_prefix", name]
        libcommon.system(cmd)

    if not os.path.exists("%s.200.3mers"%name) or not os.path.exists("%s.200.9mers"%name):
        return False
    if not os.path.exists("t000_.3mers"):
        os.symlink("%s.200.3mers"%name, "t000_.3mers")
    if not os.path.exists("t000_.9mers"):
        os.symlink("%s.200.9mers"%name, "t000_.9mers")
    return True

def prep_run(init_s):
    if not os.path.exists("model.out"):
        cmd = []
        cmd.append("%s/main/source/bin/combine_silent.%s"%(ROSETTA_HOME, EXTENSION))
        cmd.append("-in:file:s")
        cmd.extend([fn.short() for fn in init_s])
        cmd.extend(["-out:file:silent", "model.out"])
        cmd.extend(["-out:file:silent_struct_type","binary"])
        cmd.extend(["-ignore_zero_occupancy","false"])
        libcommon.system(cmd)
    #
    if not os.path.exists("ref.out"):
        cmd = []
        cmd.append("%s/main/source/bin/iterhybrid_selector.%s"%(ROSETTA_HOME, EXTENSION))
        cmd.extend(["-in:file:silent", "model.out"])
        cmd.extend(["-in:file:template_pdb", "init.pdb"])
        cmd.extend(["-cm:similarity_cut", "0.2"])
        cmd.extend(["-out:file:silent ref.out"])
        cmd.extend(["-out:nstruct", "%d"%len(init_s)])
        cmd.extend(["-out:prefix", "iter0"])
        cmd.extend(["-score:weights", "ref2015_cart"])
        cmd.extend(["-silent_read_through_errors"])
        cmd.extend(["-in:file:silent_struct_type", "binary"])
        cmd.extend(["-out:file:silent_struct_type", "binary"])
        cmd.extend(["-mute", "core", "basic"])
        libcommon.system(cmd)
    #
    if not os.path.exists("cen.cst"):
        cmd = []
        cmd.extend([fn.short() for fn in init_s])
        with open("cen.cst", 'wt') as fout:
            libcommon.asystem(module="hybrid.generate_harmonic_cst",
                              args=cmd, outfile=fout)
    if not os.path.exists("fa.cst"):
        os.symlink("cen.cst", 'fa.cst')

def main(args=None):

    # Set the up the input variables.
    if args is None:
        init_fn_s = path.Path(sys.argv[1])
        if len(sys.argv) > 2:
            n_proc = int(sys.argv[2])
        else:
            n_proc = 20
    else:
        init_fn_s = path.Path(args[0])
        if len(args) > 1:
            n_proc = int(args[1])
        else:
            n_proc = 20
    #
    hybrid_home = init_fn_s.dirname()
    hybrid_home.chdir()
    #
    init_s = []
    with init_fn_s.open() as fp:
        for line in fp:
            init_s.append(path.Path(line.strip()))


    # Prepare the input for the hybridization.
    init_s = prep_init_models(init_s)
    #
    fa_fn0 = hybrid_home.fn("input0.fa")
    fa_fn = hybrid_home.fn("input.fa")
    if (not fa_fn0.status()) or (not fa_fn.status()):

        out_s = libcommon.asystem(module="hybrid.exec_pdb_seq",
                                  args=[init_s[0].short()],
                                  stdout=True).split("\n")[:-1]
        with fa_fn0.open("wt") as fout:
            fout.write(">init\n")
            for line in out_s:
                if not line.startswith(">"):
                    fout.write("%s\n"%line)
        with fa_fn.open("wt") as fout:
            fout.write(">init\n")
            for line in out_s[1:]:
                if line.startswith(">"):
                    fout.write("/\n")
                else:
                    fout.write("%s"%line)
            fout.write("\n")


    # Search for fragments.
    if not build_frag(fa_fn0, n_proc):
        sys.exit("Failed to build fragments\n")


    # Actually run the hybridization using a modified IterationMaster.py
    # script.
    prep_run(init_s)
    #
    with open("NODEFILE", 'wt') as fout:
        fout.write("%d/localhost"%n_proc)
    #
    cmd = []
    cmd.append("python")
    cmd.append(os.path.join(SCRIPT_HOME, "IterationMaster.py"))
    cmd.extend(["-iha", "60"])
    cmd.extend(["-nodefile", "NODEFILE"])
    cmd.append("-simple")
    cmd.append("-niter")
    cmd.append(str(N_ITERATIONS))
    with open("iter_hybrid.log", 'wt') as fout:
        libcommon.system(cmd, outfile=fout, errfile=fout)

if __name__ == '__main__':
    main()
