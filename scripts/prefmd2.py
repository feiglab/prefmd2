#!/usr/bin/env python

"""
Main script of PREFMD2. Basic usage:
    python prefmd2.py -t job_title -i initial_model.pdb
"""

import os
import sys
import shutil
import re
import argparse

import path
import libcommon
import init_refine
import hybrid
import locPREFMD
import define_topology
import equil
import prod
import score
import average
import scwrl
import qa


#------------------------------------------------------------------------------
# Commandline arguments.                                                      -
#------------------------------------------------------------------------------

arg = argparse.ArgumentParser(prog='PREFMD2')
arg.add_argument('-t', dest='title',  help='job title')
arg.add_argument('-i', '--input', dest='input_pdb', help='input PDB file')
arg.add_argument('-d', '--dir', dest='work_dir',
                 help='working directory (default=./)')
# arg.add_argument('--keep', dest='keep', action='store_true', default=False,
#                  help='set temporary file mode (default=False)')
arg.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                 default=False, help='set verbose mode')
arg.add_argument('--hybrid', dest='use_hybrid', action='store_true',
                 help='use the multiple initial models mode')
arg.add_argument('--extensive', dest='use_extensive', action='store_true',
                 help='use extensive sampling')
# arg.add_argument('--membrane', dest='is_membrane_protein',
#                  action='store_true')
# arg.add_argument('--ligand', dest='has_ligand', action='store_true')
# arg.add_argument('--oligomer', dest='is_oligomer', action='store_true')
arg.add_argument('--cpus', dest='n_cpus', type=int,
                 help=('number of CPUs to use for tasks running on CPUs'
                       ' (default=4)'))
arg.add_argument('--gpus', dest='gpu_ids', type=libcommon.parse_gpu_arg,
                 help=('ids of the GPUs to use in the job (default=0).'
                       ' Examples: 0 (only GPU 0 will be used), 0:1 (use GPU 0'
                       ' and 1), 0:1:3 (use GPU 0, 1 and 3)'))
arg.add_argument('--force', action='store_true', default=False)
arg.add_argument('--check', action='store_true', default=False,
                 help=('only perform a check on the PREFMD2 dependencies and'
                       ' exit'),)
arg.add_argument('--stage', type=str, default="all",
                 help=('stage of the refinement pipeline ("all" will perform'
                       ' the whole pipeline)'),
                 choices=["all", "init", "hybrid", "locprefmd", "topology",
                          "equilibration", "production", "scoring",
                          "averaging", "modeling", "evaluation"])
arg.add_argument('-j', '--json', dest='json_fp',
                 help='JSON file of an already initialized job')
arg.add_argument('--label', type=str, default='new_job_label',
                 help=('optional label to identify jobs (useful for keeping'
                       ' track of multiple jobs)'))


#------------------------------------------------------------------------------
# Parse and checks the commandline and the environment.                       -
#------------------------------------------------------------------------------

script_dirpath = os.path.dirname(os.path.realpath(__file__))


arg = arg.parse_args()
libcommon.update_verbosity_env(arg.verbose)

if os.getenv("PREFMD2_TEST") == "1":
    print("* Test mode is active")

if arg.check:
    libcommon.check_dependencies(check_hybrid=arg.use_hybrid)
    print("- All dependencies are satisfied and the PREFMD2 environment is"
          " ready.")
    sys.exit(0)


# --title and --input must be provided when starting a new job.
if arg.stage in ("all", "init"):

    for argname, argdest in (("title", "title"), ("input", "input_pdb")):
        if getattr(arg, argdest) is None:
            raise ValueError("A --%s for the refinement job must be provided"
                             " when --stage is 'all' or 'init' (that is, when"
                             " beginning a new refinement job)." % argname)

    if re.search("[^a-zA-Z\d_\-\.]", arg.title):
        raise ValueError("The --title for the refinement job must contain"
                         " only a-z, 0-9 and '-', '_', '.' characters.")

    # Set the default CPUs and GPUs.
    arg.gpu_ids = arg.gpu_ids if arg.gpu_ids is not None else ["0"]
    arg.n_cpus = arg.n_cpus if arg.n_cpus is not None else 8

    # Set the working directory.
    arg.work_dir = arg.work_dir if arg.work_dir is not None else "./"

    # Set the PDB file of the protein to refine.
    if arg.input_pdb is not None:
        if not os.path.isfile(arg.input_pdb):
            raise FileNotFoundError(" Input PDB file %s was not found."
                                    " Cannot start a job." % arg.input_pdb)
        arg.input_pdb = path.Path(arg.input_pdb)


# Some arguments can be provided only when starting a new job, not when
# resuming a previous one.
else:
    argerror_msg = ("The --%s argument can only be provided when --stage is"
                    " 'all' or 'init' (that is, when beginning a new"
                    " refinement job).")

    # Check non-flag arguments.
    for argname, argdest in [("title", "title"), ("input", "input_pdb"),
                             ("dir", "work_dir"),
                             ("cpus", "n_cpus"), ("gpus", "gpu_ids")]:
        if getattr(arg, argdest) is not None:
            raise ValueError(argerror_msg % argname)

    # Check flags.
    for argname, argdest in [("hybrid", "use_hybrid"),
                             ("extensive", "use_extensive"),
                             # ("membrane", "is_membrane_protein"),
                             # ("ligand", "has_ligand"),
                             # ("oligomer", "is_oligomer"),
                             # ("keep", "keep"),
                            ]:
        if getattr(arg, argdest):
            raise ValueError(argerror_msg % argname)


    # When resuming a previous job, a JSON file must be provided.
    if arg.json_fp is None:
        raise ValueError("When resuming a previous job, a --json argument that"
                         " points at the JSON file (job.json) in its directory"
                         " must be specified.")


#------------------------------------------------------------------------------
# Init.                                                                       -
#------------------------------------------------------------------------------

if arg.stage in ("all", "init"):  # Start a new job.

    # Make some checks for PREFMD2 dependencies.
    libcommon.check_dependencies(check_hybrid=arg.use_hybrid)

    # Checks if an output directory is already present.
    project_dp = os.path.join(arg.work_dir, arg.title)
    if os.path.isdir(project_dp):
        if not arg.force:
            raise FileExistsError("A refinement job directory named '%s'"
                " already exists at %s. If you want to overwrite it, use the"
                " --force flag." % (arg.title, arg.work_dir))
        else:
            shutil.rmtree(project_dp)


    # Check for an input PDB, clean it and set up a new refinement job
    # directory.
    print("\n# Starting the 'init' stage.")

    job = init_refine.prep(arg)

    print("- Completed the 'init' stage.")

    if arg.stage == "init":
        sys.exit(0)


else:  # Resume a previous new job.

    print("\n# Resuming a previous job...")

    if not os.path.isfile(arg.json_fp):
        raise FileNotFoundError("%s file was not found."
                                " Cannot resume a job." % arg.json_fp)

    # Make some checks for PREFMD2 dependencies.
    libcommon.check_dependencies(check_hybrid=job.use_hybrid)

    # Actually load a previous job.
    input_json = path.Path(arg.json_fp)
    job = libcommon.Job.from_json(input_json)

    print("- Completed resuming a previous job.")


#------------------------------------------------------------------------------
# Hybridize.                                                                  -
#------------------------------------------------------------------------------

if arg.stage in ("all", "hybrid"):

    if job.use_hybrid:

        print("\n# Starting the 'hybrid' stage.")

        # Prepare the task.
        if not job.hybrid_init:
            hybrid.prep(job=job, input_pdb=job.init_pdb[0])
        # Run the task.
        hybrid.run(job=job)
        # Get the output.
        hybrid_out = libcommon.get_outputs(job, method='hybrid',
                                           expand='model_s',
                                           force_complete=True)
        job.init_pdb.extend(hybrid_out[0][0])
        job.init_pdb = job.init_pdb[:5]
        job.to_json()

        print("- Completed the 'hybrid' stage.")

    else:
        if arg.stage == "hybrid":
            raise ValueError("Please initialize a job with the --hybrid flag"
                             " if you want to run the script with --stage set"
                             " as 'hybrid'.")

    if arg.stage == "hybrid":
        sys.exit(0)

# Number of initial models. When not performing hybridization is 1.
n_init = len(job.init_pdb)


#------------------------------------------------------------------------------
# locPREFMD.                                                                  -
#------------------------------------------------------------------------------

if arg.stage in ("all", "locprefmd"):

    print("\n# Starting the 'locprefmd' stage.")

    # Prepare the task.
    if not job.locprefmd_init:
        locPREFMD.prep(job=job, input_pdbs=job.init_pdb, stage="locprefmd")

    # Run the task.
    locPREFMD.run(job=job, stage="locprefmd")

    print("- Completed the 'locprefmd' stage.")

    if arg.stage == "locprefmd":
        sys.exit(0)

locPREFMD_out = libcommon.get_outputs(job, method="locPREFMD",
                                      stage="locprefmd",
                                      force_complete=True)[:n_init]


#------------------------------------------------------------------------------
# Define topology.                                                            -
#------------------------------------------------------------------------------

if arg.stage in ("all", "topology"):

    if job.topology_init:
        raise ValueError("Topology was already defined. Can not define it"
                         " again.")

    print("\n# Starting the 'topology' stage.")

    define_topology.prep(job, locPREFMD_out[0][0])

    print("- Completed the 'topology' stage.")

    if arg.stage == "topology":
        sys.exit(0)

if not job.topology_init:
    raise ValueError("Topology was not defined. Can not continue.")


#------------------------------------------------------------------------------
# Equilibration phase.                                                        -
#------------------------------------------------------------------------------

if arg.stage in ("all", "equilibration"):

    print("\n# Starting the 'equilibration' stage.")

    # Prepare the task.
    if not job.equilibration_init:

         # Soluble protein.
        equil_fn = "equil.json" if os.getenv("PREFMD2_TEST") != "1" else "equil.test.json"

        equil.prep(job=job, equil_index=0,
                   input_pdb=[out[0] for out in locPREFMD_out],
                   input_json=path.Path("%s/%s"%(libcommon.DEFAULT_HOME,
                                                 equil_fn)))

    # Actually run the task.
    equil.run(job=job)

    print("- Completed the 'equilibration' stage.")

    if arg.stage == "equilibration":
        sys.exit(0)


#------------------------------------------------------------------------------
# Production phase.                                                           -
#------------------------------------------------------------------------------

if arg.stage in ("all", "production"):

    print("\n# Starting the 'production' stage.")

    # Prepare the production task.
    if not job.production_init:

        # Soluble proteins.
        if job.use_extensive:
            if n_init == 1:
                prod_json = "prod_ext.json"
            else:
                prod_json = "prod_hybrid_ext.json"
            n_traj = 10
        else:
            prod_json = "prod.json" if os.getenv("PREFMD2_TEST") != "1" else "prod.test.json"
            n_traj = 5

        prod_input = path.Path("%s/%s"%(libcommon.DEFAULT_HOME, prod_json))

        for i in range(n_init):
            prod.prep(job=job, prod_index=i, input_equil=i,
                      input_json=prod_input, n_replica=n_traj)

    # Run the production run.
    prod.run(job=job)

    print("- Completed the 'production' stage.")

    if arg.stage == "production":
        sys.exit(0)

prod_out = libcommon.get_outputs(job, method="prod", force_complete=True)


#------------------------------------------------------------------------------
# Scoring.                                                                    -
#------------------------------------------------------------------------------

if arg.stage in ("all", "scoring"):

    print("\n# Starting the 'scoring' stage.")

    # Prepare the scoring task.
    if not job.scoring_init:
        score.prep(job=job, input_dcd=[out[0] for out in prod_out])

    # Actually perform the scoring.
    score.run(job=job)

    print("- Completed the 'scoring' stage.")

    if arg.stage == "scoring":
        sys.exit(0)


#------------------------------------------------------------------------------
# Averaging.                                                                  -
#------------------------------------------------------------------------------

if arg.stage in ("all", "averaging"):

    print("\n# Starting the 'averaging' stage")

    # Single initial model mode.
    if n_init == 1:

        if not job.averaging_init:
            average.prep(job=job, output_prefix='%s.prod_0'%job.title,
                         input_prod=[0],
                         input_json=path.Path("%s/average.json"%(
                                              libcommon.DEFAULT_HOME)),
                         rule='score')
            average.prep(job=job, output_prefix='%s.prod_0.cluster'%job.title,
                         input_prod=[0],
                         input_json=path.Path("%s/average.json"%(
                                              libcommon.DEFAULT_HOME)),
                         rule='cluster')

        average.run(job=job, n_init=1)

    # Multiple initial models mode.
    else:

        average.prep(job=job, output_prefix=job.title,
                     input_prod=[i for i in range(n_init)],
                     input_json=path.Path("%s/average.json"%(
                                          libcommon.DEFAULT_HOME)),
                     rule='casp12', idx=0)
        for i in range(n_init):
            average.prep(job=job, output_prefix='%s.prod_%d'%(job.title, i),
                         input_prod=[i],
                         input_json=path.Path("%s/average.json"%(
                                              libcommon.DEFAULT_HOME)),
                         rule='score', idx=i+1)

        average.run(job=job, n_init=n_init)

    print("- Completed the 'averaging' stage.")

    if arg.stage == "averaging":
        sys.exit(0)

average_out = []
for output in libcommon.get_outputs(job, 'average', expand='pdb_s',
                                    force_complete=True):
    average_out.extend(output[0])
average_out = average_out[:libcommon.N_MODEL]


#------------------------------------------------------------------------------
# Build the final models.                                                     -
#------------------------------------------------------------------------------

if arg.stage in ("all", "modeling"):

    print("\n# Starting the 'modeling' stage.")

    # Copies the files of the averaged models.
    print("- Copying averaged models files...")
    job.work_home.chdir()
    model_home = job.work_home.subdir("model", build=True)
    prep_s = []
    for i, out in enumerate(average_out):
        prep_fn = model_home.fn("prep_%d.pdb"%(i+1))
        shutil.copy(out.short(), prep_fn.short())
        prep_s.append(prep_fn)


    # Runs scwrl.
    print("- Repacking side chains with scrwl4...")
    if not job.scwrl_init:
        scwrl.prep(job, prep_s)
    scwrl.run(job)
    scwrl_out = libcommon.get_outputs(job, 'scwrl', force_complete=True)


    # Runs locPREFMD.
    print("- Refining with locPREFMD...")
    if not job.locprefmd_mod_init:
        locPREFMD.prep(job=job, input_pdbs=[out[0] for out in scwrl_out],
                       stage="modeling")
    locPREFMD.run(job=job, stage="modeling")

else:
    job.work_home.chdir()
    model_home = job.work_home.subdir("model", build=True)


# Copy the final models in the 'final' directory.
locPREFMD_out = libcommon.get_outputs(job, "locPREFMD", stage="modeling",
                                      force_complete=True)

model_s = []
for i,out in enumerate(locPREFMD_out):
    model_fn = model_home.fn("model_%d.pdb"%(i+1))
    if arg.stage in ("all", "modeling"):
        shutil.copy(out[0].short(), model_fn.short())

    model_s.append(model_fn)


if arg.stage in ("all", "modeling"):
    print("- Completed the 'modeling' stage.")

    if arg.stage == "modeling":
        sys.exit(0)


#------------------------------------------------------------------------------
# Quality estimation.                                                         -
#------------------------------------------------------------------------------

if arg.stage in ("all", "evaluation"):

    print("\n# Starting the 'evaluation' stage.")

    if not job.qa_init:
        qa_fn = "qa.json" if os.getenv("PREFMD2_TEST") != "1" else "qa.test.json"
        qa.prep(job=job, input_pdb=model_s,
                input_json=path.Path("%s/%s"%(libcommon.DEFAULT_HOME, qa_fn)))
    qa.run(job)

    qa_out = libcommon.get_outputs(job, 'qa')

    final_home = job.work_home.subdir("final", build=True)
    for i,out in enumerate(qa_out):
        pdb_fn = final_home.fn("model_%d.pdb"%(i+1))
        shutil.copy(out[0].short(), pdb_fn.short())

    print("- Completed the 'evaluation' stage.")
    print("- Final models written inside %s/" % final_home.path())
