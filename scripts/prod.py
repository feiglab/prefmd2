"""
Module for initializing and starting the production MD runs in the refinement
pipeline. Its functions will be called in the main locprefmd2.py script or may
be executed from the subprocess_launcher.py script.
Launches functions from the exec_prod.py script, which actually contain the
code to run MD production.
"""


import os
import json
import argparse

import path
import libcommon


METHOD = 'prod'


def prep(job, prod_index, input_equil, input_json, n_replica):
    """
    Set up the production tasks.
    """
    job.prod_home = job.work_home.subdir("prod", build=True)
    job.prod_home.chdir()
    #
    iter_home = job.prod_home.subdir("%d"%prod_index, build=True)
    options = {}
    options['input_equil'] = input_equil
    options['input_json'] = input_json.path()
    options['n_replica'] = n_replica
    with iter_home.fn("input.json").open("wt") as fout:
        fout.write(json.dumps(options, indent=2))
    #
    for i in range(n_replica):
        run_home = iter_home.subdir("%d"%i, build=True)
        #
        input_s = [run_home, input_equil, input_json]
        output_s = [run_home.fn("solute.dcd")]
        log_s = [run_home.fn("%s.%s.status" % (METHOD, i))]
        job.add_task(method=METHOD,
                     input_arg=input_s, output_arg=output_s, log_arg=log_s,
                     use_gpu=True, idx=i, n_proc=job.n_cpus)
    #
    job.production_init = True
    job.to_json()


def run(job):
    """
    Run all production tasks.
    """

    task_s = libcommon.get_tasks_to_run(job, method=METHOD)
    libcommon.print_tasks_to_run_message(task_s, "production")

    if len(job.gpu_ids) == 1 or len(task_s) == 1 or \
       os.getenv("PREFMD2_OPENMM_PLATFORM") == "CPU":
        run_tasks(job, task_s)
    else:
        libcommon.run_parallel_tasks(job, task_s, METHOD)


def run_tasks(job, task_s):
    """
    Run a subset of production tasks.
    Executes the exec_prod.py script, where the code for the actual MD
    production run is.
    """

    for index, task in task_s:

        print("- Running production task %s..." % (index+1))

        run_home = task['input'][0]
        input_equil_index  = task['input'][1]
        input_json = task['input'][2]
        #
        equil_home = job.get_task("equil")[input_equil_index][1]['input'][0]
        #
        output_s = task['output']
        log_p = task['log'][0]
        #
        with input_json.open() as fp:
            options = json.load(fp)
        #
        run_home.build()
        run_home.chdir()
        run_name = 'r%s'%(run_home.split("/")[-1])
        #
        if 'restraint' in options:
            options['restraint']['reference'] = equil_home.fn("%s.orient.pdb"%job.title).short()
        options['input'] = {}
        options['input']['psf'] = equil_home.fn("%s.psf"%job.title).short()
        options['input']['pdb'] = equil_home.fn("%s.equil.pdb"%job.title).short()
        options['input']['n_atom'] = job.n_atom
        if options['restart']:
            options['input']['restart'] = equil_home.fn("%s.equil.restart.pkl"%job.title).short()
        if job["has_ligand"]:
            # options['ligand_json'] = job.ligand_json.short()
            raise NotImplementedError("Ligands.")
        #
        options['openmm']['platform'] = os.getenv("PREFMD2_OPENMM_PLATFORM", "CUDA")
        options['gpu_id'] = task['_resource']['gpu_id']
        #
        options['ff']['toppar'] = libcommon.complete_data_paths(
                                      options['ff']['toppar'])
        options['ff']['cgenff'] = libcommon.complete_data_paths(
                                      options['ff']['cgenff'])
        #
        run_json = run_home.fn("input.json")
        with run_json.open("wt") as fout:
            fout.write(json.dumps(options, indent=2))
        #
        cmd = [run_name]
        cmd.extend(["--input", run_json.short()])
        # if job.verbose:  cmd.append('--verbose')
        # if job.keep_tmp: cmd.append('--keep')
        #
        libcommon.asystem(module="exec_prod", args=cmd)

        with log_p.open("w") as lfh:
            lfh.write("DONE\n")
