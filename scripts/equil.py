"""
Module for initializing and starting the equilibration MD runs in the
refinement pipeline. Its functions will be called in the main locprefmd2.py
script or may be executed from the subprocess_launcher.py script.
Launches functions from the exec_equil.py script, which actually contain the
code to run MD equilibration.
"""

import os
import json

import libcommon


METHOD = 'equil'


def prep(job, equil_index, input_pdb, input_json):
    """
    Set up the equilibration tasks.
    """

    #
    OUTs = ['%s.psf'%job.title, '%s.orient.pdb'%job.title,
            '%s.equil.restart.pkl'%job.title, '%s.equil.pdb'%job.title]
    #
    job.equil_home = job.work_home.subdir("equil", build=True)
    job.equil_home.chdir()
    #
    for idx, fn in enumerate(input_pdb):
        run_home = job.equil_home.subdir("%d"%equil_index, build=True)
        #
        input_s = [run_home, fn, input_json]
        output_s = [run_home.fn(X) for X in OUTs]
        log_s = [run_home.fn("%s.%s.status" % (METHOD, idx))]

        job.add_task(method=METHOD,
                     input_arg=input_s, output_arg=output_s, log_arg=log_s,
                     use_gpu=True, idx=idx, n_proc=job.n_cpus)
        equil_index += 1
    #
    job.equilibration_init = True
    job.to_json()


# def prep_membrane(job, equil_index, input_pdb, input_psf, input_crd, input_json):
#     if len(job.get_task(METHOD, not_status='DONE')) > 0:
#         return
#     #
#     OUTs = ['%s.psf'%job.title, '%s.orient.pdb'%job.title, '%s.equil.restart.pkl'%job.title, '%s.equil.pdb'%job.title]
#     #
#     job.equil_home = job.work_home.subdir("equil", build=True)
#     job.equil_home.chdir()
#     #
#     for pdb_fn, psf_fn, crd_fn in zip(input_pdb, input_psf, input_crd):
#         run_home = job.equil_home.subdir("%d"%equil_index)
#         #
#         input_s = [run_home, pdb_fn, psf_fn, crd_fn, input_json]
#         output_s = [run_home.fn(X) for X in OUTs]
#         status = True
#         for output in output_s:
#             if not output.status():
#                 status = False ; break
#         #
#         if status:
#             job.add_task(METHOD, input_s, output_s, use_gpu=True, n_proc=12, status='DONE')
#         else:
#             job.add_task(METHOD, input_s, output_s, use_gpu=True, n_proc=12)
#         equil_index += 1
#     #
#     job.to_json()


def run(job):
    """
    Run all equilibration tasks.
    """

    task_s = libcommon.get_tasks_to_run(job, method=METHOD)
    libcommon.print_tasks_to_run_message(task_s, "equilibration")

    if len(job.gpu_ids) == 1 or len(task_s) == 1 or \
       os.getenv("PREFMD2_OPENMM_PLATFORM") == "CPU":
        run_tasks(job, task_s)
    else:
        libcommon.run_parallel_tasks(job, task_s, METHOD)


def run_tasks(job, task_s):
    """
    Run a subset of equilibration tasks.
    Executes the exec_equil.py script, where the code for the actual MD
    production run is.
    """
    #
    for index, task in task_s:
        #
        print("- Running equilibration task %s..." % (index+1))
        #
        run_home = task['input'][0]
        input_pdb  = task['input'][1]
        input_json = task['input'][-1]
        #
        if job["is_membrane_protein"]:
            # is_membrane = True
            # EXEC = EXEC_MEMBRANE
            # input_psf = task['input'][2]
            # input_crd = task['input'][3]
            raise NotImplementedError("Membrane protein.")
        else:
            is_membrane = False
        #
        output_s = task['output']
        log_p = task['log'][0]
        #
        with input_json.open() as fp:
            options = json.load(fp)
        options['ssbond'] = []
        for line in job.ssbond:
            chain_1 = line[15]
            chain_2 = line[29]
            if chain_1 == ' ' and chain_2 == ' ':
                line = '%sA%sA%s'%(line[:15], line[16:29], line[30:])
            options['ssbond'].append(line)
        #
        options['input_pdb'] = input_pdb
        options['input_json'] = input_json
        if job["has_ligand"]:
            # options['ligand_json'] = job.ligand_json
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
        run_home.build()
        run_home.chdir()
        #
        equil_json = run_home.fn("input.json")
        with equil_json.open("wt") as fout:
            fout.write(json.dumps(options, indent=2,
                       default=libcommon.JSONserialize))
        #
        if not is_membrane:
            cmd = [job.title, input_pdb.short()]
        else:
            # cmd = [job.title, input_pdb.short(), input_psf.short(), input_crd.short()]
            raise NotImplementedError("Membrane protein.")
        cmd.extend(['--input', equil_json.short()])
        # if job.verbose:  cmd.append('--verbose')
        # if job.keep_tmp: cmd.append('--keep')
        #
        libcommon.asystem(module="exec_equil", args=cmd)

        with log_p.open("w") as lfh:
            lfh.write("DONE\n")
