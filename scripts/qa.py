"""
Module for initializing and starting the model quality assessment tasks in the
refinement pipeline. Its functions will be called in the main locprefmd2.py
script or may be executed from the subprocess_launcher.py script.
Launches functions from the exec_local_qa.py script, which actually contain the
code to run the assessment.
"""

import os
import json

import libcommon


METHOD = 'qa'


def get_ssbond(pdb_fn):
    ssbond = []
    with pdb_fn.open() as fp:
        for line in fp:
            if line.startswith("SSBOND"):
                ssbond.append(line.rstrip())
    return ssbond


def prep(job, input_pdb, input_json):
    #
    job.qa_home = job.work_home.subdir("qa", build=True)
    job.qa_home.chdir()
    #
    for i, fn in enumerate(input_pdb):
        name = fn.name()
        run_home = job.qa_home.subdir(name, build=True)
        #
        input_s = [run_home, fn, input_json]
        output_s = [run_home.fn('%s.qa.pdb'%name)]
        log_s = [run_home.fn("%s.%s.status"%(METHOD, i))]
        #
        job.add_task(method=METHOD,
                     input_arg=input_s, output_arg=output_s, log_arg=log_s,
                     use_gpu=True, idx=i, n_proc=job.n_cpus)
    #
    job.qa_init = True
    job.to_json()


def run(job):

    task_s = libcommon.get_tasks_to_run(job, method=METHOD)
    libcommon.print_tasks_to_run_message(task_s, "quality assessment")

    if len(job.gpu_ids) == 1 or len(task_s) == 1 or \
       os.getenv("PREFMD2_OPENMM_PLATFORM") == "CPU":
        run_tasks(job, task_s)
    else:
        libcommon.run_parallel_tasks(job, task_s, METHOD)


def run_tasks(job, task_s):

    for index, task in task_s:
        #
        print("- Running quality assessment task %s..." % (index+1))
        #
        run_home = task['input'][0]
        input_pdb  = task['input'][1]
        input_json = task['input'][2]
        #
        output_s = task['output']
        log_p = task['log'][0]
        #
        with input_json.open() as fp:
            options = json.load(fp)
        options['ssbond'] = []
        for line in get_ssbond(input_pdb):
            chain_1 = line[15]
            chain_2 = line[29]
            if chain_1 == ' ' and chain_2 == ' ':
                line = '%sA%sA%s'%(line[:15], line[16:29], line[30:])
            options['ssbond'].append(line)
        #
        run_home.build()
        run_home.chdir()
        #
        options['input_pdb'] = input_pdb.short()
        options['input_json'] = input_json.short()
        #
        options['openmm']['platform'] = os.getenv("PREFMD2_OPENMM_PLATFORM", "CUDA")
        options['gpu_id'] = task['_resource']['gpu_id']
        #
        options['ff']['toppar'] = libcommon.get_c36m_paths()
        options['ff']['cgenff'] = None
        #
        run_json = run_home.fn("input.json")
        with run_json.open("wt") as fout:
            fout.write(json.dumps(options, indent=2,
                       default=libcommon.JSONserialize))
        #
        cmd = [input_pdb.name()]
        cmd.extend(['--input', run_json.short()])
        # if job.verbose:  cmd.append('--verbose')
        # if job.keep_tmp: cmd.append('--keep')
        #
        libcommon.asystem(module="exec_local_qa", args=cmd)

        with log_p.open("w") as lfh:
            lfh.write("DONE\n")
