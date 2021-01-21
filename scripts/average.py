"""
Module for initializing and starting the averaging procedure in the
refinement pipeline. Its functions will be called in the main locprefmd2.py
script or may be executed from the subprocess_launcher.py script.
Calls functions from the exec_average.py script, which actually contains the
code to run the averaging process.
"""

import os
import json
import numpy as np

import path
import libcommon

METHOD = 'average'

PARAM = {}
PARAM['score'] = ("RWplus", 25.0)
PARAM['casp12'] = ("RWplus", 0.5, 225., 45.)
PARAM['cluster'] = (2.0, 20, 5) # rmsd_cutoff, subsample, max_run_md
PARAM['msm'] = ('RWplus', 25.0)

def prep(job, output_prefix, input_prod, input_json, rule='score', idx=0):
    #
    prod_s = job.get_task("prod")
    prod_s = [prod for _, prod in prod_s if int(prod['input'][0].split("/")[-2]) in input_prod]
    #
    score_s = job.get_task("score")
    #
    job.average_home = job.work_home.subdir("average", build=True)
    job.average_home.chdir()
    #
    if rule == 'score':
        input_s = [output_prefix, (rule, PARAM[rule]), input_json, [], []]
    elif rule == 'casp12':
        input_s = [output_prefix, (rule, PARAM[rule]), input_json, [], [], []]
    elif rule == 'cluster':
        input_s = [output_prefix, (rule, PARAM[rule]), input_json, []]
    #

    for prod in prod_s:

        if rule in ['score', 'casp12']:
            score = None
            for _,s in score_s:
                if s['input'][0] == prod['output'][0]:
                    score = s
                    break
            if score is None:
                raise Exception("Could not get the score.")
        #
        input_s[3].append(prod['output'][0])
        if rule in ['score', 'casp12']:
            input_s[4].append(score['output'][0])
        if rule in ['casp12']:
            input_s[5].append(score['output'][1])

    #
    output_s = [job.average_home.fn("%s.pdb_s"%output_prefix)]
    log_s = [job.average_home.fn("%s.status" % (METHOD))]
    job.add_task(method=METHOD,
                 input_arg=input_s, output_arg=output_s, log_arg=log_s,
                 use_gpu=True, idx=idx, n_proc=job.n_cpus)
    #
    job.averaging_init = True
    job.to_json()


# def prep_from_msm(job, output_prefix, msm_fn, input_json):
#     job.average_home = job.work_home.subdir("average", build=True)
#     job.average_home.chdir()
#     #
#     rule = 'msm'
#     input_s = [output_prefix, (rule, PARAM[rule]), input_json, [], [], msm_fn]
#     #
#     msm = np.load(msm_fn.short())
#     traj_fn_s = [path.Path(fn) for fn in msm['traj_s']]
#     #
#     score_s = job.get_task("score")
#     for traj_fn in traj_fn_s:
#         score = None
#         for _,s in score_s:
#             if s['input'][0] == traj_fn:
#                 score = s
#                 break
#         if score is None:
#             return
#         if score['resource'][0] != 'DONE':
#             return
#         if not score['output'][0].status():
#             return
#         #
#         input_s[3].append(traj_fn)
#         input_s[4].append(score['output'][0])
#     #
#     output_s = [job.average_home.fn("%s.pdb_s"%output_prefix)]
#     job.add_task(METHOD, input_s, output_s, use_gpu=True, n_proc=1)
#     #
#     job.averaging_init = True
#     job.to_json()


def run(job, n_init):

    task_s = libcommon.get_tasks_to_run(job, method=METHOD)
    libcommon.print_tasks_to_run_message(task_s, "averaging")

    if len(job.gpu_ids) == 1 or n_init == 1 or \
       os.getenv("PREFMD2_OPENMM_PLATFORM") == "CPU":
        run_tasks(job, task_s)
    else:
        libcommon.run_parallel_tasks(job, task_s, METHOD)


def run_tasks(job, task_s):

    for index, task in task_s:
        #
        print("- Running averaging task %s..." % (index+1))
        #
        input_s = task['input']
        output_prefix = input_s[0]
        rule = input_s[1]
        input_json = input_s[2]
        input_dcd_s = input_s[3]
        output_pdb = task['output'][0]
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
        options['rule'] = rule
        #
        options['openmm']['platform'] = os.getenv("PREFMD2_OPENMM_PLATFORM", "CUDA")
        options['gpu_id'] = task['_resource']['gpu_id']
        #
        #
        options['ff']['toppar'] = libcommon.get_c36m_paths()
        options['ff']['cgenff'] = None
        #
        job.average_home.chdir()
        #
        input_json = job.average_home.fn("%s.json"%output_prefix)
        with input_json.open("wt") as fout:
            fout.write(json.dumps(options, indent=2,
                                  default=libcommon.JSONserialize))
        #
        cmd = [output_prefix, job.top_fn.short()]
        cmd.extend(['--input', input_json.short()])
        if rule[0] == 'msm':
            cmd.extend(['--msm', input_s[5].short()])
        cmd.append('--dcd')
        cmd.extend([fn.short() for fn in input_dcd_s])
        if rule[0] in ['score', 'casp12', 'msm']:
            cmd.append('--score')
            cmd.extend([fn.short() for fn in input_s[4]])
        if rule[0] in ['casp12']:
            cmd.append('--qual')
            cmd.extend([fn.short() for fn in input_s[5]])
        #
        libcommon.asystem(module="exec_average", args=cmd)

        with log_p.open("w") as lfh:
            lfh.write("DONE\n")
