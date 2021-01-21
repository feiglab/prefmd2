"""
Module for initializing and starting the scoring task in the refinement
pipeline. Its functions will be called in the main locprefmd2.py script.
Extract structures from the MD trajectories using the exec_pdb_extract.py
script, scores with a statistical potential (RWplus) the frames using the
exec_calc_statpot.py script and finally assess the level of structural
divergence between the initial model and the frames using the TMscore program
in the exec_calc_tmscore.py script.
"""

import os
import libcommon


METHOD = 'score'


def prep(job, input_dcd):
    #
    for idx, dcd_fn in enumerate(input_dcd):
        run_home = dcd_fn.dirname()
        input_s = [dcd_fn]
        output_s = [run_home.fn("statpot.dat"), run_home.fn("qual_init.dat")]
        log_s = [run_home.fn("%s.%s.status" % (METHOD, idx))]
        job.add_task(method=METHOD,
                     input_arg=input_s, output_arg=output_s, log_arg=log_s,
                     use_gpu=False, n_proc=job.n_cpus)
    #
    job.scoring_init = True
    job.to_json()


def run(job, ids=None):

    task_s = libcommon.get_tasks_to_run(job, method=METHOD)
    libcommon.print_tasks_to_run_message(task_s, "scoring")

    #
    for index, task in task_s:
        #
        print("- Running scoring task %s..." % (index+1))
        # Set up the input parameters of the task.
        input_dcd = task['input'][0]
        run_home = input_dcd.dirname()
        output_score = task['output'][0]
        output_qual = task['output'][1]
        log_p = task['log'][0]
        #
        if len(task['etc']) > 0:    # meta
            # index_fn = task['etc'][0]
            # top_fn = task['etc'][1]
            raise NotImplementedError("'meta' scoring.")
        else:
            index_fn = None
            top_fn = job.top_fn
        #
        run_home.chdir()
        #
        pdblist = run_home.fn("pdb_s")

        # Extract PDB files from the MD trajectories.
        with pdblist.open("wt") as fout:
            cmd = [top_fn.short()]
            if index_fn is not None:
                cmd.extend(['--topIndex', index_fn.short()])
            cmd.extend(['--dir', 'ens'])
            cmd.extend(['--name', 'sample'])
            cmd.append("--structured")
            cmd.extend(['--dcd', input_dcd.short()])
            libcommon.asystem(module="exec_pdb_extract", args=cmd,
                              outfile=fout)

        # Compute statistical potential scores.
        if os.getenv("PREFMD2_PYTHON_MPROCS") is None:
            n_jobs = task['_resource']['n_proc']
        else:
            n_jobs = int(os.getenv("PREFMD2_PYTHON_MPROCS"))
        cmd = []
        cmd.extend(['-j', '%d' % n_jobs])
        cmd.extend(['-l', pdblist.short()])
        cmd.extend(['-o', output_score.short()])
        cmd.append('--rwplus')
        # cmd.append('--dfire')
        libcommon.asystem(module="exec_calc_statpot", args=cmd)

        # Compute structural similarity with the initial model with TMscore.
        cmd = []
        cmd.extend(['-j', '%d' % n_jobs])
        cmd.extend(['-r', job.init_pdb[0].short()])
        cmd.extend(['-l', pdblist.short()])
        cmd.extend(['-o', output_qual.short()])
        libcommon.asystem(module="exec_calc_tmscore", args=cmd)

        with log_p.open("w") as lfh:
            lfh.write("DONE\n")
