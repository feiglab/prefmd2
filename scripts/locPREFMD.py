"""
Module with functions for running locPREFMD tasks in the refinement pipeline.
"""

import os

import libcommon


METHOD = 'locPREFMD'
LOCPREFMD = os.getenv("LOCPREFMD")
if LOCPREFMD is None:
    LOCPREFMD = ""
EXEC =  os.path.join(LOCPREFMD, 'scripts', 'locprefmd.sh')


def prep(job, input_pdbs, stage):
    #
    for idx, fn in enumerate(input_pdbs):
        out = fn.dirname().fn("%s.prefmd.pdb"%(fn.name()))
        log = job.init_home.fn("%s.%s.%s.status" % (METHOD, stage, idx))
        #
        job.add_task(method=METHOD, stage=stage,
                     input_arg=[fn], output_arg=[out], log_arg=[log],
                     use_gpu=False, n_proc=job.n_cpus)
    #
    if stage == "locprefmd":
        job.locprefmd_init = True
    elif stage == "modeling":
        job.locprefmd_mod_init = True
    else:
        raise KeyError("Unknown stage: %s" % stage)
    job.to_json()


def run(job, stage, ids=None):

    task_s = libcommon.get_tasks_to_run(job, method=METHOD, stage=stage)

    # Begin to run locPREFMD on input models.
    for _, task in task_s:

        input_pdb = task['input'][0]
        run_home = input_pdb.dirname()
        output_pdb = task['output'][0]
        log_p = task['log'][0]

        print("- Running locPREFMD for local refinement on %s ..." % (
              os.path.basename(input_pdb.path())))

        # Actually executes locPREFMD script.
        run_home.chdir()
        with output_pdb.open("wt") as fout:
            cmd = [EXEC, input_pdb.short(),
                   "--cpus", str(task['_resource']['n_proc'])]
            libcommon.system(cmd, outfile=fout, errfile=libcommon.get_stderr())

        # Writes a file to record the status of the job.
        with log_p.open("w") as lfh:
            lfh.write("DONE\n")
