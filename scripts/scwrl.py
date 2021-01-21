import os
import libcommon

METHOD = 'scwrl'
EXEC = 'scwrl4'

def prep(job, input_pdb):
    #
    for i, pdb_fn in enumerate(input_pdb):
        out_fn = pdb_fn.dirname().fn("%s.scwrl.pdb"%(pdb_fn.name()))
        #
        input_s = [pdb_fn]
        output_s = [out_fn]
        log_s = [pdb_fn.dirname().fn("%s.%s.status"%(METHOD, i))]
        job.add_task(method=METHOD,
                     input_arg=input_s, output_arg=output_s, log_arg=log_s,
                     use_gpu=False, n_proc=1)
    #
    job.scwrl_init = True
    job.to_json()

def run(job, ids=None):

    task_s = libcommon.get_tasks_to_run(job, method=METHOD)

    #
    for index, task in task_s:
        input_pdb = task['input'][0]
        run_home = input_pdb.dirname()
        output_pdb = task['output'][0]
        log_p = task['log'][0]
        #
        run_home.chdir()
        #
        cmd = [EXEC]
        cmd.extend(['-i', input_pdb.short()])
        cmd.extend(['-o', output_pdb.short()])
        #
        libcommon.system(cmd, stdout=True)
        if not os.path.isfile(output_pdb.short()):
            raise ChildProcessError("Scwrl4 could not produce an output"
                                    " model.")
        #
        with log_p.open("w") as lfh:
            lfh.write("DONE\n")
