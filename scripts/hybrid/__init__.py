import libcommon


METHOD = 'hybrid'


def prep(job, input_pdb):
    #
    job.hybrid_home = job.work_home.subdir("hybrid", build=True)
    out = job.hybrid_home.fn("model_s")
    log_s = [job.hybrid_home.fn("%s.status" % METHOD)]
    #
    job.add_task(method=METHOD,
                 input_arg=[job.title, input_pdb, job.hybrid_home],
                 output_arg=[out], log_arg=log_s,
                 use_gpu=False, n_proc=job.n_cpus)
    #
    job.hybrid_init = True
    job.to_json()


def run(job):
    """
    Run an hybridization task. This function will call the main() function of
    the hybrid.exec_hybrid module.
    """

    task_s = libcommon.get_tasks_to_run(job, method=METHOD)
    libcommon.print_tasks_to_run_message(task_s, "hybridization")

    #
    for index, task in task_s:
        #
        print("- Running hybridization task %s..." % (index+1))
        #
        title = task['input'][0]
        input_pdb = task['input'][1]
        run_home = task['input'][2]
        output_list = task['output'][0]
        log_p = task['log'][0]
        #
        n_proc = task['_resource']['n_proc']
        #
        run_home.chdir()
        cmd = [title, input_pdb.short(), '-j', '%d'%n_proc]
        libcommon.asystem(module="hybrid.exec_hybrid", args=cmd)

        with log_p.open("w") as lfh:
            lfh.write("DONE\n")
