"""
Variables, functions and classes shared among several modules of the package.
"""

import os
import time
import io
import sys
import json
import argparse
import signal
import subprocess as sp
from importlib import import_module

import path

assert sys.version_info.major == 3


#------------------------------------------------------------------------------
# Variables and paths commonly used in PREFMD2.                               -
#------------------------------------------------------------------------------

WORK_HOME = os.getenv("PREFMD2_HOME")
if WORK_HOME is None:
    WORK_HOME = ""

BIN_HOME = os.path.join(WORK_HOME, 'bin')
EXEC_HOME = os.path.join(WORK_HOME, "scripts")
DEFAULT_HOME = os.path.join(WORK_HOME, 'default')
DATA_HOME = os.path.join(WORK_HOME, 'data')

def complete_data_paths(params):
    return [os.path.join(DATA_HOME, p) for p in params]

def update_verbosity_env(verbose):
    if verbose:
        os.environ["PREFMD2_VERBOSE"] = "1"
    else:
        os.environ["PREFMD2_VERBOSE"] = "0"

def get_verbosity():
    if "PREFMD2_VERBOSE" in os.environ:
        return os.environ["PREFMD2_VERBOSE"] != "0"
    else:
        return False

def get_stderr():
    if get_verbosity():
        return None
    else:
        return "/dev/null"

def print_tasks_to_run_message(task_s, task_name):
    print("- A total of %s %s tasks will be run." % (len(task_s), task_name))


SUBMIT_MAX_PROC = 64

N_MODEL = 5
MAX_ERROR = 20

# Used in the hybridization phase to exclude templates.
TBM_EXCLUDE = None  # '%s/exclude.casp13' % DEFAULT_HOME


#------------------------------------------------------------------------------
# Executes external processes and modules.                                    -
#------------------------------------------------------------------------------

class GracefulExit(Exception):
    pass

def listen_signal():
    def gracefulExit(*arg):
        raise GracefulExit()
    signal.signal(signal.SIGTERM, gracefulExit)


def asystem(module, args, func="main",
            stdout=False, stdin=None, outfile=None, errfile=None):

    # Executes in the current Python process by importing the module specified
    # in the 'module' argument.
    if get_verbosity():
        print("- Executing module %s with parameters: %s" % (module, args))

    # Prepare the standard streams.
    redirect_stdout = stdout or outfile is not None
    if redirect_stdout:  # Redirect stdout to a StringIO object.
        original_stdout = sys.stdout
        sys.stdout = io.StringIO()
    if errfile is not None:  # Redirect the stderr to a file handle.
        original_stderr = sys.stderr
        sys.stderr = errfile

    # Actually executes the command.
    try:
        getattr(import_module(module), func)(args)
    except Exception as e:
        raise e  # The finally block is run before this.
    finally:
        if redirect_stdout:  # Revert stdout back to the original one.
            sys.stdout.seek(0)
            output = sys.stdout.read()
            sys.stdout = original_stdout
            if outfile is not None:  # Write the output to a file handle.
                outfile.write(output)
        else:
            output = None
        if errfile is not None:  # Revert back stderr to its origin.
            sys.stderr = original_stderr

    return output


def system(cmd, stdout=False, stdin=None, outfile=None, errfile=None):
    if type(cmd) == type(""):
        cmd = cmd.strip().split()
    if get_verbosity():
        print("- Executing external process: " + " ".join(cmd))  # + '\n')
    #
    if stdout or (outfile is not None):
        STDOUT = sp.PIPE
    else:
        STDOUT = None
    if errfile is None:
        STDERR = None
    elif errfile == '/dev/null':
        STDERR = sp.DEVNULL
    else:
        STDERR = errfile
    #
    listen_signal()
    #
    proc_output = None
    try:
        proc = sp.Popen(cmd, stdin=stdin, stdout=STDOUT, stderr=STDERR)
        #print (proc.pid, ' '.join(cmd))
        proc_output = proc.communicate()
    except (GracefulExit, KeyboardInterrupt) as error:
        #print ("killing... %d"%proc.pid)
        proc.terminate()
        sys.exit()
    if proc.returncode != 0:
        raise Exception("Error encountered in child process (%s): %s" % (
                        proc.returncode, proc.errors))
    #
    if (proc_output is not None) and (STDOUT is not None):
        out = proc_output[0].decode("utf8")
        if outfile is not None:
            outfile.write(out)
        return out

def JSONserialize(X):
    if isinstance(X, path.Dir):
        return 'Dir %s'%X.short()
    elif isinstance(X, path.Path):
        return 'Path %s'%X.short()
    else:
        return X

def JSONdeserialize(X):
    if isinstance(X, list):
        out = []
        for Xi in X:
            out.append(JSONdeserialize(Xi))
        return out
    elif isinstance(X, dict):
        out = {}
        for Xi in X:
            out[Xi] = JSONdeserialize(X[Xi])
        return out
    else:
        if isinstance(X, str):
            x = X.split()
            if x[0] == 'Dir':
                X = path.Dir(x[1])
            elif x[0] == 'Path':
                X = path.Path(x[1])
        return X


#------------------------------------------------------------------------------
# Classes and functions to manage refinement jobs and tasks.                  -
#------------------------------------------------------------------------------

class Job(dict):

    def __init__(self, work_dir='.', title='', build=False):

        super().__init__()
        self.title = title
        if build:
            self.work_home = path.Dir("%s/%s"%(work_dir, title), build=True)
            self.json_job = self.work_home.fn("job.json")
        self.task = {}
        self.hybrid_init = False
        self.locprefmd_init = False
        self.topology_init = False
        self.equilibration_init = False
        self.production_init = False
        self.scoring_init = False
        self.averaging_init = False
        self.scwrl_init = False
        self.locprefmd_mod_init = False
        self.qa_init = False
        self.gpu_ids = []

    def __repr__(self):
        return self.work_home.path()

    def has(self, key):
        return hasattr(self, key)

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, item):
        setattr(self, key, item)

    @classmethod
    def from_json(cls, fn):
        with fn.open() as fp:
            X = json.load(fp)
        #
        cwd = os.getcwd()
        work_home = fn.dirname()
        work_home.chdir()
        #
        job = cls()
        for key in X:
            job[key] = JSONdeserialize(X[key])
        #
        os.chdir(cwd)
        return job

    def to_json(self):
        for method in self.task:
            if len(self.task[method]) == 0:
                del (self.task[method])
        #
        cwd = os.getcwd()
        self.work_home.chdir()
        with self.json_job.open('wt') as fout:
            fout.write(json.dumps(self.__dict__, indent=2, default=JSONserialize))
        os.chdir(cwd)

    def add_task(self, method, input_arg, output_arg, log_arg, stage=None,
                 n_proc=1, use_gpu=False, idx=None, *arg):

        if method not in self.task:
            self.task[method] = []

        #
        # status allocated_resource use_gpu n_proc
        if idx is not None:
            gpu_id = idx % len(self.gpu_ids)
        else:
            gpu_id = None

        self.task[method].append({
                 '_resource': {'use_gpu': use_gpu,
                               'gpu_id': gpu_id,
                               'n_proc': min(SUBMIT_MAX_PROC, n_proc)},
                 'input': input_arg,
                 'output': output_arg,
                 'log': log_arg,
                 'stage': stage,
                 'idx': idx,
                 'etc': arg,
                 })
        # Write a log file.
        for log_arg_p in log_arg:
            with log_arg_p.open("w") as lfh:
                lfh.write("WAIT\n")

    def get_task(self, method, stage=None, completed=None):

        if method not in self.task:
            raise KeyError("No '%s' task was prepared. Can not proceed." %
                           method)

        if stage is not None:
            sel_tasks = [t for t in self.task[method] if t["stage"] == stage]
            if not sel_tasks:
                raise KeyError("No '%s' task (for stage '%s') was prepared."
                               " Can not proceed."% (method, stage))

        #
        task_s = []
        for i, task in enumerate(self.task[method]):
            if stage is not None and task["stage"] != stage:
                continue
            if completed is not None:
                if completed:  # Remove uncompleted jobs.
                    if not self._check_logs_status(task):
                        continue
                else:  # Remove completed jobs.
                    if self._check_logs_status(task):
                        continue
            task_s.append((i, task))
        return task_s

    def _check_logs_status(self, task):
        return all([self._check_log_status(lp) for lp in task["log"]])

    def _check_log_status(self, lp):
        completed = False
        with lp.open("r") as lfh:
            completed = lfh.readline().rstrip() == "DONE"
        return completed


def get_tasks_to_run(job, method, stage=None):

    task_s = job.get_task(method=method, stage=stage, completed=False)

    if len(task_s) == 0:
        if stage is None:
            message = ("All tasks '%s' have been completed. Can not run them"
                       " again." % (method))
        else:
            message = ("All tasks '%s' for (stage '%s') have been completed."
                       " Can not run them again." % (method, stage))
        raise ValueError(message)

    return task_s


def get_outputs(job, method, stage=None, force_complete=False, expand=None):
    task_s = job.get_task(method, stage=stage,
                          completed=True if force_complete else None)
    task_a = job.get_task(method, stage=stage)

    if force_complete:
        if len(task_s) < len(task_a):
            raise ValueError("Not all tasks '%s' (for stage '%s') were"
                             " completed. Can not proceed." % (method, stage))

    out_s = []
    for index, task in task_s:
        if expand is None:
            out_s.append(task['output'])
        else:
            output_expanded = []
            for out in task['output']:
                if out.endswith(expand):
                    with out.open() as fp:
                        _out = []
                        for line in fp:
                            if line.startswith("#"): continue
                            _out.append(path.Path(line.strip()))
                        output_expanded.append(_out)
                else:
                    output_expanded.append(out)
            out_s.append(output_expanded)
    return out_s


#------------------------------------------------------------------------------
# Code to run parallel processes with each process using a different GPU.     -
#------------------------------------------------------------------------------

def parse_gpu_arg(gpu_arg):

    if gpu_arg is None:
        return None

    gpus_ids = []
    sl = gpu_arg.rstrip(":").split(":")
    for si in sl:
        try:
            _id = int(si)
            if si in gpus_ids:
                raise argparse.ArgumentTypeError("Invalid argument. Each"
                    " GPU id can be used only once.")
            gpus_ids.append(si)
        except ValueError:
            raise argparse.ArgumentTypeError("Invalid argument. Must"
                " supply an integer or a list of integers separated by :"
                " characters.")
    return gpus_ids


def run_parallel_tasks(job, task_s, module, polling_time=0.1):
    """
    Run on multiple GPUs.
    """

    # Get the the tasks to execut on each GPU.
    gpu_batches = get_gpu_batches(job, task_s)

    print("- Preparing to run in parallel %s tasks distributed on %s GPUs" % (
          len(task_s), len(job.gpu_ids)))

    # Initializes the processes.
    procs = []  # Will contain a list of Popen objects.
    stdout_logs = []  # Will contain a list of file handles for logs.

    logs_dp = os.path.join(os.path.dirname(job.json_job.path()),
                           "parallel_logs")
    if not os.path.isdir(logs_dp):  # Build the directory for the log files.
        os.mkdir(logs_dp)

    for gpu_id, batch_task_s in gpu_batches:

        log_fp = os.path.join(logs_dp, "%s.module.gpu%s.log" % (module,
                                                                gpu_id))

        print("- Running %s tasks on GPU %s..." % (len(batch_task_s), gpu_id))
        if get_verbosity():
            print("- Writing output to %s" % log_fp)

        stdout_fh = open(log_fp, "w")
        # Actually run each subprocess.
        proc = sp.Popen([os.path.join(EXEC_HOME, "subprocess_launcher.py"),
                         "-j", job.json_job.path(),
                         "--module", module,
                         "--gpu", str(gpu_id)],
                        stdout=stdout_fh)
        procs.append(proc)
        stdout_logs.append(stdout_fh)

    # Start to poll the running subprocesses.
    working = True
    completed_ids = set()
    while working:
        if polling_time is not None:
            time.sleep(polling_time)
        completed = 0
        for idx, proc_i in enumerate(procs):
            status_i = proc_i.poll()
            if status_i is not None:
                if status_i == 0:  # A subprocess has sucessfully completed.
                    gpu_idx = gpu_batches[idx][0]
                    if not gpu_idx in completed_ids:
                        print("- Completed all tasks for GPU %s." % gpu_idx)
                        completed_ids.add(gpu_idx)
                    completed += 1
                else:  # There was an error in some subprocess.
                    for proc_j in procs:  # Stop all subprocesses.
                        if proc_j.poll() is None:
                            proc_j.kill()
                    working = False
        if completed == len(procs):
            working = False

    # Close the log file handles.
    for stdout_fh in stdout_logs:
        if not stdout_fh.closed:
            stdout_fh.close()

    # Check for errors.
    error_procs = [proc for proc in procs if proc.poll() != 0]
    if error_procs:
        for proc in error_procs:
            print("Returnstatus=%s in proc: %s" % (proc.poll(), proc.args))
        raise ChildProcessError("Stopping beacause there was an error in a"
                                " subprocess.")


def get_gpu_batches(job, task_s):

    if len(job.gpu_ids) == 1:
        return task_s
    else:
        gpus_ids_dict = {}
        for i, t in task_s:
            gpu_id = str(t['_resource']['gpu_id'])
            if gpu_id not in job.gpu_ids:
                raise KeyError("GPU %s was not assigned when initializing the"
                               " refinement job." % gpu_id)
            if not gpu_id in gpus_ids_dict:
                gpus_ids_dict[gpu_id] = [(i, t)]
            else:
                gpus_ids_dict[gpu_id].append((i, t))
        job_gpu_ids = job.gpu_ids[:len(task_s)]
        return [(gi, gpus_ids_dict[str(gi)]) for gi in job_gpu_ids]


def filter_tasks_by_gpu(task_s, gpu_id):
    return [(i,t) for (i,t) in task_s \
            if str(t['_resource']['gpu_id']) == gpu_id]


def setup_external_mode():
    arg = argparse.ArgumentParser(prog='subprocess_launcher')
    arg.add_argument('-j', '--json', dest='json_fp', required=True,
                     help='JSON file of an already initialized job')
    arg.add_argument('-m', '--module', required=True,
                     help='Name of the module to use')
    arg.add_argument('--gpu', dest='gpu', default="0", type=str,
                     help=('id of the GPU to use in the job. e.g.: 0 (only GPU'
                           ' 0 will be used).'))
    arg = arg.parse_args()

    input_json = path.Path(arg.json_fp)
    job = Job.from_json(input_json)

    task_s = get_tasks_to_run(job, method=arg.module)
    task_s = filter_tasks_by_gpu(task_s, arg.gpu)
    return job, arg.module, task_s


#------------------------------------------------------------------------------
# Check for dependencies.                                                     -
#------------------------------------------------------------------------------

def check_exe(exe_name):
    for path in os.getenv("PATH").split(os.pathsep):
        exe_file = os.path.join(path, exe_name)
        if os.path.isfile(exe_file):
            return exe_file
    return None


def check_dependencies(check_hybrid=False):

    verbose = get_verbosity()

    if verbose:
        print("\n# Checking dependencies...")

    #----------------------------
    # Check basic dependencies. -
    #----------------------------

    # Check PREFMD2 directory.
    WORK_HOME = os.getenv("PREFMD2_HOME")
    if WORK_HOME is None:
        raise EnvironmentError("$PREFMD2_HOME environmental variable is not"
            " defined. Please make sure to define it and set it to the top of"
            " the PREFMD2 directory structure.")

    # Check CHARMM.
    CHARMMEXEC = os.getenv("CHARMMEXEC")
    if verbose:
        print("* Checking CHARMM")
        print("    - Checking $CHARMMEXEC:", CHARMMEXEC)
    if CHARMMEXEC is None:
        raise EnvironmentError("$CHARMMEXEC environmental variable is not"
            " defined. Please make sure to install CHARMM and define"
            " $CHARMMEXEC as the path of the executable file of CHARMM.")

    # Check MMTSB.
    MMTSBDIR = os.getenv("MMTSBDIR")
    CHARMMDATA = os.getenv("CHARMMDATA")
    if verbose:
        print("* Checking MMTSB")
        print("    - Checking $MMTSBDIR:", MMTSBDIR)
    if MMTSBDIR is None:
        raise EnvironmentError("$MMTSBDIR environmental variable is not"
            " defined. Please make sure to obtain the MMTSB toolset"
            " (https://github.com/mmtsb/toolset/) and define $MMTSBDIR as"
            " path to the top of the MMTSB tree.")
    if verbose:
        print("    - Checking MMTSB scripts...")
    if check_exe("convpdb.pl") is None:
        raise FileNotFoundError("convpdb.pl from MMTSB is not present in your"
            " $PATH. Please make sure to obtain the MMTSB toolset"
            " (https://github.com/mmtsb/toolset/), define $MMTSBDIR and add"
            " $MMTSBDIR/perl and $MMSTBDIR/bin to your $PATH.")
    if verbose:
        print("    - Checking $CHARMMDATA:", CHARMMDATA)
    if CHARMMDATA is None:
        raise EnvironmentError("$CHARMMDATA environmental variable is not"
            " defined. Please make sure to obtain the MMTSB toolset"
            " (https://github.com/mmtsb/toolset/) and define $CHARMMDATA as"
            " $MMTSBDIR/data/charmm.")

    # Check locPREFMD.
    LOCPREFMD = os.getenv("LOCPREFMD")
    MOLPROBITY = os.getenv("MOLPROBITY")
    if verbose:
        print("* Checking locPREFMD")
        print("    - Checking $LOCPREFMD:", LOCPREFMD)
    if LOCPREFMD is None:
        raise EnvironmentError("$LOCPREFMD environmental variable is not"
            " defined. Please make sure to obtain locPREFMD"
            " (https://github.com/mmtsb/toolset/) and define $LOCPREFMD.")
    if verbose:
        print("    - Checking locPREFMD scripts...")
    if os.path.join(LOCPREFMD, "scripts", "locprefmd.sh") is None:
        raise EnvironmentError("locprefmd.sh from locPREFMD was not identified"
            " on your system. Please make sure to obtain locPREFMD"
            " (https://github.com/mmtsb/toolset).")
    if verbose:
        print("    - Checking $MOLPROBITY;", MOLPROBITY)
    if MOLPROBITY is None:
        raise EnvironmentError("$MOLPROBITY environmental variable is not"
            " defined. Please make sure to obtain MolProbity"
            " (https://github.com/rlabduke/MolProbity) and define $MOLPROBITY"
            " as the path to the top of the MolProbity tree.")

    # Check mdconv.
    if verbose:
        print("* Checking mdconv")
        print("    - Checking the mdconv executable...")
    if check_exe("mdconv") is None:
        raise FileNotFoundError("mdconv executable is not present in your"
            " $PATH. Please make sure to obtain mdconv"
            " (https://github.com/feiglab/mdconv) and add the directory with"
            " its executable file to your $PATH.")

    # Check the scoring tools.
    RWPLUS_HOME = os.getenv("RWPLUS_HOME")
    if verbose:
        print("* Checking the scoring tools")
        print("    - Checking $RWPLUS_HOME:", RWPLUS_HOME)
    if RWPLUS_HOME is None:
        raise EnvironmentError("$RWPLUS_HOME environmental variable is not"
            " defined. Please make sure to obtain the RWplus program"
            " (https://zhanglab.ccmb.med.umich.edu/RW/) and define"
            " $RWPLUS_HOME as the path to the RWplus home directory.")
    if verbose:
        print("    - Checking the TMscore executable...")
    if check_exe("TMscore") is None:
        raise EnvironmentError("mdconv executable is not present in your"
            " $PATH. Please make sure to obtain the TMscore program"
            " (https://zhanglab.ccmb.med.umich.edu/TM-score/) and add the"
            " directory with its executable file to your $PATH.")

    # Check scwrl4.
    if verbose:
        print("* Checking Scwrl4")
        print("    - Checking the scwrl4 executable...")
    if check_exe("scwrl4") is None:
        raise EnvironmentError("scwrl4 executable is not present in your"
            " $PATH. Please make sure to obtain the Scwrl4 program"
            " (http://dunbrack.fccc.edu/SCWRL3.php/) and add the directory"
            " with its executable file to your $PATH.")

    if not check_hybrid:
        return

    #-----------------------------
    # Check hybridization tools. -
    #-----------------------------

    # Check for HHsuite executables and databases.
    HHSUITE_SEQ_DB = os.getenv("HHSUITE_SEQ_DB")
    HHSUITE_PDB_DB = os.getenv("HHSUITE_PDB_DB")
    if verbose:
        print("* Checking HHsuite")
        print("    - Checking $HHSUITE_SEQ_DB:", HHSUITE_SEQ_DB)
    if HHSUITE_SEQ_DB is None:
        raise EnvironmentError("$HHSUITE_SEQ_DB environmental variable is not"
            " defined. Please make sure to obtain a HHsuite sequence database"
            " (http://wwwuser.gwdg.de/~compbiol/uniclust/2020_06/) and define"
            " $HHSUITE_SEQ_DB.")
    if verbose:
        print("    - Checking $HHSUITE_PDB_DB:", HHSUITE_PDB_DB)
    if HHSUITE_PDB_DB is None:
        raise EnvironmentError("$HHSUITE_PDB_DB environmental variable is not"
        " defined. Please make sure to obtain a HHsuite PDB database("
        " http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/"
        ") and define $HHSUITE_PDB_DB.")
    if verbose:
        print("    - Checking HHsuite executables...")
    for exec_name in ("hhblits", "hhsearch", "hhalign", "ffindex_get"):
        if check_exe(exec_name) is None:
            raise EnvironmentError("{} executable is not present in your"
                " $PATH. Please make sure to obtain the HHsuite programs"
                " (https://github.com/soedinglab/hh-suite) and add the"
                " directory with its executable files to your $PATH.".format(
                exec_name))

    # MODELLER.
    if verbose:
        print("* Checking MODELLER")
        print("    - Checking the modeller library...")
    try:
        import modeller
    except ImportError:
        raise ImportError("MODELLER is not installed as a Python library."
            " Please make sure to obtain MODELLER"
            " (https://salilab.org/modeller/download_installation.html).")

    ROSETTA_HOME = os.getenv("ROSETTA_HOME")
    if verbose:
        print("* Checking Rosetta")
        print("    - Checking $ROSETTA_HOME:", ROSETTA_HOME)
    if ROSETTA_HOME is None:
        raise EnvironmentError("$ROSETTA_HOME environmental variable is not"
            " defined. Please make sure to install Rosetta and define"
            " $ROSETTA_HOME as the path of the home directory of Rosetta.")
