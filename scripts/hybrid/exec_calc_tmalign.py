"""
#!/usr/bin/env python

import os
import sys
import path
import subprocess as sp
import argparse
import tempfile

EXEC=os.path.abspath(__file__)

def distribute_jobs(ref_fn, pdb_fn_s, n_proc=1):
    ref_fn = path.Path(ref_fn)
    job_s = [[] for i in range(n_proc)]
    for i,pdb_fn in enumerate(pdb_fn_s):
        job_s[i%n_proc].append(path.Path(pdb_fn))
    #
    dir_s = []
    for i in range(n_proc):
        dir = tempfile.mkdtemp(prefix='tmalign.')
        #
        with open("%s/cmd"%dir, 'wt') as fout:
            fout.write("cd %s\n"%dir)
            fout.write('%s '%EXEC)
            fout.write("--ref %s "%ref_fn)
            fout.write("--pdb ")
            fout.write(" ".join(['%s'%pdb_fn for pdb_fn in job_s[i]]))
            fout.write(" > out\n")
        dir_s.append(dir)
    #
    return dir_s

def run_single(ref_fn, pdb_fn):
    cmd = []
    cmd.append("TMalign")
    cmd.append(pdb_fn)
    cmd.append(ref_fn)
    #
    tm_s = []
    output = sp.check_output(cmd)
    if sys.version_info.major == 3:
        output = output.decode("utf-8")
    output = output.split("\n")[:-1]
    for line in output:
        if line.startswith("TM-sco"):
            tm_s.append(float(line.strip().split()[1]))
    return tm_s

def run_jobs(dir_s, n_proc=1):
    script_fn = tempfile.mkstemp(prefix='tmalign.', suffix='.sh')[1]
    with open(script_fn, 'wt') as fout:
        for dir in dir_s:
            fout.write("sh %s/cmd\n"%dir)
    #
    pwd = os.getcwd()
    os.chdir(os.path.dirname(script_fn))
    #
    cmd = []
    cmd.append('parallel')
    cmd.append('-j')
    cmd.append('%d'%n_proc)
    cmd.append('--workdir')
    cmd.append('.')
    cmd.append("::::")
    cmd.append(script_fn)
    #
    with open(os.devnull, 'w') as stdout:
        sp.call(cmd, stdout=stdout)
    #
    os.remove(script_fn)
    #
    os.chdir(pwd)

def merge_output(dir_s, n_proc=1):
    output_s = []
    for dir in dir_s:
        output = []
        with open("%s/out"%dir) as fp:
            for line in fp:
                output.append(line.strip().split())
        output_s.append(output)
    n_output = sum([len(output) for output in output_s])
    #
    score = []
    for i in range(n_output):
        out = output_s[i%n_proc][int(i/n_proc)]
        score.append([float(o) for o in out[:-1]])
    return score

def clear_jobs(dir_s):
    for dir in dir_s:
        sp.call(['rm', '-rf', dir])

def main():
    arg = argparse.ArgumentParser(prog='calc_tmalign')
    arg.add_argument('-j', '--cpu', dest='n_proc', help='number of proc.', default=1, type=int)
    arg.add_argument('-r', '--ref', dest='ref_fn', help='reference', required=True)
    arg.add_argument('-p', '--pdb', dest='pdb_fn_s', help='PDB files', default=[], nargs='*')
    arg.add_argument('-l', '--list', dest='pdb_list', help='PDB file lists', default=[], nargs='*')
    if len(sys.argv) == 1:
        return arg.print_help()
    arg = arg.parse_args()
    #
    pdb_fn_s = arg.pdb_fn_s
    for fn in arg.pdb_list:
        with open(fn) as fp:
            for line in fp:
                if not line.startswith("#"):
                    pdb_fn_s.append(line.strip())
    arg.n_proc = min(arg.n_proc, len(pdb_fn_s))
    #
    if arg.n_proc > 1:
        dir_s = distribute_jobs(arg.ref_fn, pdb_fn_s, n_proc=arg.n_proc)
        run_jobs(dir_s, n_proc=arg.n_proc)
        score = merge_output(dir_s, n_proc=arg.n_proc)
        clear_jobs(dir_s)
    else:
        score = []
        for pdb_fn in pdb_fn_s:
            try:
                score.append(run_single(arg.ref_fn, pdb_fn))
            except:
                score.append((0.0, 0.0))
    #
    for i,pdb_fn in enumerate(pdb_fn_s):
        sys.stdout.write("%6.4f %6.4f  "%tuple(score[i]))
        sys.stdout.write("%s\n"%pdb_fn)

if __name__=='__main__':
    main()
"""

#!/usr/bin/env python

import os
import sys
import argparse

import libscore


def main(args=None):
    arg = argparse.ArgumentParser(prog='calc_tmscore')
    arg.add_argument('-r', '--ref', dest='ref_fn', help='reference PDB file',
                     required=True)
    arg.add_argument('-l', '--list', dest='pdb_list', help='PDB file lists',
                     required=True)
    arg.add_argument('-o', '--out', dest='out_fp', help='output filepath',
                     default="qual_init")
    arg.add_argument('-j', '--cpu', dest='n_proc', help='number of proc.',
                     default=1, type=int)
    arg = arg.parse_args(args)

    libscore.run_scoring_tmalign(reference_fp=arg.ref_fn,
                                 input_fpl=arg.pdb_list,
                                 out_fp=arg.out_fp, n_jobs=arg.n_proc)

if __name__=='__main__':
    main()
