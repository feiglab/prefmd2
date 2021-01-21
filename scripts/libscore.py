"""
Module for functions to score the structures extracted from the production MD
trajectories.
"""

import os
import subprocess
import multiprocessing


#------------------------------------------------------------------------------
# Executables paths.                                                          -
#------------------------------------------------------------------------------

RWPLUS_HOME = os.getenv("RWPLUS_HOME")
DDFIRE_EXEC = os.getenv("DDFIRE_EXEC")  # Actually not needed for the default
                                        # pipeline.


#------------------------------------------------------------------------------
# Common.                                                                     -
#------------------------------------------------------------------------------

def _get_fpl_from_file(input_fpl):
    fp_l = []  # List of the paths for input PDB files.
    with open(input_fpl, "r") as i_fh:
        fp_l = [l.rstrip() for l in i_fh]
    return fp_l

def _write_output_file(out_fp, fp_l, scores):
    with open(out_fp, "w") as o_fh:
        for f, s in zip(fp_l, scores):
            o_fh.write("{}\n".format(s))

def _check_filepaths(input_fpl):
    for pdb_filepath in input_fpl:
        if not os.path.isfile(pdb_filepath):
            raise FileNotFoundError("File {} not found.".format(pdb_filepath))


#------------------------------------------------------------------------------
# RW plus potential.                                                          -
#------------------------------------------------------------------------------

def _convert_rw_line(line):
    return float(line.split("=")[1].replace("kcal/mol", ""))

def _rw_single_file(pdb_filepath):
    args = ["./calRWplus", pdb_filepath]
    try:
        out = subprocess.check_output(args)
        out = out.decode("utf-8")
        rw_score = _convert_rw_line(out)
    except subprocess.CalledProcessError as e:
        try:
            rw_score = _convert_rw_line(e.output)
        except IndexError:
            rw_score = 0.0
    return rw_score

def rwplus(pdb_filepaths, rwplus_home=RWPLUS_HOME, n_jobs=1):
    """
    Return the RWplus scores for a list of 'pdb_filepaths'.
    """
    pdb_filepaths_abs = [os.path.abspath(fp) for fp in pdb_filepaths]
    original_dirpath = os.getcwd()
    os.chdir(rwplus_home)

    rw_scores = []
    if n_jobs == 1:
        for fp in pdb_filepaths_abs:
            rw_scores.append(_rw_single_file(fp))
    else:
        pool = multiprocessing.Pool(n_jobs)
        rw_scores = pool.map(_rw_single_file, pdb_filepaths_abs)
        pool.close()

    os.chdir(original_dirpath)
    return rw_scores

def score_rwplus(input_fpl, out_fp=None, n_jobs=1):
    """
    Arguments:
      input_fpl: file with the list of input PDBs.
    """
    fp_l = _get_fpl_from_file(input_fpl)
    _check_filepaths(fp_l)

    scores = rwplus(fp_l, n_jobs=n_jobs)
    if out_fp is not None:
        _write_output_file(out_fp, fp_l, scores)
    return scores


#------------------------------------------------------------------------------
# dDFIRE potential.                                                           -
#------------------------------------------------------------------------------

def ddfire(pdb_filepaths, ddfire_filepath=DDFIRE_EXEC, get_components=True):
    """
    Return the dDFIRE scores for a list of 'pdb_filepaths'.
    """
    args = [ddfire_filepath] + pdb_filepaths
    out = subprocess.check_output('export DATADIR="%s/"; %s' % (
                                  os.path.dirname(ddfire_filepath),
                                  " ".join(args)), shell=True)

    ddfire_scores = []
    dfire_scores = []
    ddfire_ag1 = []
    ddfire_ag2 = []
    ddfire_ag3 = []

    out = out.decode("utf-8")

    for l in out.split("\n")[0:-1]:
        fields = l.split()
        ddfire_scores.append(float(fields[-5]))
        dfire_scores.append(float(fields[-4]))
        ddfire_ag1.append(float(fields[-3]))
        ddfire_ag2.append(float(fields[-2]))
        ddfire_ag3.append(float(fields[-1]))

    if not get_components:
        return ddfire_scores
    else:
        return ddfire_scores, dfire_scores, ddfire_ag1, ddfire_ag2, ddfire_ag3

def score_ddfire(input_fpl, out_fp=None, n_jobs=1, batch_size=50):
    """
    Arguments:
      input_fpl: file with the list of input PDBs.
    """

    fp_l = _get_fpl_from_file(input_fpl)
    _check_filepaths(fp_l)

    batches = []
    for i in range(0, len(fp_l), batch_size):
        batches.append(fp_l[i:i+batch_size])

    ddfire_scores = []
    dfire_scores = []
    if n_jobs == 1:
        for batch in batches:
            _ddfire_s, _dfire_s, _, _, _ = ddfire(batch, get_components=True)
            ddfire_scores.extend(_ddfire_s)
            dfire_scores.extend(_dfire_s)
    else:
        pool = multiprocessing.Pool(n_jobs)
        mp_results = pool.map(ddfire, batches)
        pool.close()
        for batch_results in mp_results:
            ddfire_scores.extend(batch_results[0])
            dfire_scores.extend(batch_results[1])

    if out_fp is not None:
        _write_output_file(out_fp + ".ddfire", fp_l, ddfire_scores)
        _write_output_file(out_fp + ".dfire", fp_l, dfire_scores)
    return ddfire_scores, dfire_scores


#------------------------------------------------------------------------------
# Compute structural similarity with TMscore.                                 -
#------------------------------------------------------------------------------

def compare_with_tmscore(reference_filepath, target_filepath,
                         tmscore_path="TMscore"):
    """
    Use the TMscore program to compare a target PDB file to a reference PDB
    file.
    """
    args = [tmscore_path, reference_filepath, target_filepath]
    out_str = subprocess.check_output(args)
    out_str = out_str.decode("utf-8")

    results_dict = {}
    for l in out_str.split("\n"):
        if l.startswith("TM-score"):
            results_dict["TM-score"] = float(l[14:20])
        elif l.startswith("GDT-TS-score"):
            results_dict["GDT-TS"] = float(l[14:20])
        elif l.startswith("GDT-HA-score"):
            results_dict["GDT-HA"] = float(l[14:20])
        elif l.startswith("RMSD"):
            results_dict["RMSD"] = float(l.split("=")[1])

    return results_dict

def _compare_with_tmscore(args):
    return compare_with_tmscore(args[0], args[1])

def score_tmscore(reference_fp, input_fpl, out_fp=None, n_jobs=1):
    """
    Arguments:
      reference_fp: reference PDB filepath.
      input_fpl: file with the list of input PDBs to compare with
          reference_fp.
    """
    fp_l = _get_fpl_from_file(input_fpl)
    _check_filepaths(fp_l)

    tmscore_l = []
    rmsd_l = []
    gdtts_l = []
    gdtha_l = []

    if n_jobs == 1:
        for fp in fp_l:
            r = compare_with_tmscore(reference_fp, fp)
            tmscore_l.append(r["TM-score"])
            rmsd_l.append(r["RMSD"])
            gdtts_l.append(r["GDT-TS"])
            gdtha_l.append(r["GDT-HA"])

    else:
        pool = multiprocessing.Pool(n_jobs)
        rl = pool.map(_compare_with_tmscore,
                      [(fp, reference_fp) for fp in fp_l])
        pool.close()
        tmscore_l = [r["TM-score"] for r in rl]
        rmsd_l = [r["RMSD"] for r in rl]
        gdtts_l = [r["GDT-TS"] for r in rl]
        gdtha_l = [r["GDT-HA"] for r in rl]

    if out_fp is not None:
        _write_output_file(out_fp + ".tmscore", fp_l, tmscore_l)
        _write_output_file(out_fp + ".rmsd", fp_l, rmsd_l)
        _write_output_file(out_fp + ".gdtts", fp_l, gdtts_l)
        _write_output_file(out_fp + ".gdtha", fp_l, gdtha_l)

    return tmscore_l, rmsd_l, gdtts_l, gdtha_l


#------------------------------------------------------------------------------
# Compute structural similarity with TMalign.                                 -
#------------------------------------------------------------------------------

def compare_with_tmalign(reference_filepath, target_filepath,
                         tmalign_path="TMalign"):
    """
    Use the TMalign program to compare a target PDB file to a reference PDB
    file.
    """
    args = [tmalign_path, reference_filepath, target_filepath]
    out_str = subprocess.check_output(args)
    out_str = out_str.decode("utf-8")
    results_dict = {}
    for i, l in enumerate(out_str.split("\n")):
        if l.startswith("TM-score"):
            if "TM_1" in results_dict:  # Other chain TM-score.
                results_dict["TM_2"] = float(l[10:17])
            else:  # Reference TM-score.
                results_dict["TM_1"] = float(l[10:17])
    return results_dict

def _compare_with_tmalign(args):
    return compare_with_tmalign(args[0], args[1])

def score_tmalign(reference_fp, input_fpl, out_fp=None, n_jobs=1):
    """
    Arguments:
      reference_fp: reference PDB filepath.
      input_fpl: file with the list of input PDBs to compare with
          reference_fp.
    """
    fp_l = _get_fpl_from_file(input_fpl)
    _check_filepaths(fp_l)

    tmscore_l = []
    rmsd_l = []
    gdtts_l = []
    gdtha_l = []

    if n_jobs == 1:
        for fp in fp_l:
            r = compare_with_tmalign(reference_fp, fp)
            tmscore_l.append(r["TM_1"])
    else:
        pool = multiprocessing.Pool(n_jobs)
        rl = pool.map(_compare_with_tmalign,
                      [(fp, reference_fp) for fp in fp_l])
        pool.close()
        tmscore_l = [r["TM_1"] for r in rl]

    if out_fp is not None:
        _write_output_file(out_fp + ".tmalign", fp_l, tmscore_l)

    return tmscore_l


#------------------------------------------------------------------------------
# Functions actually called in the refinement pipeline.                       -
#------------------------------------------------------------------------------

def run_scoring_statpot(input_fpl, out_fp, n_jobs=1,
                        use_rwplus=True, use_ddfire=True, verbose=False):
    results = []
    fp_l = _get_fpl_from_file(input_fpl)
    if use_rwplus:
        if verbose:
            print("- Computing rwplus.")
        rw_l = score_rwplus(input_fpl, n_jobs=n_jobs)
        results.append(rw_l)
    if use_ddfire:
        if verbose:
            print("- Computing ddfire.")
        dd_l, d_l = score_ddfire(input_fpl, n_jobs=n_jobs)
        results.append(dd_l)
        results.append(d_l)
    results.append(fp_l)

    with open(out_fp, "w") as o_fh:
        # RWplus, dDFIRE, DFIRE, file
        for t in zip(*results):
            o_fh.write(" ".join([str(ti) for ti in t]) + "\n")


def run_scoring_tmscore(reference_fp, input_fpl, out_fp, n_jobs=1,
                        verbose=False):

    fp_l = _get_fpl_from_file(input_fpl)
    if verbose:
        print("- Comparing to initial structure with TMscore.")
    tm_l, rmsd_l, ts_l, ha_l = score_tmscore(reference_fp=reference_fp,
                                             input_fpl=input_fpl,
                                             n_jobs=n_jobs)
    with open(out_fp, "w") as o_fh:
        for rmsd, tm, ts, ha, p in zip(rmsd_l, tm_l, ts_l, ha_l, fp_l):
            o_fh.write("{} {} {} {} {}\n".format(rmsd, tm, ts, ha, p))

def run_scoring_tmalign(reference_fp, input_fpl, out_fp, n_jobs=1,
                        verbose=False):

    fp_l = _get_fpl_from_file(input_fpl)
    if verbose:
        print("- Comparing to initial structure with TMalign.")
    tm_l = score_tmalign(reference_fp=reference_fp,
                         input_fpl=input_fpl,
                         n_jobs=n_jobs)
    with open(out_fp, "w") as o_fh:
        for tm, p in zip(tm_l, fp_l):
            o_fh.write("{} {}\n".format(tm, p))
