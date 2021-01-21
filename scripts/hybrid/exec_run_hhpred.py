#!/usr/bin/env python

"""
Module to search templates with programs from HHsuite (first run hhblits on a
sequence databases and then run hhsearch on the PDB70 database).
"""

import os
import sys
import argparse
import subprocess as sp

import path
import libcommon
from hybrid.libhhsuite import read_a3m, parse_hhr, read_sequence

errfile = None  # "/dev/null"

def run_hhblits(fa_fn, db, prefix, n_proc, verbose):
    a3m_fn = '%s.a3m'%prefix
    hhr_fn = '%s.hhblits.hhr'%prefix
    if os.path.exists(a3m_fn):
        return a3m_fn
    #
    if os.path.exists('%s_a3m.ffdata'%db):
        db_fn = db
    elif db in ['uniclust30', 'uc30']:
        db_fn = HH_sequence_database
    else:
        sys.stderr.write("ERROR: Cannot identify the HHblits database\n")
        sys.exit()
    #
    cmd = ['hhblits']
    cmd.extend(["-i","%s"%fa_fn])
    cmd.extend(["-o", "%s"%hhr_fn])
    cmd.extend(["-oa3m", "%s"%a3m_fn])
    cmd.extend(["-d", "%s"%db_fn])
    cmd.extend(["-cpu", "%d"%n_proc])
    libcommon.system(cmd, stdout=True, errfile=errfile)
    return a3m_fn

def run_psipred(a3m_fn, verbose):
    has_psipred = False
    with open(a3m_fn, 'r') as fp:
        for line in fp:
            if line.startswith(">ss_pred"):
                has_psipred = True
                break
    if has_psipred:
        return a3m_fn
    #
    cmd = ['addss.pl', a3m_fn]
    libcommon.system(cmd, stdout=True, errfile=errfile)
    return a3m_fn

def run_hhsearch(a3m_fn, db, prefix, use_global, use_viterbi, n_proc, verbose):
    out_fn = prefix
    if use_viterbi:
        out_fn += '.vit'
    else:
        out_fn += '.mac'
    if use_global:
        out_fn += '.global'
    else:
        out_fn += '.local'
    if os.path.exists(out_fn):
        return out_fn
    #
    if db is 'pdb70':
        db_fn = os.getenv("HHSUITE_PDB_DB")
    else:
        if os.path.exists('%s_a3m.ffdata'%db):
            db_fn = db
        else:
            sys.stderr.write("ERROR: Cannot identify the HHsearch database\n")
            sys.exit()
    #
    cmd = ['hhsearch']
    cmd.extend(["-i", "%s"%a3m_fn])
    cmd.extend(["-d", "%s"%db_fn])
    cmd.extend(["-o", "%s"%out_fn])
    cmd.extend(["-z", "100", "-b", "100"])
    cmd.extend(["-cpu", "%d"%n_proc])
    if use_viterbi:
        cmd.append("-norealign")
    if use_global:
        cmd.append("-glob")
    libcommon.system(cmd, stdout=True, errfile=errfile)
    return out_fn

def run_hhalign(q_fn, t_fn, prefix, use_global, use_viterbi, use_alt, verbose):
    out_fn = prefix
    if use_viterbi:
        out_fn += '.vit'
    else:
        out_fn += '.mac'
    if use_global:
        out_fn += '.global'
    else:
        out_fn += '.local'
    if os.path.exists(out_fn):
        return out_fn
    #
    cmd = ["hhalign"]
    cmd.extend(["-i", "%s"%q_fn])
    cmd.extend(["-t", "%s"%t_fn])
    cmd.extend(["-o", "%s"%out_fn])
    if use_viterbi:
        cmd.append("-norealign")
    if use_global:
        cmd.append("-glob")
    if use_alt:
        cmd.extend(["-alt", "10"])
    libcommon.system(cmd, stdout=True, errfile=errfile)
    return out_fn

def select_template(hh_s, Evalue_cutoff=1000.0, SeqID_cutoff=0.0):
    templ_id_s = []
    for hh in hh_s:
        if 10.0**(-hh.Evalue) > Evalue_cutoff:
            continue
        if hh.seq_id < SeqID_cutoff:
            continue
        templ_id_s.append(hh.pdb_id)
    return templ_id_s

def build_model(out_prefix, glo_fn, fa_fn, pdb70_db, n_proc, allow_fix_build, verbose, alt=0):
    out_fn = '%s.model.pdb'%out_prefix
    if os.path.exists(out_fn):
        return out_fn
    #
    if pdb70_db is 'pdb70':
        pdb70_db = os.getenv("HHSUITE_PDB_DB")
    #
    cmd = []
    cmd.append(glo_fn)
    cmd.extend(["-s", "%s"%fa_fn])
    if (os.path.exists("%s_pdb.ffdata"%pdb70_db) and \
        os.path.exists("%s_pdb.ffindex"%pdb70_db)):
        cmd.extend(["-d", "%s"%pdb70_db])
    cmd.extend(["-t", "%s"%out_prefix])
    cmd.extend(["-i", "%d"%(alt+1)])
    cmd.extend(["-o", "%s.pir"%out_prefix])
    libcommon.asystem(module="hybrid.exec_hh_msa", args=cmd,
                      stdout=True, errfile=errfile)
    #
    cmd = []
    cmd.append("%s.pir"%out_prefix)
    cmd.extend(["-n", "%d"%n_proc])
    cmd.extend(["-c", "%d"%n_proc])
    if allow_fix_build:
        cmd.append("--fix")
    libcommon.asystem(module="hybrid.exec_run_modeller", args=cmd,
                      stdout=True, errfile=errfile)
    #
    if os.path.exists(out_fn):
        return out_fn
    else:
        return None

def trim_unaligned(model, hh):
    aligned = []
    for i in range(hh.q_res[0]-1):
        aligned.append(False)
    for i in range(len(hh.q_align)):
        if hh.q_align[i] == '-':
            continue
        if hh.t_align[i] == '-':
            aligned.append(False)
        else:
            aligned.append(True)
    #
    pdb = []
    with open(model) as fp:
        for line in fp:
            if not line.startswith("ATOM"):
                pdb.append(line)
                continue
            resNo = int(line[22:26])
            if resNo < hh.q_res[0] or resNo > hh.q_res[1]:
                continue
            if aligned[resNo-1]:
                pdb.append(line)
    #
    with open(model, 'wt') as fout:
        fout.writelines(pdb)

def get_A3M(templ_id, pdb70_db, verbose):
    a3m_fn = '%s.a3m'%templ_id
    if os.path.exists(a3m_fn):
        return a3m_fn
    #
    if pdb70_db is 'pdb70':
        pdb70_db = os.getenv("HHSUITE_PDB_DB")
    #
    cmd = []
    cmd.append("ffindex_get")
    cmd.append("%s_a3m.ffdata"%pdb70_db)
    cmd.append("%s_a3m.ffindex"%pdb70_db)
    cmd.append("%s%s"%(templ_id[:4].upper(), templ_id[4:]))
    try:
        with open("%s.a3m"%templ_id, 'wt') as fout:
            libcommon.system(cmd, outfile=fout, stdout=True, errfile=errfile)
        if os.path.exists(a3m_fn):
            return a3m_fn
    except:
        return None
    return None

def write_FASTA_file(name, a3m_info):
    fa_fn = '%s.fa'%name
    if os.path.exists(fa_fn):
        return fa_fn
    #
    seq = a3m_info[name]
    seq = seq.replace("-","")
    seq = seq.upper()
    #
    with open(fa_fn, 'wt') as fout:
        fout.write(">%s\n"%name)
        fout.write("%s\n"%seq)
    return fa_fn

def main(args=None):
    arg = argparse.ArgumentParser(prog='exec_run_hhpred')
    arg.add_argument(dest='fa_fn', metavar='SEQ', \
            help='Sequence file in FASTA format')
    arg.add_argument('-d', '--database', dest='db', metavar='DB',\
            default='uniclust30',\
            help='Database for the HHblit run (default: uniclust30)')
    arg.add_argument('--hhdb', dest='hhdb', metavar='HHDB',\
            default='pdb70',\
            help='Database for the HHsearch run (default: pdb70)')
    arg.add_argument('-glob', '--global', dest='use_global', \
            default=False, action='store_true', \
            help='Use global alignment method (default: local)')
    arg.add_argument('-vit', '--viterbi', dest='use_viterbi', \
            default=False, action='store_true', \
            help='Use Viterbi alignment method (default: MAC)')
    arg.add_argument('-c', '--cpu', dest='n_proc',\
            default=8, type=int,\
            help='Number of processors to be used (default: 8)')
    arg.add_argument('-o', '--output', dest='prefix', default='', \
            help='Output prefix (default: prefix of the input file)')
    arg.add_argument('-p', '--protocol', dest='protocol',\
            help='Use protocol (default: none), (options: simple/build/build.partial/build.inverse)')
    arg.add_argument('-n', '--n_templ', dest='n_templ', default=1, type=int,\
            help='Use top N template (default: 1)')
    arg.add_argument('--evalue', dest='Evalue_cutoff', \
            help='E-value cutoff for template selection (default=1e-3)', default=0.001)
    arg.add_argument('--seq_id', dest='SeqID_cutoff', \
            help='SeqID cutoff for template selection (default=20%%)', default=20.0)
    arg.add_argument('--pdb', dest='pdb_fn0', default=None,\
            help="The template PDB file for inverse.build mode")
    arg.add_argument('--n_model', dest='n_model', default=1, type=int,\
            help='Number of models to be built (default: n_proc)')
    arg.add_argument('--alt', dest='n_alt', default=1, type=int, \
            help='Number of alt models to be built (default: 1)')
    arg.add_argument('--include', dest='include', default=None, \
            help='Template list file')
    arg.add_argument('--exclude', dest='exclude', default=None, \
            help='Template exclusion list file')
    arg.add_argument('--force', dest='forced_build', default=False, \
            help='Forced building', action='store_true')
    arg.add_argument('--allow_fix_build', dest='allow_fix_build', default=False, \
            help='Use fixed build mode', action='store_true')
    arg.add_argument('-v', '--verbose', dest='verbose',\
            default=False, action='store_true', \
            help='Use verbose mode (default: False)')
    arg = arg.parse_args(args)

    #
    if arg.prefix == '':
        arg.prefix = path.name(arg.fa_fn)
    if arg.protocol == 'inverse.build' and arg.pdb_fn0 is None:
        sys.stderr.write("Error: a template PDB file is required for inverse.build mode.\n")
        return
    #
    sequence = read_sequence(arg.fa_fn)

    # First run HHblits.
    a3m_fn = run_hhblits(arg.fa_fn, arg.db, arg.prefix, arg.n_proc, arg.verbose)
    # a3m_fn = run_psipred(a3m_fn, arg.verbose)
    #

    # Then run HHsearch.
    if arg.protocol == '':
        out_fn = run_hhsearch(a3m_fn, arg.hhdb, arg.prefix, arg.use_global, arg.use_viterbi,\
                              arg.n_proc, arg.verbose)
    elif arg.protocol == 'simple':
        # vit.local
        vit_loc = run_hhsearch(a3m_fn, arg.hhdb, arg.prefix, False, True, arg.n_proc, arg.verbose)
        # mac.global
        mac_glo = run_hhsearch(a3m_fn, arg.hhdb, arg.prefix, True, False, arg.n_proc, arg.verbose)
    elif arg.protocol == 'build' or arg.protocol == 'build.partial':
        # vit.local
        vit_loc = run_hhsearch(a3m_fn, arg.hhdb, arg.prefix, False, True, arg.n_proc, arg.verbose)
        # mac.global
        mac_glo = run_hhsearch(a3m_fn, arg.hhdb, arg.prefix, True, False, arg.n_proc, arg.verbose)
        #
        exclude = []
        if arg.exclude is not None:
            with open(arg.exclude) as fp:
                for line in fp:
                    x = line.strip() ; pdb_id = '%s%s'%(x[:4].lower(), x[4:])
                    exclude.append(pdb_id)
        #
        vit_loc_s = parse_hhr(vit_loc, exclude_s=exclude)
        for hh in vit_loc_s:
            hh.extend_with_full_sequence(sequence)
        #
        if arg.include is None:
            templ_id_s = select_template(vit_loc_s, Evalue_cutoff=arg.Evalue_cutoff, SeqID_cutoff=arg.SeqID_cutoff)[:arg.n_templ]
        else:
            templ_id_s = select_template(vit_loc_s)
            include = [line.strip() for line in open(arg.include)]
            templ_id_from_include = []
            for templ_id in include:
                alt_id = '%s%s'%(templ_id[:4].lower(), templ_id[4:])
                if templ_id in templ_id_s:
                    templ_id_from_include.append(templ_id)
                elif alt_id in templ_id_s:
                    templ_id_from_include.append(alt_id)
                elif arg.forced_build:
                    templ_id_from_include.append(templ_id)
            templ_id_s = templ_id_from_include
        with open("%s.templ_s"%arg.prefix, 'wt') as fout:
            for templ_id in templ_id_s:
                fout.write("%s\n"%templ_id)
        #
        summary = []
        for templ_id in templ_id_s:
            out_prefix = '%s-%s'%(arg.prefix, templ_id)
            #
            templ_a3m_fn = get_A3M(templ_id, arg.hhdb, arg.verbose)
            if templ_a3m_fn is None: continue
            #
            align_fn = run_hhalign(a3m_fn, templ_a3m_fn, out_prefix, \
                                   (arg.protocol == 'build'), False, (arg.n_alt > 1), arg.verbose)
            if not os.path.exists(align_fn): continue
            hh_s = parse_hhr(align_fn)[:arg.n_alt]
            out = vit_loc_s[vit_loc_s.index(templ_id)]
            for ih,hh in enumerate(hh_s):
                if ih > 0:
                    if hh.score < hh_s[0].score * 0.5: break
                    out_prefix_ih = '%s.%d'%(out_prefix, ih)
                else:
                    out_prefix_ih = out_prefix
                model = build_model(out_prefix_ih, align_fn, arg.fa_fn, arg.hhdb, arg.n_proc, arg.allow_fix_build, arg.verbose, alt=ih)
                if (arg.protocol == 'build.partial') and (model is not None):
                    trim_unaligned(model, hh)
                hh.extend_with_full_sequence(sequence)
                summary.append(' '.join(["%-8.1e"%(10.0**(-out.Evalue)), '%5.1f'%hh.seq_id, '%5.1f'%hh.seq_cov, hh.pdb_id, out_prefix_ih])+'\n')
        with open("%s.templ_s.summary"%arg.prefix, 'wt') as fout:
            fout.writelines(summary)

    elif arg.protocol == 'build.inverse':
        if arg.exclude is None:
            exclude = []
        else:
            exclude = [line.strip() for line in open(arg.exclude)]
        #
        template_indexs, template_ids = select_template(hhr_fn, hhr_fn, exclude)
        a3m_info, a3m_name = read_a3m(a3m_fn)
        #
        selected = []
        for templ in template_ids:
            if templ in a3m_name:
                selected.append(templ)
        for templ in a3m_name:
            if not templ.startswith("ss_") and templ not in selected:
                selected.append(templ)
        selected = selected[:arg.n_templ]
        #
        fout = open("%s.templ_s"%arg.prefix, 'wt')
        for templ in selected:
            out_prefix = '%s-%s'%(arg.prefix, template_id)
            fa_fn = write_FASTA_file(templ, a3m_info)
            #
            a3m, hhr = run_hhblits(fa_fn, arg.db, templ, arg.n_proc, arg.verbose)
            a3m = run_psipred(a3m, arg.verbose)
            #
            # mac.global
            ali_fn = run_hhalign(a3m, a3m_fn, templ, True, False, False, arg.verbose)
            build_model(ali_fn, fa_fn, out_prefix, arg.hhdb, templ, arg.n_proc, arg.verbose)
            #
            fout.write("%s\n"%templ)
        fout.close()

if __name__ == '__main__':
    main()
