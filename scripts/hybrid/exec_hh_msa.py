#!/usr/bin/env python

from __future__ import print_function

import os
import sys
import argparse

from path import Path, name
from hybrid.libhhsuite import parse_hhr, read_sequence

class MultipleSequenceAlign:
    def __init__(self, query_name, query_seq):
        self.n_templ = 0
        self.query_name = query_name
        self.query = query_seq
        self.l_query = len(self.query)
        self.templ_info = []
        self.templ = []
    def __repr__(self):
        return self.write_in_FASTA()
    def write_in_FASTA(self):
        self.trim_common_gaps()
        wrt = []
        wrt.append(">%s\n"%self.query_name)
        wrt.append("%s\n"%self.query)
#
        for i in range(self.n_templ):
            wrt.append(">%s\n"%self.templ_info[i].name)
            wrt.append("%s\n"%self.templ[i])
        return ''.join(wrt)
    def write_in_PIR(self):
        self.trim_common_gaps()
        wrt = []
        wrt.append(">P1;%s\n"%self.query_name)
        wrt.append("sequence:%s:1:A:%d:A::::\n"%(self.query_name, self.l_query))
        wrt.append("%s*\n"%self.query)
        #
        for i in range(self.n_templ):
            info = self.templ_info[i]
            wrt.append(">P1;%s\n"%info.name)
            if info.t_res is not None:
                wrt.append("structure:%s:%d:%s:%d:%s::::\n"%(info.pdb_id, \
                    info.t_res[0], info.chain_id, info.t_res[1], info.chain_id))
            else:
                wrt.append("structure:%s::%s::%s::::\n"%(info.pdb_id, \
                    info.chain_id, info.chain_id))
            wrt.append("%s*\n"%self.templ[i])
        return ''.join(wrt)
    def append(self, templ):
        if not templ.valid: return
        #
        self.n_templ += 1
        self.templ_info.append(templ)
        #
        q = []
        ts = [[] for i in range(self.n_templ)]
        k = 0
        for i in range(len(templ.q_align)):
            while(templ.q_align[i] != '-' and self.query[k] == '-'):
                q.append(self.query[k])
                for j in range(self.n_templ-1):
                    ts[j].append(self.templ[j][k])
                ts[-1].append('-')
                k += 1

            if templ.q_align[i] != '-':
                q.append(templ.q_align[i])
                for j in range(self.n_templ-1):
                    ts[j].append(self.templ[j][k])
                ts[-1].append(templ.t_align[i])
                k += 1
            else:
                q.append(templ.q_align[i])
                for j in range(self.n_templ-1):
                    ts[j].append('-')
                ts[-1].append(templ.t_align[i])
        self.query = ''.join(q)
        self.templ = [''.join(t) for t in ts]
    def trim_common_gaps(self):
        q = []
        t = [[] for i in range(self.n_templ)]
        for i in range(len(self.query)):
            status = True
            if self.query[i] == '-':
                status = False
                for j in range(self.n_templ):
                    if self.templ[j][i] != '-':
                        status = True
                        break

            if status:
                q.append(self.query[i])
                for j in range(self.n_templ):
                    t[j].append(self.templ[j][i])
        self.query = ''.join(q)
        self.templ = [''.join(t[i]) for i in range(self.n_templ)]


def main(args=None):
    arg = argparse.ArgumentParser(prog='hh_msa',\
            description='converts hhr file to MSA file in PIR format for selected templates')
    arg.add_argument(dest='hhr_fn', metavar='HHR', \
            help='HHsearch log file')
    arg.add_argument('-t', '--title', dest='title', default=None, metavar='NAME', \
            help='Prefix of the query')
    arg.add_argument('-s', '--sequence', dest='fa_fn', required=True, metavar='SEQ',\
            help='Sequence file in FASTA format')
    arg.add_argument('-i', '--index', dest='index', nargs='*', default=[], metavar='INDEX', type=int, \
            help='Index for templates')
    arg.add_argument('-n', '--top', dest='n_top', default=1, metavar='TOP', type=int, \
            help='Top N templates (default, n=1)')
    arg.add_argument('--id', dest='seq_id', default=None, metavar='SEQID', type=float, \
            help='Sequence identity cutoff')
    arg.add_argument('-e', '--evalue', dest='E_cut', default=None, metavar='EVALUE', type=float, \
            help='Sequence identity cutoff')
    arg.add_argument('-d', '--database', dest='pdb70_db', default=None, metavar='DB', type=str, \
            help='PDB70 database path')
    arg.add_argument('-o', '--output', dest='output', default=None, metavar='OUTPUT', \
            help='Output file name')
    arg.add_argument('--noPDB', dest='noPDB', default=False, action='store_true',\
            help='Use sequence alignment only')
    arg.add_argument('--format', dest='format', default='PIR', metavar='FORMAT',\
            help='Output file format (PIR/FASTA)')
    arg.add_argument('--partial', default=False, action='store_true', \
            help='Do not extend to full sequence.')
    arg = arg.parse_args(args)
    if arg.title is None:
        arg.title = name(arg.fa_fn)
    #
    if arg.pdb70_db is not None:
        if not os.path.exists("%s_pdb.ffindex"%arg.pdb70_db):
            arg.pdb70_db = None
        elif not os.path.exists("%s_pdb.ffdata"%arg.pdb70_db):
            arg.pdb70_db = None
    #
    sequence = read_sequence(arg.fa_fn)
    hhresult = parse_hhr(arg.hhr_fn)
    for hh in hhresult:
        hh.extend_with_full_sequence(sequence)
    #
    selected = []
    if len(arg.index) != 0:
        for i in arg.index:
            selected.append(hhresult[i-1])
    elif arg.seq_id is not None:
        for hh in hhresult:
            if hh.seq_id >= arg.seq_id:
                selected.append(hh)
    elif arg.E_cut is not None:
        E_cut = -log10(arg.E_cut)
        for hh in hhresult:
            if hh.Evalue >= E_cut:
                selected.append(hh)
    else:
        for i in range(arg.n_top):
            selected.append(hhresult[i])
    #
    if len(selected) == 0:
        selected.append(hhresult[0])
    else:
        selected.sort(key=lambda hh:hh.seq_id, reverse=True)
    #
    msa = MultipleSequenceAlign(arg.title, sequence)
    for templ in selected:
        if not arg.noPDB:
            templ.prepare_pdb(pdb70_db=arg.pdb70_db)
        msa.append(templ)
    #
    if arg.output is None:
        if arg.format == 'FASTA':
            sys.stdout.write(msa.write_in_FASTA())
        else:
            sys.stdout.write(msa.write_in_PIR())
    else:
        with open(arg.output, 'wt') as fout:
            if arg.format == 'FASTA':
                fout.write(msa.write_in_FASTA())
            else:
                fout.write(msa.write_in_PIR())

if __name__ == '__main__':
    main()
