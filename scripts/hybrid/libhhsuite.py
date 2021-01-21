#!/usr/bin/env python

import os
import sys
import subprocess as sp
from math import log10

from path import Path
from hybrid.libpdb import PDB, Sequence
from hybrid import libsalign

PDB_HOME = None

class HHresult:
    def __init__(self, name, pdb_id, pdb70=True):
        self.valid = True
        self.name = name
        self.pdb_id = pdb_id
        self.is_pdb70 = pdb70
        if pdb70:
            self.chain_id = pdb_id[-1]
        else:
            self.chain_id = ' '
        self._pdb_fn0 = None
        self._pdb_fn = None
        self.q_res = (9999, -9999)
        self.t_res = (9999, -9999)
        self.q_align = ''
        self.t_align = ''
        self.prob = 0.0
        self.Evalue = 100
        self.score = 0.0
        self.seq_id = 0.0
        self.seq_cov = 0.0
    def __repr__(self):
        return self.write_in_FASTA()
    def __eq__(self, pdb_id):
        return self.pdb_id == pdb_id
    def write_in_FASTA(self, query_name="Query"):
        wrt = []
        wrt.append(">%s\n"%query_name)
        wrt.append("%s\n"%self.q_align)
        wrt.append(">%s\n"%self.name)
        wrt.append("%s\n"%self.t_align)
        return ''.join(wrt)
    def write_in_PIR(self, query_name="Query"):
        wrt = []
        wrt.append(">P1;%s\n"%query_name)
        wrt.append("sequence:%s:%4d:A:%4d:A::::\n"%(query_name, self.q_res[0], self.q_res[1]))
        wrt.append("%s*\n"%self.q_align)
        #
        wrt.append(">P1;%s\n"%self.name)
        if self.t_res is None:
            wrt.append("structure:%s::%s::%s::::\n"%(self.pdb_id, self.chain_id, self.chain_id))
        else:
            wrt.append("structure:%s:%4d:%s:%4d:%s::::\n"%(self.pdb_id,\
                    self.t_res[0], self.chain_id, self.t_res[1], self.chain_id))
        wrt.append("%s*\n"%self.t_align)
        return ''.join(wrt)
    def extend_with_full_sequence(self, seq):
        ext = (seq[:self.q_res[0]-1], seq[self.q_res[1]:])
        self.q_res = (1, len(seq))
        self.q_align = '%s%s%s'%(ext[0], self.q_align, ext[1])
        self.t_align = '%s%s%s'%('-'*len(ext[0]), self.t_align, '-'*len(ext[1]))
        #
        match = 0 ; covered = 0
        for ia in range(len(self.q_align)):
            if self.q_align[ia] != '-' and self.q_align[ia] == self.t_align[ia]:
                match += 1
            if self.q_align[ia] != '-':
                if self.t_align[ia] != '-':
                    covered += 1
        self.seq_id = float(match)/len(seq)*100.0
        self.seq_cov = float(covered)/len(seq)*100.0
    def pdb_fn(self, pdb70_db=None):
        if self._pdb_fn is not None:
            return self._pdb_fn
        if self._pdb_fn0 is not None:
            return self._pdb_fn0
        #
        if os.path.exists("%s.pdb"%self.pdb_id):
            self._pdb_fn0 = Path("%s.pdb"%self.pdb_id)
        elif os.path.exists("../%s.pdb"%self.pdb_id):
            self._pdb_fn0 = Path("../%s.pdb"%self.pdb_id)
        elif os.path.exists("%s.pdb"%self.pdb_id[:4]):
            self._pdb_fn0 = Path("%s.pdb"%self.pdb_id[:4])
        elif os.path.exists("../%s.pdb"%self.pdb_id[:4]):
            self._pdb_fn0 = Path("../%s.pdb"%self.pdb_id[:4])
        if self._pdb_fn0 is not None and self._pdb_fn0.status(size=1024):
            return self._pdb_fn0
        #
        if pdb70_db is not None:
            out = sp.check_output("ffindex_get %s %s %s 2>/dev/null"%(\
                                  '%s_pdb.ffdata'%pdb70_db,\
                                  '%s_pdb.ffindex'%pdb70_db,\
                                  self.pdb_id), shell=True)
            self._pdb_fn0 = Path("%s.pdb"%self.pdb_id)
            with self._pdb_fn0.open('wt') as fout:
                fout.write(out)
            if self._pdb_fn0.status(size=1024):
                return self._pdb_fn0
        #
        if os.path.exists("%s/pdb%s.ent"%(PDB_HOME, self.pdb_id[:4])):
            self._pdb_fn0 = Path("%s/pdb%s.ent"%(PDB_HOME, self.pdb_id[:4]))
        else:
            sp.call("exec_pdb_get.py -f %s &> /dev/null"%self.pdb_id, shell=True)
            self._pdb_fn0 = Path("%s.pdb"%self.pdb_id)
        #
        if self._pdb_fn0.status(size=1024):
            return self._pdb_fn0
        else:
            sys.stderr.write("Error: failed to get PDB file for %s\n"%(self.pdb_id))
            self.valid = False
            return None
    def prepare_pdb(self, pdb70_db=None):
        pdb_fn = self.pdb_fn(pdb70_db=pdb70_db)
        if pdb_fn is None:
            self.valid = False
            return
        protein, chainID_s = libsalign.read_pdb_chunk(pdb_fn)
        #
        if self.chain_id not in chainID_s:
            if self.is_pdb70:
                sys.stderr.write("Error: failed to find the template chain in the PDB file, %s\n"%self.pdb_fn())
                self.valid = False
                return
            else:
                self.chain_id = chainID_s[0]
        #
        pdb = PDB(pdb_fn, read_het=False)
        resNo0 = []
        for res in pdb[0].get_residues():
            if res.chainID() != self.chain_id:
                continue
            if not res.check_bb():
                continue
            resNo0.append(res.resNo())
        seq0 = pdb.get_sequence(seqres=True)[self.chain_id].fasta
        protein = protein[self.chain_id]
        aligned_sequence = self.t_align.replace("-","")
        t_align = self.t_align.replace("X","-")
        #
        status, align0 = libsalign.align_structure_sequence(seq0, protein)
        if not status:
            self.t_res = None
            sys.stderr.write("Warning: failed to identify aligned region (1) for %s\n"%self.name)
            libsalign.align_structure_sequence(seq0, protein, verbose=True)
            pdb_fn = Path("%s.pdb"%self.pdb_id)
            with pdb_fn.open("wt") as fout:
                fout.writelines(Sequence.write_SEQRES(self.chain_id, self.t_align.replace("-","")))
                fout.writelines(pdb[0].write_domain(chain_id=self.chain_id))
            return None
        #
        status, align1 = libsalign.align_sequence_sequence(seq0, t_align)
        if not status:
            self.t_res = None
            sys.stderr.write("Warning: failed to identify aligned region (2) for %s\n"%self.name)
            libsalign.align_sequence_sequence(seq0, t_align, verbose=True)
            pdb_fn = Path("%s.pdb"%self.pdb_id)
            with pdb_fn.open("wt") as fout:
                fout.writelines(Sequence.write_SEQRES(self.chain_id, self.t_align.replace("-","")))
                fout.writelines(pdb[0].write_domain(chain_id=self.chain_id))
            return None
        #
        resNo = []
        k = 0
        for i in range(len(align0.fasta)):
            if align0.align[i] != '-':
                resNo.append(resNo0[k])
                k += 1
            else:
                resNo.append(None)
        #
        aligned_resNo = resNo[align1.align.index(aligned_sequence[0]):\
                              align1.align.rindex(aligned_sequence[-1])+1]
        k = 0
        t_align = []
        self.t_res = (9999, -9999)
        res_range = []
        for i in range(len(self.t_align)):
            if self.t_align[i] == '-' or k == len(aligned_resNo):
                t_align.append("-")
            elif aligned_resNo[k] is None:  # disorder
                t_align.append("-")
                k += 1
            else:
                t_align.append(self.t_align[i])
                self.t_res = (min(aligned_resNo[k], self.t_res[0]),\
                              max(aligned_resNo[k], self.t_res[1]))
                res_range.append(aligned_resNo[k])
                k += 1
        if self.t_res == (9999, -9999):
            return None
        self.t_align = ''.join(t_align)
        #
        self.pdb_id = '%s.%04d_%04d'%(self.pdb_id, self.t_res[0], self.t_res[1])
        pdb_fn = Path("%s.pdb"%self.pdb_id)
        with pdb_fn.open("wt") as fout:
            fout.writelines(Sequence.write_SEQRES(self.chain_id, self.t_align.replace("-","")))
            fout.writelines(pdb[0].write_domain(res_range=res_range, chain_id=self.chain_id))
        return pdb_fn

def parse_hhr(hhr_fn, exclude_s=[]):
    result = []
    with open(hhr_fn) as fp:
        hh_id = ''
        exclude = False
        for line in fp:
            if line.startswith("Query "):
                query_name = line.strip().split()[1][:14]
            elif line.startswith(">"):
                hh_id = line.strip().split()[0][1:15]
                if '|' not in hh_id:
                    pdb_id = '%s_%s'%(hh_id[:4].lower(), hh_id[-1])
                else:
                    pdb_id = hh_id.split("|")[1]
                if (hh_id[:4].lower() in exclude_s) or (hh_id[:4].upper() in exclude_s) or \
                   (hh_id in exclude_s) or (pdb_id in exclude_s):
                    exclude = True
                    continue
                exclude = False
                if os.path.exists("%s.pdb"%hh_id):
                    pdb_id = hh_id
                if pdb_id not in result:
                    hh = HHresult(pdb_id, pdb_id, pdb70=(pdb_id[4] == '_'))
                else:
                    name = '%s.%d'%(pdb_id, result.count(pdb_id))
                    hh = HHresult(name, pdb_id, pdb70=(pdb_id[4] == '_'))
                result.append(hh)
            elif exclude:
                pass
            elif line.startswith("Probab"):
                x = line.strip().split()
                hh.prob = float(x[0].split("=")[1])
                hh.Evalue = -log10(float(x[1].split("=")[1]))
                hh.score = float(x[2].split("=")[1])
                hh.seq_id = float(x[4].split("=")[1][:-1])
            elif line.startswith("Q %s"%query_name):
                x = line.strip().split()
                hh.q_align += x[3]
                hh.q_res = (min(hh.q_res[0], int(x[2])), max(hh.q_res[1], int(x[4])))
            elif line.startswith("T %s"%hh_id):
                x = line.strip().split()
                hh.t_align += x[3]
                hh.t_res = (min(hh.t_res[0], int(x[2])), max(hh.t_res[1], int(x[4])))
    return result

def read_sequence(fa_fn):
    seq = []
    with open(fa_fn) as fp:
        for line in fp:
            if not line.startswith(">"):
                seq.append(line.strip())
    return ''.join(seq)

def read_a3m(a3m_fn):
    a3m = {}
    name_s = []
    with open(a3m_fn) as fp:
        for line in fp:
            if line.startswith(">"):
                name = line.strip().split()[0][1:]
                if '|' in name:
                    name = name.split("|")[1]
                a3m[name] = []
                name_s.append(name)
            else:
                a3m[name].append(line.strip())
    for name in a3m:
        a3m[name] = ''.join(a3m[name])
    return a3m, name_s
