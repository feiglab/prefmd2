#!/usr/bin/env python

import os
import sys
import copy
from string  import ascii_uppercase as ins_code
#
from hybrid.seqName import to_one_letter, convert_to_stdres_line, stdres

BACKBONE_ATOMS = ['N','CA','C','O']

class ProteinResidue:
    def __init__(self, resNo, resName):
        if resNo.strip().replace("-",'').isdigit():
            self.resNo = int(resNo.strip())
            self.altNo = 0
            self.is_inserted = False
        else:
            self.resNo = int(resNo.strip()[:-1])
            self.altNo = ins_code.index(resNo.strip()[-1])+1
            self.is_inserted = True
        self.resNo0 = resNo
        if len(resName) == 3:
            self.resName = to_one_letter(resName)
        else:
            self.resName = resName
    def __eq__(self,othr):
        return self.resNo == othr.resNo and self.altNo == othr.altNo and \
                self.resName == othr.resName
    def __repr__(self):
        return '%s %s'%(self.resNo0, self.resName)
    def is_prev_residue(self, prev):
        if self.altNo == 0 and prev.altNo == 0:
            # in case of without residue insertion code
            if self.resNo == prev.resNo+1:
                return True
            elif (prev.resNo < 0) and (self.resNo > 0) and \
                    (self.resNo == prev.resNo+2):
                return True
            else:
                return False
        else:
            # in case of with residue insertion code
            #  ie, 10-10A-10B-11 -> all continuous
            if self.resNo == prev.resNo:
                if self.altNo == prev.altNo+1:
                    return True
                else:
                    return False
            else:
                if (self.altNo == 0) and (prev.altNo != 0):
                    if self.resNo == prev.resNo+1:
                        return True
                    elif (prev.resNo < 0) and (self.resNo > 0) and \
                            (self.resNo == prev.resNo+2):
                        return True
                    else:
                        return False
                else:
                    return False
    def compare(self, ref):
        # if   self is previous than ref -> -1
        # elif self is the same with ref ->  0
        # else                           ->  1
        if self.resNo < ref.resNo:
            return -1
        elif self.resNo > ref.resNo:
            return 1
        else:
            if self.altNo < ref.altNo:
                return -1
            elif self.altNo > ref.altNo:
                return 1
            else:
                return 0
    def __lt__(self, ref):
        return self.compare(ref)
    def __cmp__(self, ref):
        return self.compare(ref)

class ProteinSegment:
    def __init__(self, segNo):
        self.segNo = segNo
        self.residue_s = []
        self.aligned = []
    def n_residue(self):
        return len(self.residue_s)
    def get_sequence(self):
        fa = []
        for residue in self.residue_s:
            fa.append(residue.resName)
        return ''.join(fa)
    def trim_sequence(self, trim):
        if trim[0] != 0:
            self.residue_s = self.residue_s[trim[0]:]
        if trim[1] != 0:
            self.residue_s = self.residue_s[:-trim[1]]
    def __cmp__(self, other):
        if self.n_residue() != other.n_residue():
            return -cmp(self.n_residue(), other.n_residue())
        else:
            return cmp(self.segNo, other.segNo)
    def __lt__(self, other):
        if sys.version_info.major == 2:
            return self.__cmp__(other)
        if self.n_residue() != other.n_residue():
            return self.n_residue().__gt__(other.n_residue())
        else:
            return self.segNo.__lt__(other.segNo)

class ProteinStructure:
    def __init__(self):
        self.segment_s = []
    def n_segment(self):
        return len(self.segment_s)
    def add_residue(self, resNo, resName, force_to_add_segment=False):
        residue = ProteinResidue(resNo, resName)
        #
        # initialize
        if len(self.segment_s) == 0:
            self.segment_s.append(ProteinSegment(0))
        elif residue in self.segment_s[-1].residue_s:
            return
        # add another segment if residues are discontinous
        elif force_to_add_segment:
            n_seg = len(self.segment_s)
            self.segment_s.append(ProteinSegment(n_seg))
        elif not residue.is_prev_residue(self.segment_s[-1].residue_s[-1]) and self.segment_s[-1].n_residue() > 0:
            n_seg = len(self.segment_s)
            self.segment_s.append(ProteinSegment(n_seg))
        # add the current residue to the recent segment
        self.segment_s[-1].residue_s.append(residue)
    def get_sequence(self):
        fa = []
        for segment in self.segment_s:
            fa.append(segment.get_sequence())
        return ''.join(fa)
    def get_sorted_segment(self):
        segment_s = copy.deepcopy(self.segment_s)
        segment_s.sort()
        return segment_s

class SegmentAlign:
    def __init__(self, fasta):
        self.status = True
        self.fasta = fasta
        self.l_ali = len(fasta)
        self.align = '-'*(self.l_ali)
        self.aligned_segment = {}
        self.aligned_residue = [None for i in range(self.l_ali)]
    def __repr__(self):
        wrt = []
        wrt.append(">Sequence")
        wrt.append("%s"%self.fasta)
        wrt.append(">Structure")
        wrt.append("%s"%self.align)
        return '\n'.join(wrt)
    def place_segment(self, i, seq, segment):
        l_seq = segment.n_residue()
        self.align = '%s%s%s'%(self.align[:i], seq, self.align[i+l_seq:])
        self.aligned_segment[segment.segNo] = (i,i+l_seq)
        for k in range(l_seq):
            self.aligned_residue[i+k] = segment.residue_s[k]
    def similar_segment(self, seq, resStart, resEnd):
        l_fasta = resEnd-resStart
        l_seq = len(seq)
        l_cand = l_fasta-l_seq+1
        #
        n_diff = []
        for i in range(l_cand):
            fasta_i = self.fasta[resStart+i:resStart+i+l_seq]
            i_diff = []
            for j in range(l_seq):
                if fasta_i[j] != seq[j]:
                    i_diff.append(j)
            n_diff.append((resStart+i, len(i_diff), i_diff))
        n_diff.sort(key=lambda x:x[1])
        if len(n_diff) == 0:
            return 0, range(l_seq)
        return n_diff[0][0], n_diff[0][2]
    def find_segment(self, seq, resStart, resEnd):
        l_f = self.fasta.find(seq, resStart, resEnd)
        r_f = self.fasta.rfind(seq, resStart, resEnd)
        if l_f == r_f:
            if l_f != -1:            # found match
                return 1, [(l_f,0,0)]
            else:                   # failed to find
                n_tol = int(len(seq)/5)
                for i in range(1,n_tol):
                    # seq[:-i]
                    l_fr = self.fasta.find(seq[:-i], resStart, resEnd)
                    r_fr = self.fasta.rfind(seq[:-i], resStart, resEnd)
                    if l_fr != -1:
                        return 1, [(l_fr, 0, i)]
                    #
                    # seq[i:]
                    l_fl = self.fasta.find(seq[i:], resStart, resEnd)
                    r_fl = self.fasta.rfind(seq[i:], resStart, resEnd)
                    if l_fl != -1:
                        return 1, [(l_fl, i, 0)]
                #
                return 1, [(-1,0,0)]

        n_fold = 0
        l_seq = len(seq)
        l_res = resEnd-resStart+1
        l_pos = l_res-l_seq
        found = []
        for k in range(l_pos):
            i = self.fasta.find(seq, resStart+k, resEnd)
            if i != -1 and (i,0,0) not in found:
                n_fold += 1
                found.append((i,0,0))
        return n_fold, found
    def get_segment_boundary(self, segment):
        prev_segment = list(self.aligned_segment.keys())
        prev_segment.append(segment.segNo)
        prev_segment.sort()
        i_seg = prev_segment.index(segment.segNo)
        if i_seg == 0:
            seg_prev = None
        else:
            seg_prev = prev_segment[i_seg-1]
        if i_seg == len(prev_segment)-1:
            seg_next = None
        else:
            seg_next = prev_segment[i_seg+1]
        #
        if seg_prev != None:
            resStart = self.aligned_segment[seg_prev][1]
        else:
            resStart = 0
        if seg_next != None:
            resEnd = self.aligned_segment[seg_next][0]
        else:
            resEnd = self.l_ali
        return resStart, resEnd
    def align_segment(self, segment, verbose=False):
        resStart, resEnd = self.get_segment_boundary(segment)
        #
        seq = segment.get_sequence()
        n_fold,found = self.find_segment(seq, resStart, resEnd)
        #
        if n_fold == 1:    # found the unique match
            i = found[0][0]
        elif n_fold == 0:  # failed to find
            i = -1
        else:              # found more than two matches -> repeat
            return -1, found

        if i != -1:
            if found[0][1] != 0 or found[0][2] != 0:
                segment.trim_sequence(found[0][1:])
            self.place_segment(i, segment.get_sequence(), segment)
            if verbose:
                sys.stderr.write("INFO: %s\n"%self.fasta)
                sys.stderr.write("INFO: %s\n"%self.align)
                sys.stderr.write("INFO: %s%s\n"%(' '*i,'*'*(segment.n_residue())))
                sys.stderr.write("INFO: %s%s\n\n"%(' '*resStart,'^'*(resEnd-resStart)))
            return 1, found
        elif verbose:
            k,diff = self.similar_segment(seq, resStart, resEnd)
            siml = [' '*k]
            for j in range(max(diff)+1):
                if j in diff:
                    siml.append("x")
                else:
                    siml.append(" ")
            siml = ''.join(siml)
            sys.stderr.write("FAIL: %s\n"%seq)
            sys.stderr.write("INFO: %s\n"%self.fasta)
            sys.stderr.write("INFO: %s\n"%self.align)
            sys.stderr.write("SIML: %s%s\n"%(' '*k,seq))
            sys.stderr.write("SIML: %s\n"%siml)
            sys.stderr.write("INFO: %s%s\n\n"%(' '*resStart,'^'*(resEnd-resStart)))
        return 0, found

def read_pdb_chunk(pdb_fn, res_start=None, res_end=None, required_atoms=BACKBONE_ATOMS):
    protein = {}
    chainID_s = []
    #
    with open('%s'%pdb_fn) as fp:
        if res_start == None or res_start[0].strip() == '':
            read = True
            res_start_no = None
        else:
            read = False
            res_start_no = ProteinResidue(res_start[0], 'X')
        if res_end == None or res_end[0].strip() == '':
            res_end_no = None
        else:
            res_end_no = ProteinResidue(res_end[0], 'X')

        status_backbone = {}

        for line in fp:
            if line.startswith("END"):
                break
            if line[:4] not in ['ATOM','HETA']:
                continue
            if line[16] not in [' ','A']:
                continue
            atmName = line[12:16].strip()
            if atmName not in required_atoms:
                continue
            #
            line = convert_to_stdres_line(line)
            if line.startswith("HETA"): continue
            #
            resName = line[17:20]
            if resName not in stdres: continue
            #
            chainID = line[21]
            resNo = line[22:27]
            #
            residue = ProteinResidue(resNo, resName)
            #
            if (not read):
                if (chainID.strip() == res_start[1] or res_start[1].strip() == ''):
                    if residue.compare(res_start_no) != -1:
                        read = True
            if (not read):
                continue
            #
            if res_end_no != None:
                if (chainID.strip() == res_end[1] or res_end[1].strip() == ''):
                    if residue.compare(res_end_no) > 0:
                        break
            #
            # Filtering residues without complete backbone atoms
            key = '%s.%s'%(resNo.strip(),chainID)
            if key not in status_backbone.keys():
                status_backbone[key] = [False for _ in range(len(required_atoms))]
            status_backbone[key][required_atoms.index(atmName)] = True
            if resName == 'GLY' and atmName == 'CA' and 'CB' in required_atoms:
                status_backbone[key][required_atoms.index("CB")] = True
            if False in status_backbone[key]:
                continue
            #
            if chainID not in chainID_s:
                chainID_s.append(chainID)
                protein[chainID] = ProteinStructure()
            prot = protein[chainID]
            #
            prot.add_residue(resNo, resName)
            #
    return protein, chainID_s

def read_sequence_chunk(sequence, structure=''):
    if structure == '':
        prot = ProteinStructure()
        l_seq = len(sequence)
        for i in range(l_seq):
            resName = sequence[i]
            if resName in ['-','!']: continue
            resNo = '%d'%(i+1)
            prot.add_residue(resNo, resName)
        return prot
    else:
        return read_sequence_chunk_with_structure(sequence, structure)

def read_sequence_chunk_with_structure(seq, str):
    prot = ProteinStructure()
    #
    seq0 = seq.replace("-","")
    st, ali = align_structure_sequence(seq0, str)
    #
    chain_break_s = []
    segNo_s = ali.aligned_segment.keys()
    segNo_s.sort()
    for segNo in segNo_s[:-1]:
        chain_break_s.append(ali.aligned_segment[segNo][1]-1)
    #
    updated_seq = []
    k = -1
    for i in range(len(seq)):
        if seq[i] == '-':
            updated_seq.append('-')
        else:
            k += 1
            updated_seq.append(seq[i])
            if k in chain_break_s:
                updated_seq.append('!')
    updated_seq = ''.join(updated_seq)
    return read_sequence_chunk(updated_seq), updated_seq

def report_segment_align(ali, fp):
    import math
    n_seq_in_row = 100
    l_ali = ali.l_ali
    l_row = int(math.ceil(l_ali/n_seq_in_row))+1
    fp.write("REMARK  REFERENCE SEQUENCE\n")
    for i in range(l_row):
        fp.write("REMARK  %s\n"%ali.fasta[i*n_seq_in_row:(i+1)*n_seq_in_row])
    fp.write("REMARK  TARGET SEQUENCE\n")
    for i in range(l_row):
        fp.write("REMARK  %s\n"%ali.align[i*n_seq_in_row:(i+1)*n_seq_in_row])
    #
    segNo_s = ali.aligned_segment.keys()
    segNo_s.sort()
    fp.write("REMARK  ALIGNED SEGMENTS\n")
    for segNo in segNo_s:
        fp.write("REMARK  SEGMENT %3d %4d %4d\n"%(segNo+1,\
                ali.aligned_segment[segNo][0],\
                ali.aligned_segment[segNo][1]))
    #
    for i in range(l_ali):
        res = ali.aligned_residue[i]
        if res != None:
            if res.is_inserted:
                fp.write("ALIGN   %s %s %4d %5s\n"%\
                        (ali.fasta[i], ali.align[i], i+1, res.resNo0))
            else:
                fp.write("ALIGN   %s %s %4d %4s\n"%\
                        (ali.fasta[i], ali.align[i], i+1, res.resNo0))
        else:
            fp.write("ALIGN   %s %s %4d %4s\n"%\
                    (ali.fasta[i], ali.align[i], i+1, '-'))

def update_pdb_resNo(pdb_fn, align_s):
    status = True
    #
    resNo_s = {}
    for chainID in align_s.keys():
        ali = align_s[chainID]
        for i in range(ali.l_ali):
            res = ali.aligned_residue[i]
            if res != None:
                resNo_s[(chainID, res.resNo0)] = i+1
    #
    pdb = []
    with file('%s'%pdb_fn) as fp:
        for line in fp:
            if not line.startswith("ATOM"):
                pdb.append(line)
                continue
            chain = line[21]
            resNo = line[22:27]
            key = (chain,resNo)
            if key in resNo_s:
                resNo_new = resNo_s[key]
                line = '%s%4d %s'%(line[:22], resNo_new, line[27:])
                pdb.append(line)
            else:
                sys.stdout.write('\nERROR  %s'%line)
                status = False
                return status, pdb
    return status, pdb

def build_align(align, pdb, verbose=False):
    align0 = copy.deepcopy(align)
    #
    is_repeat = False
    sorted_segment = pdb.get_sorted_segment()
    for segment in sorted_segment:
        st, found = align.align_segment(segment, verbose=verbose)
        if st == 0:
            return False, align
        else:
            segment.aligned = found
            if st == -1:
                is_repeat = True
    if not is_repeat:
        return True, align
    #
    n_state = []
    for seg in sorted_segment:
        n_state.append(len(seg.aligned))
    i_state = [0 for i in range(len(n_state))]
    while(check_state(n_state, i_state)):
        success = True
        align = copy.deepcopy(align0)
        for k, segment in enumerate(sorted_segment):
            st, found = align.align_segment(segment, verbose=verbose)
            if st == 0: # candidate disappeared
                i_state[k-1] += 1
                success = False
                break
            elif st == -1:
                align.place_segment(segment.aligned[i_state[k]][0], \
                        segment.get_sequence(), segment)
        if success: break
    if success:
        return True, align
    else:
        return False, align

def align_structure_sequence(fasta, pdb, verbose=False):
    align = SegmentAlign(fasta)
    return build_align(align, pdb, verbose=verbose)

def align_sequence_sequence(fasta, seq, structure=None, verbose=False):
    align = SegmentAlign(fasta)
    if structure == None:
        pdb = read_sequence_chunk(seq)
        align.sequence = seq
    else:
        pdb, updated_seq = read_sequence_chunk(seq, structure=structure)
        align.sequence = updated_seq
    return build_align(align, pdb, verbose=verbose)

def check_state(n_state, i_state):
    for k in range(len(n_state)):
        if i_state[k] >= n_state[k]:
            return False
    return True
