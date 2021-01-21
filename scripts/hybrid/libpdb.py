import os
import sys
import tempfile
import importlib
import numpy as np
from string import digits

import path
from hybrid.seqName import (to_one_letter, to_three_letter,
                            convert_to_stdres_line)


PDBfmt = '%-6s%5d %-4s %3s %1s%5s   %8.3f%8.3f%8.3f  %4.2f%6.2f%s\n'


class Sequence(object):
    def __init__(self, title, chain_id='', res_no=None):
        self.title = title
        self.chain_id = chain_id
        self.fasta = ''
        self.res_no = res_no
    @classmethod
    def read(cls, fa_fn):
        seq = [] ; title = ''
        with open(fa_fn) as fp:
            for line in fp:
                if line.startswith(">"):
                    title = line.strip()[1:]
                else:
                    seq.append(line.strip())
        s = cls(title)
        s.fasta = ''.join(seq)
        return s
    def write(self):
        wrt = []
        if self.res_no != None:
            if self.chain_id != '':
                wrt.append(">%s_%s_%s_%s\n"%(self.title, self.chain_id, self.res_no[0], self.res_no[1]))
            else:
                wrt.append(">%s_%s_%s\n"%(self.title, self.res_no[0], self.res_no[1]))
        else:
            if self.chain_id.strip() != '':
                wrt.append(">%s_%s\n"%(self.title, self.chain_id))
            else:
                wrt.append(">%s\n"%self.title)
        wrt.append("%s\n"%self.fasta)
        return ''.join(wrt)
    @staticmethod
    def write_SEQRES(chain_id, sequence):
        n_in_line = 13
        #
        wrt = []
        for k,resName in enumerate(sequence):
            if k%n_in_line == 0:
                wrt.append("SEQRES%4d %s %4d "%(int(k/n_in_line)+1, chain_id, len(sequence)))
            wrt.append(" %3s"%to_three_letter(resName))
            if (k+1)%n_in_line == 0:
                wrt.append("\n")
        if len(wrt) != 0 and wrt[-1] != '\n':
            wrt.append("\n")
        return wrt

    # TODO: remove.
    # def trim_exp_tag(self):
    #     seq0 = self.fasta
    #     #
    #     if self.res_no is None:
    #         res_no = [1, len(seq0)]
    #     else:
    #         res_no = [int(self.res_no[0]), int(self.res_no[1])]
    #     #
    #     found = False ; seq = seq0
    #     for lib in read_EXP_TAG():
    #         l_seq = len(lib)
    #         if seq0[:l_seq] == lib:    # N-ter
    #             found = True
    #             seq = seq0[l_seq:]
    #             res_no[0] += l_seq
    #         elif seq0[-l_seq:] == lib: # C-ter
    #             found = True
    #             seq = seq0[:-l_seq]
    #             res_no[1] -= l_seq
    #         if found: break
    #     #
    #     res_no = ('%04d'%res_no[0], '%04d'%res_no[1])
    #     out = Sequence(self.title, chain_id=self.chain_id, res_no=res_no)
    #     out.fasta = seq
    #     return out

    # def trim_signal_peptide(self, EXEC='/home/huhlim/apps/signalp/bin/signalp', options=""):
    #     fa_fn = tempfile.mkstemp(prefix="signalP")[1]
    #     with open(fa_fn, 'wt') as fout:
    #         fout.write(self.write())
    #     #
    #     if self.res_no is None:
    #         res_no = [1, len(seq0)]
    #     else:
    #         res_no = [int(self.res_no[0]), int(self.res_no[1])]
    #     #
    #     EXEC = path.Path(EXEC)
    #     os.environ['PATH'] = os.environ['PATH'] + ':%s'%EXEC.dirname()
    #     #
    #     cmd = []
    #     cmd.append('signalp')
    #     cmd.extend(['-fasta', fa_fn])
    #     cmd.extend(['-stdout', '-verbose=false'])
    #     cmd.extend(options.split())
    #     #
    #     sp = importlib.import_module("subprocess")
    #     output = sp.check_output(cmd)
    #     if sys.version_info.major == 3:
    #         output = output.decode("utf8")
    #     os.remove(fa_fn)
    #     #
    #     signalP = None
    #     for line in output.split("\n"):
    #         if 'CS pos: ' in line:
    #             signalP = int(line.split("CS pos: ")[1].split()[0].split("-")[0])
    #     #
    #     if signalP is None:
    #         return self
    #     else:
    #         res_no = ('%04d'%(res_no[0]+signalP), '%04d'%res_no[1])
    #         out = Sequence(self.title, chain_id=self.chain_id, res_no=res_no)
    #         out.fasta = self.fasta[signalP:]
    #         return out

    def __len__(self):
        return len(self.fasta)
    def __repr__(self):
        return self.fasta
    def __str__(self):
        return self.__repr__()
    def __eq__(self, othr):
        return self.fasta == othr.fasta
    def __getitem__(self, key):
        return self.fasta[key]

# def read_EXP_TAG():
#     lib_s = []
#     with open("/home/huhlim/lib/EXP_TAG_s") as fp:
#         for line in fp:
#             lib_s.append(line.strip())
#     return lib_s

class PDB:
    def __init__(self, pdb_fn, title=None, \
            read=True, read_het=True, to_std=True, model_index=[], res_range=[]):

        if isinstance(pdb_fn, PDB):
            self.pdb_fn = pdb_fn.pdb_fn
        else:
            self.pdb_fn = path.Path(pdb_fn)
        #
        self.title = '.'.join(self.pdb_fn.split("/")[-1].split(".")[:-1])
        self.model_s = []
        self.use_model = False
        if read:
            self.read(model_index=model_index, res_range=res_range, \
                    read_het=read_het, to_std=to_std)
    def __repr__(self):
        return self.pdb_fn.__repr__()
    def __len__(self):
        return len(self.model_s)
    def __getitem__(self, i):
        return self.model_s[i]
    def read(self, read_het=True, to_std=True, model_index=[], res_range=[]):
        read_model = True
        self.model_s = []
        i_model = 0
        resNo_prev = None
        chain_prev = None
        #
        model = Model()
        self.model_s.append(model)
        #
        if 'seqres' in dir(self):
            self.seqres = {}
            self.seqres_lines = {}
        #
        with self.pdb_fn.open('r') as fp:
            lines = fp.readlines()
        #
        for line in lines:
            if line.startswith("MODEL"):
                self.use_model = True
                model_no = int(line.strip().split()[1])
                if len(model) != 0 and model.use:
                    i_model += 1
                if len(model_index) != 0 and (model_no not in model_index):
                    read_model = False
                    continue
                else:
                    read_model = True
                if len(model) != 0 and model.use:
                    model = Model(model_no=model_no)
                    self.model_s.append(model)
                else:
                    model.model_no = model_no
                model.append(PDBline(line))
            elif line.startswith("SEQRES"):
                if 'seqres' not in dir(self):
                    self.seqres = {}
                    self.seqres_lines = {}

                x = line.strip().split()
                chain_id = x[2]
                seq = ''.join([to_one_letter(xi) for xi in x[4:]])
                seq = seq.replace("X","")
                #
                if chain_id not in self.seqres:
                    self.seqres[chain_id] = Sequence(self.title, chain_id=chain_id)
                    self.seqres_lines[chain_id] = []
                self.seqres[chain_id].fasta = '%s%s'%(self.seqres[chain_id].fasta, seq)
                self.seqres_lines[chain_id].append(line)
            elif line.startswith("ATOM") or line.startswith("HETATM"):
                _std = True
                if not read_model: continue
                if line[16] not in ['A', ' ']:
                    continue
                if line.startswith("HETATM") and to_std:
                    line = convert_to_stdres_line(line)
                    _std = False
                if line.startswith("HETATM") and (not read_het):
                    continue
                #
                model.use = True
                #
                resNo = line[22:27]
                chain_id = line[21]
                #
                if resNo != resNo_prev or chain_id != chain_prev:
                    resNo_prev = resNo
                    chain_prev = chain_id
                    #
                    residue = Residue(line)
                    residue._std = _std
                    model.append(residue)
                residue.append(line)
                #
            elif read_model:
                model.append(PDBline(line))
        #
        for i_model in range(len(self)):
            self.model_s[i_model].title = self.title
            self.model_s[i_model].convert_to_nparray()
    def write(self, exclude_remark=False, exclude_symm=False, exclude_missing_bb=False,\
            model_index=[], exclude_nucl=False,\
            remark_s=[], exclude_unread=True, chain_id=None):
        if len(model_index) == 0:
            model_index = list(range(len(self)))
        #
        wrt = []
        wrt.extend(remark_s)
        #
        if 'seqres' in dir(self):
            for chainID in self.get_chainID():
                if chain_id != None and chainID != chain_id:
                    continue
                if chainID not in self.seqres_lines:
                    continue
                wrt.extend(self.seqres_lines[chainID])
        #
        model_no = 0
        for j in model_index:
            if self.use_model:
                model_no += 1
                wrt.append("MODEL %4d\n"%(model_no))
            wrt.extend(self.model_s[j].write(exclude_remark=exclude_remark, exclude_symm=exclude_symm,\
                                             exclude_missing_bb=exclude_missing_bb, chain_id=chain_id))
            if self.use_model:
                wrt.append("ENDMDL\n")
        wrt.append("END\n")
        return wrt
    def rewrite(self, exclude_remark=False, exclude_symm=False, model_index=[], backup=False, \
                exclude_missing_bb=False, remark_s=[], report_topN=0, exclude_unread=True, chain_id=None):
        import os
        if backup: os.system("cp %s %s.bak"%(self.pdb_fn, self.pdb_fn))
        fout = file("%s"%self.pdb_fn, 'wt')
        fout.writelines(self.write(exclude_remark=exclude_remark, exclude_symm=exclude_symm, \
                        exclude_missing_bb=exclude_missing_bb, model_index=model_index, report_topN=report_topN,\
                        exclude_unread=exclude_unread, chain_id=chain_id))
        fout.close()
    def write_domain(self, res_range=[], model_index=[], chain_id=None):
        wrt = []
        for i in range(len(self)):
            if len(model_index) != 0 and (i not in model_index):
                continue
            if self.use_model:
                wrt.append("MODEL %4d\n"%(i+1))
            wrt.extend(self.model_s[i].write_domain(res_range=res_range, chain_id=chain_id))
            if self.use_model:
                wrt.append("ENDMDL\n")
        wrt.append("END\n")
        return wrt
    def get_chainID(self, model_no=0):
        return self[model_no].get_chainID()
    def get_residue_number_and_chainID(self, model_no=0):
        return self[model_no].get_residue_number_and_chainID()
    def get_residue_number(self, model_no=0):
        return self[model_no].get_residue_number()
    def get_residue_name(self, model_no=0):
        return self[model_no].get_residue_name()
    def write_fasta_file(self, fa_fn, model_no=0, title=None, seqres=False):
        self.fa_fn = path.Path(fa_fn)
        if seqres and ('seqres' in dir(self)):
            chain_id_s = list(self.seqres.keys())
            chain_id_s.sort()
            fout = self.fa_fn.open('wt')
            for chain_id in chain_id_s:
                fout.write('%s'%self.seqres[chain_id].write())
            fout.close()
        else:
            self[model_no].write_fasta_file(fa_fn, title=title)
    def get_sequence(self, model_no=0, seqres=False):
        if seqres and ('seqres' in dir(self)):
            return self.seqres
        else:
            if 'fasta' not in dir(self):
                self.fasta = self[model_no].get_sequence()
            return self.fasta
    def write_seqres(self, seqres=True):
        n_in_line = 13
        self.seqres_lines = {}
        for chainID in self.get_chainID():
            self.seqres_lines[chainID] = []
            #
            seq = self.get_sequence(seqres=seqres)[chainID]
            for k,resName in enumerate(seq):
                if k%n_in_line == 0:
                    self.seqres_lines[chainID].append("SEQRES%4d %s %4d "%(int(k/n_in_line)+1, chainID, len(seq)))
                self.seqres_lines[chainID].append(' %3s'%to_three_letter(resName))
                if (k+1)%n_in_line == 0:
                    self.seqres_lines[chainID].append("\n")
            if self.seqres_lines[chainID][-1] != '\n':
                self.seqres_lines[chainID].append("\n")
        return self.seqres_lines
    def get_phipsi(self, res_range=[]):
        self.rama_s = []
        for model in self:
            self.rama_s.append(model.get_phipsi(res_range=res_range))
        return self.rama_s
    def prep_fix(self, sc=[]):
        pdb = self.model_s[0]
        resNo_s = pdb.get_residue_number()
        resName_s = pdb.get_residue_name()
        if len(sc) == 0:
            sc = resNo_s
        wrt = []
        for i_res in sc:
            if i_res not in resNo_s: continue
            i = resNo_s.index(i_res)
            if resName_s[i] in ['GLY','ALA','PRO']:
                continue
            wrt.append("%4d USC\n"%i_res)
        return wrt
    def extend(self, pdb):
        self.model_s.extend(pdb.model_s)

class Model:
    def __init__(self, model_no=0):
        self.use = False
        self.model_no = model_no
        self.lines = []
        self.res_index = {}
    def append(self, X):
        self.lines.append(X)
        if X.isResidue():
            self.res_index[(X.chainID(), X.resNo_char())] = len(self.lines)-1
    def index(self, key):
        return self.res_index[key]
    def __getitem__(self, i):
        return self.lines[i]
    def __len__(self):
        return len(self.get_residues())
    def write(self, exclude_remark=False, exclude_symm=False, exclude_missing_bb=False,\
            exclude_nucl=False, exclude_SSbond=False, remark_s=[], chain_id=None):
        wrt = []
        wrt.extend(remark_s)
        for line in self.lines:
            if line.isResidue():
                if chain_id != None and chain_id != line.chainID():
                    continue
                if exclude_nucl and line.isAtom() and \
                        (line.resName().strip() in ['DA','DC','DG','DT','DU','A','C','G','T','U']):
                    continue
                if exclude_missing_bb and (not line.check_bb()):
                    continue
                wrt.append('%s'%line)
            else:
                if line.startswith("MODEL"):
                    continue
                elif line.startswith("END"):
                    continue
                elif line.startswith("TER"):
                    if len(wrt) != 0 and (not wrt[-1].startswith("TER")):
                        wrt.append("TER\n")
                    continue
                if line.startswith('REMARK 350') and (not exclude_symm):
                    wrt.append('%s'%line)
                    continue
                elif line.startswith('SSBOND') and (not exclude_SSbond):
                    wrt.append('%s'%line)
                    continue
                elif exclude_remark:
                    continue
                wrt.append('%s'%line)
        return wrt
    def __repr__(self):
        return ''.join(self.write())
    def write_domain(self, res_range=[], chain_id=None):
        return ['%s'%line for line in \
                self.get_residue_lines(res_range=res_range, chain_id=chain_id)]
    def convert_to_nparray(self):
        i_atm = 0
        for residue in self:
            if residue.isResidue():
                i_atm = residue.convert_to_nparray(i_atm)
    def get_chainID(self):
        chainID_s = []
        for line in self.get_residues():
            chainID = line.chainID()
            if chainID not in chainID_s:
                chainID_s.append(chainID)
        return chainID_s
    def split_chain(self):
        pdb_s = {}
        for chainID in self.get_chainID():
            pdb_s[chainID] = []
        for residue in self.get_residues():
            pdb_s[residue.chainID()].append(residue)
        return pdb_s
    def split_chain_line(self):
        pdb_s = {}
        for chainID in self.get_chainID():
            pdb_s[chainID] = []
        for residue in self.get_residues():
            pdb_s[residue.chainID()].append('%s'%residue)
        return pdb_s
    def get_residue_number_and_chainID(self):
        chain_s = {}
        for res in self.get_residues():
            chainID = res.chainID()
            resNo = res.resNo()
            if chainID not in chain_s:
                chain_s[chainID] = []
            chain_s[chainID].append(resNo)
        return chain_s
    def get_residue_number(self):
        resNo_s = []
        for line in self.get_residues():
            resNo_s.append(line.resNo())
        return resNo_s
    def get_residue_name(self):
        sequence_s = []
        for line in self.get_residues():
            sequence_s.append(line.resName())
        return sequence_s
    def write_fasta_file(self, fa_fn, title=None):
        self.fa_fn = fa_fn
        seq = self.get_sequence(title=title)
        fout = file('%s'%fa_fn, 'wt')
        for chain_id in self.get_chainID():
            fout.write('%s'%seq[chain_id].write())
        fout.close()
    def get_sequence(self, title=None):
        if title == None: title = self.title
        sequence_s = {}
        for line in self.get_residues():
            if line.isHetatm(): continue
            chainID = line.chainID()
            if chainID not in list(sequence_s.keys()):
                sequence_s[chainID] = Sequence(title, chain_id=chainID)
            sequence_s[chainID].fasta = '%s%s'%\
                    (sequence_s[chainID].fasta, to_one_letter(line.resName()))
        return sequence_s
    def get_atoms(self, atmName=None, update=False):
        if ('_atoms' not in dir(self) or update):
            if '_atoms' not in dir(self):
                self._atoms = {}
                self._atoms[atmName] = []
        elif atmName in self._atoms:
            return self._atoms[atmName]
        else:
            self._atoms[atmName] = []
        for line in self.lines:
            if not line.isAtom():
                continue
            self._atoms[atmName].append(Atom(line[atmName]))
        return self._atoms[atmName]
    def get_atom_lines(self, atmName=None):
        if '_atom_lines' not in dir(self):
            self._atom_lines = {}
            self._atom_lines[atmName] = []
        elif atmName in self._atom_lines:
            return self._atom_lines[atmName]
        else:
            self._atom_lines[atmName] = []
        for line in self.lines:
            if not line.isResidue():
                continue
            self._atom_lines[atmName].append(line[atmName])
        return self._atom_lines[atmName]
    def get_residues(self, res_range=[], chain_id=None):
        lines = []
        for line in self.lines:
            if not line.isResidue():
                continue
            if len(res_range) != 0 and \
                    (line.resNo() not in res_range) and \
                    (line.resNo_char() not in res_range):
                continue
            if chain_id is not None and line.chainID() != chain_id:
                continue
            lines.append(line)
        return lines
    def get_residue_lines(self, res_range=[], chain_id=None):
        lines = []
        for line in self.get_residues(res_range=res_range, chain_id=chain_id):
            lines.append('%s'%line)
        return lines
    def R(self, atmName=None, atmIndex=None, res_range=[]):
        _R = []
        if atmName != None:
            for residue in self.get_residues(res_range=res_range):
                _R.append([residue.R(atmName=atmName)])
        else:
            for residue in self.get_residues(res_range=res_range):
                _R.append(residue._R)
        return np.concatenate(_R)
    def tr(self, T, R):
        for line in self.lines:
            if line.isResidue():
                line.tr(T,R)

class Residue:
    def __init__(self, line):
        self._diso = False
        self._std = True
        self._header = line[:6]
        self._resName = line[17:20]
        self._resNo = line[22:27]
        self._chainID = line[21]
        #
        self._R = []
        self._i_atm = []
        self._atmName = []
        self._occ = []
        self._bfac = []
        self._tail = []
    def __len__(self):
        return len(self._R)
    def append(self, line):
        atmName = line[12:16].strip()
        if len(atmName) == 4 and atmName[0] in digits:
            atmName = '%s%s'%(atmName[1:], atmName[0])
        self._atmName.append(atmName)
        self._i_atm.append(int(line[6:11]))
        self._R.append((float(line[30:38]),\
                        float(line[38:46]),\
                        float(line[46:54])))
        try:
            self._occ.append(float(line[54:60]))
        except:
            self._occ.append(1.0)
        try:
            self._bfac.append(float(line[60:66]))
        except:
            self._bfac.append(0.0)
        try:
            self._tail.append(line[66:].rstrip())
        except:
            self._tail.append("")
    def convert_to_nparray(self, i_atm):
        self._R = np.array(self._R)
        return self.update_i_atm(i_atm)
    def isResidue(self):
        return True
    def isAtom(self):
        return self._header[:4] == 'ATOM'
    def isHetatm(self):
        return self._header == 'HETATM'
    def exists(self, atmName):
        return atmName in self._atmName
    def isStd(self):
        return self._std
    def check_bb(self):
        stat = [False, False, False, False]
        if 'N' in self._atmName:
            stat[0] = True
        if 'CA' in self._atmName:
            stat[1] = True
        if 'C' in self._atmName:
            stat[2] = True
        if 'O' in self._atmName:
            stat[3] = True
        #
        if False in stat:
            return False
        else:
            return True
    def write(self):
        wrt = []
        if not self._diso:
            for i in range(len(self._R)):
                wrt.append(self[i])
        return ''.join(wrt)
    def __repr__(self):
        return self.write()
    def __getitem__(self, i):
        if isinstance(i, str):
            i = self.atmIndex(i)
        if len(self._atmName[i]) == 4:
            if self._atmName[i][-1] in digits:
                atmName = '%s%s'%(self._atmName[i][-1],self._atmName[i][:3])
            else:
                atmName = '%s'%self._atmName[i]
        else:
            atmName = ' %s'%self._atmName[i]
        line = PDBfmt%(self._header, self._i_atm[i], atmName,\
                       self._resName, self._chainID, self._resNo,\
                       self._R[i][0], self._R[i][1], self._R[i][2],\
                       self._occ[i], self._bfac[i], self._tail[i])
        return line
    def resName(self):
        return self._resName
    def resNo(self):
        return int(self._resNo[:4])
    def resNo_char(self):
        return self._resNo
    def chainID(self):
        return self._chainID
    def atmName(self):
        return self._atmName
    def bfac(self):
        return self._bfac
    def i_atm(self, atmName=None, atmIndex=None):
        if atmName != None:
            return self._i_atm[self.atmIndex(atmName)]
        elif atmIndex != None:
            return self._i_atm[atmIndex]
        else:
            return self._i_atm
    def R(self, atmName=None, atmIndex=None):
        if atmName != None:
            return self._R[self.atmIndex(atmName)]
        elif atmIndex != None:
            return self._R[atmIndex]
        else:
            return self._R
    def atmIndex(self, atmName):
        return self._atmName.index(atmName)
    def get_backbone(self):
        return [self._atmName.index("N"), self._atmName.index("CA"),\
                self._atmName.index("C"), self._atmName.index("O")]
    def get_heavy(self):
        heavy = []
        for i,atm in enumerate(self._atmName):
            if atm[0] != 'H':
                heavy.append(i)
        return heavy
    def get_CB(self):
        if self._resName == 'GLY':
            return [self._atmName.index("CA")]
        else:
            return [self._atmName.index("CB")]
    def update_i_atm(self, i_atm):
        self._i_atm = list(range(i_atm+1, i_atm+len(self._i_atm)+1))
        return self._i_atm[-1]
    def tr(self, T, R):
        for i in range(len(self)):
            self._R[i] = R.dot(self._R[i]) + T

class Atom(Residue):
    def __init__(self, line):
        self._header = line[:6]
        self._resName = line[17:20]
        self._resNo = line[22:27]
        self._chainID = line[21]
        #
        atmName = line[12:16].strip()
        if len(atmName) == 4 and atmName[0] in digits:
            atmName = '%s%s'%(atmName[1:], atmName[0])
        self._atmName = atmName
        self._i_atm = int(line[6:11])
        self._R = np.array((float(line[30:38]),\
                            float(line[38:46]),\
                            float(line[46:54])))
        self._occ = float(line[54:60])
        self._bfac = float(line[60:66])

    def __repr__(self):
        if len(self._atmName) == 4:
            if self._atmName[-1] in digits:
                atmName = '%s%s'%(self._atmName[-1],self._atmName[:3])
            else:
                atmName = '%s'%self._atmName
        else:
            atmName = ' %s'%self._atmName
        line = PDBfmt%(self._header, self._i_atm, atmName,\
                   self._resName, self._chainID, self._resNo,\
                   self._R[0], self._R[1], self._R[2],\
                   self._occ, self._bfac)
        return line
    def R(self):
        return self._R
    def i_atm(self):
        return self._i_atm
    def atmName(self):
        return self._atmName
    def occ(self):
        return self._occ
    def bfac(self):
        return self._bfac

class PDBline:
    def __init__(self, line):
        self.line = line
    def __repr__(self):
        return self.line
    def isResidue(self):
        return False
    def isAtom(self):
        return False
    def isHetatm(self):
        return False
    def startswith(self, key):
        return self.line.startswith(key)
