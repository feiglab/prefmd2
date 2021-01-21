#!/usr/bin/env python

import os
import sys
import argparse
import subprocess as sp

from path import name
from hybrid.seqName import to_one_letter

PDBDIR = "none"  # '/blue/pdb'

def read_SEQRES(pdb_fn):
    seq = {}
    chain_s = []
    with open(pdb_fn) as fp:
        for line in fp:
            if not line.startswith("SEQRES"):
                continue
            chain_id = line[10:13].strip()
            if chain_id not in seq:
                seq[chain_id] = []
                chain_s.append(chain_id)
            seq[chain_id].extend(line[18:].strip().split())
    return seq, chain_s

def read_ATOMS(pdb_fn):
    seq = {}
    chain_s = []
    with open(pdb_fn) as fp:
        for line in fp:
            if not (line.startswith("ATOM") or line.startswith("HETA")):
                if line.startswith("END"):
                    break
                continue
            if line[12:16].strip() != 'CA':
                continue
            if line[16] not in [' ', 'A']:
                continue
            chain_id = line[21].strip()
            if chain_id not in seq:
                seq[chain_id] = []
                chain_s.append(chain_id)
            seq[chain_id].append(line[17:20].strip())
    return seq, chain_s

def main(args=None):
    arg = argparse.ArgumentParser(prog='exec_pdb_seq')
    arg.add_argument(dest='pdb_fn', metavar='PDB', help='input PDB file')
    arg.add_argument('-o', '--output', dest='fa_fn', \
                     default=None,\
                     help='output FA file')
    arg.add_argument('-s', '--seqres', dest='seqres',\
                     default=False, action='store_true',\
                     help='use SEQRES info')
    arg.add_argument('--std', dest='std',\
                     default=False, action='store_true',\
                     help='No conversion to STDRES name')
    arg.add_argument('-i', '--ignore', dest='ignore',\
                     default=False, action='store_true',\
                     help='ignore Non-standard residues')
    arg.add_argument('-p', '--pdbid', dest='pdbid',\
                     default=False, action='store_true',\
                     help='PDB ID as an input')
    #
    arg = arg.parse_args(args)
    #
    if arg.pdbid:
        pdb_fn_local = '%s/pdb%s.ent'%(PDBDIR, arg.pdb_fn[:4])
        if os.path.exists(pdb_fn_local):
            arg.pdb_fn = pdb_fn_local
        else:
            sp.call(["wget", "-q", "https://files.rcsb.org/download/%s.pdb"%(arg.pdb_fn.lower()[:4])],
                     stderr=sp.DEVNULL)
            arg.pdb_fn = '%s.pdb'%(arg.pdb_fn.lower()[:4])
    #
    if arg.seqres:
        seq,chain_s = read_SEQRES(arg.pdb_fn)
    else:
        seq,chain_s = read_ATOMS(arg.pdb_fn)
    #
    for chain_id in chain_s:
        seq[chain_id] = ''.join([to_one_letter(aa, std=not arg.std) for aa in seq[chain_id]])
        if arg.ignore:
            seq[chain_id] = seq[chain_id].replace("X",'')
    #
    title = name(arg.pdb_fn)
    output = []
    for chain_id in chain_s:
        if chain_id != '' and '_' not in title and len(chain_s) != 1:
            output.append(">%s_%s\n"%(title, chain_id))
        else:
            output.append(">%s\n"%(title))
        #
        output.append("%s\n"%seq[chain_id])
    #
    if arg.fa_fn is None:
        sys.stdout.writelines(output)
    else:
        with open(arg.fa_fn, 'wt') as fout:
            fout.writelines(output)

if __name__=='__main__':
    main()
