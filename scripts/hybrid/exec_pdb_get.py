#!/usr/bin/env python

import os
import sys
import argparse
import subprocess as sp
from string import ascii_uppercase, ascii_lowercase, digits

from hybrid.seqName import convert_to_stdres_line

chainID_s = ascii_uppercase+ascii_lowercase+digits
PDBDIR = "none"  # '/blue/pdb'


def parse_pdb(pdb_fn, chain_s, convert_to_std=False):
    out_s = {}
    for chain in chain_s:
        out_s[chain] = []
    #
    n_disu = 0
    with open(pdb_fn) as fp:
        for line in fp:
            if line.startswith("SEQRES"):
                chain_id = line.strip().split()[2]
                if chain_id in chain_s:
                    out_s[chain_id].append(line)
            elif line.startswith("SSBOND"):
                x = line.strip().split()
                chain_1 = x[3] ; chain_2 = x[6]
                if (chain_1 not in chain_s) or (chain_2 not in chain_s):
                    continue
                n_disu += 1
                line = 'SSBOND%4d'%n_disu + line[10:]
                out_s[chain_1].append(line)
                if chain_2 != chain_1:
                    out_s[chain_2].append(line)
            elif line.startswith("ATOM") or line.startswith("HETA"):
                chain_id = line[21]
                if chain_id in chain_s:
                    if convert_to_std:
                        line = convert_to_stdres_line(line)
                    out_s[chain_id].append(line)
            elif line.startswith("END"):
                break
    #
    for chain_id in chain_s:
        if len(out_s[chain_id]) != 0:
            with open("%s_%s.pdb"%(pdb_fn[:4], chain_id), 'wt') as fout:
                fout.writelines(out_s[chain_id])

def main(args=None):
    arg = argparse.ArgumentParser(prog='exec_pdb_get')
    arg.add_argument(dest='pdb_id_s', metavar='PDB ID', help='input PDB ID', nargs='+')
    arg.add_argument('-f', '--force', dest='force_download',\
                     default=False, action='store_true',\
                     help='forced download')
    arg.add_argument('-m', '--mmcif', dest='mmCIF',\
                     default=False, action='store_true',\
                     help='download mmCIF format PDB file')
    arg.add_argument('-o', '--original', dest='leave_original',\
                     default=False, action='store_true',\
                     help='leave the original PDB file after parsing')
    arg.add_argument('--std', dest='convert_to_std', \
                     default=False, action='store_true', \
                     help='convert to std resname')
    arg = arg.parse_args(args)
    #
    for pdb_id in arg.pdb_id_s:
        id = pdb_id[:4].lower()
        chain_s = [chain_id for chain_id in pdb_id[4:] if chain_id in chainID_s]
        #
        if arg.mmCIF:
            try:
                sp.call(['wget', '-q', 'https://files.rcsb.org/download/%s.cif'%id], stderr=sp.DEVNULL)
            except:
                sp.call(['wget', '-q', 'https://files.rcsb.org/download/%s.cif'%id])
            continue
        #
        pdb_fn_local = '%s/pdb%s.ent'%(PDBDIR, id)
        if os.path.exists(pdb_fn_local) and not arg.force_download:
            sp.call("cp %s %s.pdb"%(pdb_fn_local, id), shell=True)
        else:
            try:
                sp.call(['wget', '-q', 'https://files.rcsb.org/download/%s.pdb'%id], stderr=sp.DEVNULL)
            except:
                sp.call(['wget', '-q', 'https://files.rcsb.org/download/%s.pdb'%id])
        #
        if not os.path.exists("%s.pdb"%id):
            sys.stderr.write("Error: failed to get %s\n"%id)
            continue
        #
        if len(chain_s) != 0:
            parse_pdb('%s.pdb'%id, chain_s, convert_to_std=arg.convert_to_std)
            # if not arg.leave_original:
            #     os.remove("%s.pdb"%id)

if __name__=='__main__':
    main()
