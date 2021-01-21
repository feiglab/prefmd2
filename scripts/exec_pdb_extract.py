#!/usr/bin/env python

import os
import sys
import argparse
import numpy as np
import mdtraj

def main(args=None):
    arg = argparse.ArgumentParser(prog='extract_pdb')
    arg.add_argument(dest='ref_fn')
    arg.add_argument('--dcd', dest='dcd_fn_s', nargs='*', default=[])
    arg.add_argument('--dir', dest='out_dir', default='.')
    arg.add_argument('--name', dest='out_prefix', default='sample')
    arg.add_argument('--select', dest='selection', default=None)
    arg.add_argument('--topIndex', dest='topIndex_fn', default=None)
    arg.add_argument('--stride', dest='stride', default=None, type=int)
    arg.add_argument('--frame', dest='frame', default=None, type=int, nargs='*')
    arg.add_argument('--structured', dest='structured', default=False, action='store_true')
    arg.add_argument('--split_chain', dest='split_chain', default=False, action='store_true')
    arg.add_argument('--center', dest='center', default=False, action='store_true')
    arg.add_argument('--init_number', dest='init_number', default=0, type=int)
    if len(sys.argv) == 1:
        return arg.print_help()
    arg = arg.parse_args(args)
    #
    top = mdtraj.load(arg.ref_fn)
    if arg.topIndex_fn is not None:
        topIndex = np.load(arg.topIndex_fn)['select']
        top = top.atom_slice(topIndex)
    if arg.selection is None:
        atomIndex = None
    else:
        atomIndex = top.topology.select(arg.selection)
    if arg.split_chain:
        n_chain = top.topology.n_chains
    #
    if not os.path.exists(arg.out_dir):
        os.makedirs(arg.out_dir)
    #
    model_no = arg.init_number-1
    for dcd_fn in arg.dcd_fn_s:
        pdb_s = mdtraj.load(dcd_fn, top=top, stride=arg.stride, atom_indices=atomIndex)
        if arg.frame is not None:
            pdb_s = [pdb for i,pdb in enumerate(pdb_s) if i in arg.frame]
        for pdb in pdb_s:
            model_no += 1
            if arg.structured:
                out_dir = '%s/%d'%(arg.out_dir, int(model_no/100))
                if not os.path.exists(out_dir):
                    os.makedirs(out_dir)
                if arg.split_chain:
                    out_fn_s = ['%s/%s.%d.%d.pdb'%(out_dir, arg.out_prefix, model_no, chain_no)\
                            for chain_no in range(n_chain)]
                else:
                    out_fn = '%s/%s.%d.pdb'%(out_dir, arg.out_prefix, model_no)
            else:
                if arg.split_chain:
                    out_fn_s = ['%s/%s.%d.%d.pdb'%(arg.out_dir, arg.out_prefix, model_no, chain_no)\
                            for chain_no in range(n_chain)]
                else:
                    out_fn = '%s/%s.%d.pdb'%(arg.out_dir, arg.out_prefix, model_no)
            if arg.center:
                pdb.center_coordinates(mass_weighted=True)
            if arg.split_chain:
                for i in range(n_chain):
                    sys.stdout.write("%s\n"%out_fn_s[i])
                    atomIndex_byChain = pdb.top.select("chainid %d"%i)
                    pdb_per_chain = pdb.atom_slice(atomIndex_byChain)
                    pdb_per_chain.save_pdb(out_fn_s[i])
            else:
                sys.stdout.write("%s\n"%out_fn)
                pdb.save_pdb(out_fn)

if __name__ == '__main__':
    main()
