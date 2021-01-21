#!/usr/bin/env python

import os
import sys
import mdtraj
from mdtraj.core.element import sodium, chlorine
import numpy as np
from multiprocessing import Pool
from string import ascii_uppercase
import argparse

from libcommon import DEFAULT_HOME

import path
from libmd import update_residue_name

CHAINs = ascii_uppercase
AVOGADRO = 6.0221420e23

water_pdb_fn = '%s/water.pdb'%DEFAULT_HOME
water_box = 1.8662

def calc_dist(arg):
    # prot.shape = (n_atom, 3)
    # wat.shape = (n_wat, 3, 3)
    prot = arg[0] ; wat = arg[1]
    dr = prot[None,None,:,:] - wat[:,:,None,:]
    d = np.sqrt(np.sum(dr**2, axis=-1))
    return d

def place_water(pdb, water, water_cutoff):
    boxsize = pdb.unitcell_lengths[0]
    xyz = pdb.xyz[0].copy()
    for i in range(3):
        xyz[xyz[:,i] > boxsize[i]] -= boxsize[i]
        xyz[xyz[:,i] < 0.] += boxsize[i]
    #
    water_s = []
    n_box = np.ceil(boxsize / water_box).astype(int)
    for i in range(n_box[0]):
        for j in range(n_box[1]):
            for k in range(n_box[2]):
                wat = water + np.array([i,j,k], dtype=np.float32) * water_box
                for w in range(3):
                    wat = wat[wat[:,0,w] < boxsize[w]]
                water_s.append(wat)
    #
    proc = Pool(min(12, len(water_s)))
    dist_s = proc.map(calc_dist, [(xyz, wat) for wat in water_s])
    dist_s = np.concatenate(dist_s)
    proc.close()
    #
    placed = np.concatenate(water_s)[np.all(dist_s > water_cutoff, axis=(1,2))]
    np.random.shuffle(placed)
    #
    return placed

def place_ions(pdb, ion_conc, net_charge):
    for res in ['ASP', 'GLU']:
        net_charge -= pdb.top.select("resname %s and name CA"%res).shape[0]
    for res in ['LYS', 'ARG']:
        net_charge += pdb.top.select("resname %s and name CA"%res).shape[0]
    unitCell = pdb.unitcell_lengths[0]
    volume = unitCell[0]*unitCell[1]*unitCell[2]
    #
    n_ion = [0, 0]
    if net_charge < 0:
        n_ion[0] = abs(net_charge)
    else:
        n_ion[1] = abs(net_charge)
    conc = int(np.round(ion_conc / (1/AVOGADRO) * volume * 1e-24))
    n_ion[0] += conc ; n_ion[1] += conc
    return n_ion

def xyz_to_pdb(pdb0, water0, placed, n_ion, n_water):
    top = pdb0.top.copy()
    xyz = pdb0.xyz[0]
    #
    segNo = 0 ; resNo_prev = None
    for i,chain in enumerate(top.chains):
        for residue in chain.residues:
            if residue.segment_id == '':
                residue.segment_id='P%03d'%segNo
            #
            chain_break = False
            for atom in residue.atoms:
                if atom.name in ['OXT', 'OT2', 'OT1']:
                    chain_break = True
            if chain_break:
                segNo += 1
    #
    ia = 0
    chain = top.add_chain() ; resNo = 0
    #
    segName = 'SOD0'
    for i in range(n_ion[0]):
        if chain.n_residues >= 9999:
            chain = top.add_chain()
            resNo = 0
        #
        resNo += 1
        residue = top.add_residue('SOD', chain, resSeq=resNo, segment_id=segName)
        atom = top.add_atom('SOD', sodium, residue)
    xyz = np.append(xyz, placed[ia:ia+n_ion[0],0], axis=0)
    ia += n_ion[0]
    #
    segName = 'CLA0'
    for i in range(n_ion[1]):
        if chain.n_residues >= 9999:
            chain = top.add_chain()
            resNo = 0
        #
        resNo += 1
        residue = top.add_residue('CLA', chain, resSeq=resNo, segment_id=segName)
        atom = top.add_atom('CLA', chlorine, residue)
    xyz = np.append(xyz, placed[ia:ia+n_ion[1],0], axis=0)
    ia += n_ion[1]
    #
    water_residue = water0.top.residue(0)
    water_residue.name = 'TIP3'
    water_residue.atom(0).name = 'OH2'
    water_residue.atom(1).name = 'H1'
    water_residue.atom(2).name = 'H2'
    #
    water_bond = []
    for b in water0.top.bonds:
        if b[0].index >= water_residue.n_atoms:
            continue
        if b[1].index >= water_residue.n_atoms:
            continue
        water_bond.append((b[0].index, b[1].index))
    #
    water_index = 0
    segName = 'W%03d'%water_index
    for i in range(n_water):
        if chain.n_residues >= 9999:
            chain = top.add_chain()
            resNo = 0
            water_index += 1
            segName = 'W%03d'%water_index
        #
        atom_s = {}
        resNo += 1
        residue = top.add_residue(water_residue.name, chain, resSeq=resNo, segment_id=segName)
        for a in water_residue.atoms:
            atom = top.add_atom(a.name, a.element, residue)
            atom_s[a.index] = atom
        #
        for b in water_bond:
            atom_1 = atom_s[b[0]]
            atom_2 = atom_s[b[1]]
            top.add_bond(atom_1, atom_2)
    xyz = np.append(xyz, np.concatenate(placed[ia:]), axis=0)
    #
    pdb = mdtraj.Trajectory(xyz[None,:], top, \
            unitcell_lengths=pdb0.unitcell_lengths,\
            unitcell_angles=pdb0.unitcell_angles)
    #
    return pdb

def main(args=None):
    arg = argparse.ArgumentParser(prog='solvate')
    arg.add_argument(dest="in_pdb")
    arg.add_argument(dest="out_pdb")
    arg.add_argument(dest='ion_conc', type=float)
    arg.add_argument('--charge', dest='net_charge', type=int, default=0)
    arg.add_argument('--cutoff', dest='water_cutoff', type=float, default=2.)
    #
    if len(sys.argv) == 1:
        arg.print_help()
        return
    arg = arg.parse_args(args)
    #
    pdb0 = mdtraj.load(arg.in_pdb)
    update_residue_name(path.Path(arg.in_pdb), pdb0)
    #
    water0 = mdtraj.load(water_pdb_fn)
    water = water0.xyz[0].reshape((-1, 3, 3)) + water_box/2.0
    #
    placed = place_water(pdb0, water, arg.water_cutoff*0.1)
    n_ion = place_ions(pdb0, arg.ion_conc, arg.net_charge)
    n_water = placed.shape[0] - sum(n_ion)
    #
    pdb = xyz_to_pdb(pdb0, water0, placed, n_ion, n_water)
    pdb.save(arg.out_pdb)

if __name__ == '__main__':
    main()
