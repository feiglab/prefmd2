import os
import shutil
import sys

import warnings
warnings.filterwarnings("ignore")

import numpy as np
from sklearn.decomposition import PCA

import mdtraj
from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *

import path
import libcommon

from libquat import Quaternion


def solvate_pdb(output_prefix, pdb, options, verbose):
    orient_fn = path.Path('%s.orient.pdb'%output_prefix)
    if not orient_fn.status():
        xyz = pdb.xyz[0]
        xyz -= xyz.mean(axis=0)
        #
        for i in range(2):
            pca = PCA(n_components=(i+1))
            pca.fit(xyz)
            #
            axis_0 = pca.components_[i]
            axis_1 = np.zeros(3, dtype=float)
            axis_1[i] = 1.
            if np.all(axis_0 == axis_1):
                continue
            #
            axis_r = np.cross(axis_0, axis_1)
            angle_r = np.arccos(np.dot(axis_0, axis_1))
            #
            q = Quaternion.from_axis_and_angle(axis_r, angle_r)
            #
            xyz = np.dot(xyz, q.rotate().T)
        #
        solvateBox = options['md'].get("solvateBox", 'rectangular')
        if solvateBox == 'cubic':
            system_size = np.ones(3)*max(np.max(xyz, axis=0) - np.min(xyz, axis=0))
        else:
            system_size = np.max(xyz, axis=0) - np.min(xyz, axis=0)
        system_size += 2.0*(options['md']['solvate'] * 0.1)
        translate = (system_size - (np.max(xyz, axis=0) + np.min(xyz, axis=0))) * 0.5
        xyz += translate
        pdb.xyz[0] = xyz
        #
        pdb.unitcell_vectors = (system_size * np.eye(3))[None,:]
        pdb.save(orient_fn.short())
    #
    solv_fn = path.Path('%s.solvate.pdb'%output_prefix)
    if not solv_fn.status():
        cmd = []
        cmd.append(orient_fn.short())
        cmd.append(solv_fn.short())
        cmd.append("%8.5f"%options['md']['ion_conc'])
        libcommon.asystem(module="exec_solvate", args=cmd)
        #
        cmd = []
        cmd.append(solv_fn.short())
        libcommon.asystem(module="exec_update_water_name", args=cmd)
        #
        if 'use_modified_CMAP' in options['ff'] and options['ff']['use_modified_CMAP']:
            cmd = [solv_fn.short()]
            output = libcommon.asystem(module="exec_resName_modified_CMAP",
                                       args=cmd, stdout=True)
        else:
            with solv_fn.open() as fp:
                output = fp.read()
        #
        with solv_fn.open('wt') as fout:
            if 'ssbond' in options:
                for line in options['ssbond']:
                    fout.write("%s\n"%line)
            fout.write(output)
    #
    return orient_fn, solv_fn

def generate_PSF(output_prefix, solv_fn, options, verbose, fix_crd=True):
    psf_fn = path.Path("%s.psf"%output_prefix)
    crd_fn = path.Path("%s.crd"%output_prefix)
    if psf_fn.status() and crd_fn.status():
        return psf_fn, crd_fn
    #
    cmd = []
    cmd.append(solv_fn.short())
    cmd.extend(['-psf', psf_fn.short()])
    cmd.extend(['-crd', crd_fn.short()])
    cmd.append("--toppar")
    cmd.extend(options['ff']['toppar'])
    if 'blocked' in options['ff'] and options['ff']['blocked']:
        cmd.append("--blocked")
        if 'terminal' in options['ff']:
            cmd.append("--terminal")
            cmd.extend(options['ff']['terminal'])
    if 'patch' in options['ff']:
        cmd.append("--patch")
        cmd.extend(options['ff']['patch'])
    if 'ligand' in options and len(options['ligand']['str_fn_s']) > 0:
        cmd.extend(options['ff']['cgenff'])
        cmd.extend([fn.short() for fn in options['ligand']['str_fn_s']])
    libcommon.asystem(module="exec_genPSF", args=cmd)
    #**************************************************************************
    # Fix the crd file, since OpenMM has troubles with crd files with too many
    # atoms.
    if fix_crd:
        # Backups and reads the original crd file.
        shutil.copy(crd_fn.path(), crd_fn.path() + ".bak")
        with open(crd_fn.path(), "r") as i_fh:
            crd_lines = i_fh.readlines()
        # Make a modified copy of the crd file.
        with open(crd_fn.path(), "w") as o_fh:
            found_crds = False
            n_atoms = None
            for l in crd_lines:
                if l.startswith("*"):
                    o_fh.write(l)
                    continue
                if not found_crds:  # Get the number of atoms.
                    n_atoms_string = l.rstrip()
                    found_crds = True
                    o_fh.write(l)
                else:  # Coordinate lines.
                    # Add some spaces between the atom number and residue
                    # number. The OpenMM parser just uses .split().
                    o_fh.write(l[:len(n_atoms_string)] + "   " + \
                                l[len(n_atoms_string):])
    #**************************************************************************
    return psf_fn, crd_fn

def update_residue_name(pdb_fn, pdb):
    resName_s = ['HIS', 'HSP', 'HSD', 'HSE']
    residue_s = {} ; chain_s = []
    with pdb_fn.open() as fp:
        for line in fp:
            if not line.startswith("ATOM") and not line.startswith("HETA"):
                continue
            resName = line[17:20]
            if resName not in resName_s:
                continue
            atmName = line[12:16].strip()
            if atmName != 'CA':
                continue
            chain = line[21].strip()
            if chain not in chain_s:
                chain_s.append(chain)
            resNo = line[22:27].strip()
            chain_index = chain_s.index(chain)
            try:
                residue_s[(chain_index, int(resNo))] = resName
            except:
                residue_s[(chain_index, resNo)] = resName
    #
    for residue in pdb.top.residues:
        chain_index = residue.chain.index
        resNo = residue.resSeq
        if (chain_index, resNo) in residue_s:
            resName = residue_s[(chain_index, resNo)]
            residue.name = resName

def construct_restraint(psf, pdb, force_const):
    rsr = CustomExternalForce("k0*d^2 ; d=periodicdistance(x,y,z, x0,y0,z0)")
    rsr.addPerParticleParameter("x0")
    rsr.addPerParticleParameter("y0")
    rsr.addPerParticleParameter("z0")
    rsr.addPerParticleParameter('k0')
    #
    calphaIndex = []
    for i,atom in enumerate(psf.topology.atoms()):
        if atom.name == 'CA':
            calphaIndex.append(i)
    #
    k = -1
    for i,atom in enumerate(pdb.top.atoms):
        if atom.name != 'CA': continue
        #
        k += 1
        mass = atom.element.mass
        param = pdb.xyz[0,i].tolist()
        param.append(force_const * mass * kilocalories_per_mole/angstroms**2)
        rsr.addParticle(calphaIndex[k], param)
    return rsr

def construct_membrane_restraint(psf, pdb, force_const):
    rsr = CustomExternalForce("k0*d^2 ; d=periodicdistance(x,y,z, x0,y0,z0)")
    rsr.addPerParticleParameter("x0")
    rsr.addPerParticleParameter("y0")
    rsr.addPerParticleParameter("z0")
    rsr.addPerParticleParameter('k0')
    #
    heavyIndex = []
    for i,atom in enumerate(psf.topology.atoms()):
        #if atom.element.mass > 4.0*amu:
        if atom.name == 'P':
            heavyIndex.append(i)
    #
    k = -1
    for i,atom in enumerate(pdb.top.atoms):
        #if atom.element.mass < 4.0:
        if atom.name != 'P':
            continue
        #
        k += 1
        mass = atom.element.mass
        param = pdb.xyz[0,i].tolist()
        param.append(force_const * mass * kilocalories_per_mole/angstroms**2)
        rsr.addParticle(heavyIndex[k], param)
    return rsr

def construct_ligand_restraint(pair_s):
    bond = CustomBondForce("k * (r-r0)^2")
    bond.addPerBondParameter('k')
    bond.addPerBondParameter('r0')
    #
    for pair in pair_s:
        bond.addBond(pair[0], pair[1], \
                (pair[2]*kilocalories_per_mole/angstroms**2, pair[3]*nanometers))
    return bond
