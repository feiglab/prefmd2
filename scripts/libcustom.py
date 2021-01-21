#!/usr/bin/env python

from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *

class Custom:
    def __init__(self, r_type, n_atom, n_param, i_atm, prm):
        self.r_type = r_type
        self.n_atom = n_atom
        self.n_param = n_param
        self.i_atm = i_atm
        self.prm = prm

def construct_custom_restraint(ref, custom_s):
    rsr = []
    if len(custom_s) == 0:
        return rsr
    #
    crd = ref.xyz[0]
    mass = [atom.element.mass for atom in ref.topology.atoms]
    #
    bond = CustomBondForce("k * (r-r0)^2")
    bond.addPerBondParameter('k')
    bond.addPerBondParameter('r0')
    #
    bond_flat = CustomBondForce("k * (max(abs(r-r0)-rflat, 0.0))^2")
    bond_flat.addPerBondParameter('k')
    bond_flat.addPerBondParameter('r0')
    bond_flat.addPerBondParameter('rflat')
    #
    #pos = CustomExternalForce("k0*dsq ; dsq=((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    pos = CustomExternalForce("k0*d^2 ; d=periodicdistance(x,y,z, x0,y0,z0)")
    pos.addPerParticleParameter("x0")
    pos.addPerParticleParameter("y0")
    pos.addPerParticleParameter("z0")
    pos.addPerParticleParameter('k0')
    #
    #pos_flat = CustomExternalForce("k0*(max(d-flat, 0.0))^2 ; d=sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    pos_flat = CustomExternalForce("k0*(max(d-flat, 0.0))^2 ; d=periodicdistance(x,y,z, x0,y0,z0)")
    pos_flat.addPerParticleParameter('k0')
    pos_flat.addPerParticleParameter("x0")
    pos_flat.addPerParticleParameter("y0")
    pos_flat.addPerParticleParameter("z0")
    pos_flat.addPerParticleParameter('flat')
    #
    for custom in custom_s:
        if custom.r_type == 'bond':
            p = (custom.prm[0]*kilocalories_per_mole/angstroms**2,\
                 custom.prm[1]*angstroms)
            bond.addBond(custom.i_atm[0], custom.i_atm[1], p)
        elif custom.r_type == 'bond_flat':
            p = (custom.prm[0]*kilocalories_per_mole/angstroms**2,\
                 custom.prm[1]*angstroms, \
                 custom.prm[2]*angstroms)
            bond_flat.addBond(custom.i_atm[0], custom.i_atm[1], p)
        elif custom.r_type == 'angle':
            pass
        elif custom.r_type == 'torsion':
            pass
        elif custom.r_type == 'position':
            i_atm = custom.i_atm[0]
            p = crd[i_atm].tolist()
            p.append(custom.prm[0]*mass[i_atm]*kilocalories_per_mole/angstroms**2)
            pos.addParticle(i_atm, p)
        elif custom.r_type == 'position_flat':
            i_atm = custom.i_atm[0]
            p = [custom.prm[0]*mass[i_atm]*kilocalories_per_mole/angstroms**2]
            p.extend(crd[i_atm].tolist())
            p.append(custom.prm[1]*angstroms)
            pos_flat.addParticle(i_atm, p)
    #
    if bond.getNumBonds() > 0:
        rsr.append(bond)
    if bond_flat.getNumBonds() > 0:
        rsr.append(bond_flat)
    if pos.getNumParticles() > 0:
        rsr.append(pos)
    if pos_flat.getNumParticles() > 0:
        rsr.append(pos_flat)
    return rsr

def read_custom_restraint(custom_file):
    custom_restraints = []
    if custom_file is None:
        return None, custom_restraints

    with open('%s'%custom_file) as fp:
        ref_fn = None
        for line in fp:
            if line.startswith("#"):
                continue
            elif line.startswith("REFERENCE"):
                ref_fn = line.strip().split()[1]
                continue
            x = line.strip().split()
            r_type = x[0]
            n_atom = int(x[1])
            n_param = int(x[2])
            i_atm = [int(xi)-1 for xi in x[3:3+n_atom]]
            prm = [float(xi) for xi in x[3+n_atom:]]
            custom = Custom(r_type, n_atom, n_param, i_atm, prm)
            custom_restraints.append(custom)
    return ref_fn, custom_restraints
