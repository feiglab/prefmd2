import os
import sys

import warnings
warnings.filterwarnings("ignore")

import numpy as np

from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *

# from libcommon import *

ALPHA = 1.57
MIN_DENSITY = 1.0e-6
np.set_printoptions(linewidth=1000, edgeitems=10000, suppress=True)

def gaussian(m, s, X0):
    return np.exp(-0.5*((X0[None,:]-m[:,None])/s)**2) / (np.sqrt(2.*np.pi) * s)
    
def get_dist_distr(frame_s, method='CA', d_cut=10.0, d_max=20.0, d_bin=0.5, \
        density_cut=0.15, sequence_separation=3):
    if method == 'CA':
        index = frame_s.top.select("name CA")
    elif method == 'CB':
        index = frame_s.top.select("(resn != GLY and name CA) or (resn == GLY and name CA)")
    R = frame_s.xyz[:,index] ; l_seq = R.shape[1]
    #
    dist_bin = np.arange(0., d_max+0.001*d_bin, d_bin) / 10.    # nm
    cnt = np.where(dist_bin <= (d_cut/10.))
    #
    pair_s = []
    distr_s = []
    for i in range(l_seq):
        for j in range(i+sequence_separation, l_seq):
            dR = R[:,j] - R[:,i]
            d = np.sqrt(np.sum(dR**2, axis=-1))
            #
            distr = gaussian(d, (d_bin/10.), dist_bin) * (d_bin/10.)
            distr = distr.sum(axis=0)
            distr /= np.float32(d.shape[0])
            distr[-1] = 1.0 - distr[:-1].sum()
            distr = np.maximum(MIN_DENSITY, distr)
            #
            prob = distr[cnt].sum()
            if prob > density_cut:
                pair_s.append((i,j))
                distr_s.append(distr)
    return dist_bin, pair_s, distr_s

def construct_restraint_contact(psf, distr, force_const, method='CA'):
    atomIndex = []
    for i,atom in enumerate(psf.topology.atoms()):
        if atom.name == 'CA':
            atomIndex.append(i)
    #
    dist_cntr, pair_s, distr_s = distr
    #
    n_pair = len(pair_s)
    n_dist_bin = len(dist_cntr)
    potential_s = -np.log(distr_s) + ALPHA * np.log(np.maximum(0.2, dist_cntr)/dist_cntr[-1])
    potential_s -= np.min(potential_s, axis=-1)[:,None]
    potential_s = potential_s.T.ravel() * force_const
    #potential_s = potential_s.ravel() * force_const
    #
    bond = CustomCompoundBondForce(2, \
            "potential(pair, d) ; d = max(d_min, min(d_max, distance(p1,p2)))")
    bond.addGlobalParameter("d_min", dist_cntr[0])
    bond.addGlobalParameter("d_max", dist_cntr[-1])
    bond.addPerBondParameter("pair")
    bond.addTabulatedFunction("potential", \
            Continuous2DFunction(n_pair, n_dist_bin, potential_s, 0., float(n_pair-1), dist_cntr[0], dist_cntr[-1]))
    #
    for i,pair in enumerate(pair_s):
        bond.addBond([atomIndex[p] for p in pair], [float(i)])
    return bond
