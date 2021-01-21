import os
import sys
import json
import mdtraj
import numpy as np
from itertools import product
from mdtraj.core.element import hydrogen

import path
import libcommon  # from libcommon import *

PARAM_DIST_SAME  = 0.2 # nm
PARAM_DIST_COORD = 0.3 # nm
PARAM_DIST_BSITE = 0.5 # nm
PARAM_DIST_ALIGN = 0.8 # nm

def read_ligand_json(fn):
    if isinstance(fn, str):
        fn = path.Path(fn)
    #
    cwd = os.getcwd()
    work_home = fn.dirname()
    work_home.chdir()
    #
    with fn.open() as fp:
        info = json.load(fp)
        for key in info:
            info[key] = libcommon.JSONdeserialize(info[key])
    #
    os.chdir(cwd)
    return info

def get_ligand_info(job, wait_after_run, sleep=30):
    ligand_home = job.work_home.subdir("ligand", build=True)
    job.ligand_json = ligand_home.fn("ligand.json")
    if job.ligand_json.status():
        return True
    #
    with open("%s/STDRES" % libcommon.DEFAULT_HOME) as fp:
        STDRES = []
        for line in fp:
            STDRES.append(line.strip())
    with open("%s/METALs" % libcommon.DEFAULT_HOME) as fp:
        METALs = []
        for line in fp:
            METALs.append(line.strip())
    #
    str_fn_s = []
    while True:
        status = True
        #
        pdb_fn = ligand_home.glob("*.pdb")
        if len(pdb_fn) == 0:
            status = False ; break
        ligand_pdb_fn = pdb_fn[0]
        #
        ligand_s = read_ligand_pdb(ligand_pdb_fn, STDRES, METALs)
        for ligName in ligand_s:
            is_metal = ligand_s[ligName]
            if is_metal: continue
            #
            str_fn = ligand_home.fn("%s.str"%ligName)
            if not str_fn.status():
                status = False ; break
            str_fn_s.append(str_fn)
        #
        if wait_after_run:
            sys.stderr.write("waiting for ligand info... \n")
            time.sleep(sleep)
        else:
            break
    #
    if not status: return status
    #
    info = {}
    #
    info['pdb_fn'] = ligand_pdb_fn
    info['ligand_s'] = ligand_s # as dict
    info['str_fn_s'] = str_fn_s
    #
    set_bsite(info, job.top_fn, ligand_pdb_fn)
    #
    cwd = os.getcwd()
    ligand_home.chdir()
    with job.ligand_json.open("wt") as fout:
        fout.write(json.dumps(info, indent=2, default=JSONserialize))
    os.chdir(cwd)

    return status

def read_ligand_pdb(fn, STDRES, METALs):
    ligand_s = {}
    with fn.open() as fp:
        for line in fp:
            if (not line.startswith("ATOM")) and (not line.startswith("HETATM")):
                continue
            resName = line[17:21].strip()
            if resName not in STDRES:
                if resName not in ligand_s:
                    ligand_s[resName] = (resName in METALs)
    return ligand_s

def get_aligned_residues(ref, pdb):
    dr = pdb[None,:] - ref[:,None]
    dmtx = np.sqrt(np.sum(dr**2, -1))
    aligned = np.argmin(dmtx, -1)
    status = (np.min(dmtx, -1) > PARAM_DIST_SAME)
    aligned[status] = -1
    return aligned

def get_bsite_geometry(dist_s, atomIndex, cbIndex, caIndex):
    d_sc = dist_s[atomIndex]
    if cbIndex.shape[0] != 0:
        d_cb = dist_s[cbIndex]
    d_ca = dist_s[caIndex]
    #
    geom_sc = []
    for i,j in zip(*np.where(d_sc < PARAM_DIST_COORD)):
        i_atm = atomIndex[i]
        d = d_sc[i,j]
        geom_sc.append((i_atm, j, d))
    if len(geom_sc) == 0:
        d = np.min(d_sc)
        i,j = np.where(d_sc == d)
        i_atm = atomIndex[i[0]]
        j = j[0]
        geom_sc.append((i_atm, j, d))
    #
    geom_ca = [(caIndex[0], np.argmin(d_ca), np.min(d_ca))]
    if cbIndex.shape[0] != 0:
        geom_cb = [(cbIndex[0], np.argmin(d_cb), np.min(d_cb))]
    else:
        geom_cb = None
    return geom_sc, geom_cb, geom_ca

def set_bsite(info, top_fn, ligand_pdb_fn):
    ligand_pdb = mdtraj.load(ligand_pdb_fn.short())
    ligand_pdb = ligand_pdb.atom_slice(ligand_pdb.top.select("element != H"))
    #
    ligandIndex = ligand_pdb.top.select("resname %s"%(" ".join(info['ligand_s'].keys())))
    proteinIndex = ligand_pdb.top.select("protein")
    proteinIndex = proteinIndex[~np.isin(proteinIndex, ligandIndex)]
    #
    protein = ligand_pdb.atom_slice(proteinIndex)
    ligand = ligand_pdb.atom_slice(ligandIndex)
    ligandResidueNumber_s = np.unique([r.index for r in ligand.top.residues]).tolist()
    #
    atomPair_s = [atomPair for atomPair in product(proteinIndex, ligandIndex)]
    dist_s = mdtraj.compute_distances(ligand_pdb, atomPair_s)[0].reshape((len(proteinIndex), len(ligandIndex)))
    #
    bsiteAtom_s = proteinIndex[np.where(dist_s < PARAM_DIST_BSITE)[0]]
    bsiteResidue_s = np.unique([ligand_pdb.top.atom(i_atm).residue.index for i_atm in bsiteAtom_s])
    bsiteCalphaIndex = ligand_pdb.top.select("name CA")[bsiteResidue_s]
    bsiteCa_s = ligand_pdb.xyz[0,bsiteCalphaIndex]
    bsiteGeom_s = []
    for i_res in bsiteResidue_s:
        bsiteGeom_s.append(get_bsite_geometry(dist_s, \
                ligand_pdb.top.select("resi %d"%i_res), \
                ligand_pdb.top.select("resi %d and name CB"%i_res), \
                ligand_pdb.top.select("resi %d and name CA"%i_res)))
    #
    alignAtom_s = proteinIndex[np.where(dist_s < PARAM_DIST_ALIGN)[0]]
    alignResidue_s = np.unique([ligand_pdb.top.atom(i_atm).residue.index for i_atm in alignAtom_s])
    alignCalphaIndex = ligand_pdb.top.select("name CA")[alignResidue_s]
    alignCa_s = ligand_pdb.xyz[0,alignCalphaIndex]
    #
    top = mdtraj.load(top_fn.short())
    topCalphaIndex = top.top.select("name CA")
    topCa_s = top.xyz[0,topCalphaIndex]
    #
    bsiteAligned = get_aligned_residues(bsiteCa_s, topCa_s)
    bsiteAlignedCalphaIndex = topCalphaIndex[bsiteAligned[bsiteAligned != -1]]
    bsiteAlignedResidue_s = [top.top.atom(i_atm).residue for i_atm in bsiteAlignedCalphaIndex]
    #
    info['bsite_geometry'] = [] ; j = -1
    for i,align in enumerate(bsiteAligned):
        if align == -1: continue
        j += 1
        #
        bsiteResName = ligand_pdb.top.residue(bsiteResidue_s[i]).name
        bsiteAlignedResName = bsiteAlignedResidue_s[j].name
        #
        if bsiteAlignedResName == bsiteResName:  # use geom_sc
            geom = bsiteGeom_s[i][0]
        elif bsiteAlignedResName == 'GLY' or bsiteResName == 'GLY':  # use geom_ca
            geom = bsiteGeom_s[i][2]
        else:   # use either geom_ca or geom_cb depending on distance
            if bsiteGeom_s[i][1][0][2] < bsiteGeom_s[i][2][0][2]:
                geom = bsiteGeom_s[i][1]
            else:
                geom = bsiteGeom_s[i][2]
        #
        for i_prot, i_lig, d in geom:
            protResidueIndex = bsiteAlignedResidue_s[j].index
            protAtmName = protein.top.atom(i_prot).name
            ligAtom = ligand.top.atom(i_lig)
            ligResidueIndex = ligandResidueNumber_s.index(ligAtom.residue.index)
            ligAtmName = ligAtom.name
            info['bsite_geometry'].append((protResidueIndex, protAtmName, ligResidueIndex, ligAtmName, float(d)))
    #
    alignAligned = get_aligned_residues(alignCa_s, topCa_s)
    alignAlignedCalphaIndex = topCalphaIndex[alignAligned[alignAligned != -1]]
    alignAlignedResidue_s = [top.top.atom(i_atm).residue for i_atm in alignAlignedCalphaIndex]
    #
    info['bsite_align'] = [[], []] ; j = -1
    for i,align in enumerate(alignAligned):
        if align == -1: continue
        j += 1
        #
        alignResidueIndex = ligand_pdb.top.residue(alignResidue_s[i]).index
        alignAlignedResidueIndex = alignAlignedResidue_s[j].index
        info['bsite_align'][0].append(alignResidueIndex)
        info['bsite_align'][1].append(alignAlignedResidueIndex)

def add_missing_hydrogen(ligand, str_fn_s):
    def read_str(str_fn):
        atom_s = []
        bond_s = {}
        with str_fn.open() as fp:
            for line in fp:
                if line.startswith("ATOM "):
                    atom_s.append(line.strip().split()[1])
                elif line.startswith("BOND "):
                    x = line.strip().split()
                    if x[1] not in bond_s:
                        bond_s[x[1]] = []
                    bond_s[x[1]].append(x[2])
                    if x[2] not in bond_s:
                        bond_s[x[2]] = []
                    bond_s[x[2]].append(x[1])
        return atom_s, bond_s
    #
    str_s = {}
    for str_fn in str_fn_s:
        ligName = str_fn.name()
        str_s[ligName] = read_str(str_fn)
    #
    xyz = []
    for residue in ligand.top.residues:
        resName = residue.name
        if resName not in str_s:
            for atom in residue.atoms:
                xyz.append(ligand.xyz[0,atom.index])
            continue
        #
        atom_s = str_s[resName][0]
        status = [False for _ in atom_s]
        for atom in residue.atoms:
            if atom.name in atom_s:
                status[atom_s.index(atom.name)] = True
                xyz.append(ligand.xyz[0,atom.index])
        #
        for atmName,placed in zip(atom_s, status):
            if placed: continue
            #
            if not atmName.startswith("H"):
                sys.exit("ERROR: missing heavy atom (%s) in the input ligand (%s).\n"%(atmName, resName))
            bonded_atmName = str_s[resName][1][atmName][0]
            bonded_atom = [atom for atom in residue.atoms_by_name(bonded_atmName)][0]
            r0 = ligand.xyz[0,bonded_atom.index]
            dr = 2.0 * np.random.random(3) - 1.0
            dr /= np.sqrt(np.sum(dr**2))    # unit-vec
            r = r0 + dr * 0.1
            #
            xyz.append(r)
            atom = ligand.top.add_atom(atmName, hydrogen, residue)
        #
        for atmName in str_s[resName][1]:
            curr_atom = [atom for atom in residue.atoms_by_name(atmName)][0]
            for bonded_atmName in str_s[resName][1][atmName]:
                bonded_atom = [atom for atom in residue.atoms_by_name(bonded_atmName)][0]
                if curr_atom.index < bonded_atom.index:
                    ligand.top.add_bond(curr_atom, bonded_atom)
    #
    xyz = np.array(xyz)
    new = mdtraj.Trajectory(xyz[None,:], ligand.top)
    return new

def add_ligand(info, pdb_fn, out_fn):
    ligand_pdb = mdtraj.load(info['pdb_fn'].short())
    calphaIndex = ligand_pdb.top.select("name CA")
    ligand_align = calphaIndex[info['bsite_align'][0]]
    #
    target_pdb = mdtraj.load(pdb_fn.short())
    calphaIndex = target_pdb.top.select("name CA")
    target_align = calphaIndex[info['bsite_align'][1]]
    #
    ligand_pdb.superpose(target_pdb, atom_indices=ligand_align, ref_atom_indices=target_align)
    ligandIndex = ligand_pdb.top.select("resname %s"%(" ".join(info['ligand_s'].keys())))
    ligand = ligand_pdb.atom_slice(ligandIndex)
    #
    ligand = add_missing_hydrogen(ligand, info['str_fn_s'])
    #
    for i,residue in enumerate(ligand.top.residues):
        residue.resSeq = i+1
        residue.segment_id = 'LIG'
    #
    segNo = 0 ; resNo_prev = None
    for i,chain in enumerate(target_pdb.top.chains):
        for residue in chain.residues:
            residue.segment_id='P%03d'%segNo
            #
            chain_break = False
            for atom in residue.atoms:
                if atom.name in ['OXT', 'OT2', 'OT1']:
                    chain_break = True
            if chain_break:
                segNo += 1
    #
    top = target_pdb.top.join(ligand.top)
    xyz = np.concatenate([target_pdb.xyz[0], ligand.xyz[0]])[None,:]
    #
    pdb = mdtraj.Trajectory(xyz, top)
    pdb.save(out_fn.short())

def update_ligand_name(pdb_fn, ligand_s):
    het_s = {het[:3].strip(): het.strip() for het in ligand_s if len(het.strip()) > 3}
    if len(het_s) == 0:
        return
    #
    wrt = []
    with pdb_fn.open() as fp:
        for line in fp:
            if (not line.startswith("ATOM")) and (not line.startswith("HETATM")):
                wrt.append(line)
                continue
            #
            resName3 = line[17:21].strip()
            if resName3 not in het_s:
                wrt.append(line)
                continue
            #
            resName = het_s[resName3]
            line = '%s%4s%s'%(line[:17], resName, line[21:])
            wrt.append(line)
    with pdb_fn.open('wt') as fout:
        fout.writelines(wrt)

def get_ligand_restratint(pdb, psf_fn, info):
    def read_psf(psf_fn):
        psf = {}
        with psf_fn.open() as fp:
            read = False
            for line in fp:
                if '!NATOM' in line:
                    read = True ; continue
                elif line.strip() == '':
                    read = False
                if not read: continue
                #
                x = line.strip().split()
                i_atm = int(x[0])
                segName = x[1]
                resNo = int(x[2])
                resName = x[3]
                atmName = x[4]
                #
                if segName not in psf:
                    psf[segName] = {}
                psf[segName][(resNo, atmName)] = i_atm
        return psf
    #
    if isinstance(psf_fn, str):
        psf_fn = path.Path(psf_fn)
    #
    ligandAtom = np.array([atom.residue.segment_id == 'LIG' for atom in pdb.top.atoms], dtype=bool)
    ligand = pdb.atom_slice(np.where(ligandAtom)[0])
    #
    calphaIndex = pdb.top.select("name CA")
    restraint_pair = []
    for geom in info['bsite_geometry']:
        resIndex, atmName, ligIndex, ligAtmName, d = geom
        #
        residue = pdb.top.atom(calphaIndex[resIndex]).residue
        lig = ligand.top.residue(ligIndex)
        #
        pair = [(residue.segment_id, residue.resSeq, atmName),\
                (lig.segment_id, lig.resSeq, ligAtmName), d]
        restraint_pair.append(pair)
    #
    psf = read_psf(psf_fn)
    #
    FMT = 'bond 2 2  %6d %6d     %8.5f %6.3f\n'
    #
    restraint_s = []
    for prot, lig, d in restraint_pair:
        i_atm = psf[prot[0]][prot[1:]] - 1
        i_lig = psf[lig[0]][lig[1:]] - 1
        #
        force_const = 10.0 / (1.0 + max(0.0, (10*(d-PARAM_DIST_COORD))**2))
        restraint = (i_atm, i_lig,  force_const, d)
        restraint_s.append(restraint)
    return restraint_s
