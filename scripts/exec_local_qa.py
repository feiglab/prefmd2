#!/usr/bin/env python

import os
import sys
import json
import argparse

import warnings
warnings.filterwarnings("ignore")

import numpy as np
import scipy.spatial.distance

import mdtraj
from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *

import path
import libcommon
import libgpu
from libcustom import construct_custom_restraint, read_custom_restraint
from libmd import solvate_pdb, generate_PSF, construct_restraint

SPEED = '%s/exec_check_md_speed.py'%libcommon.EXEC_HOME

def equil_md(output_prefix, solv_fn, psf_fn, crd_fn, options, verbose):
    psf = CharmmPsfFile(psf_fn.short())
    crd = CharmmCrdFile(crd_fn.short())
    #
    pdb = mdtraj.load(solv_fn.short())
    #
    box = np.array(crd.positions.value_in_unit(nanometers), dtype=float)
    boxsize = np.max(box, 0) - np.min(box, 0)
    psf.setBox(*boxsize)
    #
    ff = CharmmParameterSet(*options['ff']['toppar'])
    platform = Platform.getPlatformByName(options['openmm']['platform'])
    properties = {}
    properties = libgpu.update_gpu(properties, options)
    #
    sys = psf.createSystem(ff, \
                           nonbondedMethod=PME, \
                           switchDistance=0.8*nanometers, \
                           nonbondedCutoff=1.0*nanometers, \
                           constraints=HBonds)
    #
    sys.addForce(construct_restraint(psf, pdb, 0.5))
    #
    if 'custom' in options['ff'] and options['ff']['custom'] is not None:
        custom_restrains = read_custom_restraint(options['ff']['custom'])
        custom_s = construct_custom_restraint(pdb, custom_restraints[1])
        for custom in custom_s:
            sys.addForce(custom)
    #
    steps_left = options['md']['equil'][0]
    steps_heat = options['md']['heat'][0]
    temp = options['md']['heat'][1]
    temp_incr = options['md']['heat'][2]
    #
    i = 0
    while (temp < options['md']['dyntemp']):
        integrator = LangevinIntegrator(temp*kelvin, \
                                        options['md']['langfbeta']/picosecond, \
                                        options['md']['dyntstep']*picosecond)
        #
        simulation = Simulation(psf.topology, sys, integrator, platform, properties)
        simulation.context.setPositions(crd.positions)
        if i == 0:
            simulation.minimizeEnergy(maxIterations=500)
            simulation.context.setVelocitiesToTemperature(temp*kelvin)
        else:
            with open(chk_fn, 'rb') as fp:
                simulation.context.loadCheckpoint(fp.read())
        simulation.reporters.append(StateDataReporter('%s.heat.%03d.log'%(output_prefix, temp), 500, step=True, \
            time=True, kineticEnergy=True, potentialEnergy=True, temperature=True, progress=True, \
            remainingTime=True, speed=True, volume=True, totalSteps=steps_heat, separator='\t'))
        #
        simulation.step(steps_heat)
        #
        chk_fn = '%s.heat.restart'%(output_prefix)
        with open(chk_fn, 'wb') as fout:
            fout.write(simulation.context.createCheckpoint())
        simulation = None   # have to free CUDA-related variables
        temp += temp_incr ; i += 1
        steps_left -= steps_heat
    #
    sys.addForce(MonteCarloBarostat(1.0*bar, options['md']['dyntemp']*kelvin))
    integrator = LangevinIntegrator(options['md']['dyntemp']*kelvin, \
                                    options['md']['langfbeta']/picosecond, \
                                    options['md']['dyntstep']*picosecond)
    #
    simulation = Simulation(psf.topology, sys, integrator, platform, properties)
    simulation.context.setPositions(crd.positions)
    with open(chk_fn, 'rb') as fp:
        simulation.context.loadCheckpoint(fp.read())

    simulation.reporters.append(StateDataReporter('%s.equil.log'%output_prefix, 2500, step=True, \
        time=True, kineticEnergy=True, potentialEnergy=True, temperature=True, progress=True, \
        remainingTime=True, speed=True, volume=True, totalSteps=steps_left, separator='\t'))
        #
    simulation.step(steps_left)
    #
    chk_fn = '%s.equil.restart'%(output_prefix)

    state = simulation.context.getState(getPositions=True,\
                                        getVelocities=True,\
                                        getForces=True,\
                                        getEnergy=True, \
                                        enforcePeriodicBox=True)

    return state, crd

def prod_md(output_prefix, solv_fn, psf_fn, crd, restart_state, options, verbose):
    psf = CharmmPsfFile(psf_fn.short())
    box = restart_state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(nanometer)
    boxsize = np.diag(box)
    psf.setBox(*boxsize)
    #
    pdb = mdtraj.load(solv_fn.short())
    #
    ff = CharmmParameterSet(*options['ff']['toppar'])
    platform = Platform.getPlatformByName(options['openmm']['platform'])
    properties = {}
    properties = libgpu.update_gpu(properties, options)
    #
    sys = psf.createSystem(ff, \
                           nonbondedMethod=PME, \
                           switchDistance=0.8*nanometers, \
                           nonbondedCutoff=1.0*nanometers, \
                           constraints=HBonds)
    #
    if options['md']['force_const'] > 0.0:
        sys.addForce(construct_restraint(psf, pdb, options['md']['force_const']))
    #
    if 'custom' in options['ff'] and options['ff']['custom'] is not None:
        custom_restrains = read_custom_restraint(options['ff']['custom'])
        custom_s = construct_custom_restraint(pdb, custom_restraints[1])
        for custom in custom_s:
            sys.addForce(custom)
    #
    n_prod = options['md']['prod'][0]
    n_step = options['md']['prod'][1]
    dynoutfrq = options['md']['prod'][2]
    #
    proteinIndex = pdb.top.select("protein")
    is_protein = (proteinIndex.shape[0] > 0)
    #
    dcd_fn_s = []
    for i in range(n_prod):
        dcd_fn = path.Path('%s.prod.r%d.dcd'%(output_prefix, i))
        if dcd_fn.status():
            dcd_fn_s.append(dcd_fn)
            continue
        #
        integrator = LangevinIntegrator(options['md']['dyntemp']*kelvin, \
                                        options['md']['langfbeta']/picosecond, \
                                        options['md']['dyntstep']*picosecond)
        #
        simulation = Simulation(psf.topology, sys, integrator, platform, properties)
        #simulation.context.setPositions(pdb.openmm_positions(0))
        simulation.context.setPositions(crd.positions)
        simulation.context.setState(restart_state)
        simulation.context.setTime(0.0)
        #
        if is_protein:
            simulation.reporters.append(mdtraj.reporters.DCDReporter(dcd_fn.short(), dynoutfrq, atomSubset=proteinIndex))
        else:
            simulation.reporters.append(DCDReporter(dcd_fn.short(), dynoutfrq))

        try:
            simulation.step(n_step)
        except:
            pass
        #
        simulation = None
        dcd_fn_s.append(dcd_fn)
    return dcd_fn_s

def get_error_estimation(output_prefix, pdb_fn, dcd_fn_s):
    pdb = mdtraj.load(pdb_fn.short())
    #
    proteinIndex = pdb.top.select("protein")
    is_protein = (proteinIndex.shape[0] > 0)
    if is_protein:
        pdb = pdb.atom_slice(proteinIndex)
    calphaIndex = pdb.top.select("name CA")
    #
    rmsf_s = []
    for dcd_fn in dcd_fn_s:
        traj = mdtraj.load(dcd_fn.short(), top=pdb)
        if len(traj) < 2:
            continue
        traj.superpose(pdb, atom_indices=calphaIndex)
        #
        xyz = traj.xyz[:,calphaIndex]
        #
        rmsf = []
        for i in range(xyz.shape[1]):
            pd = scipy.spatial.distance.pdist(xyz[:,i,:])
            rmsf.append(10.0*np.sqrt(np.mean(pd**2)))
        rmsf_s.append(rmsf)
    rmsf_s = np.array(rmsf_s)
    qa = np.mean(rmsf_s, axis=0)
    #
    bfactor = []
    for i,residue in enumerate(pdb.top.residues):
        bfactor.append(np.ones(residue.n_atoms)*qa[i])
    bfactor = np.concatenate(bfactor)
    #
    out_fn = path.Path("%s.qa.pdb"%output_prefix)
    pdb.save(out_fn.short(), bfactors=bfactor)
    return out_fn

def run(input_pdb, output_prefix, options, verbose, nonstd):
    tempfile_s = []
    #
    pdb = mdtraj.load(input_pdb.short())
    orient_fn, solv_fn = solvate_pdb(output_prefix, pdb, options, verbose)
    tempfile_s.extend([orient_fn, solv_fn])
    #
    psf_fn, crd_fn = generate_PSF(output_prefix, solv_fn, options, verbose)
    tempfile_s.append(crd_fn)
    #
    restart_state, crd = equil_md(output_prefix, solv_fn, psf_fn, crd_fn, options, verbose)
    #
    dcd_fn_s = prod_md(output_prefix, solv_fn, psf_fn, crd, restart_state, options, verbose)
    #
    qa = get_error_estimation(output_prefix, solv_fn, dcd_fn_s)

def check_speed(output_prefix, input_json, options, verbose, nonstd):
    if 'time_limit' not in options['md'] or options['md']['time_limit'] < 0.0:
        return options
    #
    cmd = [output_prefix]
    cmd.extend(["--input", input_json])
    libcommon.asystem(module="exec_check_md_speed", args=cmd)
    #
    with open("SPEED") as fp:
        speed = float(fp.read().strip())    # steps/s
    #
    time_limit = options['md']['time_limit']
    max_steps = time_limit * speed
    #
    n_heat_steps = (1+ int((options['md']['dyntemp']-options['md']['heat'][1])/options['md']['heat'][2]))
    heat_steps = options['md']['heat'][0] * n_heat_steps
    equil_steps = options['md']['equil'][0]
    prod_steps = options['md']['prod'][0] * options['md']['prod'][1]
    needed_steps = equil_steps + prod_steps

    if needed_steps < max_steps:
        return options
    #
    dynoutfrq = options['md']['prod'][2]
    if equil_steps < 0.4*max_steps:
        dynsteps = int((max_steps-equil_steps)/dynoutfrq/options['md']['prod'][0]) * dynoutfrq
        options['md']['prod'][1] = dynsteps
    else:
        equil_steps = int(0.4*max_steps/dynoutfrq) * dynoutfrq
        options['md']['equil'][0] = equil_steps
        if heat_steps > equil_steps*0.8:
            options['md']['heat'][0] = int(equil_steps*0.8/n_heat_steps/100)*100
        #
        dynsteps = int((max_steps-equil_steps)/dynoutfrq/options['md']['prod'][0]) * dynoutfrq
        options['md']['prod'][1] = dynsteps
    #
    with open(input_json, 'wt') as fout:
        fout.write(json.dumps(options, indent=2))
    return options

def main(args=None):
    arg = argparse.ArgumentParser(prog='local_qa')
    arg.add_argument(dest='output_prefix')
    arg.add_argument('--input', dest='input_json', required=True)
    arg.add_argument('--non_standard', dest='non_standard', default=False, action='store_true')
    arg.add_argument('-v', '--verbose', default=False, action='store_true')
    arg.add_argument('--keep', dest='keep', action='store_true', default=False,\
            help='set temporary file mode (default=False)')
    #
    if len(sys.argv) == 1:
        arg.print_help()
        return
    arg = arg.parse_args(args)
    #
    with open(arg.input_json) as fp:
        options = json.load(fp)
    #
    input_pdb = path.Path(options['input_pdb'])
    #
    options = check_speed(arg.output_prefix, arg.input_json, options, arg.verbose, arg.non_standard)
    #
    run(input_pdb, arg.output_prefix, options, arg.verbose, arg.non_standard)

if __name__=='__main__':
    main()
