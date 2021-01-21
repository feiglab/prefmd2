#!/usr/bin/env python

import os
import sys
import time
import json
import argparse

import mdtraj
import numpy as np

import path
import libcommon


EXEC = '%s/exec_prod_runner.py'%libcommon.EXEC_HOME
SPEED = '%s/exec_check_md_speed.py'%libcommon.EXEC_HOME

SPEED_SCALE = 0.85  # slow-down because of restraints

def build_restraint_Cartesian(ref):
    calphaIndex = ref.top.select("name CA")
    return calphaIndex + 1

def build_restraint_distance(ref, distance_cutoff=10.0, sequence_separation=4):
    calphaIndex = ref.top.select("name CA")
    chainIndex = np.array([ref.top.atom(i).residue.chain.index for i in calphaIndex])
    xyz = ref.xyz[0][calphaIndex]
    dr = xyz[:,None] - xyz[None,:]
    dist = np.sqrt(np.sum(dr**2, axis=-1)) * 10.0
    #
    rsr_s = []
    for i,ca_i in enumerate(calphaIndex):
        chain_i = chainIndex[i]
        for j,ca_j in enumerate(calphaIndex):
            if j <= i: continue
            chain_j = chainIndex[j]
            if chain_i == chain_j and (j-i < sequence_separation):
                continue
            d = dist[i,j]
            if d < distance_cutoff:
                rsr_s.append((ca_i+1, ca_j+1, d))
    return rsr_s

def build_restraint(output_prefix, options, verbose):
    ref_fn = path.Path(options['restraint']['reference'])
    ref = mdtraj.load(ref_fn.short())
    #
    rsr_fn_s = []
    #
    if options['restraint']['mode'] == 'Cartesian':
        wrt = []
        wrt.append("REFERENCE   %s\n"%ref_fn.path())
        #
        rsr_C = build_restraint_Cartesian(ref)
        par = options['restraint']['Cartesian']
        n_param = len(par)
        if n_param == 1:
            format = 'position 1 1  %5d  %8.5f\n'
            for ca in rsr_C:
                wrt.append(format%(ca, par[0]))
        else:
            format = 'position_flat 1 2  %5d  %8.5f %8.5f\n'
            for ca in rsr_C:
                wrt.append(format%(ca, par[0], par[1]))
        with open("%s.rsr"%output_prefix, 'wt') as fout:
            fout.writelines(wrt)
        for _ in range(options['md']['iter']):
            rsr_fn_s.append("%s.rsr"%output_prefix)

    elif options['restraint']['mode'] == 'distance':
        wrt = []
        wrt.append("REFERENCE   %s\n"%ref_fn.path())
        #
        rsr_d = build_restraint_distance(ref)
        par = options['restraint']['distance']
        n_param = len(par)
        if n_param == 1:
            format = 'bond 2 2  %5d %5d  %8.5f %8.5f\n'
            for ca_i, ca_j, d in rsr_d:
                wrt.append(format%(ca_i, ca_j, par[0], d))
        else:
            format = 'bond_flat 2 3  %5d %5d  %8.5f %8.5f %8.5f\n'
            for ca_i, ca_j, d in rsr_d:
                wrt.append(format%(ca_i, ca_j, par[0], d, par[1]))
        with open("%s.rsr"%output_prefix, 'wt') as fout:
            fout.writelines(wrt)
        for _ in range(options['md']['iter']):
            rsr_fn_s.append("%s.rsr"%output_prefix)

    elif options['restraint']['mode'] == 'dual':
        rsr_C = build_restraint_Cartesian(ref)
        rsr_d = build_restraint_distance(ref)
        #
        par_C = options['restraint']['Cartesian']
        n_param_C = len(par_C)
        if n_param_C == 1:
            format_C = 'position 1 1  %5d  %8.5f\n'
        else:
            format_C = 'position_flat 1 2  %5d  %8.5f %8.5f\n'
        par_d = options['restraint']['distance']
        n_param_d = len(par_d)
        if n_param_d == 1:
            format_d = 'bond 2 2  %5d %5d  %8.5f %8.5f\n'
        else:
            format_d = 'bond_flat 2 3  %5d %5d  %8.5f %8.5f %8.5f\n'
        #
        w_d = 0.5 ; w_C = 1.0 - w_d
        #
        wrt = []
        wrt.append("REFERENCE   %s\n"%ref_fn.path())
        #
        if w_C > 0.0:
            if n_param_C == 1:
                for ca in rsr_C:
                    wrt.append(format_C%(ca, par_C[0]*w_C))
            else:
                for ca in rsr_C:
                    wrt.append(format_C%(ca, par_C[0]*w_C, par_C[1]))
        if w_d > 0.0:
            if n_param_d == 1:
                for ca_i, ca_j, d in rsr_d:
                    wrt.append(format_d%(ca_i, ca_j, par_d[0]*w_d, d))
            else:
                for ca_i, ca_j, d in rsr_d:
                    wrt.append(format_d%(ca_i, ca_j, par_d[0]*w_d, d, par_d[1]))
        with open("%s.rsr"%(output_prefix), 'wt') as fout:
            fout.writelines(wrt)
        for _ in range(options['md']['iter']):
            rsr_fn_s.append("%s.rsr"%output_prefix)

    elif options['restraint']['mode'] == 'dual.anneal':
        rsr_C = build_restraint_Cartesian(ref)
        rsr_d = build_restraint_distance(ref)
        #
        par_C = options['restraint']['Cartesian']
        n_param_C = len(par_C)
        if n_param_C == 1:
            format_C = 'position 1 1  %5d  %8.5f\n'
        else:
            format_C = 'position_flat 1 2  %5d  %8.5f %8.5f\n'
        par_d = options['restraint']['distance']
        n_param_d = len(par_d)
        if n_param_d == 1:
            format_d = 'bond 2 2  %5d %5d  %8.5f %8.5f\n'
        else:
            format_d = 'bond_flat 2 3  %5d %5d  %8.5f %8.5f %8.5f\n'
        #
        for k in range(options['md']['iter']):
            w_d = float(k)/float(options['md']['iter']-1)
            w_C = 1.0 - w_d
            #
            wrt = []
            wrt.append("REFERENCE   %s\n"%ref_fn.path())
            #
            if w_C > 0.0:
                if n_param_C == 1:
                    for ca in rsr_C:
                        wrt.append(format_C%(ca, par_C[0]*w_C))
                else:
                    for ca in rsr_C:
                        wrt.append(format_C%(ca, par_C[0]*w_C, par_C[1]))
            if w_d > 0.0:
                if n_param_d == 1:
                    for ca_i, ca_j, d in rsr_d:
                        wrt.append(format_d%(ca_i, ca_j, par_d[0]*w_d, d))
                else:
                    for ca_i, ca_j, d in rsr_d:
                        wrt.append(format_d%(ca_i, ca_j, par_d[0]*w_d, d, par_d[1]))
            with open("%s.%d.rsr"%(output_prefix, k), 'wt') as fout:
                fout.writelines(wrt)
            rsr_fn_s.append("%s.%d.rsr"%(output_prefix, k))
    return rsr_fn_s


def run(output_prefix, input_json, options, verbose, requeue=False, generate_solute_file=True):
    run_home = path.Dir(".")
    #
    DONE = run_home.fn("solute.dcd")
    if DONE.status():
        return
    #
    if 'restraint' in options:
        rsr_s = build_restraint(output_prefix, options, verbose)
        assert len(rsr_s) == options['md']['iter']
    #
    init_pdb_fn = run_home.fn("%s.init.pdb"%output_prefix)
    if not init_pdb_fn.status():
        os.symlink(options['input']['pdb'], init_pdb_fn.short())
    init_pdb = mdtraj.load(init_pdb_fn.short())
    boxsize = init_pdb.unitcell_lengths[0]
    boxsize = ' '.join(['%8.3f'%(width*10.0) for width in boxsize])
    #
    k_iter = 0
    n_error = 0
    out_dcd_fn_s = []
    if 'restart' in options:
        restart_fn = path.Path(options['input']['restart'])
    else:
        restart_fn = None

    n_frame_per_iter = int(options['md']['dynsteps'] / options['md']['dynoutfrq'])
    if 'n_atom' in options['input']:
        dcdOut_topology = init_pdb.atom_slice(np.arange(options['input']['n_atom']))
    else:
        dcdOut_topology = init_pdb

    while (k_iter < options['md']['iter']):
        out_dcd_fn = run_home.fn("%s.%d.dcd"%(output_prefix, k_iter))
        out_chk_fn = run_home.fn("%s.%d.restart.pkl"%(output_prefix, k_iter))
        if out_dcd_fn.status(size=1000) and out_chk_fn.status(size=1000):
            k_iter += 1
            out_dcd_fn_s.append(out_dcd_fn)
            continue
        elif out_dcd_fn.status(size=1000) and (not out_chk_fn.status(size=1000)):
            out_dcd_fn.remove()
        #
        if k_iter > 0:
            restart_fn = run_home.fn("%s.%d.restart.pkl"%(output_prefix, k_iter-1))
        #
        cmd = []
        cmd.extend(['--input', input_json.short()])
        cmd.extend(['--pdb', init_pdb_fn.short()])
        cmd.extend(['--dcd', out_dcd_fn.short()])
        cmd.extend(['--chk', out_chk_fn.short()])
        if 'restraint' in options:
            cmd.extend(['--custom', rsr_s[k_iter]])
        if restart_fn is not None:
            cmd.extend(['--restart', restart_fn.short()])
        cmd.extend(['--log', run_home.fn("%s.%d.log"%(output_prefix, k_iter)).short()])
        if 'use_hmr' in options['ff'] and options['ff']['use_hmr']:
            cmd.append("--hmr")
        if 'use_barostat' in options['md']:
            cmd.append("--barostat")
            cmd.append(options['md']['use_barostat'])
        #
        t_begin = time.time()
        with open("%s.err"%output_prefix, 'wt') as ferr:
            libcommon.asystem(module="exec_prod_runner", args=cmd,
                              errfile=ferr)
        dt_per_step = time.time() - t_begin
        #
        try:
            n_frame = mdtraj.load(out_dcd_fn.short(), top=dcdOut_topology, atom_indices=[0]).n_frames
        except:
            n_frame = 0
        n_frame_status = (n_frame == n_frame_per_iter)
        #
        if out_chk_fn.status(size=1000) and n_frame_status:
            k_iter += 1
            out_dcd_fn_s.append(out_dcd_fn)
            restart_fn = out_chk_fn
        else:
            n_error += 1
            libcommon.system("mv %s.err %s.err.%d"%(output_prefix, output_prefix, n_error))
            if n_error >= libcommon.MAX_ERROR:
                with run_home.fn("ERROR").open("wt") as fout:
                    fout.write("#")
                sys.exit("ERROR: failed to run %s\n"%run_home)
    #
    if generate_solute_file:
        libcommon.system("mdconv -out %s -atoms 1:%d -unwrap -box %s %s"%\
                (DONE.short(), options['input']['n_atom'], boxsize, \
                ' '.join([dcd_fn.short() for dcd_fn in out_dcd_fn_s])))

def check_speed(output_prefix, input_json, options, verbose):
    if 'time_limit' in options['md'] and options['md']['time_limit'] > 0.0:
        cmd = [output_prefix]
        cmd.extend(["--input", input_json.short()])
        libcommon.asystem(module="exec_check_md_speed", args=cmd)
        #
        with open("SPEED") as fp:
            speed = float(fp.read().strip())    # steps/s
            speed *= SPEED_SCALE
        #
        time_limit = options['md']['time_limit']
        max_steps = time_limit * speed / options['md']['iter']
        #
        if max_steps < options['md']['dynsteps']:
            dynoutfrq = options['md']['dynoutfrq']
            dynsteps = int(max_steps/dynoutfrq) *dynoutfrq
            #
            options['md']['dynsteps'] = dynsteps
            with input_json.open("wt") as fout:
                fout.write(json.dumps(options, indent=2))
    return options

def main(args=None):
    arg = argparse.ArgumentParser(prog='prod')
    arg.add_argument(dest='output_prefix')
    arg.add_argument('--input', dest='input_json', required=True)
    arg.add_argument('--toppar', dest='toppar', nargs='*', default=None)
    arg.add_argument('--custom', dest='custom_file', default=None)
    arg.add_argument('-v', '--verbose', default=False, action='store_true')
    arg.add_argument('--keep', dest='keep', action='store_true', default=False,\
            help='set temporary file mode (default=False)')
    arg.add_argument('--requeue', dest='requeue', action='store_true', default=False)
    arg.add_argument('--no-solute', dest='generate_solute_file', action='store_false', default=True)
    #
    if len(sys.argv) == 1:
        arg.print_help()
        return
    arg = arg.parse_args(args)
    #
    arg.input_json = path.Path(arg.input_json)
    #
    with arg.input_json.open() as fp:
        options = json.load(fp)
    if arg.toppar is not None:
        options['ff']['toppar']= arg.toppar
    if arg.custom_file is not None:
        options['ff']['custom'] = arg.custom_file

    options = check_speed(arg.output_prefix, arg.input_json, options, arg.verbose)
    run(arg.output_prefix, arg.input_json, options, arg.verbose, requeue=arg.requeue, \
        generate_solute_file=arg.generate_solute_file)

if __name__=='__main__':
    main()
