"""
Module for initializing a refinement job and processing the file of the
initial protein 3D model.
Makes use of the convpdb.pl script from MMSTSB.
"""

import os
import numpy as np

import path
import libcommon


def prep(arg):

    #
    assert arg.input_pdb is not None
    #
    print("- Preparing the working directory...")
    #
    work_home = path.Dir("%s/%s"%(arg.work_dir, arg.title))
    json_job = work_home.fn("job.json")
    #
    def is_continuous_domain(pdb_fn):
        Rs = [] ; R = [] ; Rs.append(R)
        with pdb_fn.open() as fp:
            for line in fp:
                if not line.startswith("ATOM"):
                    if line.startswith("TER"):
                        R = [] ; Rs.append(R)
                    continue
                atmName = line[12:16].strip()
                if atmName == 'CA':
                    R.append([line[30:38], line[38:46], line[46:54]])
        for R in Rs:
            R = np.array(R, dtype=float)
            b = np.sqrt(np.sum((R[1:] - R[:-1])**2, -1))
            b_max = np.max(b)
            if b_max > 5.0:
                return False
        return True
    #
    job = libcommon.Job(arg.work_dir, arg.title, build=True)
    job.run_type = 'refine'
    if arg.use_hybrid and is_continuous_domain(arg.input_pdb):
        job.use_hybrid = True
    else:
        job.use_hybrid = False
    job.use_extensive = bool(arg.use_extensive)
    # if arg.is_membrane_protein:
    #     job.is_membrane_protein = True
    job.is_membrane_protein = False
    # if arg.has_ligand:
    #     job.has_ligand = True
    job.has_ligand = False
    # if arg.is_oligomer:
    #     job.is_oligomer = True
    #     if job.use_hybrid:
    #         job.use_hybrid = False
    job.is_oligomer = False
    job.n_cpus = arg.n_cpus
    job.gpu_ids = arg.gpu_ids
    #
    job.init_home = job.work_home.subdir("init", build=True)
    # job.verbose = arg.verbose
    # job.keep_tmp = arg.keep
    #
    out = job.init_home.fn("init.pdb")
    cmd = ['convpdb.pl', '-out', 'generic', arg.input_pdb.short()]
    output = libcommon.system(cmd, stdout=True)
    with out.open("wt") as fout:
        fout.write(output)

    has_hydrogens, pdb_lines_renum = delete_hydrogens(out.path())
    if has_hydrogens:
        # Structures with hydrogens will cause problems in later stages.
        print("- Removing hydrogens from the initial model...")
        with out.open("wt") as fout:
            fout.writelines(pdb_lines_renum)

    job.init_pdb = [out]
    job.to_json()
    return job


def delete_hydrogens(in_filepath, out_filepath=None):

    atom_lines = []
    found_h = False
    with open(in_filepath, "r") as i_fh:
        for line in i_fh:
            if line.startswith(("ATOM", "HETATM")):
                atm_type = line[11:16].replace(" ", "")
                if atm_type[0] == "H":
                    found_h = True
                    continue
                resname = line[17:20]
                if resname in ("HSP", "HSD"):
                    _line = line[:17] + "HIS" + line[20:]
                else:
                    _line = line
                atom_lines.append(_line)

    if found_h:
        atom_lines_renum = []
        for i, line in enumerate(atom_lines):
            # print("---")
            # print(line)
            atm_serial = int(line[5:11].replace(" ", ""))
            line_renum = line[:5] + str(i+1).rjust(6) + line[11:]
            # print(line_renum)
            atom_lines_renum.append(line_renum)
    else:
        atom_lines_renum = atom_lines

    if out_filepath is not None:
        with open(out_filepath, "w") as o_fh:
            o_fh.writelines(atom_lines_renum)

    return found_h, atom_lines_renum
