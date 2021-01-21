"""
Functions to set up the topology files for the refinement job.
"""


def prep(job, input_pdb):

    print("- Defining topology...")

    job.top_fn = job.init_home.fn("solute.pdb")
    n_atom = 0
    if (not job.top_fn.status()) or (not job.has("n_atom") or not job.has("ssbond")):
        pdb = [] ; ssbond = []
        with input_pdb.open() as fp:
            for line in fp:
                if line.startswith("ATOM"):
                    n_atom += 1
                    pdb.append(line)
                elif line.startswith("SSBOND"):
                    pdb.append(line)
                    ssbond.append(line.rstrip())
        with job.top_fn.open("wt") as fout:
            fout.writelines(pdb)
            fout.write("TER\n")
            fout.write("END\n")
    else:
        return

    job.n_atom = n_atom
    job.ssbond = ssbond
    job.topology_init = True
    #
    job.to_json()


# from libligand import get_ligand_info
#
# def get_membrane_topology(job, n_init, wait_after_run, sleep=30):
#     membrane_home = job.work_home.subdir("membrane", build=True)
#     #
#     while True:
#         job.membrane_pdb = []
#         job.membrane_psf = []
#         job.membrane_crd = []
#         #
#         status = True
#         for i in range(n_init):
#             m_home = membrane_home.subdir("%d"%i, build=True)
#             pdb_fn = m_home.glob("*.pdb")
#             psf_fn = m_home.glob("*.psf")
#             crd_fn = m_home.glob("*.crd")
#             if len(pdb_fn) == 0:
#                 status = False ; break
#             if len(psf_fn) == 0:
#                 status = False ; break
#             if len(crd_fn) == 0:
#                 status = False ; break
#             #
#             job.membrane_pdb.append(pdb_fn[0])
#             job.membrane_psf.append(psf_fn[0])
#             job.membrane_crd.append(crd_fn[0])
#         #
#         if wait_after_run and (not status):
#             sys.stderr.write("waiting for CHARMM-GUI membrane topology... \n")
#             time.sleep(sleep)
#         else:
#             break
#         #
#     if status: job.to_json()
#     return status
#
# def get_oligomer_topology(job, wait_after_run, sleep=30):
#     oligomer_home = job.work_home.subdir("oligomer", build=True)
#     #
#     while True:
#         oligomer_pdb = oligomer_home.glob("*.pdb")
#         #
#         status = (len(oligomer_pdb) != 0)
#         #
#         if wait_after_run and (not status):
#             sys.stderr.write("waiting for Oligomer topology... \n")
#             time.sleep(sleep)
#         else:
#             break
#         #
#     if status:
#         job.oligomer_pdb = oligomer_pdb
#         job.to_json()
#     return status
