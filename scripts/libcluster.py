import os
import copy
import mdtraj
import numpy as np
from sklearn.cluster import DBSCAN

from libmd_contact import get_dist_distr

def get_clusters(frame_s, rmsd_cutoff=2.0, subsample=0, get_distr=False):
    if subsample > 0:
        subsampled = np.zeros(len(frame_s), dtype=bool)
        subsampled[::subsample] = True
        notsampled = np.logical_not(subsampled)
        frame_s0 = copy.deepcopy(frame_s)
        frame_s = frame_s0[subsampled]
        frame_x = frame_s0[notsampled]
    #
    dmtx = np.array([mdtraj.rmsd(frame_s, frame) for frame in frame_s], dtype=float)
    if os.getenv("PREFMD2_PYTHON_MPROCS") is None:
        n_jobs = -1
    else:
        n_jobs = int(os.getenv("PREFMD2_PYTHON_MPROCS"))
    dbscan = DBSCAN(eps=rmsd_cutoff*0.1, min_samples=1, metric='precomputed',
                    n_jobs=n_jobs)
    labels = dbscan.fit_predict(dmtx)
    n_cluster = labels.max()+1
    #
    if subsample > 0:
        rmsd_to_center = []
        for i in range(n_cluster):
            c = (labels == i)
            index = np.ix_(c,c)
            d = dmtx[index].sum(axis=0)
            #
            member = frame_s[c]
            center = member[np.argmin(d)]
            rmsd_to_center.append(mdtraj.rmsd(frame_x, center))
        rmsd_to_center = np.array(rmsd_to_center, dtype=float)
        label_x = np.argmin(rmsd_to_center, axis=0)
    #
    cluster_s = []
    for i in range(n_cluster):
        c = (labels == i)
        index = np.ix_(c,c)
        d = dmtx[index].sum(axis=0)
        #
        member = frame_s[c]
        center = member[np.argmin(d)]
        if subsample > 0:
            x = (label_x == i)
            member = mdtraj.join([member, frame_x[x]])
        n_member = len(member)
        #
        superposed = member.superpose(center)
        xyz = superposed.xyz.mean(axis=0)
        averaged = copy.deepcopy(center)
        averaged.xyz = xyz
        if get_distr:
            distr = get_dist_distr(superposed)
            cluster_s.append((n_member, center, averaged, distr))
        else:
            cluster_s.append((n_member, center, averaged))
    cluster_s.sort(key=lambda x: x[0], reverse=True)
    return cluster_s

#def main():
#    arg = argparse.ArgumentParser(prog='cluster_models')
#    arg.add_argument(dest='output_prefix')
#    arg.add_argument(dest='init_pdb')
#    arg.add_argument('--dcd', dest='dcd_fn_s', required=True, nargs='*')
#    arg.add_argument('--rmsd', dest='rmsd_cutoff', type=float, default=2.0)
#    arg.add_argument('--subsample', dest='subsample', type=int, default=0)
#    #
#    if len(sys.argv) == 1:
#        arg.print_help()
#        return
#    arg = arg.parse_args()
#    #
#    run(arg)
#
#if __name__ == '__main__':
#    main()
