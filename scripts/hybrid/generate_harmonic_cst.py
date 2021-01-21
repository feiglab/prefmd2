#!/usr/bin/env python

import os
import sys
import numpy as np

from hybrid import seqName

SEQUENCE_SEPARATION = 9
PARAM_DEVIATION_CUTOFF = 2.0

def read_pdb(pdb_fn):
    pdb = {}
    with open(pdb_fn) as fp:
        for line in fp:
            if not line.startswith("ATOM"): continue
            #
            resName = seqName.to_std(line[17:20])
            if resName not in seqName.stdres: continue
            #
            atmName = line[12:16].strip()
            if resName != 'GLY' and atmName != 'CB':
                continue
            elif resName == 'GLY' and atmName != 'CA':
                continue
            #
            resNo = int(line[22:26])
            chain_id = ''
            #
            key = (chain_id, resNo, resName, atmName)
            if key in pdb: continue
            #
            pdb[key] = np.array([line[30:38], line[38:46], line[46:54]], dtype=float)
    return pdb

def dist(r_1, r_2):
    return np.sqrt(np.sum((r_2-r_1)**2))

def get_gremlin_d_and_width(aa_1, aa_2):
    aa_s = ["GLY","ALA","SER","VAL","CYS","THR","PRO","ASP","ASN","ILE","LEU","GLU","GLN","MET","HIS","LYS","PHE","TYR","ARG","TRP"]
    i_1 = aa_s.index(aa_1) ; i_2 = aa_s.index(aa_2)
    #
    dw = [[4.467,5.201,5.510,5.671,5.777,5.619,6.140,6.135,6.321,6.413,6.554,7.036,7.297,7.383,7.472,8.216, 7.966, 9.098, 9.166, 8.966],\
          [0,    5.381,5.829,5.854,6.057,5.982,6.412,6.388,6.766,6.587,6.707,7.124,7.583,7.605,7.591,8.327, 8.162, 9.121, 9.365, 9.252],\
          [0,    0,    6.190,6.567,6.590,6.450,6.937,6.760,7.081,7.142,7.394,7.483,7.807,8.010,8.051,8.792, 8.694, 9.594, 9.753, 9.770],\
          [0.017,0,    0,    6.759,6.941,6.791,7.063,6.972,7.219,7.441,7.633,7.404,8.008,8.335,8.179,8.077, 9.057, 9.442, 9.513,10.021],\
          [0.269,0.262,0,    0,    6.426,6.801,7.157,6.985,7.205,7.476,7.685,7.449,7.962,8.265,8.422,8.494, 9.026, 9.362, 9.460, 9.752],\
          [0.153,0.291,0.292,0,    0,    6.676,7.062,6.971,7.159,7.442,7.642,7.628,8.055,8.397,8.221,8.715, 9.030, 9.813, 9.764, 9.980],\
          [0.107,0.312,0.205,0.145,0,    0,    7.288,7.321,7.497,7.554,7.751,7.938,8.308,8.247,8.537,9.198, 8.895, 9.965,10.266, 9.719],\
          [0.129,0.394,0.240,0.173,0.178,0,    0,    8.001,7.672,7.472,7.696,8.945,8.601,8.401,8.634,9.306, 9.111, 9.979,10.123, 9.867],\
          [0.120,0.378,0.214,0.138,0.181,0.188,0,    0,    7.682,7.631,7.889,8.485,8.502,8.550,8.672,9.319, 9.168,10.039,10.135, 9.976],\
          [0.245,0.399,0.321,0.298,0.259,0.320,0.339,0,    0,    8.096,8.342,7.949,8.302,8.874,8.523,8.329, 9.602, 9.719, 9.746,10.470],\
          [0.193,0.289,0.323,0.287,0.299,0.307,0.416,0.392,0,    0,    8.522,8.077,8.480,9.122,8.676,8.479, 9.900, 9.889, 9.852,10.707],\
          [0.169,0.349,0.305,0.232,0.240,0.262,0.334,0.337,0.249,0,    0,    9.863,9.328,8.870,9.454,9.842, 9.403,10.544,10.713,10.303],\
          [0.179,0.214,0.342,0.242,0.295,0.259,0.336,0.341,0.341,0.321,0,    0,    9.074,9.102,9.391,9.667, 9.506,10.534,10.610,10.429],\
          [0.125,0.250,0.287,0.179,0.206,0.190,0.317,0.348,0.279,0.261,0.198,0,    0,    9.530,9.396,9.096,10.253,10.400,10.250,11.110],\
          [0.249,0.340,0.446,0.510,0.538,0.409,0.475,0.354,0.423,0.453,0.475,0.389,0,    0,   10.607,9.582, 9.602,10.843,10.879,10.661],\
          [0.216,0.356,0.408,0.359,0.347,0.378,0.410,0.357,0.373,0.406,0.411,0.450,0.436,0,    0,   10.662, 9.344,10.627,11.322,10.136],\
          [0.255,0.394,0.369,0.295,0.439,0.292,0.388,0.361,0.310,0.327,0.318,0.511,0.498,0.457,0,    0,    10.903,10.999,10.577,11.758],\
          [0.206,0.380,0.435,0.383,0.203,0.417,0.457,0.325,0.289,0.379,0.401,0.443,0.401,0.342,0.333,0,     0,    11.536,11.615,11.807],\
          [0.358,0.550,0.445,0.634,0.521,0.464,0.550,0.343,0.398,0.582,0.591,0.434,0.521,0.611,0.714,0.738, 0,     0,    12.050,11.355],\
          [0.219,0.260,0.394,0.246,0.286,0.264,0.425,0.351,0.393,0.347,0.260,0.512,0.451,0.377,0.542,0.441, 0.460, 0,     0,    12.806],\
	  [0.267,0.443,0.467,0.535,0.585,0.430,0.506,0.676,0.586,0.589,0.611,0.469,0.547,0.661,0.554,0.704, 0.767, 0.855, 0,     0,   ],\
	  [0.334,0.485,0.483,0.514,0.491,0.477,0.506,0.327,0.372,0.557,0.578,0.363,0.535,0.641,0.595,0.648, 0.738, 0.822, 0.704, 0,   ],\
	  [0.239,0.290,0.497,0.271,0.417,0.315,0.462,0.475,0.458,0.397,0.331,0.493,0.490,0.397,0.458,0.470, 0.447, 0.684, 0.889, 0.473]]
    #
    if i_1 >= i_2:
        return dw[i_2][i_1], dw[i_1+2][i_2]
    else:
        return dw[i_1][i_2], dw[i_2+2][i_1]

def get_contact_list(pdb):
    contact_s = []
    d0_s = []
    dij_s = []
    #
    key_s = sorted(pdb.keys()) ; n_res = len(key_s)
    for i in range(n_res-1):
        res_1 = key_s[i]
        for j in range(i+1, n_res):
            res_2 = key_s[j]
            #
            if res_1[0] == res_2[0] and abs(res_2[1]-res_1[1]) < SEQUENCE_SEPARATION:
                continue
            #
            d0, w = get_gremlin_d_and_width(res_1[2], res_2[2])
            dij = dist(pdb[res_1], pdb[res_2])
            if dij < d0 + 1.0 + 2.0*w:
                contact_s.append((res_1, res_2))
                d0_s.append(d0)
                dij_s.append(dij)
    return contact_s, d0_s, dij_s

def report_cst(pool, contact_s, d0_s, dij_s):
    n_pool = len(pool)
    for i in range(len(contact_s)):
        res_1, res_2 = contact_s[i]
        #
        dist_s = []
        dev_s = []
        for k in range(n_pool):
            if res_1 not in pool[k]: continue
            if res_2 not in pool[k]: continue
            dk = dist(pool[k][res_1], pool[k][res_2])
            dist_s.append(dk-dij_s[i])
            dev = dk - d0_s[i]
            dev_s.append(dev)
        dev_s = np.array(dev_s, dtype=float)
        w = np.mean(1.0/(1.0+dev_s**2))
        #
        dist_s = np.sqrt(np.mean(np.array(dist_s, dtype=float)**2))
        #
        wrt = []
        wrt.append("AtomPair %3s %4d %3s %4d "%(res_1[3], res_1[1], res_2[3], res_2[1]))
        #
        if dist_s > PARAM_DEVIATION_CUTOFF:
            wrt.append("SCALARWEIGHTEDFUNC %8.5f "%w)
            wrt.append("BOUNDED 0 %8.3f 1.0 0.5\n"%(d0_s[i] + 2.0))
        else:
            wrt.append("SCALARWEIGHTEDFUNC %8.5f "%0.25)
            wrt.append("FLAT_HARMONIC %8.3f 1.0 1.0\n"%dij_s[i])
        #
        sys.stdout.write("".join(wrt))

def main(args=None):
    if args is None:
        if len(sys.argv) == 1:
            sys.exit("USAGE: %s [initial model] [other models]"%__file__)
            return
        pdb_fn_s = sys.argv[1:]
    else:
        pdb_fn_s = args[0:]
    #
    #
    pdb = read_pdb(pdb_fn_s[0])
    contact_s, d0_s, dij_s = get_contact_list(pdb)
    #
    pool = [pdb] + [read_pdb(fn) for fn in pdb_fn_s[1:]]
    report_cst(pool, contact_s, d0_s, dij_s)

if __name__ == '__main__':
    main()
