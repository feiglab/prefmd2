#!/usr/bin/env python

import sys

def read_resNo(pdb_fn):
    ssbond = []
    resNo_s = []
    #
    resNo_prev = None
    with open(pdb_fn) as fp:
        for line in fp:
            if not (line.startswith("ATOM")):
                if line.startswith("TER"):
                    resNo_prev = None
                elif line.startswith("SSBOND"):
                    ssbond.append(line)
                continue
            #resNo = line[22:26]
            resNo = line[21:26]
            if resNo != resNo_prev:
                resNo_prev = resNo
                resNo_s.append(resNo)
    return resNo_s, ssbond

def update_resNo(pdb_fn, resNo_s, ssbond):
    k = -1
    for line in ssbond:
        sys.stdout.write(line)
    with open(pdb_fn) as fp:
        resNo_prev = None
        for line in fp:
            if not (line.startswith("ATOM") or line.startswith("HETA")):
                if line.startswith("TER"):
                    resNo_prev = None
                    sys.stdout.write("TER\n")
                elif line.startswith("SSBOND"):
                    pass
                elif line.startswith("ENDMDL"):
                    sys.stdout.write(line)
                elif line.startswith("END"):
                    break
                else:
                    sys.stdout.write(line)
                continue
            #resNo = line[22:26]
            resNo = line[21:26]
            if resNo != resNo_prev:
                k += 1
                if k < len(resNo_s):
                    resNo_prev = resNo
                    resNo_new = resNo_s[k]
                else:
                    break
            #new = 'ATOM  %s %4s%s'%(line[6:21], resNo_new, line[26:])
            new = 'ATOM  %s%4s%s'%(line[6:21], resNo_new, line[26:])
            sys.stdout.write(new)
    #sys.stdout.write("TER\n")
    sys.stdout.write("END\n")

def main(args=None):
    if args is None:
        if len(sys.argv) < 3:
            sys.stderr.write("USAGE: %s [REFERENCE] [PDB]\n"%__file__)
            return
        ref_fn = sys.argv[1]
        pdb_fn = sys.argv[2]
    else:
        ref_fn = args[0]
        pdb_fn = args[1]
    #
    resNo_s, ssbond = read_resNo(ref_fn)
    update_resNo(pdb_fn, resNo_s, ssbond)

if __name__=='__main__':
    main()
