#!/usr/bin/env python

import sys
import argparse

resName_s = ('ALA','CYS','CYM','ASP','GLU','PHE',\
             'GLY','HIS','HSE','HSD','HSP','ILE','LYS','LEU',\
             'MET','ASN','PRO','GLN','ARG',\
             'SER','THR','VAL','TRP','TYR')

def main(args=None):
    arg = argparse.ArgumentParser(prog='cmap_resname_convert')
    arg.add_argument(dest='fp', metavar='FILE', nargs='?', default=sys.stdin,
                     type=argparse.FileType('r'))
    arg.add_argument('--cmap', dest='use_cmap', default=False, action='store_true')
    arg.add_argument('--res', dest='residue_s', default=None)
    #
    arg = arg.parse_args(args)
    #
    if arg.residue_s is not None:
        residue_s = []
        for r in arg.residue_s.split(","):
            if '-' in r:
                x = r.split("-")
                residue_s.extend(range(int(x[0]), int(x[1])+1))
            else:
                residue_s.append(int(r))
        arg.residue_s = residue_s
    #
    for line in arg.fp:
        if not (line.startswith("ATOM") or line.startswith("HETA")):
            sys.stdout.write(line)
            continue
        #
        resNo = int(line[22:26])
        if (arg.residue_s is not None) and (resNo not in arg.residue_s):
            sys.stdout.write(line)
            continue
        #
        resName = line[17:20]
        if resName not in resName_s:
            sys.stdout.write(line)
            continue
        #
        if arg.use_cmap:
            sys.stdout.write("ATOM  %s %s"%(line[6:20], line[21:]))
        else:
            sys.stdout.write("ATOM  %s9%s"%(line[6:20], line[21:]))

if __name__=='__main__':
    try:
        main()
    except IOError:
        sys.exit()
