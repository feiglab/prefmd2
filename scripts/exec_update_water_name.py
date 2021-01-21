#!/usr/bin/env python

import sys

def main(args=None):
    if args is not None:
        pdb_fn = args[-1]
    else:
        pdb_fn = sys.argv[1]
    #
    pdb = []
    with open(pdb_fn) as fp:
        for line in fp:
            if not line.startswith("ATOM"):
                pdb.append(line)
                continue
            #
            resName = line[17:20].strip()
            if resName not in ['TIP', 'HOH', 'WAT']:
                pdb.append(line)
                continue
            #
            atmName = line[12:16].strip()
            if atmName == 'O':
                atmName = ' OH2'
            else:
                atmName = ' %-3s'%atmName
            line = ''.join([line[:12], atmName, ' TIP3', line[21:]])
            pdb.append(line)
    with open(pdb_fn, 'wt') as fout:
        fout.writelines(pdb)

if __name__ == '__main__':
    main()
