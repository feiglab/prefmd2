#!/usr/bin/env python

import os
import sys
import argparse
import subprocess as sp
from tempfile import mkdtemp

AAs = ('ALA','CYS','ASP','GLU','PHE','GLY','HIS','HIP','HID','HIE','HSE','HSD','HSP','ILE','LYS','LEU',\
       'MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR')
WATERs = ('TIP')
atmName_update = {"O":" OT1", "OXT":" OT2"}

def read_pdb(fn):
    segName_s = []
    seg_s = {}
    disu_s = []
    his_s = {}
    with open(fn) as fp:
        for line in fp:
            if (not line.startswith("ATOM")) and (not line.startswith("HETA")):
                if line.startswith("SSBOND"):
                    chain_1 = line[15] ; resNo_1 = int(line[17:21])
                    chain_2 = line[29] ; resNo_2 = int(line[31:35])
                    disu_s.append([[chain_1, resNo_1], [chain_2, resNo_2]])
                continue
            segName = line[72:76].strip()
            if segName not in segName_s:
                is_aa = (line[17:20].strip() in AAs)
                is_wat = (line[17:20].strip() in WATERs)
                segName_s.append(segName)
                if is_aa:
                    seg_s[segName] = [0, []]
                elif is_wat:
                    seg_s[segName] = [2, []]
                else:
                    seg_s[segName] = [1, []]
            resName = line[17:20].strip()
            atmName = line[12:16].strip()
            if resName in ['HIS']:
                if atmName == 'HD1':    # HSD
                    resNo = line[22:26]
                    his_s[(segName, resNo)] = 'HSD'
                elif atmName == 'HE2':  # HSE
                    resNo = line[22:26]
                    his_s[(segName, resNo)] = 'HSE'
                #resName_new = resName_update[resName]
                #line = line[:17] + resName_new + line[20:]
            if (resName in AAs) and (atmName in atmName_update):
                atmName_new = atmName_update[atmName]
                line = line[:12] + atmName_new + line[16:]
            if line.startswith("HETA"):
                line = 'ATOM  %s'%line[6:]

            seg_s[segName][1].append(line)
    #
    for segName in seg_s:
        if seg_s[segName][0] != 0: continue
        for i,line in enumerate(seg_s[segName][1]):
            resName = line[17:20].strip()
            if resName != 'HIS': continue
            resNo = line[22:26]
            if (segName, resNo) in his_s:
                hisName = his_s[(segName, resNo)]
            else:
                hisName = 'HSD'
            line = line[:17] + hisName + line[20:]
            seg_s[segName][1][i] = line
    #
    def find_seg(cys, seg_s):
        seg = None
        for segName in seg_s:
            if seg_s[segName][0] != 0: continue
            for line in seg_s[segName][1]:
                chain_id = line[21]
                resNo = int(line[22:26])
                if chain_id == cys[0] and resNo == cys[1]:
                    return segName
        return seg
    #
    disu_out = []
    for i,disu in enumerate(disu_s):
        seg_0 = find_seg(disu[0], seg_s)
        seg_1 = find_seg(disu[1], seg_s)
        if seg_0 is not None and seg_1 is not None:
            disu_out.append((seg_0, disu[0][1], seg_1, disu[1][1]))
    return segName_s, seg_s, disu_out

def parse_patch(args):
    patch_s = []
    for arg in args:
        patch,res = arg.split(":")
        segName,resNo = res.split(".")
        patch_s.append((patch, segName, resNo))
    return patch_s

def write_top_cmd(toppar):
    rtf_s = [] ; prm_s = [] ; str_s = []
    for fn in toppar:
        if fn.endswith("rtf"):
            rtf_s.append(fn)
        elif fn.endswith("prm"):
            prm_s.append(fn)
        else:
            str_s.append(fn)
    #
    cmd = []
    cmd.append("bomlev -1\n")
    cmd.append("!\n")
    for i,fn in enumerate(rtf_s):
        cmd.append('open unit 10 read form name "%s"\n'%fn)
        cmd.append('read rtf card unit 10')
        if i > 0:
            cmd.append(' append')
        if len(str_s) > 0:
            cmd.append(' flex')
        cmd.append("\n")
        cmd.append('close unit 10\n')
    cmd.append("!\n")
    for i,fn in enumerate(prm_s):
        cmd.append('open unit 10 read form name "%s"\n'%fn)
        cmd.append('read para card unit 10')
        if i > 0:
            cmd.append(' append')
        if len(str_s) > 0:
            cmd.append(' flex')
        cmd.append("\n")
        cmd.append('close unit 10\n')
    cmd.append("!\n")
    cmd.append("faster on\n")
    cmd.append("!\n")
    for fn in str_s:
        cmd.append('open unit 78 read form name "%s"\n'%fn)
        cmd.append("stream unit 78\n")
    cmd.append("!\n")
    return cmd

def split_seg(segName_s, seg_s):
    tmpdir = mkdtemp(prefix='genPSF.')
    tmpdir = os.path.abspath(tmpdir)
    #
    for segName in segName_s:
        out = '%s/%s'%(tmpdir, segName)
        with open(out, 'wt') as fout:
            fout.writelines(seg_s[segName][1])
            fout.write("TER\n")
            fout.write("END\n")
        seg_s[segName].append(out)
    return tmpdir

def write_pdb_cmd(tmpdir, segName_s, seg_s, disu_s, patch_s, write_crd=False, blocked=False, terminal=['ACE', 'CT3']):
    cmd = []
    for segName in segName_s:
        if seg_s[segName][0] != 0: continue
        #
        fn = seg_s[segName][2]
        cmd.append('open unit 10 read form name "%s"\n'%fn)
        cmd.append("read sequ pdb unit 10\n")
        if blocked:
            cmd.append("generate firs %s last %s"%tuple(terminal) + " %s setup warn\n"%segName)
        else:
            cmd.append("generate %s setup warn\n"%segName)
        cmd.append("close unit 10\n")
    for disu in disu_s:
        cmd.append("patch DISU %s %d %s %d\n"%(disu))
    for patch in patch_s:
        cmd.append("patch %s %s %s setup\n"%patch)
    cmd.append("auto angle dihe\n")
    for segName in segName_s:
        if seg_s[segName][0] != 1: continue
        #
        fn = seg_s[segName][2]
        cmd.append('open unit 10 read form name "%s"\n'%fn)
        cmd.append("read sequ pdb unit 10\n")
        cmd.append("generate %s setup warn\n"%segName)
        cmd.append("close unit 10\n")
    for segName in segName_s:
        if seg_s[segName][0] != 2: continue
        #
        fn = seg_s[segName][2]
        cmd.append('open unit 10 read form name "%s"\n'%fn)
        cmd.append("read sequ pdb unit 10\n")
        cmd.append("generate %s setup noangl nodihe\n"%segName)
        cmd.append("close unit 10\n")
    cmd.append("!\n")
    #
    for segName in segName_s:
        if seg_s[segName][0] != 0: continue
        #
        fn = seg_s[segName][2]
        cmd.append('open unit 10 read form name "%s"\n'%fn)
        cmd.append("read coor pdb unit 10 resi\n")
        cmd.append("close unit 10\n")
    for segName in segName_s:
        if seg_s[segName][0] != 1: continue
        #
        fn = seg_s[segName][2]
        cmd.append('open unit 10 read form name "%s"\n'%fn)
        cmd.append("read coor pdb unit 10 resi\n")
        cmd.append("close unit 10\n")
    for segName in segName_s:
        if seg_s[segName][0] != 2: continue
        #
        fn = seg_s[segName][2]
        cmd.append('open unit 10 read form name "%s"\n'%fn)
        cmd.append("read coor pdb unit 10 resi\n")
        cmd.append("close unit 10\n")
    #
    cmd.append("!\n")
    cmd.append("bomlev -2\n")
    cmd.append("ic param\n")
    cmd.append("ic build\n")
    cmd.append("!\n")
    cmd.append('open unit 10 write form name "out.psf"\n')
    cmd.append("write psf card unit 10\n")
    cmd.append("*\n")
    cmd.append("close unit 10\n")
    cmd.append("!\n")
    #
    if write_crd:
        cmd.append('open unit 10 write form name "out.crd"\n')
        cmd.append("write coor card unit 10\n")
        cmd.append("*\n")
        cmd.append("close unit 10\n")
        cmd.append("!\n")
    #
    return cmd

def main(args=None):
    arg = argparse.ArgumentParser(prog='genPSF')
    arg.add_argument(dest='init_pdb')
    arg.add_argument('-ff',  '--toppar', dest='toppar', nargs='*')
    arg.add_argument('-psf', '--psfout', dest='psfout', default='out.psf')
    arg.add_argument('-crd', '--crdout', dest='crdout', default=None)
    arg.add_argument('-patch', '--patch', dest='patch_s', default=[], nargs='*')
    arg.add_argument('-blocked', '--blocked', dest='blocked', default=False, action='store_true')
    arg.add_argument('-terminal', '--terminal', dest='terminal', default=['ACE', 'CT3'], nargs=2)
    #
    if len(sys.argv) == 1:
        arg.print_help()
        return
    arg = arg.parse_args(args)
    #
    arg.toppar = [os.path.realpath(fn) for fn in arg.toppar]
    #
    patch_s = parse_patch(arg.patch_s)
    segName_s, seg_s, disu_s = read_pdb(arg.init_pdb)
    tmpdir = split_seg(segName_s, seg_s)
    pwd = os.getcwd()
    #
    cmd_s = []
    cmd_s.extend(write_top_cmd(arg.toppar))
    cmd_s.extend(write_pdb_cmd(tmpdir, segName_s, seg_s, disu_s, patch_s, \
            write_crd=(arg.crdout is not None), blocked=arg.blocked, terminal=arg.terminal))
    cmd_s.append("stop\n")
    #
    with open("%s/genPSF.cmd"%tmpdir, 'wt') as fout:
        fout.writelines(cmd_s)
    #
    stdin = open("%s/genPSF.cmd"%tmpdir)
    stdout = open('%s/genPSF.log'%tmpdir, 'wt')
    #
    os.chdir(tmpdir)
    sp.call([os.environ['CHARMMEXEC']], stdin=stdin, stdout=stdout)
    os.chdir(pwd)
    #
    stdin.close()
    stdout.close()
    #
    if os.path.exists("%s/out.psf"%tmpdir):
        sp.call(['mv', '%s/out.psf'%tmpdir, arg.psfout])
        if arg.crdout is not None:
            if os.path.exists("%s/out.crd"%tmpdir):
                sp.call(['mv', '%s/out.crd'%tmpdir, arg.crdout])
            else:
                sys.stdout.write("ERROR: %s\n"%tmpdir)
                return
        #
        os.system("rm -rf %s"%tmpdir)
    else:
        sys.stdout.write("ERROR: %s\n"%tmpdir)

if __name__=='__main__':
    main()
