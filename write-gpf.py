#!/usr/bin/env python

import os
import sys
import math
import subprocess
import argparse

supported_atypes = ['HD', 'C', 'A', 'N', 'NA', 'OA', 'F', 'P', 'SA', 'S',
                    'Cl', 'Br', 'I', 'Mg', 'Ca', 'Mn', 'Fe', 'Zn', 'H', 'OC']

gpf = """npts NPTS_X NPTS_Y NPTS_Z
gridfld PREFIX.maps.fld
spacing 0.375
receptor_types RECTYPES
ligand_types HD C A N NA OA F P SA S Cl Br I H
receptor REC
gridcenter CENTER_X CENTER_Y CENTER_Z
smooth 0.5
map         PREFIX.HD.map
map         PREFIX.C.map
map         PREFIX.A.map
map         PREFIX.N.map
map         PREFIX.NA.map
map         PREFIX.OA.map
map         PREFIX.F.map
map         PREFIX.P.map
map         PREFIX.SA.map
map         PREFIX.S.map
map         PREFIX.Cl.map
map         PREFIX.Br.map
map         PREFIX.I.map
map         PREFIX.H.map
elecmap     PREFIX.e.map
dsolvmap    PREFIX.d.map
dielectric -0.1465
"""

def getrectypes(fname):
    command = 'cut -c 77-79 %s | sort -u' % fname
    out = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    out = out.communicate()[0]
    rectypes = []
    for atype in out.split():
        atype = atype.strip().decode()
        if atype in supported_atypes:
            rectypes.append(atype)
    return ' '.join(rectypes)

def getbox(fname):
    with open(fname) as f:
        for line in f:
            if line.startswith('center_x'): center_x = float(line.split()[2])
            if line.startswith('center_y'): center_y = float(line.split()[2])
            if line.startswith('center_z'): center_z = float(line.split()[2])
            if line.startswith('size_x'): size_x = float(line.split()[2])
            if line.startswith('size_y'): size_y = float(line.split()[2])
            if line.startswith('size_z'): size_z = float(line.split()[2])
    npts_x = 2 * int(size_x / 0.75)
    npts_y = 2 * int(size_y / 0.75)
    npts_z = 2 * int(size_z / 0.75)
    return center_x, center_y, center_z, npts_x, npts_y, npts_z

def calcbox(fname, pad, spacing=0.375):
    x_min = float('inf')
    y_min = float('inf')
    z_min = float('inf')
    x_max = float('-inf')
    y_max = float('-inf')
    z_max = float('-inf')
    with open(fname) as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                x_max = max(x, x_max)
                y_max = max(y, y_max)
                z_max = max(z, z_max)
                x_min = min(x, x_min)
                y_min = min(y, y_min)
                z_min = min(z, z_min)
    center_x = (x_min + x_max) / 2.0
    center_y = (y_min + y_max) / 2.0
    center_z = (z_min + z_max) / 2.0
    npts_x = math.ceil((2 * pad + x_max - x_min) / spacing)  
    npts_y = math.ceil((2 * pad + y_max - y_min) / spacing)  
    npts_z = math.ceil((2 * pad + z_max - z_min) / spacing)  
    return center_x, center_y, center_z, npts_x, npts_y, npts_z


class MyParser(argparse.ArgumentParser):
    """display help message for every error"""
    def error(self, message):
        self.print_help()
        sys.stderr.write('\nERROR:\n  %s\n' % message)
        sys.exit(2)

def get_args():
    parser = MyParser()
    parser.add_argument('rec', help='receptor [.pdbqt]')
    parser.add_argument('-b', '--box', help='vina box [.config]')
    parser.add_argument('-l', '--lig', help='ligand to center box  [.pdbqt]')
    parser.add_argument('-p', '--pad', help='padding around ligand [angstrom]', default=8.0, type=float)
    parser.add_argument('--mapprefix')
    args = parser.parse_args()
    if args.mapprefix == None:
        args.mapprefix = os.path.splitext(os.path.basename(args.rec))[0]
    if (args.box == None) and (args.lig == None):
        sys.stderr.write('Use either --box or --lig\n')
        sys.exit(2)
    if args.box and args.lig:
        sys.stderr.write('Use either --box or --lig\n')
        sys.exit(2)
    args.gpf = args.mapprefix + '.gpf'
    return args

args = get_args()

if os.path.exists(args.gpf):
    sys.stderr.write('Aborting! %s already exists.\n' % args.gpf)    

rectypes = getrectypes(args.rec)
if args.box:
    center_x, center_y, center_z, npts_x, npts_y, npts_z = getbox(args.box)
else:
    center_x, center_y, center_z, npts_x, npts_y, npts_z = calcbox(args.lig, args.pad)

gpf = gpf.replace('RECTYPES',   rectypes)
gpf = gpf.replace('PREFIX',     args.mapprefix)
gpf = gpf.replace('REC',        args.rec)
gpf = gpf.replace('NPTS_X',     '%d' % npts_x)
gpf = gpf.replace('NPTS_Y',     '%d' % npts_y)
gpf = gpf.replace('NPTS_Z',     '%d' % npts_z)
gpf = gpf.replace('CENTER_X',   '%.3f' % center_x)
gpf = gpf.replace('CENTER_Y',   '%.3f' % center_y)
gpf = gpf.replace('CENTER_Z',   '%.3f' % center_z)

with open(args.gpf, 'w') as f:
    f.write(gpf)
