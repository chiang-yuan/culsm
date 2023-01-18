'''
Img2Particle

Img2Particel is a simple, in-house developed python code for conversion from
gray image to hexagonal packing particle model in LAMMPS Data format.

MIT License

Copyright (c) 6th/June/2020 Yuan Chiang
'''
import os
import argparse
import glob
import numpy as np
import math
from matplotlib import pyplot as plt
from PIL import Image
from datetime import datetime

version = '1.0'

parser = argparse.ArgumentParser(prog='Img2Particle',
                                 description='Convert images into particle models.')

parser.add_argument('infile', metavar='images', type=str,
                    help='absolute or relative address of image files')
parser.add_argument('-o', dest='outfile', metavar='data file', type=str)
parser.add_argument('-lx', dest='lx', metavar='Lx', type=float,
                    help='model size in x dimension')
parser.add_argument('-ly', dest='ly', metavar='Ly',type=float,
                    help='model size in y dimension')
parser.add_argument('-lz', dest='lz', metavar='Lz',type=float, default=1,
                    help='model size in z dimension\t(default = 1)')
parser.add_argument('-s', dest='s', metavar='lc', type=float,
                    help='lattice constant = particle spacing along (100)')
parser.add_argument('-t', dest='types', metavar='types', type=int, default=2,
                    help='number of particle types\t(default = 2)')
parser.add_argument('-n',dest='nsize', metavar='len', nargs=2, type=float, default=[0, 0],
                    help='create notch with lengths in x and y dimension\t(default = [0, 0])')
parser.add_argument('--single', action='store_true')
parser.add_argument('--verbose', action='store_true')
parser.add_argument('--version', action='version', version='%(prog)s {:s}'.format(version))

args = parser.parse_args()

infile = args.infile

def import_img(infile):
    img = Image.open(infile).convert('L')
    img.load()
    array = np.asarray(img, dtype="int32")
    return array

img = import_img(infile)
print("Load image from file: {}".format(infile))
print("\timage size: {}".format(img.shape))

# ===== Define Model

lx = args.lx
ly = args.ly
lz = args.lz
s = args.s

print("Set up model...")
print("\tlx = {:f}\tly = {:f}\tlz = {:f}".format(lx, ly, lz))
print("\tlattice constant = {:f}".format(s))

# ===== Define Lattice

# create hexgonal packing unit cell
unit = s*np.array([[1,0,0],[0,math.sqrt(3.),0],[0,0,0]], dtype=float) # distant unit

nx = int(round(lx/unit[0,0]))
ny = int(round(ly/unit[1,1]))
nz = 1

mx = 2
my = 0
mz = 0

print("Hexagonal lattice:\n{}".format(unit))
print("\tnumber of lattices and margins in x: %3d / %3d" % (nx,mx))
print("\tnumber of lattices and margins in y: %3d / %3d" % (ny,my))
print("\tnumber of lattices and margins in z: %3d / %3d" % (nz,mz))

nx = nx + 2*mx
ny = ny + 2*my
nz = nz + 2*mz

# ===== Create Atoms

dtype = dict(zip(['id', 'molecule', 'type', 'q', 'x', 'y', 'z', 'ix', 'iy', 'iz'], [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]))

natoms = int(2*nx*ny)

atoms = np.zeros((natoms,10))
atoms[:,dtype.get('id')] = np.arange(1,natoms+1)
atoms[:,dtype.get('molecule')] = np.ones((natoms,1)).reshape(-1)

id_ = 0
colx = dtype.get('x')
coly = dtype.get('y')
colz = dtype.get('z')
for i in range(nx):
    for j in range(ny):
        for k in range(nz):
            atoms[id_,colx:colz+1] = np.dot(np.transpose(unit),[[i-mx],[j-my],[k-mz]]).reshape(-1)
            id_ = id_ + 1
            atoms[id_,colx:colz+1] = np.dot(np.transpose(unit),[[i+0.5-mx],[j+0.5-my],[k-mz]]).reshape(-1)
            id_ = id_ + 1


print("Create atoms: {:d}".format(natoms))

# ===== Change Types

ntypes = args.types

rows, cols = img.shape
colt = dtype.get('type')

hist, bin_edges = np.histogram(img.reshape(-1), bins=ntypes)

print(hist, bin_edges)

for i in range(natoms):
    x = int(round(atoms[i,colx]/lx*cols))
    y = rows - 1 - int(round(atoms[i,coly]/ly*rows))
    # if x < 0 or x >= cols or y < 0 or y >= rows:
    #     atoms[i,colt] = 1
    if x < 0:
        atoms[i,colt] = ntypes+1
    elif x >= cols:
        atoms[i,colt] = ntypes+2
    else:
        for b in range(ntypes):
            if img[y,x] >= bin_edges[b]:
                atoms[i,colt] = ntypes - b

ntypes = int(max(atoms[:,colt]))

# ===== Create Notch

wn = args.nsize[0]
hn = args.nsize[1]

indices = np.logical_not(
    np.logical_and(
        np.logical_and(atoms[:,colx] < lx/2.0 + wn/2.0, atoms[:,colx] > lx/2.0 - wn/2.0),
        atoms[:,coly] < hn
    )
)

atoms = atoms[indices]
natoms = atoms.shape[0]
atoms[:,dtype.get('id')] = np.arange(1,natoms+1)
color = [str(item/(ntypes+1)) for item in atoms[:,colt]]

# ===== Single or Composite

if args.single:
    atoms = atoms[atoms[:,dtype.get('type')] == 1]
    natoms = atoms.shape[0]
    atoms[:,dtype.get('id')] = np.arange(1,natoms+1)
    color = [str(item/(ntypes+1)) for item in atoms[:,colt]]

if args.verbose:
    plt.figure(figsize=[12.8, 9.6], dpi=330, facecolor=None, edgecolor=None, frameon=True)

    plt.subplot(121)
    plt.imshow(img,cmap='gist_gray')

    plt.subplot(122)
    plt.scatter(atoms[:,colx],atoms[:,coly],s=0.1,marker="o",c=color)
    ax = plt.gca()
    ax.set_aspect(1.0)
    plt.savefig(infile+'_conversion.png')

if args.outfile is None:
    outfile = os.path.splitext(infile)[0] + '.data'
else:
    if args.outfile.lower().endswith('data'):
        outfile = args.outfile
    else:
        outfile = args.outfile + '.data'

print('Write data file: {}'.format(outfile))

with open(outfile, "w") as outfile:
    now = datetime.now().strftime("%Y/%m/%d %H:%M:%S")
    outfile.write("LAMMPS data created via Img2Particle at "+now+"\n")
    outfile.write("\n%d atoms\n" % natoms)
    outfile.write('# use create_bonds command in lammps to create springs\n')
    outfile.write("\n%d atom types\n" % ntypes)
    outfile.write("%d bond types\n" % int(math.factorial(ntypes+2-1)/(math.factorial(2)*math.factorial(ntypes-1))))

    outfile.write("\n%12.5E %12.5E xlo xhi\n" % (np.min(atoms[:,colx]),np.max(atoms[:,colx])))
    outfile.write("%12.5E %12.5E ylo yhi\n" % (np.min(atoms[:,coly]),np.max(atoms[:,coly])))
    outfile.write("%12.5E %12.5E zlo zhi\n" % (np.min(atoms[:,colz])-lz/2,np.max(atoms[:,colz])+lz/2))

    outfile.write("\nMasses\n\n")
    for i in range(ntypes):
        outfile.write("%5d\t%9.3E\n" % (i+1,1))
    '''
    outfile.write("\nAtoms # full\n\n")
    for i in range(natoms):
        outfile.write("%5d\t%d\t%d\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%d\t%d\t%d\n" % (atoms[i,dtype.get('id')],
                                                                          atoms[i,dtype.get('molecule')],
                                                                          atoms[i,dtype.get('type')],
                                                                          atoms[i,dtype.get('q')],
                                                                          atoms[i,dtype.get('x')],
                                                                          atoms[i,dtype.get('y')],
                                                                          atoms[i,dtype.get('z')],
                                                                          atoms[i,dtype.get('ix')],
                                                                          atoms[i,dtype.get('iy')],
                                                                          atoms[i,dtype.get('iz')]))
    '''
    outfile.write("\nAtoms # bond\n\n")
    for i in range(natoms):
            outfile.write("%5d\t%d\t%d\t%10.3f\t%10.3f\t%10.3f\t%d\t%d\t%d\n" % (atoms[i,dtype.get('id')],
                                                                          atoms[i,dtype.get('molecule')],
                                                                          atoms[i,dtype.get('type')],
                                                                          atoms[i,dtype.get('x')],
                                                                          atoms[i,dtype.get('y')],
                                                                          atoms[i,dtype.get('z')],
                                                                          atoms[i,dtype.get('ix')],
                                                                          atoms[i,dtype.get('iy')],
                                                                          atoms[i,dtype.get('iz')]))
