#!/usr/bin/env python

import os
import argparse
import pyscf
from pyscf import gto
import numpy as np
import subprocess

#Find the value of the environnement variable QUACK_ROOT. If not present we use the current repository
QuAcK_dir=os.environ.get('QUACK_ROOT','./')

#Create the argument parser object and gives a description of the script
parser = argparse.ArgumentParser(description='This script is the main script of QuAcK, it is used to run the calculation.\n If $QUACK_ROOT is not set, $QUACK_ROOT is replaces by the current directory.')

#Initialize all the options for the script
parser.add_argument('-b', '--basis', type=str, required=True, help='Name of the file containing the basis set in the $QUACK_ROOT/basis/ directory')
parser.add_argument('--bohr', default='Angstrom', action='store_const', const='Bohr', help='By default QuAcK assumes that the xyz files are in Angstrom. Add this argument if your xyz file is in Bohr.')
parser.add_argument('-c', '--charge', type=int, default=0, help='Total charge of the molecule. Specify negative charges with "m" instead of the minus sign, for example m1 instead of -1. Default is 0')
parser.add_argument('--cartesian', default=False, action='store_true', help='Add this option if you want to use cartesian basis functions.')
parser.add_argument('-fc', '--frozen_core', type=bool, default=False, help='Freeze core MOs. Default is false')
parser.add_argument('-m', '--multiplicity', type=int, default=1, help='Spin multiplicity. Default is 1 therefore singlet')
parser.add_argument('--working_dir', type=str, default=QuAcK_dir, help='Set a working directory to run the calculation.')
parser.add_argument('-x', '--xyz', type=str, required=True, help='Name of the file containing the nuclear coordinates in xyz format in the $QUACK_ROOT/mol/ directory without the .xyz extension')

#Parse the arguments
args = parser.parse_args()
input_basis=args.basis
unit=args.bohr
charge=args.charge
frozen_core=args.frozen_core
multiplicity=args.multiplicity
xyz=args.xyz + '.xyz'
cartesian=args.cartesian
working_dir=args.working_dir

#Read molecule
f = open(QuAcK_dir + '/mol/' + xyz,'r')
lines = f.read().splitlines()
nbAt = int(lines.pop(0))
lines.pop(0)
list_pos_atom = []
for line in lines:
    tmp = line.split()
    atom = tmp[0]
    pos = (float(tmp[1]),float(tmp[2]),float(tmp[3]))
    list_pos_atom.append([atom,pos])
f.close()

#Definition of the molecule
mol = gto.M(
    atom = list_pos_atom,
    basis = input_basis,
    charge = charge,
    spin = multiplicity - 1
)

#Fix the unit for the lengths
mol.unit=unit
#
mol.cart = cartesian

#Update mol object
mol.build()

#Accessing number of electrons
nelec=mol.nelec #Access the number of electrons
nalpha=nelec[0]
nbeta=nelec[1]

f = open(working_dir+'/input/molecule','w')
f.write('# nAt nEla nElb nCore nRyd\n')
f.write(str(mol.natm)+' '+str(nalpha)+' '+str(nbeta)+' '+str(0)+' '+str(0)+'\n')
f.write('# Znuc x  y  z\n')
for i in range(len(list_pos_atom)):
    f.write(list_pos_atom[i][0]+' '+str(list_pos_atom[i][1][0])+' '+str(list_pos_atom[i][1][1])+' '+str(list_pos_atom[i][1][2])+'\n')
f.close()

#Compute nuclear energy and put it in a file
subprocess.call(['rm', working_dir + '/int/ENuc.dat'])
f = open(working_dir+'/int/ENuc.dat','w')
f.write(mol.energy_nuc())
f.write(' ')
f.close()

#Compute 1e integrals
ovlp = mol.intor('int1e_ovlp')#Overlap matrix elements
v1e  = mol.intor('int1e_nuc') #Nuclear repulsion matrix elements
t1e  = mol.intor('int1e_kin') #Kinetic energy matrix elements
dipole = mol.intor('int1e_r') #Matrix elements of the x, y, z operators
x,y,z = dipole[0],dipole[1],dipole[2]

norb = len(ovlp)
subprocess.call(['rm', working_dir + '/int/nBas.dat'])
f = open(working_dir+'/int/nBas.dat','w')
f.write(str(norb))
f.write(' ')
f.close()


def write_matrix_to_file(matrix,size,file,cutoff=1e-15):
    f = open(file, 'a')
    for i in range(size):
        for j in range(i,size):
            if abs(matrix[i][j]) > cutoff:
                f.write(str(i+1)+' '+str(j+1)+' '+"{:.16E}".format(matrix[i][j]))
                f.write('\n')
    f.close()
    
#Write all 1 electron quantities in files
#Ov,Nuc,Kin,x,y,z
subprocess.call(['rm', working_dir + '/int/Ov.dat'])
write_matrix_to_file(ovlp,norb,working_dir+'/int/Ov.dat')
subprocess.call(['rm', working_dir + '/int/Nuc.dat'])
write_matrix_to_file(v1e,norb,working_dir+'/int/Nuc.dat')
subprocess.call(['rm', working_dir + '/int/Kin.dat'])
write_matrix_to_file(t1e,norb,working_dir+'/int/Kin.dat')
subprocess.call(['rm', working_dir + '/int/x.dat'])
write_matrix_to_file(x,norb,working_dir+'/int/x.dat')
subprocess.call(['rm', working_dir + '/int/y.dat'])
write_matrix_to_file(y,norb,working_dir+'/int/y.dat')
subprocess.call(['rm', working_dir + '/int/z.dat'])
write_matrix_to_file(z,norb,working_dir+'/int/z.dat')

#Write two-electron integrals
eri_ao = mol.intor('int2e')

def write_tensor_to_file(tensor,size,file,cutoff=1e-15):
    f = open(file, 'a')
    for i in range(size):
        for j in range(i,size):
            for k in range(i,size):
                for l in range(j,size):
                    if abs(tensor[i][k][j][l]) > cutoff:
                        f.write(str(i+1)+' '+str(j+1)+' '+str(k+1)+' '+str(l+1)+' '+"{:.16E}".format(tensor[i][k][j][l]))
                        f.write('\n')
    f.close()

subprocess.call(['rm', working_dir + '/int/ERI.dat'])
write_tensor_to_file(eri_ao,norb,working_dir+'/int/ERI.dat')

#Execute the QuAcK fortran program
subprocess.call(QuAcK_dir+'/bin/QuAcK')
