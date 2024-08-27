
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="numpy.core.getlimits")
import numpy as np
from pyscf import gto



def read_xyz(xyz):

    f = open(xyz, 'r')

    lines = f.read().splitlines()

    # nb of atoms
    nbAt = int(lines.pop(0))
    lines.pop(0)

    # list of atoms positions
    list_pos_atom = []
    for line in lines:
        tmp = line.split()
        atom = tmp[0]
        pos = (float(tmp[1]), float(tmp[2]), float(tmp[3]))
        list_pos_atom.append([atom,pos])

    f.close()

    return list_pos_atom


def mol_prop(xyz: str, input_basis: str, charge: int, multiplicity: int, unit: str, cartesian: bool):

    list_pos_atom = read_xyz(xyz)

    mol = gto.M(
        atom = list_pos_atom,
        basis = input_basis,
        charge = charge,
        spin = multiplicity - 1
    )

    mol.unit = unit
    mol.cart = cartesian

    mol.build()

    natoms = mol.natm       # nb of atoms
    nalpha = mol.nelec[0]   # nb of alpha-electrons
    nbeta = mol.nelec[1]    # nb of beta-electrons
    Enuc = mol.energy_nuc() # nuclear energy

    return natoms, nalpha, nbeta, Enuc


