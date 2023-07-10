# QuAcK: a software for emerging quantum electronic structure methods

**Contributors:**
- [Pierre-Francois Loos](https://pfloos.github.io/WEB_LOOS)
- [Enzo Monino](https://enzomonino.github.io)
- [Antoine Marie](https://antoine-marie.github.io)
- [Anthony Scemama](https://scemama.github.io)
  
# What is it?



# Installation guide
The QuAcK software can be downloaded on GitHub as a Git repository
```
git clone https://github.com/pfloos/QuAcK.git
```

Then, one must define the variable `QUACK_ROOT`. For example, 
```
export QUACK_ROOT=$HOME/Work/QuAcK
```
You must also install [PySCF](https://pyscf.org) (for example using `pip`)
```
pip install pyscf
```

PySCF is used for the computation of one- and two-electron integrals (mainly).

# Quick start

```
QuAcK 💩 % cd $QUACK_ROOT
QuAcK 💩 % python PyDuck.py -h
usage: PyDuck.py [-h] -b BASIS [--bohr] [-c CHARGE] [--cartesian] [-fc FROZEN_CORE] [-m MULTIPLICITY] [--working_dir WORKING_DIR] -x XYZ

This script is the main script of QuAcK, it is used to run the calculation. If $QUACK_ROOT is not set, $QUACK_ROOT is replaces by the current
directory.

options:
  -h, --help            show this help message and exit
  -b BASIS, --basis BASIS
                        Name of the file containing the basis set in the $QUACK_ROOT/basis/ directory
  --bohr                By default QuAcK assumes that the xyz files are in Angstrom. Add this argument if your xyz file is in Bohr.
  -c CHARGE, --charge CHARGE
                        Total charge of the molecule. Specify negative charges with "m" instead of the minus sign, for example m1 instead of -1.
                        Default is 0
  --cartesian           Add this option if you want to use cartesian basis functions.
  -fc FROZEN_CORE, --frozen_core FROZEN_CORE
                        Freeze core MOs. Default is false
  -m MULTIPLICITY, --multiplicity MULTIPLICITY
                        Number of unpaired electrons 2S. Default is 0 therefore singlet
  --working_dir WORKING_DIR
                        Set a working directory to run the calculation.
  -x XYZ, --xyz XYZ     Name of the file containing the nuclear coordinates in xyz format in the $QUACK_ROOT/mol/ directory without the .xyz
                        extension
'''
