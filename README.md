# QuAcK: a open-source software for emerging quantum electronic structure methods

<img src="logo/logo_quack.png"  width="250">

**Contributors:**
- [Pierre-Francois Loos](https://pfloos.github.io/WEB_LOOS)
- [Enzo Monino](https://enzomonino.github.io)
- [Antoine Marie](https://antoine-marie.github.io)
- [Abdallah Ammar](https://scholar.google.com/citations?user=y437T5sAAAAJ&hl=en)
- [Anthony Scemama](https://scemama.github.io)
  
# What is it?

QuAcK is a small electronic structure program written in `Fortran 90` and developed at the Laboratoire de Chimie et Physique Quantiques [LCPQ](https://www.lcpq.ups-tlse.fr) (Toulouse, France).
QuAcK is usually used for prototyping purposes and the successful ideas are usually implemented more efficiently in [Quantum Package](https://quantumpackage.github.io/qp2/). QuAcK is an excellent place to start for experienced PhD students or postdocs as the code is simple and written with a fairly well-known and straightforward language. For beginners, we suggest having a look at [qcmath](https://github.com/LCPQ/qcmath/), a [Mathematica](https://www.wolfram.com/mathematica/)-based program to help newcomers in quantum chemistry easily develop their ideas. 

QuAcK is under continuous and active development, so it is very (very) likely to contain many bugs and errors. QuAcK is a code for experts, which means that you must know what you're doing and you have to make sure you're not doing anything silly (QuAcK may allow silly things to happen on purpose!). You have been warned.

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

PySCF is used for the computation of one- and two-electron integrals (mainly) which are dumped in files and read by QuAcK.
Therefore, it is very easy to use other software to compute the integrals or to add other types of integrals.

# Quick start

```
~ 💩 % cd $QUACK_ROOT
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
```

The two most important files are:
- `$QUACK_ROOT/input/methods` that gathers the methods you want to use.
- `$QUACK_ROOT/input/options` that gathers the different options associated these methods.

These files look like this
```
QuAcK 💩 % cat input/methods 
# RHF UHF RMOM UMOM KS
  T   F   F    F    F  
# MP2* MP3 
  F   F   
# CCD pCCD DCD CCSD CCSD(T) 
  F   F    F   F    F
# drCCD rCCD crCCD lCCD
  F     F    F     F
# CIS* CIS(D) CID CISD FCI
  F    F      F   F    F
# phRPA* phRPAx* crRPA ppRPA 
  F      F       F     F 
# G0F2* evGF2* qsGF2* G0F3 evGF3
  F     F      F      F    F
# G0W0* evGW* qsGW* SRG-qsGW ufG0W0 ufGW
  T     F     F     F        F      F
# G0T0pp* evGTpp* qsGTpp* G0T0eh evGTeh qsGTeh
  F       F       F       F      F      F
# * unrestricted version available
```
and
```
QuAcK 💩 % cat input/options 
# HF: maxSCF thresh   DIIS n_diis guess_type ortho_type mix_guess level_shift stability
      512    0.0000001  T    5     1          1          F         0.0         F
# MP: reg
      F
# CC: maxSCF thresh   DIIS n_diis
      64     0.0000001  T    5
# spin: TDA singlet triplet spin_conserved spin_flip 
        F   T       F       T              T 
# GF: maxSCF thresh  DIIS n_diis lin eta renorm reg
      256    0.00001 T    5      T   0.0 0      F
# GW: maxSCF thresh  DIIS n_diis lin eta COHSEX TDA_W reg
      256     0.00001  T  5      T    0.0  F    F     F 
# GT: maxSCF thresh  DIIS n_diis lin eta TDA_T reg
      256    0.00001  T    5      T   0.1 F     F  
# ACFDT: AC Kx  XBS
         F  T   T
# BSE: BSE dBSE dTDA evDyn ppBSE BSE2
        T    T    T    F   F     F
```

For example, if you want to run a calculation on water using the cc-pvdz basis set:
```
QuAcK 💩 % python PyDuck.py -x water -b cc-pvdz
```

QuAcK runs calculations in the `QUACK_ROOT` directory which is quite unusual but it can be easily modified to run calculations elsewhere.
You just have to make sure that QuAcK reads/writes the integrals and molecular information at the right spot.

<img src="https://lcpq.github.io/PTEROSOR/img/ERC.png" width="200" />

QuAcK is supported by the [PTEROSOR](https://lcpq.github.io/PTEROSOR/) project that has received funding from the European Research Council (ERC) under the European Union’s Horizon 2020 research and innovation programme (Grant agreement No. 863481).
