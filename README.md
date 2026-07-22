
# 🦆 QuAcK: an open-source software for emerging quantum electronic structure methods

![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)
![Fortran 90](https://img.shields.io/badge/language-Fortran%2090-yellow)
![Stars](https://img.shields.io/github/stars/pfloos/QuAcK?style=social)
![Forks](https://img.shields.io/github/forks/pfloos/QuAcK?style=social)

---

<img src="logo/logo_quack.png"  width="250">

---

## 👥 Contributors

- [Pierre-Francois Loos](https://pfloos.github.io/WEB_LOOS)
- [Anthony Scemama](https://scemama.github.io)
- [Enzo Monino](https://enzomonino.github.io)
- [Antoine Marie](https://antoine-marie.github.io)
- [Abdallah Ammar](https://scholar.google.com/citations?user=y437T5sAAAAJ&hl=en)
- [Mauricio Rodriguez-Mayorga](https://scholar.google.com/citations?user=OLGOgQgAAAAJ&hl=es)
- [Loris Burth](https://github.com/lburth)

---

## 🚀 Features

- **Rapid Prototyping:** Ideal for testing and developing new quantum chemistry methods.
- **Modular Design:** Easily integrate with other tools and libraries.
- **Educational Tool:** Serves as an excellent entry point for researchers familiar with electronic structure theory.
- **Integration with PySCF:** Utilizes [PySCF](https://github.com/pyscf/pyscf) for computing one- and two-electron integrals.

---

## What is it?

**QuAcK** is a lightweight electronic structure program written in `Fortran 90`, developed at the [Laboratoire de Chimie et Physique Quantiques (LCPQ)](https://www.lcpq.ups-tlse.fr) in Toulouse, France. Originally designed as a platform for rapid prototyping of new ideas in quantum chemistry, QuAcK serves as a flexible and accessible environment for testing novel methods before integrating them more efficiently into larger-scale projects such as the [Quantum Package](https://quantumpackage.github.io/qp2/).

Thanks to its compact and transparent codebase, QuAcK is particularly well-suited for experienced PhD students and postdoctoral researchers who are already familiar with electronic structure theory and want to quickly implement or explore new concepts. Written in a clean and relatively straightforward programming language, it provides an excellent entry point for those looking to dive into method development.

For beginners in the field or those with less programming experience, we recommend starting with [qcmath](https://github.com/LCPQ/qcmath/), a symbolic and numerical quantum chemistry toolkit built in [Mathematica](https://www.wolfram.com/mathematica/). qcmath is specifically designed to help newcomers explore and develop ideas without the complexity of full-fledged numerical implementations.

> ⚠️ **Note:** QuAcK is under active and ongoing development, which means that bugs, inconsistencies, and incomplete features are to be expected. It is a tool made *by* experts *for* experts—users are expected to understand what they’re doing and to remain cautious when interpreting results. The code may allow questionable inputs or behavior *on purpose*, to encourage flexibility during prototyping—so always double-check your results and assumptions. In short: use QuAcK at your own risk—but also to your advantage, if you're ready to experiment and explore.

---

## 🛠️ Installation

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

Then, go to the `src` directory and compile
```
cd src; make
```

## ⚡ Quick Start

```
~ 💩 % cd $QUACK_ROOT
QuAcK 💩 % python PyDuck.py -h
<!-- BEGIN PYDUCK_HELP -->
Module pyopencap is not installed.
usage: PyDuck.py [-h] [--working_dir WORKING_DIR] -b BASIS -x XYZ [--bohr]
                 [--cartesian] [-c CHARGE] [-m MULTIPLICITY] [--print_2e]
                 [--formatted_2e] [--mmap_2e] [--aosym_2e] [-nc] [-dm]

This script is the main script of QuAcK, it is used to run the calculation. If
$QUACK_ROOT is not set, $QUACK_ROOT is replaced by the current directory.

options:
  -h, --help            show this help message and exit

Input files:
  --working_dir WORKING_DIR
                        Set a working directory to run the calculation.
                        Default is $QUACK_ROOT, if set, otherwise the current
                        working directory (./)
  -b, --basis BASIS     Name or path of the file containing the basis set
                        information. If the argument is not a filepath QuAcK
                        searches in WORKING_DIR/basis for a file matching the
                        argument.
  -x, --xyz XYZ         Name of the file containing the nuclear coordinates in
                        xyz format in the WORKING_DIR/mol/ directory without
                        the .xyz extension.

Input format:
  --bohr                By default QuAcK assumes that the xyz files are in
                        Angstrom. Add this argument if your xyz file is in
                        Bohr.
  --cartesian           Add this option if you want to use cartesian basis
                        functions.

Molecule:
  -c, --charge CHARGE   Total charge of the molecule. Specify negative charges
                        with "m" instead of the minus sign, for example m1
                        instead of -1. Default is 0.
  -m, --multiplicity MULTIPLICITY
                        Spin multiplicity. Default is 1 (singlet).

Electron integrals:
  --print_2e            If True, print ERIs to disk.
  --formatted_2e        Add this option if you want to print formatted ERIs.
  --mmap_2e             If True, avoid using DRAM when generating ERIs.
  --aosym_2e            If True, use 8-fold symmetry in ERIs.
  -nc, --no_cap         If true, no CAP integrals are calculated and stored.If
                        the python module pyopencap is not available, this is
                        set automatically true.

Molecular orbitals:
  -dm, --dump_molden    Dump a molden file with the molecular orbitals. If
                        this is true WriteMOs in the QuAcK options is
                        automatically set to true an the molecular orbitals
                        are dumbed.
<!-- END PYDUCK_HELP -->
```

The two most important files are:
- `$QUACK_ROOT/input/methods` that gathers the methods you want to use.
- `$QUACK_ROOT/input/options` that gathers the different options associated these methods.

Copy the files `methods.default` and `options.default` to `methods` and `options`, respectively.
```
cp $QUACK_ROOT/input/methods.default $QUACK_ROOT/input/methods
cp $QUACK_ROOT/input/options.default $QUACK_ROOT/input/options
```
You can then edit these files to run the methods you'd like (by replacing `F` with `T`) with specific options.
These files look like this
```
QuAcK 💩 % cat input/methods 
<!-- BEGIN methods -->
# RHF UHF GHF ROHF HFB cRHF cUHF eRHF scGHF
  F   F   F   F    F   F    F	 F    F
# MP2 MP3
  F   F
# CCD pCCD DCD CCSD CCSD(T)
  F   F    F   F    F
# drCCD rCCD crCCD lCCD
  F     F    F     F
# CIS CIS(D) CID CISD FCI
  F   F      F   F    F
# phRPA phRPAx crRPA ppRPA BRPA
  F     F      F     F     F   
# OORPA
  F
# G0F2 evGF2 qsGF2 G0F3 evGF3 psdG0F3 scGF2
  F    F     F     F    F     F     F
# G0W0 evGW qsGW scGW
  F    F    F    F
# G0T0pp evGTpp qsGTpp ufG0T0pp
  F      F        F       F
# G0T0eh evGTeh qsGTeh
  F      F      F
# evParquet qsParquet
  F         F
# IPEA-ADC2 IPEA-ADC3 SOSEX 2SOSEX G3W2 psdG3W2
  F         F         F     F      F    F
# ADC-GW ADC-2SOSEX ADC(3)-G3W2 ADC(3x)-G3W2 ADC-G3W2
  F      F          F           F            F
# Rtest Utest Gtest
  F     F     F
<!-- END methods -->
```
and
```
QuAcK 💩 % cat input/options 
<!-- BEGIN options -->
# HF: maxSCF thresh  DIIS guess mix shift stab search aordm readFCIDUMP MOM WriteMOs
      256    0.0000001 5    1     0.0 0.0   F    F      F     F		F   F
# MP: reg
      F
# CC: maxSCF thresh  DIIS
      64     0.00001 5
# spin: TDA singlet triplet
        F   T       T
# RPA/MP2: CVS-Alpha CVS-Beta  FC-Alpha FC-Beta
               0         0          0       0
# OORPA:  maxIter thresh  dRPA state diagHess
          256     0.00001 T    0     T    
# GF: maxSCF thresh  DIIS lin eta renorm reg linDM restart_scGF2 verbose_scGF2
      256    0.00001 5    F   0.0 0      F   F     F		 F
# GW: maxSCF thresh  DIIS lin eta TDA_W reg nfreqs shift linDM restart_scGW verbose_scGW adjust_mu_scGW change_sign_XoB
      256    0.00001 5    F   0.0 F     F   1      0.0   F     F            F            F              F
# GT: maxSCF thresh  DIIS lin eta TDA_T reg linDM
      256    0.00001 5    F   0.0 F     F   F
# ACFDT: AC Kx  XBS
         F  F   F
# BSE: phBSE phBSE2 ppBSE dBSE dTDA
       F     F      F     F    F
# HFB: temperature  sigma chem_pot_HF restart_HFB error_P
         0.05       1.0   F           F		  F
# Parquet: TDAeh TDApp max_it_1b conv_1b max_it_2b conv_2b max_diis_1b max_diis_2b lin_parquet reg_1b reg_2b reg_PA
            T     T        1      0.01     1         0.01      1            1          F       1000.0 1000.0 T
# ADC: dyson diag_approx sig_inf lin reg eta
       T     F           F       F   F   0.0
# Ensemble: weight_N  N+1[for N-1 use F]
	    0.00000   T
<!-- END options -->
```

For example, if you want to run a calculation on water using the cc-pvdz basis set:
```
QuAcK 💩 % python PyDuck.py -x water -b cc-pvdz
```

QuAcK runs calculations in the `QUACK_ROOT` directory which is quite unusual but it also use the `--working_dir` option to run calculations elsewhere.

---

## 📄 License

QuAcK is licensed under the [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html).

---

## 📫 Contact

For questions or contributions, please open an issue or submit a pull request on the [GitHub repository](https://github.com/pfloos/QuAcK).

--- 

## 💰 Funding

<img src="https://lcpq.github.io/PTEROSOR/img/ERC.png" width="200" />

QuAcK is supported by the [PTEROSOR](https://lcpq.github.io/PTEROSOR/) project that has received funding from the European Research Council (ERC) under the European Union’s Horizon 2020 research and innovation programme (Grant agreement No. 863481).

---
