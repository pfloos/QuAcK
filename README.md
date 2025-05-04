
# 🦆 QuAcK: Quantum Chemistry Prototyping Toolkit

![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)
![Fortran 90](https://img.shields.io/badge/language-Fortran%2090-yellow)
![Stars](https://img.shields.io/github/stars/pfloos/QuAcK?style=social)
![Forks](https://img.shields.io/github/forks/pfloos/QuAcK?style=social)

QuAcK is an open-source, lightweight electronic structure program written in **Fortran 90**, developed at the [Laboratoire de Chimie et Physique Quantiques (LCPQ)](https://www.lcpq.ups-tlse.fr/) in Toulouse, France. Designed primarily for rapid prototyping of new ideas in quantum chemistry, QuAcK provides a flexible environment for testing novel methods before integrating them into larger-scale projects like the [Quantum Package](https://github.com/QuantumPackage/qp2).

> ⚠️ **Note:** QuAcK is under active development. Users should be cautious and validate results, as the code may allow unconventional inputs to facilitate flexibility during prototyping.

---

## 🚀 Features

- **Rapid Prototyping:** Ideal for testing and developing new quantum chemistry methods.
- **Modular Design:** Easily integrate with other tools and libraries.
- **Educational Tool:** Serves as an excellent entry point for researchers familiar with electronic structure theory.
- **Integration with PySCF:** Utilizes [PySCF](https://github.com/pyscf/pyscf) for computing one- and two-electron integrals.

---

## 🛠️ Installation

1. **Clone the Repository:**

   ```bash
   git clone https://github.com/pfloos/QuAcK.git
   ```

2. **Set the `QUACK_ROOT` Environment Variable:**

   ```bash
   export QUACK_ROOT=$HOME/Work/QuAcK
   ```

3. **Install PySCF:**

   ```bash
   pip install pyscf
   ```

   *PySCF is used for computing one- and two-electron integrals, which are then read by QuAcK. It's also possible to use other software for integral computations.*

---

## ⚡ Quick Start

Navigate to the QuAcK directory and run the main script:

```bash
cd $QUACK_ROOT
python PyDuck.py -h
```

**Usage:**

```bash
usage: PyDuck.py [-h] -b BASIS [--bohr] [-c CHARGE] [--cartesian]
                 [--print_2e] [--formatted_2e] [--mmap_2e] [--aosym_2e]
                 [-fc FROZEN_CORE] [-m MULTIPLICITY]
                 [--working_dir WORKING_DIR] -x XYZ
```

**Options:**

- `-b, --basis BASIS`: Name of the basis set file in `$QUACK_ROOT/basis/`.
- `--bohr`: Specify if the XYZ file is in Bohr units (default is Angstrom).
- `-c, --charge CHARGE`: Total charge of the molecule (e.g., `m1` for -1).
- `-x, --xyz XYZ`: Path to the XYZ file containing molecular geometry.
- Additional options available; use `-h` for full list.

---

## 👥 Contributors

- [Pierre-François Loos](https://github.com/pfloos)
- [Anthony Scemama](https://github.com/scemama)
- [Enzo Monino](https://github.com/enzomonino)
- [Antoine Marie](https://github.com/antoine-marie)
- [Abdallah Ammar](https://scholar.google.com)
- [Mauricio Rodriguez-Mayorga](https://scholar.google.com)
- [Loris Burth](https://github.com/lorisburth)

---

## 📄 License

QuAcK is licensed under the [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html).

---

## 📫 Contact

For questions or contributions, please open an issue or submit a pull request on the [GitHub repository](https://github.com/pfloos/QuAcK).
