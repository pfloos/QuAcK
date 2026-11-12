# Calculation inputs

The [generated input reference](generated/inputs.md) shows the exact templates at the documented revision. Copy the templates into a dedicated calculation directory so each log has a matching set of inputs.

| Path in the calculation directory | Who writes it | Contents |
|---|---|---|
| `mol/NAME.xyz` | User or example setup | Geometry with an atom count and comment line. PyDuck defaults to Å |
| `basis/BASIS` | User or example setup | NWChem-format basis file. An explicit file path is also accepted |
| `input/methods` | User or example setup | Reference, method, and test flags |
| `input/options` | User or example setup | Numerical controls and method-specific choices |
| `input/molecule` | PyDuck | Electron counts and molecular data for the Fortran driver |
| `int/` | PyDuck | Nuclear repulsion, one-electron integrals, and electron-repulsion integrals |

The Fortran readers interpret the main input records by position. Keep the template order and retain the header lines. Some optional records are recognized by their headers. The methods template lists standalone `FTDA` and `FRPA` above the test flags. Older files can keep that row below the test flags or omit it. A header can describe a flag even when a particular driver has no solver for it.

## Common choices

| Setting | What to record |
|---|---|
| Charge and multiplicity | PyDuck `--charge` and `--multiplicity`. Use `--charge=-1` for a singly charged anion |
| Coordinate units | Å by default. Pass `--bohr` for bohr |
| SCF convergence | Maximum iterations, threshold, DIIS size, initial guess, and any level shift or orbital mixing |
| Response sector | TDA and the requested spin sectors. The restricted driver uses singlet/triplet labels. The unrestricted driver uses spin-conserving/spin-flip labels |
| Integral representation | Spherical or Cartesian basis, and formatted or binary two-electron integrals |
| Quasiparticle solver | Linearization, regularization, broadening, and any one-body self-consistency |
| Active orbital space | Frozen orbitals, CVS, and orbital occupations, as applicable to the method |

The shipped value is the one used when you copy the template. A reader fallback applies only where the code allows an absent record. For example, an old two-switch Parquet methods record leaves `GWParquet` disabled. The optional standalone Faddeev record has its own fallback values, shown in the [Faddeev guide](faddeev.md). Driver-specific restrictions can further limit the combinations that run.

## Reproduce a calculation

Save the source revision, geometry and its units, basis file, charge and multiplicity, input files, and output. For method development, also record the compiler, linked BLAS/LAPACK library, and whether the run used the standard or HPC driver. Use a separate working directory for every concurrent PyDuck run because it writes and replaces integral files.
