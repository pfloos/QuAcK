# Green's functions and response

For a given orbital reference, QuAcK can compute one-particle poles with Green's-function methods and two-particle excitation energies with response methods. See [method availability](generated/methods.md) before choosing a reference or a self-consistency scheme.

| Family | Approximation or observable | Settings to record |
|---|---|---|
| GF2 and GF3 | One-particle self-energy through second or third order, with additional variants in the drivers | Reference, self-consistency, regularization, linearization, and root-selection procedure |
| GW | Self-energy built from a Green's function and screened interaction | Reference, TDA for screening, regularization, and quasiparticle solver settings |
| T-matrix | Repeated particle–particle or electron–hole scattering | Channel, self-consistency, and screening or pair-space options |
| IP/EA-ADC and vertex corrections | One-particle poles with ADC or diagrammatic self-energies | Dyson or alternate formulation, diagonal approximation, and static self-energy options |
| phRPA and ppRPA | Neutral excitations and pair addition/removal poles, respectively | Channel, spin sector, TDA, and active orbital space |
| BSE | Neutral or charged two-particle response using a parent calculation | Parent self-energy, channel, static or dynamical kernel, and solver |
| Parquet | Coupled two-particle scattering channels and a one-particle self-energy | One- and two-body iteration counts, DIIS, channel TDA switches, and convergence |

## GF2 and GF3

GF2 uses a second-order self-energy. The `G0F2`, `evGF2`, `qsGF2`, and
`scGF2` flags select different levels of self-consistency. The third-order
routes are `G0F3`, `evGF3`, and `psdG0F3`. See the catalogue for driver coverage.

## Single-shot and self-consistent GW

| Input flag | Scheme | Drivers |
|---|---|---|
| `G0W0` | Single-shot GW from the input orbitals and energies | R, U, G |
| `evGW` | Eigenvalue self-consistency with fixed orbitals | R, U, G |
| `qsGW` | Quasiparticle self-consistency with an updated effective Hamiltonian and orbitals | R, U, G |
| `scGW` | Fully self-consistent GW on the imaginary axis | R and a dedicated Bogoliubov branch |

## T-matrix

The particle–particle and electron–hole T-matrix families sum repeated scattering
in their respective channels. The catalogue lists the single-shot, eigenvalue
and quasiparticle self-consistent paths, and the upfolded particle–particle path.

## ADC and vertex corrections

QuAcK includes IP/EA-ADC(2), IP/EA-ADC(3), SOSEX, 2SOSEX, and G3W2 variants.
The catalogue lists the ADC representations and identifies generalized
routes limited to the diagonal Dyson implementation.

## Parquet

`evParquet` and `qsParquet` have restricted and generalized paths.
[`GWParquet`](gw-parquet.md) is a single-shot workflow for canonical RHF.

## Experimental Faddeev

[`FTDA` and `FRPA`](faddeev.md) are standalone experimental implementations
for canonical RHF. The guide describes the pair kernels, couplings, and solver options.

## Response options

The screening TDA switch `TDA_W` and the excitation-space TDA switch `TDA` control different problems. Likewise, `phBSE` and `ppBSE` request different channels. Read the corresponding parent-method routine before using a combination for production results. Record the active orbital space and any frozen-core or CVS settings together with the parent-method options.

## Reading spectra

The quasiparticle equation can admit several roots. Save the pole energy and weight together, especially for inner-valence ionization. For an upfolded Hermitian problem, the weight comes from the squared norm in the one-particle block. A Davidson calculation targeting selected roots gives a subset of that spectrum.

## BSE availability

These entries cover the standard real-reference drivers. They indicate active calls, with dynamical corrections listed separately where checked.

| Parent | R | U | G |
|---|---|---|---|
| `G0W0` | ph and pp, both with optional dynamical correction | ph with optional dynamical correction | ph and pp. The ph dynamical call is inactive |
| `evGW`, `qsGW` | ph and pp, both with optional dynamical correction | ph with optional dynamical correction | ph. The ph dynamical call is inactive |
| `G0F2` | ph and pp | — | ph and pp |
| `evGF2` | ph and pp | — | ph |
| `qsGF2` | ph and pp | — | — |
| `G0T0pp` | ph and pp | — | ph and pp |
| `evGTpp`, `qsGTpp` | ph | — | — |

The [Fortran source index](generated/routines.md) links to the parent and response routines. The table covers the listed families and standard references. CAP, complex references, and alternate workflows require their own checks.

The `phBSE2` kernel switch is used inside the restricted GW particle–hole
BSE routines. Set `phBSE` to enter that response path.
