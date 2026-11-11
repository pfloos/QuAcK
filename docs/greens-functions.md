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

## Single-shot and self-consistent GW

`G0W0` uses fixed input orbitals and energies to form the starting propagator and screening. `evGW` updates quasiparticle energies while keeping the orbitals fixed. `qsGW` updates the effective one-particle Hamiltonian and its orbitals. The restricted `scGW` route uses an imaginary-axis self-consistent calculation. Compare like quantities across these schemes. For example, the pole weight associated with a frequency-dependent self-energy differs from the weight of an eigenvector of a static quasiparticle Hamiltonian.

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
