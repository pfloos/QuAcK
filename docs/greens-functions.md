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
| `osCOT` | One-shot correlated orbital theory based on the GW supermatrix | R |
| `scCOT` | Quasiparticle self-consistent correlated orbital theory | R |
| `STEOM-eB` | Static second-order electron–boson neutral excitation energies | R |

Enable `RHF` and put the COT methods row above the test row in `input/methods`.

```text
# osCOT scCOT
  F     T
```

The legacy `# COT` flag selects `osCOT`. The shared solver settings and
self-consistency controls in `input/options` are

```text
# COT: Davidson TDA_COT Davidson_threshold maxSCF thresh DIIS root_width SRG flow
       T        F       1.0e-10            64     1.0e-6 5    0.0        F   500.0
```

`osCOT` builds the screening from the RHF orbitals and energies and applies
one Hermitian COT update. `scCOT` starts from RHF and rebuilds the density,
full Hartree–Fock matrix, MO integrals, and screening at each iteration.
The screening and auxiliary denominators use the current COT eigenvalues.
The orbital update diagonalizes the Hermitian COT Hamiltonian in the active
space, with frozen and removed orbitals held fixed. `thresh` is the maximum
absolute element of the Hamiltonian fixed-point residual in the initial
orthonormal RHF orbital basis, in hartree. `DIIS` sets the history length,
and 0 or 1 disables extrapolation. `maxSCF` limits the number of updates.
The optional `root_width` (hartree) favors poles close to each input
quasiparticle energy. It multiplies the orbital weight by
\(1/[1+((\omega-\varepsilon)/\mathrm{root\_width})^2]\) before the usual
one-to-one assignment. The default of 0 uses unweighted orbital character.
A positive width can change the selected poles. `osCOT` always uses
unweighted orbital character.

`SRG = T` enables complementary dynamical SRG damping for both COT
methods. The physical Fock block and auxiliary energies are unchanged.
Each coupling is multiplied by
\(\sqrt{1-\exp[-2t(\varepsilon_p-D_\alpha)^2]}\), using the reference
orbital energies and auxiliary poles. This preserves a Hermitian
supermatrix and reproduces the diagonal SRG-GW self-energy at the
reference energies. The numerator is frozen in frequency and the
off-diagonal residues use the product of the two coupling factors.
`flow` is in inverse squared hartree. Zero flow with SRG enabled
removes the dynamical couplings. Set `SRG = F` for the bare GW supermatrix.

`Davidson = T` uses the current reference orbitals as seeds to follow one
pole per active occupied and virtual orbital without storing the supermatrix.
`Davidson = F` uses dense diagonalization, with a dimension limit of 4000.
`TDA_COT` selects TDA screening for COT independently of GW. Older options
files without the shared COT row default to Davidson and use `TDA_W` for
COT screening. `Davidson_threshold` is the eigenpair residual tolerance in
hartree, defaulting to `1.0e-10`. Older files with a separate `# scCOT:`
row are also accepted.

The `COT total energy` is the RHF energy functional evaluated on the
determinant formed by filling the lowest eigenvectors of the Hermitian COT
active block and keeping the frozen core doubly occupied. QuAcK reports the
kinetic, electron–nuclear, Hartree, exchange, and nuclear-repulsion terms.

COT also reports two full direct-RPA correlation energies for the same active
space. `dRPA@RHF` uses the RHF reference for `osCOT`. The corresponding
`dRPA@scCOT input` uses the reference of the final self-consistent update.
`dRPA@COT` uses
the Hermitian COT orbitals, their eigenvalues, and the integrals transformed
to that orbital basis. `TDA_COT` changes the screening used to construct COT.
Both printed correlation energies use the full direct-RPA equations.

### STEOM-eB neutral excitation energies

Enable `RHF` and set the STEOM-eB flag in `input/methods` alongside the COT row.

```text
# STEOM-eB
  T
```

Select the response approximation and spin manifolds in `input/options`.

```text
# STEOM-eB: TDA singlet triplet Davidson n_roots threshold
            F   T       T       F         5       1.0e-8
```

The default is full response for both singlets and triplets. Set `TDA = T`
to diagonalize the resonant block alone. To select one spin manifold, set
the other flag to `F`. The existing `TDA_W` option controls the screening
independently of the response approximation.

Set `Davidson = T` to compute the lowest `n_roots` in each selected spin
manifold, up to the size of the active singles space. The threshold is the
Euclidean Ritz residual norm in hartree for TDA and in hartree squared for
full response. The full-response solver applies the Hermitian operator
\(L^{\mathsf T}(A+B)L\) using triangular multiplications, where
\(A-B=LL^{\mathsf T}\). Matrix assembly, screening and the Cholesky
factorization still use dense routines. The original three-option row
selects the dense solver.

STEOM-eB uses the canonical RHF orbitals, singlet dRPA screening, and the
static second-order electron–boson Hamiltonian. It can run independently
or in the same job as COT. The response starts from the TDHF blocks,
including the bare Hartree and exchange kernel. It adds the symmetrized
on-shell GW one-body correction and the induced electron–hole terms from
the eB transformation. With the electron–boson coupling set to zero, the
response reduces to TDHF, or to CIS in TDA.
The static approximation requires nonresonant denominators, and full
response requires a stable reference matrix. QuAcK reports the affected
orbitals and screening mode when it encounters a resonance.

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
