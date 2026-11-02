# TDA Faddeev-Parquet with ADC(2) or ADC(3) external couplings

Both restricted (`R_TDA_Faddeev_self_energy`) and generalized
(`G_TDA_Faddeev_self_energy`) routines are available in `src/Parquet`.
They implement the complete-TDA physical reduction in
`faddeev_parquet_notes.tex`, using one physical 2h1p space and one physical
2p1h space. The external vertices are the bare-integral ADC(3)
`U = U^(1) + U^(2)` at the current one-particle energies by default.
An optional switch selects ADC(2) endpoints, `U = U^(1)`. The one-particle
block contains HF orbital energies. No gauge parameters are introduced.

## Run

Build QuAcK normally. In the calculation's `input` directory use:

- `tests/inp/methods.RTDAFaddeev` as `methods` for RHF, or
  `tests/inp/methods.GTDAFaddeev` for GHF;
- `tests/inp/options.TDAFaddeev` as `options`.

The self-energy is selected by the hard-coded logical switches near the
start of `src/Parquet/R_evParquet.f90` and `src/Parquet/G_evParquet.f90`,
just like the previous ADC-Parquet switch:

```fortran
logical :: do_ADC_self_energy = .false.
logical :: do_Faddeev_self_energy = .true.
logical :: do_Faddeev_ADC2 = .false.
logical :: do_principal_IP = .true.
```

`do_Faddeev_ADC2=.true.` selects ADC(2) external coupling blocks in both
the full-spectrum and matrix-free Davidson solvers, for both R and G.
The default `.false.` retains ADC(3). The switch changes only the external
vertices: the HF one-particle block and dressed Faddeev charged blocks stay
the same. It therefore does not turn the entire model into ordinary ADC(2).
ADC(2) endpoints require no MP2 amplitudes or denominators; endpoint SRG
regularization has no effect on them.

The matrix-free principal-IP mode is enabled by default in both drivers.
Set `do_principal_IP=.false.` for the full spectrum. Rebuild QuAcK after
changing these switches, then launch PyDuck or QuAcK normally.

Both `TDAeh` and `TDApp` must be `T`; the driver rejects other values.
The method uses the existing `evParquet` switch and its pair iterations,
convergence thresholds, and DIIS settings. The supplied example requests
one one-body evaluation and up to 64 pair iterations. With `max_it_2b=1`
at the HF reference the kernels are bare, giving the ordinary ADC(3)
limit. A dressed-parquet test requires further pair iterations and checking
the reported pair convergence. Fixed unconverged snapshots can differ
between R and G when DIIS uses different residual metrics; compare fixed
snapshots with `max_diis_2b=1`.

`lin_parquet` is unused in either Faddeev solver. Principal-IP mode requires
`max_it_1b=1`: the driver returns the IP without feeding a partial orbital
energy array into one-body self-consistency. Pair iterations are unchanged.
In full-spectrum mode, more than one one-body iteration feeds all dominant
poles into the existing eigenvalue-self-consistent parquet loop.

To select the previous ADC-Parquet self-energy, set
`do_ADC_self_energy=.true.` and `do_Faddeev_self_energy=.false.`.
Set both to `.false.` for the native parquet self-energy. Setting both to
`.true.` is rejected. These switches apply to the R/G `evParquet` drivers,
not `qsParquet` or `R_GW_Parquet`.

## Matrix-free principal IP

`TDA_Faddeev_IP.f90` applies the same physical upfolded matrix using pair
contractions. Both 2h1p and 2p1h branches and the full one-particle block
are retained. The charged matrix is never built. Endpoint vertices are
cached; in ADC(3) mode their MP2 denominators are evaluated once. The R solver passes
spatial ERIs directly and applies spin deltas without constructing a
four-index spin-orbital ERI tensor.

The block harmonic Davidson solver starts from the entire active HF HOMO
manifold and up to two neighboring occupied orbitals. Candidate selection
uses HOMO-subspace weight, with previous-vector overlap breaking near ties,
and excludes poles above the HF gap midpoint. Several candidates are kept
through a thick restart. Defaults in `davidson` are a 96-vector subspace,
256 iterations, and a residual tolerance of `1e-8` Ha; the diagonal
preconditioner uses a `0.05` Ha denominator floor. The final residual is
checked with a fresh matrix-vector product. Failure is reported explicitly;
there is no automatic dense fallback.

Output contains the IP, negative-energy pole, total one-particle weight Z,
HOMO-subspace weight, residual, iteration/product counts, and separate setup
and solve timings. This targets the principal HOMO-connected ionization.
It does not search every cation state: if orbital ordering changes or HOMO
character fragments strongly, compare the full spectrum. A degenerate or
split HOMO manifold is treated as a block; the smallest IP among its
converged selected roots is reported.

At fixed subspace and orbital counts, solver storage grows linearly with
the charged dimension, rather than quadratically. Pair kernels and cached
external vertices still consume memory; the existing parquet pair solves
and four-index reducible vertices are unchanged.

## Matrix and conventions

The pair interactions are captured from the actual TDA matrices used in
the last pair solve, before updating the reducible vertices:

```
Weh = Aeh - diag(e_a-e_i)
Whh = Dpp + diag(e_i+e_j)
Wee = Cpp - diag(e_a+e_b)
```

If a pp chemical-potential reference shift is active, it is removed in
extracting these interactions, and the pp mode energies are restored to
absolute energies before parquet feedback. The R routine receives separate
singlet/triplet interactions with QuAcK's normalized spatial pair bases.

`G_TDA_Faddeev_build` forms the physical matrices and selected ADC(2)/ADC(3) couplings;
`G_TDA_Faddeev_spectrum` diagonalizes the resulting matrix

```
H = [ fHF   U-^T   U+^T ]
    [ U-    Mbar-    0  ]
    [ U+      0    Mbar+]
```

This realizes `Sigma(omega) = U-^T (omega I-Mbar-)^(-1) U-
+ U+^T (omega I-Mbar+)^(-1) U+`. The projected charged matrices are
symmetric for symmetric static pair interactions. Asymmetric or nonfinite
inputs are rejected instead of silently symmetrized.

The R implementation lifts its singlet/triplet pair kernels into spin
orbitals and constructs only the alpha charged determinant block
(Ms=1/2). It retains dark quartet states, which do not couple to the
one-particle space for spin-invariant inputs. This exact representation
shares the G endpoint implementation; it does not average screened
singlet and triplet kernels. Spatial external orbitals and their weights
are returned once. The full-spectrum solver uses dense matrices; the
principal-IP solver uses this same spin representation with pair contractions.

Frozen occupied orbitals and removed virtual orbitals are excluded from
all internal pairs and external correlations. Output includes significant
poles, the active spectral-weight sum, and a dominant pole and residue for
each active orbital. Exact degeneracies are grouped when choosing a
dominant pole.

## Regularization

`reg_PA=T` applies the existing one-body SRG factor
`(1-exp(-2*reg_1b*Delta^2))/Delta` to the explicit ADC(3) second-order
endpoint denominators. Its zero-denominator limit is evaluated smoothly.
The free charged propagation and supplied pair interactions are untouched
by this endpoint regularization. `reg_PA=F` uses ordinary inverse
denominators and rejects singular ones.

The existing parent pair iteration continues to use `reg_2b` in constructing
its reducible vertices. As in that driver, `reg_PA` is not an on/off switch
for the parent's two-body regulator. Large positive SRG flows approach the
unregularized result.

## Scope and checks

This implements the static complete-TDA Faddeev model with selectable
ADC(2) or ADC(3) endpoints. Identity with the full parquet Schwinger-Dyson self-energy has
not been established. It is a separate self-energy approximation to test.

The older standalone `G_TDA_Faddeev_Parquet.f90` prototype is not called by
this implementation. It stopped at the redundant auxiliary eigensystem;
the physical reduction avoids its extra sectors and normalization issues.

Algebraic regression checks cover the matrix and external couplings
against QuAcK's existing `G_IPEA_ADC3`, R/G spin reconstruction, active
spaces and empty triplet sectors, the exact Faddeev intertwining identity
for both branches with dense dressed kernels, and SRG limits. They also
compare matrix-free products and diagonals with the dense construction,
principal roots and residues with dense eigenpairs, split HOMO manifolds,
forced thick restart, and explicit nonconvergence. ADC(2) checks verify
bare endpoint elements, unchanged charged blocks, zero-denominator safety,
and R/G full-spectrum and Davidson agreement:

```sh
python tests/check_tda_faddeev.py --workdir /tmp/new-algebraic-checks
```

Build QuAcK first; the check requires NumPy and gfortran. Set `FC` and
`LAPACK_LIBS` if necessary. The full-spectrum molecular test script requires
PySCF and a build with `do_principal_IP=.false.` in both drivers. It runs
H2O/6-31G and H2/STO-3G, with logs retained:

```sh
python tests/check_tda_faddeev_molecules.py --workdir /tmp/new-molecular-checks
```
