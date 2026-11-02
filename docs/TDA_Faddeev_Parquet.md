# TDA Faddeev-Parquet with ADC(3) external couplings

Both restricted (`R_TDA_Faddeev_self_energy`) and generalized
(`G_TDA_Faddeev_self_energy`) routines are available in `src/Parquet`.
They implement the complete-TDA physical reduction in
`faddeev_parquet_notes.tex`, using one physical 2h1p space and one physical
2p1h space. The external vertices are the bare-integral ADC(3)
`U = U^(1) + U^(2)` at the current one-particle energies. The one-particle
block contains HF orbital energies. No gauge parameters are introduced.

## Run

Build QuAcK normally. In the calculation's `input` directory use:

- `tests/inp/methods.RTDAFaddeev` as `methods` for RHF, or
  `tests/inp/methods.GTDAFaddeev` for GHF;
- `tests/inp/options.TDAFaddeev` as `options`.

Enable the new self-energy when launching QuAcK or PyDuck:

```sh
QUACK_PARQUET_SELF_ENERGY=TDA_FADDEEV python PyDuck.py [your usual arguments]
```

For an existing integral directory:

```sh
QUACK_PARQUET_SELF_ENERGY=TDA_FADDEEV /path/to/QuAcK/bin/QuAcK /path/to/calculation
```

Both `TDAeh` and `TDApp` must be `T`; the driver rejects other values.
The method uses the existing `evParquet` switch and its pair iterations,
convergence thresholds, and DIIS settings. The supplied example requests
one one-body evaluation and up to 64 pair iterations. With `max_it_2b=1`
at the HF reference the kernels are bare, giving the ordinary ADC(3)
limit. A dressed-parquet test requires further pair iterations and checking
the reported pair convergence. Fixed unconverged snapshots can differ
between R and G when DIIS uses different residual metrics; compare fixed
snapshots with `max_diis_2b=1`.

The full upfolded matrix is diagonalized, so `lin_parquet` is unused in this
path. More than one one-body iteration feeds the dominant poles back into
the existing eigenvalue-self-consistent parquet loop. For one evaluation,
the inherited driver may report one-body nonconvergence because an energy
update was intentionally not iterated.

Leaving the environment variable unset retains the existing ADC-Parquet
selection. Explicit `ADC` and `NATIVE` values select the legacy ADC and
native self-energy paths, respectively. This selector applies to the R/G
`evParquet` drivers, not `qsParquet` or `R_GW_Parquet`.

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

`G_TDA_Faddeev_build` forms the physical matrices and ADC(3) couplings;
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
are returned once. The implementation uses dense matrices and is intended
initially for small test systems.

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

This implements the static complete-TDA Faddeev model with ADC(3)
endpoints. Identity with the full parquet Schwinger-Dyson self-energy has
not been established. It is a separate self-energy approximation to test.

The older standalone `G_TDA_Faddeev_Parquet.f90` prototype is not called by
this implementation. It stopped at the redundant auxiliary eigensystem;
the physical reduction avoids its extra sectors and normalization issues.

Algebraic regression checks cover the matrix and external couplings
against QuAcK's existing `G_IPEA_ADC3`, R/G spin reconstruction, active
spaces and empty triplet sectors, the exact Faddeev intertwining identity
for both branches with dense dressed kernels, and SRG limits:

```sh
python tests/check_tda_faddeev.py --workdir /tmp/new-algebraic-checks
```

Build QuAcK first; the check requires NumPy and gfortran. Set `FC` and
`LAPACK_LIBS` if necessary. The molecular test script requires PySCF and
runs H2O/6-31G and H2/STO-3G, with logs retained:

```sh
python tests/check_tda_faddeev_molecules.py --workdir /tmp/new-molecular-checks
```
