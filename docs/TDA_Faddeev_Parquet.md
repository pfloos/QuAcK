# TDA Faddeev-Parquet with ADC(2) or ADC(3) external couplings

**Archive.** The inputs and driver calls below apply to revisions before `63ca470b`. Use the [standalone Faddeev guide](faddeev.md) for the current driver and input format.

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

Append this optional section at the end of the calculation's `input/options`:

```text
# Parquet self-energy: method(0=native,1=ADC,2=Faddeev) coupling(0=default,2=ADC2,3=ADC3,4=dressed) Davidson n_IP_roots
2 3 T 4
```

Use method `0` for native Parquet, `1` for ADC-Parquet, or `2` for
Faddeev. The Faddeev coupling codes are `2` for ADC(2), `3` for ADC(3),
and `4` for dressed. Davidson takes `T` or `F` (dense). The final integer
is the Davidson root count. Place the header directly after the preceding
option row and put the values on the next line.

Use `2 4 T 4` for the corresponding dressed calculation. Omitting the
section, or using `2 0 T 1`, retains the previous defaults of dressed
couplings for R and ADC(3) couplings for G.

The source matrices are the occupation-changing eh and pp kernels (`Bph`
and `Bpp`, including the bare interaction). They are evaluated separately
from the TDA propagator, from the same old reducible-vertex snapshot as
the charged pair kernels. They do not introduce backward TDA amplitudes.
Writing `Dph=diag(e_a-e_i)`, `Dh=diag(e_i+e_j)`, and
`De=diag(e_a+e_b)`, the unregularized amplitudes solve

```
Dph Teh + Teh (Dph+Weh) = -Bph
(Dh-Whh) Thh - Thh De = transpose(Bpp)
Dh Tee - Tee (De+Wee) = transpose(Bpp)
```

The endpoints are (ordered hh and ee pair sums)

```
U-_{ija,p} = <ij||pa> + sum_cd <cd||pa> Thh_{ij,cd}
            + sum_kc [<cj||pk> Teh_{kc,ia} - <ci||pk> Teh_{kc,ja}]
U+_{iab,p} = <ab||pi> + sum_kl <kl||pi> Tee_{kl,ab}
            + sum_kc [<kb||pc> Teh_{kc,ia} - <ka||pc> Teh_{kc,ib}]
```

The spectral denominators are `-(e_c-e_k+Omega_eh)`,
`Omega_hh-e_c-e_d`, and `e_k+e_l-Omega_ee`. The eigenvectors are used
directly, avoiding division by the transition-density reconstruction
factors. With `reg_PA=T`, the existing SRG inverse is applied to these
spectral denominators; the linear equations above then describe only the
unregularized limit. Without regularization singular denominators fail
explicitly. The source approximation reduces to the ADC(3) endpoints
through second order in the interaction, with differences starting at
third order in the endpoints (fourth order in the self-energy). Bare
pair kernels at finite interaction therefore need not give exactly ADC(3).
The charged blocks and HF one-particle block are independent of this switch.

R constructs separate singlet/triplet source amplitudes and transforms them
to the same spin-orbital determinant basis used by G. The eh source
transformation reverses the spin-flip column relative to the propagation
transformation. Pair source tensors use rectangular normalized
Clebsch-Gordan transformations. Empty triplet pair sectors are supported.

The matrix-free principal-IP mode is enabled by default in both drivers.
Set `n_IP_roots` to the number of IP eigenpairs wanted, for example `5`,
in either driver. This setting is used only in Davidson mode.
Set Davidson to `F` for the full spectrum. Launch PyDuck or QuAcK normally after editing the input.

Both `TDAeh` and `TDApp` must be `T`; the driver rejects other values.
The method uses the existing `evParquet` switch and its pair iterations,
convergence thresholds, and DIIS settings. The supplied example requests
one one-body evaluation and up to 64 pair iterations. With `max_it_2b=1`
at the HF reference the pair kernels are bare. With ADC(3) couplings
and endpoint regularization disabled this gives the ordinary ADC(3) limit.
For a bare-kernel comparison use `max_it_2b=1` in both runs. The
`dressed` couplings still use the TDA pair modes from that single iteration. A dressed-parquet test requires further pair iterations and checking
the reported pair convergence. Fixed unconverged snapshots can differ
between R and G when DIIS uses different residual metrics; compare fixed
snapshots with `max_diis_2b=1`.

`lin_parquet` is unused in either Faddeev solver. Principal-IP mode requires
`max_it_1b=1`: the driver returns the IP without feeding a partial orbital
energy array into one-body self-consistency. Pair iterations are unchanged.
In full-spectrum mode, more than one one-body iteration feeds all dominant
poles into the existing eigenvalue-self-consistent parquet loop.

These options apply to the R/G `evParquet` drivers.

## Matrix-free IP roots

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

For `n_IP_roots>1`, block Davidson targets the highest occupied active HF
orbitals, starting from their one-hole vectors (plus up to two neighboring
occupied HF guesses). Satellite configurations are not used as initial
guesses. Ritz roots are assigned one-to-one to the HF targets by maximizing
the sum of squared orbital overlaps. Root selection is therefore based on
HF-hole character, rather than on IP order or proximity to the HF gap.
The HF gap midpoint only excludes addition poles. Ordinary Ritz extraction
lets deeper ionizations follow their own energies. Thick restarts retain
the assigned roots and additional candidates carrying target character.

`n_IP_roots` cannot exceed the number of active occupied HF orbitals.
If its lower boundary cuts an HF degeneracy, the solver follows the whole
manifold and reports the requested number of poles in increasing IP order.
G counts spin partners separately; R uses the alpha charged block. The
default subspace is `max(96,4*n_targets+8)` vectors (capped by the matrix
dimension), with a 256-iteration limit and `1e-8` Ha residual tolerance.
Every reported root is checked with a fresh matrix-vector product.
Nonconvergence of any targeted root fails explicitly.

The output table gives root index, IP, pole energy, total one-particle
weight Z, and residual. It now reports the dominant HF-hole-connected
ionizations, rather than filling the root count with low-IP dark states.
If an HF hole fragments strongly, a single dominant pole may be insufficient
to describe its spectrum; compare the full spectrum in that case. Root count
1 retains the original HOMO-character selection described above.

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
and R/G full-spectrum and Davidson agreement. Multiple-root checks cover
HF-character assignments and weights against dense spectra, rejection of
intervening dark satellites, spin
degeneracies, bounded-subspace restarts, and explicit nonconvergence:

```sh
python tests/check_tda_faddeev.py --workdir /tmp/new-algebraic-checks
```

Build QuAcK first; the check requires NumPy and gfortran. Set `FC` and
`LAPACK_LIBS` if necessary. The full-spectrum molecular test script requires
PySCF and selects the dense solver in the input file. It runs
H2O/6-31G and H2/STO-3G, with logs retained:

```sh
python tests/check_tda_faddeev_molecules.py --workdir /tmp/new-molecular-checks
```

## Frozen GW energies and BSE pair kernels

The restricted `GWParquet` driver can run TDA Faddeev by setting its method
switch to `T` and appending this optional record to `input/options`.

```text
# GWParquet Faddeev: use TDA Faddeev with frozen GW energies and BSE kernels
T
```

The default is `F`, which selects the existing native GW-Parquet calculation.
With this flag enabled, `TDAeh` and `TDApp` must both be `T`. The existing
`Parquet self-energy` record selects the coupling (2, 3, or 4), dense versus
Davidson solution, and the number of IP roots. Coupling 0 uses dressed sources,
as in the restricted evParquet driver.

The one-body block uses HF energies. Frozen G0W0 energies enter the free charged
configuration energies, BSE gaps, ADC(3)-form denominators, and dressed source
equations. The pair interactions are the singlet/triplet eh and pp/hh BSE
matrices with their GW free pair energies removed. G0W0 and the BSE kernels use
the same W0 response constructed from HF. W0 can use full RPA even when the
Faddeev pair propagation uses TDA.

The optional third-order static correction retains its canonical HF density
expansion. It is added once to the external block. The external one-body diagonal defaults to HF. Set the following optional
record to `T` to use the frozen GW quasiparticle energies in that diagonal.

```text
# GWParquet Faddeev external GW: use frozen GW energies in the one-body block
T
```

This option directly replaces the HF diagonal. The Faddeev self-energy is
then applied on top of those GW energies, without a GW subtraction. Full-spectrum mode returns orbital energies and weights;
Davidson mode reports the IP roots and HF assignments in its printed tables.

Run `python tests/check_gw_faddeev.py --workdir /tmp/new-gw-ftda-checks` after
building QuAcK. It compares dense and Davidson results for all five coupling
and static-correction combinations on STO-3G water and H2, and checks the pp
chemical-potential convention and the required TDA options.
