# Single-shot GW-Parquet (restricted RHF)

Set the third parquet method switch in `input/methods`:

```text
# evParquet qsParquet GWParquet
  F         F         T
```

Old two-switch parquet records remain valid and leave GW-Parquet disabled.
Select canonical RHF. For a calculation using Faddeev with frozen GW/BSE inputs, use the [standalone Faddeev driver](faddeev.md). `R_GW_Parquet.f90` runs its own G0W0 calculation, then static
singlet/triplet ehBSE@G0W0 and ppBSE@G0W0. A separately enabled G0W0 switch
does not duplicate that G0W0 calculation. The BSE and internal parquet
energies are frozen at G0W0; no parquet iteration or ADC completion is used.
Both BSE kernels use the same HF dRPA W0, as in the existing G0W0 BSE drivers.
Their external spectral couplings use that same static GW kernel, not a
subsequent parquet update of the irreducible vertices.

The existing GW `lin`, `eta`, `TDA_W`, and `reg` options control the initial
G0W0 calculation. Its density-matrix correction and optional dynamic BSE or
ACFDT calculations are not part of this workflow. The singlet/triplet BSE
channels needed by native restricted parquet are always included.

In the existing Parquet options record:

- `TDAeh` and `TDApp` independently select TDA in the two BSE problems.
- `lin_parquet=T` gives `ePA=eGW+ZPA*(eHF+SigmaPA(eGW)-eGW)`.
- `lin_parquet=F` solves the same frozen-input Dyson equation by a damped
  root search. This does not update any input energy, kernel, or response.
- `reg_PA=T` applies SRG regularization: `reg_1b` is the flow in the native
  self-energy denominators; `reg_2b` is the flow on the common W0 screening
  modes entering both static BSE kernels and their external couplings.
  Both flows must be positive. Large flows approach the unregularized result.
- `reg_PA=F` disables these parquet SRG factors. GW broadening `eta` and
  the initial GW SRG setting remain independent.
- The parquet iteration counts, convergence tolerances, and DIIS lengths
  have no effect on this single-shot method.

The native denominator-safe self-energy is evaluated at the GW expansion
point (or the final Dyson root), keeping all internal GW energies fixed.
For this method, the derivative used in `ZPA` includes the derivative of the
SRG factors. The existing `R_Parquet_self_energy_alt_diag` entry point keeps
its historical derivative convention and calling sequence.

The output reports HF, G0W0, and final parquet energies, the native
self-energy, and its derivative-based weight. It does not report a
Galitskii-Migdal energy, which the reused native routine does not compute.

Example settings are provided in `tests/inp/methods.GWParquet` and
`tests/inp/options.GWParquet` (linearized final Dyson equation, full BSEs,
and SRG regularization). Copy these to the calculation's `input/methods`
and `input/options` to use them.

After building QuAcK, molecular regression checks can be run in a new
directory with the Python environment used for PyDuck:

```sh
python tests/check_gw_parquet.py --workdir /tmp/quack-gw-parquet-tests
```

These cover independent TDA switches, regularization, both Dyson solvers,
GW input invariance, frozen/restricted spaces, empty triplet pair spaces,
and old methods files. `tests/gw_parquet_checks.f90` additionally checks
finite-difference derivatives, legacy entry-point compatibility, and the
GW-kernel/half-vertex normalization. It links against the built QuAcK
libraries and the same BLAS/LAPACK libraries used by QuAcK.
