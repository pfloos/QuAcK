# Hartree–Fock and ground-state correlation

## Choose the orbital reference

| Reference | Input | When to use it |
|---|---|---|
| Restricted HF | `RHF` | Closed-shell determinants with shared spatial orbitals |
| Unrestricted HF | `UHF` | Collinear spin with independent alpha and beta orbitals |
| Generalized HF | `GHF` | Spin-orbital mixing beyond collinear UHF |
| Restricted open-shell HF | `ROHF` | Restricted spatial orbitals for open shells |
| Complex HF | `cRHF`, `cUHF` | Complex orbitals and dedicated response branches |
| Bogoliubov reference | `HFB` | A quasiparticle reference with anomalous density |

Specify charge and multiplicity to PyDuck. Reference selection is a separate decision in `input/methods`. For open-shell work, record the converged determinant and the occupations used, particularly when applying maximum-overlap methods. See [availability](generated/methods.md) for the post-HF paths.

## Perturbation theory and coupled cluster

MP2 and MP3 add the second- and third-order correlation contributions, respectively. Enable MP2 as well when the output should report both orders. QuAcK also exposes restricted and generalized doubles models, including ring and ladder variants. CCD and CCSD solve amplitude equations iteratively. Use the CC thresholds and iteration limits from the [input reference](generated/inputs.md) and inspect convergence before comparing correlation energies.

The input switch `CCSD(T)` requests perturbative triples on top of CCSD. For comparisons across references, start with a closed-shell molecule at the same geometry, basis, orbital space, and thresholds, then change the reference.

## Configuration interaction

Restricted CIS, CIS(D), CID, and CISD and unrestricted CIS enter the CI driver. For CIS(D), enable the CIS calculation and its doubles correction together. Select the spin blocks in the response options and check which transitions the output actually lists.
