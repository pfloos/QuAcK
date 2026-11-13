# Standalone Faddeev

**Experimental implementation.** FTDA and FRPA are research methods under development.

Enable canonical `RHF` and one of `FTDA` or `FRPA` in `input/methods`. The standard restricted driver dispatches the standalone calculation. FTDA uses forward pair propagation. FRPA includes backward pair amplitudes. The physical charged spaces contain 2h1p and 2p1h configurations.

Run the bundled small example with `python docs/examples/run.py ftda`.

## Options

The optional standalone record in `input/options` has seven positional fields. Keep its header, which the reader uses to find the row.

```text
# Faddeev: reference(0=HF,1=GW) pairs(0=bare,1=BSE,2=CCD) coupling(2=ADC2,3=ADC3,4=dressed,5=CCD) static3 external_GW Davidson n_IP_roots
0 0 3 T F T 1
```

| Field | Values and meaning |
|---|---|
| `reference` | `0` uses HF energies for internal propagation. `1` uses frozen GW energies |
| `pairs` | `0` bare, `1` BSE, `2` CCD-amplitude pair interactions |
| `coupling` | `2` ADC(2), `3` ADC(3), `4` static dressed source, `5` CCD-amplitude coupling |
| `static3` | Include the third-order static self-energy correction |
| `external_GW` | Use GW energies in the external one-particle diagonal |
| `Davidson` | `T` targets principal IP roots. `F` uses dense diagonalization |
| `n_IP_roots` | Positive root count. In Davidson mode, the count is bounded by the number of active occupied orbitals |

For FRPA with `pairs=2`, QuAcK solves CCD at the HF reference and reshapes the
converged doubles amplitudes into each pair channel. For a pair Hessian
`S = [[A, B], [B^T, C]]` and amplitude matrix `T`, it projects onto
`P = [I; T]`, giving the pair Hamiltonian `P^T S P` and metric `I - T^T T`.
The removal channel interchanges the forward and backward pair spaces.
This projection recovers the RPA pair pencil for exact RPA amplitudes.
Using CCD amplitudes defines an approximation. QuAcK checks that the projected
pair metric and Hessian are positive definite. FTDA retains its static
Hermitian closure.

For example, `0 2 3 T F T 1` selects CCD pairs, ADC(3) couplings and
Davidson for one IP root. Use `5` in the coupling column to reuse the same
CCD amplitudes for the external couplings.

BSE pairs require `reference=1`. CCD pairs are available with both FTDA and FRPA and reject coupling `4`. The reader uses the displayed defaults if the standalone record is absent.
Select one standalone flag per calculation. At this revision, the common
Faddeev driver uses the CC convergence settings for CCD amplitudes, the GW
settings for frozen GW energies, and `reg_PA`, `reg_1b`, and `reg_2b` from
the Parquet options record for regularization.

## Upfolding and spectral weights

In the Hermitian FTDA formulation, the one-particle block F is coupled
through U to a charged configuration block K. Eliminating the charged
amplitudes yields a frequency-dependent contribution to the one-particle
self-energy.

```text
Σ(ω) = Σ(static) + U† (ωI − K)⁻¹ U
```

The position of a pole and its one-particle weight answer different
questions. For a normalized eigenvector of the Hermitian upfolded matrix,
the weight is the norm squared of its one-particle component. The Davidson
example also reports the residual of the selected root in hartree.

## Frozen GW/BSE inputs

Set `reference=1` and `pairs=1` in the Faddeev record to use frozen GW energies
and BSE pair kernels. The standalone reader takes regularization controls from
the common Parquet options row.

For a reproducible comparison, report the reference energies, pair kernel, endpoint coupling, static correction, regularization, solver, and root count. When comparing different endpoints, keep the pair kernel and external one-particle block fixed.
