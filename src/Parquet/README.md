# Overview of the Parquet implementation

## Parameters controling the run

The parameters provided by the user are:
- `max_it_macro` and `max_it_micro` which set the maximum number of iterations of the macro (one-body) and micro (two-body) self-consistent cycles.
- `conv_one_body` and `conv_two_body` which set the convergence threshold of the macro (one-body) and micro (two-body) self-consistent cycles.

The hard-coded parameters are:
- `linearize` which control whether the quasiparticle equation will be linearized or not. Note that the Newton-Raphson has not been implemented yet.
- `TDA` which control whether the Tamm-Dancoff approximation is enforced for the BSE problems or not.
- `print_phLR` and `print_ppLR` control the print of eigenvalues at each diagonalization.
-

## Files and their routines

- `RParquet.f90` is the main file for the restricted Parquet calculation, it is called by `RQuack.f90`. The main task of this file is to control the self-consistent cycles.
- `R_screened_integrals.f90` gathers four subroutines, each one dedicated to the computation of screened integrals in a given channel.
- There are four files dedicated to computed effective interactions in a each channel. For example, `R_eh_singlet_Gam.f90` contains three subroutines: one for the OVOV block, one for the VOVO block and one for the full $N^4$ tensor.

## TODO list

### Check
- [ ] Comment m,s,t channels and perform ehBSE@$GW$ and ppBSE@$GW$
- [ ] Comment d,m channels and perform ehBSE@$GT$ and ppBSE@$GT$
 
### Required

- [ ] Implement diagonal self-energy
- [ ] Implement screened integrals in every channels

### Improvement

- [ ] OpenMP pp Gamma
- [ ] OpenMP eh Gamma
- [ ] DGEMM pp Gamma
- [ ] DGEMM eh Gamma

### Long-term

- [ ] Implement Newton-Raphson solution of the quasiparticle equation
- [ ] Implement Galitskii-Migdal self-energy
