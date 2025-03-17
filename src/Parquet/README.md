# Overview of the Parquet implementation

## Parameters controling the run

The parameters provided by the user are:
- `max_it_macro` and `max_it_micro` which set the maximum number of iterations of the macro (one-body) and micro (two-body) self-consistent cycles.
- `conv_one_body` and `conv_two_body` which set the convergence threshold of the macro (one-body) and micro (two-body) self-consistent cycles.
-
-

The hard-coded parameters are:
- `linearize` which control whether the quasiparticle equation will be linearized or not. Note that the Newton-Raphson has not been implemented yet.
- `TDA` which control whether the Tamm-Dancoff approximation is enforced for the BSE problems or not.
-
-

## Files and their routines
`RParquet.f90` is the main file for the restricted Parquet calculation, it is called by `RQuack.f90`. The main task of this file is to control the self-consistent cycles.

## TODO list

- [ ] Write the TODO list
