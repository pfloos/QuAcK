# Proposed methods

These proposals describe implementation milestones and first validation cases.

| Priority | Method | First milestone | Validation |
|---|---|---|---|
| 1 | Singles-corrected single-shot GW | Restricted IP/EA EOM calculation with orbital relaxation through singles | Match the QCUMBRE reference for small valence and core cases. Recover the doubles-only limit with singles disabled |
| 2 | Standalone generalized TDA Faddeev | Expose a GHF charged-sector solver through its own driver | Reproduce canonical RHF in the spin-orbital limit, then check open-shell IPs and satellites |
| 3 | Small-space FCI and exact one-particle spectra | Compute N, N−1, and N+1 energies and Dyson weights for tiny systems | Check particle-number sectors and total spectral weight. Use the results to test GF, ADC, GW, and Faddeev poles |
| 4 | Unrestricted CCSD and CCSD(T) | Implement UCCSD energies, then add perturbative triples | Recover the closed-shell limit and compare small open-shell molecules to an external reference |
| 5 | Multireference Green's functions | Establish a CASSCF/Dyall zeroth-order propagator before adding screening | Verify active-space IP/EA poles, residues, and the single-reference limit |
