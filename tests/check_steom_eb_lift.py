"""Finite Fock-space reference for the joint osCOT lift (real, two spatial orbitals).

Run with python tests/check_steom_eb_lift.py. Requires NumPy, SciPy and PySCF.
Positive EA-sector poles here denote positive electron-addition energies.
"""
from itertools import combinations_with_replacement, product
from math import factorial, prod

import numpy as np
from pyscf import ao2mo, gto, scf
from scipy.linalg import solve
from scipy.optimize import least_squares


def electron_operators():
    create = np.zeros((4, 16, 16))
    for p, det in product(range(4), range(16)):
        if not det & (1 << p):
            create[p, det | (1 << p), det] = (-1) ** (det % (1 << p)).bit_count()
    excitation = np.array([[sum(create[2*p+s] @ create[2*q+s].T for s in range(2))
                            for q in range(2)] for p in range(2)])
    return create, excitation


def integrals(atoms, charge):
    mol = gto.M(atom=atoms, charge=charge, basis="sto-3g", verbose=0)
    hf = scf.RHF(mol).run(conv_tol=1.e-14)
    assert hf.converged
    c = hf.mo_coeff
    eri = ao2mo.restore(1, ao2mo.kernel(mol, c), 2)
    gap = hf.mo_energy[1] - hf.mo_energy[0]
    omega = np.sqrt(gap * (gap + 4*eri[0, 1, 0, 1]))
    vertex = np.sqrt(2*gap/omega) * eri[:, :, 0, 1]
    return c.T @ hf.get_hcore() @ c, eri, hf.mo_energy, omega, vertex


def normal_exponential(create, lower):
    # Each factor contains either E_pq or E_ai^2/E_ia^2 and one boson.
    # Normal order once relative to occupied spin orbitals 0,1 and the vacuum.
    # Pauli exclusion kills every normal product of degree greater than four.
    quasi = [op.T if p < 2 else op for p, op in enumerate(create)]
    ops = quasi + [op.T for op in quasi]

    def string(p, q, spin):
        r, s = 2*p+spin, 2*q+spin
        return (r+4 if r < 2 else r, s if s < 2 else s+4)

    terms = []
    for j, (up, p, q) in enumerate(product((1, 0), range(2), range(2))):
        for spin in range(2):
            terms.append((j, up, string(p, q, spin)))
    for up, p, spin in product((1, 0), range(2), range(2)):
        q = 1-p
        terms.append((8+2*(1-up)+p, up, string(p,q,spin)+string(p,q,1-spin)))
    coefficients = {}
    for degree in range(1, 5):
        for ids in combinations_with_replacement(range(len(terms)), degree):
            group = [terms[i] for i in ids]
            labels = sum((t[2] for t in group), ())
            if len(set(labels)) != len(labels):
                continue
            swaps = sum(x > y for j,x in enumerate(labels) for y in labels[j+1:])
            factor = (-1)**swaps / prod(factorial(ids.count(i)) for i in set(ids))
            electron = np.eye(16)
            for label in sorted(labels):
                electron = electron @ ops[label]
            up = sum(t[1] for t in group)
            boson = np.linalg.matrix_power(lower.T, up) @ np.linalg.matrix_power(lower, degree-up)
            key = tuple(t[0] for t in group)
            coefficients[key] = coefficients.get(key, 0) + factor*np.kron(electron, boson)
    identity = np.eye(16*len(lower))

    def transform(amplitudes):
        return identity + sum(np.prod(amplitudes[list(key)])*mat for key,mat in coefficients.items())
    return transform


def audit(data, cap, fit, stationary=False):
    h, eri, eps, omega, v = data
    cr, E = electron_operators()
    lower = np.diag(np.sqrt(np.arange(1, cap+1)), 1)
    unit = np.eye(cap+1)
    electronic = np.einsum("pq,pqij->ij", h, E)
    for p,q,r,s in product(range(2), repeat=4):
        electronic += 0.5*eri[p,q,r,s]*(E[p,q] @ E[r,s] - (q==r)*E[p,s])
    free = np.kron(sum(eps[q]*E[q,q] for q in range(2)), unit) + np.kron(np.eye(16), omega*lower.T@lower)
    H = (np.kron(electronic, unit) + np.kron(np.eye(16), omega*lower.T@lower)
         + np.kron(np.einsum("pq,pqij->ij",v,E)-2*v[0,0]*np.eye(16),lower+lower.T))
    fields = np.array([np.kron(cr[0], unit), np.kron(cr[2], unit),
                       np.kron(cr[0], lower), np.kron(cr[2], lower.T)])
    ref = 3*len(unit)

    def project(columns):
        return (fields[:,:,ref] @ columns[:,:,ref].T
                + fields[:,ref,:] @ columns[:,ref,:].T)

    def charged_graph(strength):
        gw = np.block([[np.diag(eps), strength*v],
                       [strength*v.T, np.diag(eps+[-omega,omega])]])
        poles, vectors = np.linalg.eigh(gw)
        picked = np.argmax(vectors[:2]**2, axis=1)
        assert len(set(picked)) == 2
        roots = vectors[:,picked]
        return gw, poles[picked], solve(roots[:2].T, roots[2:].T).T

    gw, poles, _ = charged_graph(1.)
    np.testing.assert_allclose(project(fields), np.eye(4), atol=1.e-14)
    np.testing.assert_allclose(project(H@fields-fields@H), gw, atol=1.e-10, rtol=0)

    def spectrum(number):
        det = np.array([d for d in range(16) if d.bit_count()==number])
        indices = (det[:,None]*len(unit)+np.arange(len(unit))).ravel()
        return np.linalg.eigvalsh(H[np.ix_(indices,indices)])

    ground = spectrum(2)[0]
    charged = [ground-spectrum(1), spectrum(3)-ground]
    distances = np.array([min(abs(pole-levels)) for pole,levels in zip(poles,charged)])
    if not fit:
        return distances

    normal_U = normal_exponential(cr,lower)
    neutral_dets = [d for d in range(16) if d.bit_count()==2]
    neutral = (np.array(neutral_dets)[:,None]*len(unit)+np.arange(len(unit))).ravel()
    column = np.zeros(len(H))

    def transform(x):
        U = normal_U(x)
        if stationary:
            U[:,ref] = column
        return U

    amplitudes = np.zeros(12)
    # Fix a reproducible continuation path for this underdetermined diagnostic.
    for strength in (.02,.05,.1,.2,.4,.6,.8,1.):
        current = free + strength*(H-free)
        if stationary:
            # Diagnostic only: impose the exact model ground state as the reference column.
            state = np.linalg.eigh(current[np.ix_(neutral,neutral)])[1][:,0]
            column[neutral] = state/state[np.flatnonzero(neutral==ref)[0]]
        _, _, graph = charged_graph(strength)

        def residual(x):
            U = transform(x)
            G = solve(U,current@U)
            dressed = project(U@fields[:2]@solve(U,np.eye(len(U))))
            charged_block = project(G@fields[:2]-fields[:2]@G)[2:]
            return np.r_[(dressed[2:]-graph@dressed[:2]).ravel(), charged_block.ravel()]

        result = least_squares(residual, amplitudes, ftol=1.e-14, xtol=1.e-14, gtol=1.e-13)
        amplitudes = result.x
        assert np.linalg.norm(result.fun) < 2.e-8, (strength,result.message,result.fun)

    U = transform(amplitudes)
    G = solve(U,H@U)
    def diagnostics(U,G):
        dressed = project(U@fields[:2]@solve(U,np.eye(len(U))))
        graph_error = np.linalg.norm(solve(dressed[:2].T,dressed[2:].T).T-graph)
        decoupling_error = np.linalg.norm(project(G@fields[:2]-fields[:2]@G)[2:])
        return graph_error, decoupling_error

    graph_error, decoupling_error = diagnostics(U,G)
    ref_error = np.linalg.norm(G[:,ref]-G[ref,ref]*np.eye(len(G))[:,ref])
    print(f"  graph {graph_error:.3e}, charged residual {decoupling_error:.3e} Eh,"
          f" reference residual {ref_error:.7f} Eh")
    assert max(graph_error,decoupling_error) < 2.e-8
    if stationary:
        assert ref_error < 1.e-9
        # Look inside the same generator ansatz for a direction which changes
        # the neutral singles element while preserving the charged constraints.
        singlet = np.kron(E[1,0]@np.eye(16)[:,3]/np.sqrt(2),unit[:,0])
        def singlet_energy(x):
            u = transform(x)
            g = solve(u,H@u)
            return singlet@g@singlet-g[ref,ref]
        boson = np.eye(len(H))[:,ref+1]
        def with_neutral_condition(x):
            u = transform(x)
            g = solve(u,H@u)
            return np.r_[residual(x),boson@g@singlet]
        neutral_fit = least_squares(with_neutral_condition,amplitudes,ftol=1.e-14,xtol=1.e-14,gtol=1.e-13)
        assert np.linalg.norm(neutral_fit.fun) < 2.e-8
        print(f"  residual including neutral singles-to-boson condition {np.linalg.norm(neutral_fit.fun):.3e}")

        amplitudes = neutral_fit.x
        step = 1.e-6
        axes = np.eye(len(amplitudes))
        gradient = np.array([(singlet_energy(amplitudes+step*d)-singlet_energy(amplitudes-step*d))/(2*step)
                             for d in axes])
        jac = np.array([(with_neutral_condition(amplitudes+step*d)-with_neutral_condition(amplitudes-step*d))/(2*step)
                        for d in axes]).T
        _, singular, vh = np.linalg.svd(jac,full_matrices=True)
        null = vh[np.count_nonzero(singular > 1.e-8):]
        direction = null.T@(null@gradient)
        if np.linalg.norm(direction) > 1.e-6:
            direction /= np.linalg.norm(direction)
            shift = .01
            def second_solution(x):
                return np.r_[with_neutral_condition(x), direction@(x-amplitudes)-shift]
            alt = least_squares(second_solution,amplitudes+shift*direction,
                                ftol=1.e-14, xtol=1.e-14, gtol=1.e-13)
            assert np.linalg.norm(alt.fun) < 2.e-8
            gaps = [float(singlet_energy(x)) for x in (amplitudes,alt.x)]
            print("  singlet energies for two constrained generators (Eh)", gaps)
            return distances, abs(gaps[0]-gaps[1])
        return distances, 0.
    return distances


if __name__ == "__main__":
    for atoms,charge in [("H 0 0 0; H 0 0 .74",0),("He 0 0 0; H 0 0 .8",1)]:
        print(atoms)
        data = integrals(atoms,charge)
        audit(data,4,True)
        _, ambiguity = audit(data,4,True,stationary=True)
        if charge == 1:
            assert ambiguity > 1.e-6
        gaps = audit(data,8,False)
        np.testing.assert_allclose(gaps,audit(data,12,False),atol=1.e-9,rtol=0)
        assert np.all(gaps > 1.e-4)
        print("  nearest full removal/addition pole distances (Eh)",gaps)
