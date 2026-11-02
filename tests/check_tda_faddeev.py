#!/usr/bin/env python3
"""Algebraic TDA Faddeev tests (NumPy + gfortran; build QuAcK first).

python tests/check_tda_faddeev.py --workdir /tmp/new-faddeev-checks
Set FC/LAPACK_LIBS if needed. No molecular integrals are required.
"""
import argparse
import ctypes as ct
import itertools as it
import os
from pathlib import Path
import shlex
import subprocess
import sys

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workdir', type=Path, required=True)
    work = parser.parse_args().workdir.resolve()
    work.mkdir(parents=True, exist_ok=False)
    root = Path(__file__).resolve().parents[1]
    # Independent reference: extract the unmodified matrix-building statements
    # from QuAcK's existing ADC(3), not a copy of the new implementation.
    adc = (root / 'src/ADC/G_IPEA_ADC3.f90').read_text()
    body = adc.split('  call wall_time(start_timing)', 1)[1].split('   call wall_time(end_timing)', 1)[0]
    reference = '''subroutine reference_adc(nOrb,nC,nO,nV,nR,nH,ERI,eHF,H)
implicit none
integer,intent(in) :: nOrb,nC,nO,nV,nR,nH
double precision,intent(in) :: ERI(nOrb,nOrb,nOrb,nOrb),eHF(nOrb)
double precision,intent(out) :: H(nH,nH)
integer :: p,i,j,k,l,a,b,c,d,ija,klc,iab,kcd,n2h1p
double precision,external :: Kronecker_delta
n2h1p=(nO-nC)*(nO-nC-1)*(nOrb-nR-nO)/2
H=0d0
''' + body + '\nend subroutine\n'
    (work / 'reference.f90').write_text(reference)
    shared = work / ('checks.dylib' if sys.platform == 'darwin' else 'checks.so')
    command = shlex.split(os.environ.get('FC', 'gfortran')) + [
        '-dynamiclib' if sys.platform == 'darwin' else '-shared', '-fPIC',
        '-fcheck=all', '-fbacktrace', '-ffree-line-length-none', '-I' + str(root / 'include'),
        str(root / 'src/Parquet/G_TDA_Faddeev_self_energy.f90'),
        str(root / 'src/Parquet/R_TDA_Faddeev_self_energy.f90'), str(work / 'reference.f90'),
        str(root / 'lib/utils.a'), '-o', str(shared)]
    command += shlex.split(os.environ.get('LAPACK_LIBS', '-lblas -llapack'))
    subprocess.run(command, check=True)
    lib = ct.CDLL(str(shared))

    def call(name, *args):
        keep = [ct.c_int(int(a)) if isinstance(a, (int, bool, np.integer)) else
                ct.c_double(a) if isinstance(a, float) else a for a in args]
        getattr(lib, name + '_')(*[a.ctypes.data_as(ct.c_void_p) if isinstance(a, np.ndarray)
                                  else ct.byref(a) for a in keep])

    def arr(a):
        return np.array(a, dtype=np.float64, order='F')

    def pairs(indices):
        return list(it.combinations(indices, 2))

    def kernels(v, occupied, virtual):
        ph = list(it.product(occupied, virtual))
        hh, ee = pairs(occupied), pairs(virtual)
        return [arr([[v[b, i, j, a] - v[b, i, a, j] for j, b in ph] for i, a in ph]),
                arr([[v[i, j, k, l] - v[i, j, l, k] for k, l in hh] for i, j in hh]).reshape((len(hh),)*2, order='F'),
                arr([[v[a, b, c, d] - v[a, b, d, c] for c, d in ee] for a, b in ee]).reshape((len(ee),)*2, order='F')]

    def build(e, v, o, nc=0, nr=0, reg=False, flow=1., block=0, wk=None):
        n = len(e); occ = range(nc, o); vir = range(o, n-nr)
        hh, ee = pairs(occ), pairs(vir)
        minus = [(i, j, a) for i, j in hh for a in vir]
        plus = [(i, a, b) for i in occ for a, b in ee]
        sp = lambda p: 1 - 2*(p % 2)
        if block:
            minus = [q for q in minus if sp(q[0])+sp(q[1])-sp(q[2]) == 1]
            plus = [q for q in plus if sp(q[1])+sp(q[2])-sp(q[0]) == 1]
        nh = (n//2 if block else n) + len(minus) + len(plus)
        h = np.zeros((nh, nh), order='F')
        wk = kernels(v, occ, vir) if wk is None else wk
        call('g_tda_faddeev_build', reg, flow, block, n, nc, o, nr, len(occ)*len(vir),
             len(hh), len(ee), nh, e, e, v, *wk, h)
        return h, minus, plus

    rng = np.random.default_rng(82391)
    n = 4; o = 2
    factors = rng.normal(size=(5, n, n))*.04
    factors = factors + factors.transpose(0, 2, 1)
    v = arr(np.einsum('xpr,xqs->pqrs', factors, factors))
    e = arr([-1.1, -.7, .25, .65])
    for nc, nr in [(0, 0), (1, 0), (0, 1), (1, 1)]:
        h, _, _ = build(e, v, o, nc, nr)
        ref = np.zeros_like(h, order='F')
        call('reference_adc', n, nc, o, n-o, nr, len(h), v, e, ref)
        # Existing ADC couples frozen external orbitals, unlike the new path.
        ref[:nc, n:] = 0.; ref[n:, :nc] = 0.
        if nr:
            ref[n-nr:n, n:] = 0.; ref[n:, n-nr:n] = 0.
        np.testing.assert_allclose(h, ref, atol=1e-13, rtol=1e-13)
    print('PASS bare-kernel matrices and ADC(3) endpoints vs existing G_IPEA_ADC3, including active spaces')

    # Spin-adapted bare pair kernels. Verify the R lift against direct Coulomb
    # spin expansion, including equal spatial pairs and all three triplets.
    ng = 2*n; og = 2*o
    vg = np.zeros((ng,)*4, order='F')
    for p, q, r, s in it.product(range(ng), repeat=4):
        if p % 2 == r % 2 and q % 2 == s % 2:
            vg[p, q, r, s] = v[p//2, q//2, r//2, s//2]
    eg = arr(np.repeat(e, 2))

    def r_kernels(nc, nr):
        occ = range(nc, o); vir = range(o, n-nr)
        ph = list(it.product(occ, vir))
        s = arr([[2*v[b, i, j, a]-v[b, j, a, i] for j, b in ph] for i, a in ph])
        t = arr([[-v[b, j, a, i] for j, b in ph] for i, a in ph])
        result = [s, t]
        for indices in [occ, vir]:
            spairs = list(it.combinations_with_replacement(indices, 2))
            tpairs = pairs(indices)
            result.append(arr([[(v[i, j, k, l]+v[i, j, l, k]) / np.sqrt((1+(i==j))*(1+(k==l)))
                                for k, l in spairs] for i, j in spairs]))
            result.append(arr([[v[i, j, k, l]-v[i, j, l, k] for k, l in tpairs]
                               for i, j in tpairs]).reshape((len(tpairs),)*2, order='F'))
        return result

    for nc, nr in [(0, 0), (1, 0), (0, 1), (1, 1)]:
        rk = r_kernels(nc, nr)
        gk = kernels(vg, range(2*nc, og), range(og, ng-2*nr))
        lifted = [np.zeros_like(a, order='F') for a in gk]
        dims = [len(rk[i]) for i in [0, 2, 3, 4, 5]]
        call('r_tda_faddeev_lift', n, nc, o, nr, *dims, *rk, *lifted)
        for got, expected in zip(lifted, gk):
            np.testing.assert_allclose(got, expected, atol=1e-14)
        h, minus, plus = build(eg, vg, og, 2*nc, 2*nr, wk=lifted)
        ha, ma, pa = build(eg, vg, og, 2*nc, 2*nr, block=1, wk=lifted)
        selected = list(range(0, ng, 2)) + [ng+minus.index(q) for q in ma] + [ng+len(minus)+plus.index(q) for q in pa]
        np.testing.assert_allclose(ha, h[np.ix_(selected, selected)], atol=1e-14)
        outside = sorted(set(range(len(h))) - set(selected))
        np.testing.assert_allclose(h[np.ix_(selected, outside)], 0., atol=1e-14)
        er = np.zeros(n, order='F'); zr = er.copy(order='F')
        call('r_tda_faddeev_self_energy', False, 1., n, nc, o, nr, *dims, e, e, v, *rk, er, zr)
        eig, u = np.linalg.eigh(h)
        # Full spectral measure, independent of potentially degenerate pole
        # eigenvectors: compare the R dominant-pole data to the G resolvent.
        for p in range(nc, n-nr):
            idx = np.argmin(abs(eig-er[p]))
            group = abs(eig-eig[idx]) < 1e-8
            np.testing.assert_allclose(zr[p], np.sum(u[2*p, group]**2), atol=1e-11)
    print('PASS R singlet/triplet lift, spin conservation, R/G poles and weights, empty triplet sectors')

    # Independent ordered-amplitude Faddeev construction: derive V_r on the
    # full ordered tensor space, diagonalize its pair Hamiltonians, and test
    # P M_F = Mbar P for dense, dressed kernels in both branches.
    occ = list(range(o)); vir = list(range(o, n))
    ph = list(it.product(occ, vir)); hh = pairs(occ); ee = pairs(vir)
    wk = []
    for size in [len(ph), len(hh), len(ee)]:
        w = rng.normal(size=(size, size))*.03
        wk.append(arr((w+w.T)/2))
    h, minus, plus = build(e, v, o, wk=wk)
    for removal in [True, False]:
        ordered = ([(i,j,a) for i in occ for j in occ for a in vir] if removal else
                   [(i,a,b) for i in occ for a in vir for b in vir])
        configs = minus if removal else plus
        pair = hh if removal else ee
        third = wk[1] if removal else wk[2]
        spectators = vir if removal else occ
        eh_spectators = occ if removal else vir
        sign = -1 if removal else 1
        free = np.diag([(e[i]+e[j]-e[a]) if removal else (e[j]+e[a]-e[i]) for i,j,a in ordered])
        vs = [np.zeros_like(free) for _ in range(3)]

        def antisym_pair(u, w):
            if u == w:
                return 0, 0
            return pair.index(tuple(sorted((u,w)))), 1 if u < w else -1

        for r, x in enumerate(ordered):
            for s, y in enumerate(ordered):
                if removal:
                    i,j,a = x; k,l,c = y
                    if j == l: vs[0][r,s] = -wk[0][ph.index((i,a)), ph.index((k,c))]
                    if i == k: vs[1][r,s] = -wk[0][ph.index((j,a)), ph.index((l,c))]
                    z,sz = antisym_pair(i,j); t,st = antisym_pair(k,l)
                    if a == c and sz*st: vs[2][r,s] = -.5*sz*st*third[z,t]
                else:
                    i,a,b = x; k,c,d = y
                    if b == d: vs[0][r,s] = wk[0][ph.index((i,a)), ph.index((k,c))]
                    if a == c: vs[1][r,s] = wk[0][ph.index((i,b)), ph.index((k,d))]
                    z,sz = antisym_pair(a,b); t,st = antisym_pair(c,d)
                    if i == k and sz*st: vs[2][r,s] = .5*sz*st*third[z,t]
        om,xeh = np.linalg.eigh(np.diag([e[a]-e[i] for i,a in ph])+wk[0])
        op,xpp = np.linalg.eigh(np.diag([e[i]+e[j] for i,j in pair])+sign*third)
        sizes = [len(eh_spectators)*len(ph)]*2+[len(spectators)*len(pair)]
        sm = [np.zeros((len(ordered), size)) for size in sizes]
        for r, x in enumerate(ordered):
            if removal:
                i,j,a = x; p1=(i,a); p2=(j,a); s1=j; s2=i; s3=a; u,w=i,j
            else:
                i,a,b = x; p1=(i,a); p2=(i,b); s1=b; s2=a; s3=i; u,w=a,b
            for mu in range(len(ph)):
                sm[0][r,eh_spectators.index(s1)*len(ph)+mu]=xeh[ph.index(p1),mu]
                sm[1][r,eh_spectators.index(s2)*len(ph)+mu]=xeh[ph.index(p2),mu]
            z,sz = antisym_pair(u,w)
            if sz:
                for mu in range(len(pair)):
                    sm[2][r,spectators.index(s3)*len(pair)+mu]=.5*sz*xpp[z,mu]
        left = [sm[0].T, sm[1].T, 2*sm[2].T]
        mf = np.block([[left[r] @ ((free+vs[r]) if r==s else vs[r]) @ sm[s]
                        for s in range(3)] for r in range(3)])
        amap = np.zeros((len(configs), len(ordered)))
        for r, x in enumerate(configs):
            exchange = (x[1],x[0],x[2]) if removal else (x[0],x[2],x[1])
            amap[r,ordered.index(x)]=1.; amap[r,ordered.index(exchange)]=-1.
        proj = amap @ np.hstack(sm)
        offset = n if removal else n+len(minus)
        mbar = h[offset:offset+len(configs), offset:offset+len(configs)]
        np.testing.assert_allclose(proj @ mf, mbar @ proj, atol=2e-13)
        # A pseudoinverse is used only as an independent test of the identity,
        # not in the implementation or in its derivation.
        np.testing.assert_allclose(amap @ (free+sum(vs)) @ (.5*amap.T), mbar, atol=2e-13)
    print('PASS exact TDA Faddeev intertwining and ordered-amplitude reduction, both branches')

    # Regulator endpoints stay finite at zero denominators; internal kernels
    # and free propagation must remain untouched.
    zero = arr(np.zeros(n))
    hz, _, _ = build(zero, v, o, reg=True, flow=.7)
    assert np.isfinite(hz).all()
    hr, _, _ = build(e, v, o, reg=True, flow=1e6)
    hu, _, _ = build(e, v, o)
    np.testing.assert_allclose(hr, hu, atol=1e-13)
    print('PASS SRG zero-denominator limit and unregularized large-flow limit')
    print('All algebraic TDA Faddeev checks passed.')


if __name__ == '__main__':
    main()
