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
    # Exercise bounded-subspace restart and explicit nonconvergence through
    # the same operator/solver, using a small test-only Fortran interface.
    restart_wrapper = '''
subroutine controlled_ip(nOrb,nO,nS,nOO,nVV,maxspace,maxiter,e,v,Weh,Whh,Wee,omega,z,hw,residual,info)
use TDA_Faddeev_IP_operator
implicit none
integer,intent(in) :: nOrb,nO,nS,nOO,nVV,maxspace,maxiter
double precision,intent(in) :: e(nOrb),v(nOrb,nOrb,nOrb,nOrb),Weh(nS,nS),Whh(nOO,nOO),Wee(nVV,nVV)
double precision,intent(out) :: omega,z,hw,residual
integer,intent(out) :: info
integer :: iterations
type(Faddeev_operator) :: op
call initialize(op,.false.,.false.,1d0,0,nOrb,nOrb,0,nO,0,nS,nOO,nVV,e,e,v,Weh,Whh,Wee, &
 .false.,spread(spread(0d0,1,nS),2,nS),spread(spread(0d0,1,nOO),2,nVV),spread(spread(0d0,1,nOO),2,nVV))
call davidson(op,nO,omega,z,hw,residual,iterations,info,maxspace,maxiter)
end subroutine
'''
    restart_wrapper += """
subroutine controlled_roots(nOrb,nO,nS,nOO,nVV,nroots,maxspace,maxiter,e,v,Weh,Whh,Wee,omega,z,residual,info)
use TDA_Faddeev_IP_operator
implicit none
integer,intent(in) :: nOrb,nO,nS,nOO,nVV,nroots,maxspace,maxiter
double precision,intent(in) :: e(nOrb),v(nOrb,nOrb,nOrb,nOrb),Weh(nS,nS),Whh(nOO,nOO),Wee(nVV,nVV)
double precision,intent(out) :: omega(nroots),z(nroots),residual(nroots)
integer,intent(out) :: info
integer :: iterations
type(Faddeev_operator) :: op
call initialize(op,.false.,.false.,1d0,0,nOrb,nOrb,0,nO,0,nS,nOO,nVV,e,e,v,Weh,Whh,Wee, &
 .false.,spread(spread(0d0,1,nS),2,nS),spread(spread(0d0,1,nOO),2,nVV),spread(spread(0d0,1,nOO),2,nVV))
call davidson_roots(op,nO,nroots,omega,z,residual,iterations,info,maxspace,maxiter)
end subroutine
"""
    (work / 'controlled.f90').write_text(restart_wrapper)
    shared = work / ('checks.dylib' if sys.platform == 'darwin' else 'checks.so')
    command = shlex.split(os.environ.get('FC', 'gfortran')) + [
        '-dynamiclib' if sys.platform == 'darwin' else '-shared', '-fPIC',
        '-fcheck=all', '-fbacktrace', '-ffree-line-length-none', '-I' + str(root / 'include'),
        str(root / 'src/Parquet/G_TDA_Faddeev_self_energy.f90'),
        str(root / 'src/Parquet/TDA_Faddeev_IP.f90'),
        str(root / 'src/Parquet/R_TDA_Faddeev_self_energy.f90'), str(work / 'reference.f90'), str(work / 'controlled.f90'),
        str(root / 'lib/utils.a'), '-o', str(shared)]
    command += shlex.split(os.environ.get('LAPACK_LIBS', '-lblas -llapack'))
    subprocess.run(command, check=True, cwd=work)
    lib = ct.CDLL(str(shared))

    def call(name, *args, dressed=False, sources=None):
        if name.startswith(('g_tda_faddeev_', 'r_tda_faddeev_')) and name not in (
                'r_tda_faddeev_lift', 'r_tda_faddeev_lift_sources'):
            offset = int(name.endswith('_roots'))
            if name.startswith('r_'):
                n, nc, o, nr = args[3+offset:7+offset]
                no, nv = 2*(o-nc), 2*(n-nr-o)
                ns, nh, ne = no*nv, no*(no-1)//2, nv*(nv-1)//2
            else:
                pos = (8 if name.endswith('_build') else 9) + offset
                ns, nh, ne = args[pos:pos+3]
            if sources is None:
                sources = [np.zeros(shape, order='F') for shape in [(ns,ns),(nh,ne),(nh,ne)]]
            args += (dressed, *sources)
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

    def build(e, v, o, nc=0, nr=0, reg=False, flow=1., block=0, wk=None, adc2=False, sources=None):
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
        call('g_tda_faddeev_build', adc2, reg, flow, block, n, nc, o, nr, len(occ)*len(vir),
             len(hh), len(ee), nh, e, e, v, *wk, h, dressed=sources is not None, sources=sources)
        return h, minus, plus

    def check_action(e, v, o, h, wk, nc=0, nr=0, block=0, spatial=None, reg=False, flow=1., adc2=False, sources=None):
        n = len(e); no = o-nc; nv = n-nr-o
        eri = v if spatial is None else spatial
        vectors = arr(rng.normal(size=(len(h), 5)))
        got = np.zeros_like(vectors, order='F'); diag = arr(np.zeros(len(h)))
        call('g_tda_faddeev_ip_action', adc2, reg, flow, block, len(eri), n, nc, o, nr,
             no*nv, no*(no-1)//2, nv*(nv-1)//2, len(h), 5, e, e, eri, *wk, vectors, got, diag, dressed=sources is not None, sources=sources)
        np.testing.assert_allclose(got, h @ vectors, atol=3e-13, rtol=1e-12)
        np.testing.assert_allclose(diag, np.diag(h), atol=2e-13)

    def solve_ip(e, v, o, wk, nc=0, nr=0, block=0, spatial=None, reg=False, flow=1., adc2=False, sources=None):
        n = len(e); no = o-nc; nv = n-nr-o
        eri = v if spatial is None else spatial
        omega, z, hw, residual = [arr(np.zeros(1)) for _ in range(4)]
        info = np.zeros(1, dtype=np.int32)
        call('g_tda_faddeev_ip', adc2, reg, flow, block, len(eri), n, nc, o, nr,
             no*nv, no*(no-1)//2, nv*(nv-1)//2, e, e, eri, *wk, omega, z, hw, residual, info, dressed=sources is not None, sources=sources)
        assert info[0] == 0, (omega, residual, info)
        assert residual[0] < 1e-8
        return omega[0], z[0], hw[0]

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
        wk = kernels(v, range(nc, o), range(o, n-nr))
        check_action(e, v, o, h, wk, nc, nr)
        omega,z,hw = solve_ip(e, v, o, wk, nc, nr)
        values, vec = np.linalg.eigh(h)
        principal = np.argmax(vec[o-1,:]**2)
        np.testing.assert_allclose(omega, values[principal], atol=1e-10)
        np.testing.assert_allclose(z, sum(vec[:n,principal]**2), atol=1e-9)
    print('PASS bare-kernel matrices and ADC(3) endpoints vs existing G_IPEA_ADC3, including active spaces')

    # ADC(2) changes only the external vertices. Check each element directly
    # against bare antisymmetrized integrals, including frozen external lines.
    def check_adc2(e, v, o, nc=0, nr=0, block=0, spatial=None, wk=None):
        h3, _, _ = build(e, v, o, nc, nr, block=block, wk=wk)
        h2, minus, plus = build(e, v, o, nc, nr, block=block, wk=wk, adc2=True)
        nphys = len(e)//2 if block else len(e)
        expected = np.zeros((len(minus)+len(plus), nphys))
        for q in range(nphys):
            p = 2*q if block else q
            if p < nc or p >= len(e)-nr:
                continue
            for r, (i,j,a) in enumerate(minus):
                expected[r,q] = v[i,j,p,a]-v[i,j,a,p]
            for r, (i,a,b) in enumerate(plus):
                expected[len(minus)+r,q] = v[a,b,p,i]-v[a,b,i,p]
        np.testing.assert_allclose(h2[nphys:,:nphys], expected, atol=1e-14)
        np.testing.assert_allclose(h2[:nphys,nphys:], expected.T, atol=1e-14)
        np.testing.assert_array_equal(h2[nphys:,nphys:], h3[nphys:,nphys:])
        np.testing.assert_array_equal(h2[:nphys,:nphys], h3[:nphys,:nphys])
        wk = kernels(v, range(nc,o), range(o,len(e)-nr)) if wk is None else wk
        check_action(e,v,o,h2,wk,nc,nr,block,spatial,adc2=True)
        omega,z,hw = solve_ip(e,v,o,wk,nc,nr,block,spatial,adc2=True)
        vals,vec = np.linalg.eigh(h2)
        homo = o//2-1 if block else o-1
        target = np.argmax(vec[homo,:]**2)
        np.testing.assert_allclose(omega,vals[target],atol=1e-10)
        np.testing.assert_allclose(z,sum(vec[:nphys,target]**2),atol=1e-9)
        return omega,z

    for nc,nr in [(0,0),(1,0),(0,1),(1,1)]:
        check_adc2(e,v,o,nc,nr)
    # No MP2 denominators may be evaluated in ADC(2), even if all are zero.
    hzero,_,_ = build(arr(np.zeros(n)),v,o,adc2=True)
    hzreg,_,_ = build(arr(np.zeros(n)),v,o,reg=True,flow=.7,adc2=True)
    np.testing.assert_array_equal(hzero,hzreg)
    assert np.isfinite(hzero).all()
    print('PASS ADC(2) bare endpoints, unchanged charged blocks, active spaces, Davidson and zero denominators')


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
        check_action(eg, vg, og, ha, lifted, 2*nc, 2*nr, block=1, spatial=v)
        ip,zprincipal,hwp = solve_ip(eg, vg, og, lifted, 2*nc, 2*nr, block=1, spatial=v)
        er = np.zeros(n, order='F'); zr = er.copy(order='F')
        call('r_tda_faddeev_self_energy', False, False, 1., n, nc, o, nr, *dims, e, e, v, *rk, er, zr)
        np.testing.assert_allclose(ip, er[o-1], atol=1e-10)
        np.testing.assert_allclose(zprincipal, zr[o-1], atol=1e-9)
        omega,z,hw = solve_ip(eg, vg, og, lifted, 2*nc, 2*nr)
        np.testing.assert_allclose(omega, ip, atol=1e-10)
        np.testing.assert_allclose(z, zprincipal, atol=1e-9)
        eig, u = np.linalg.eigh(h)
        # Full spectral measure, independent of potentially degenerate pole
        # eigenvectors: compare the R dominant-pole data to the G resolvent.
        for p in range(nc, n-nr):
            idx = np.argmin(abs(eig-er[p]))
            group = abs(eig-eig[idx]) < 1e-8
            np.testing.assert_allclose(zr[p], np.sum(u[2*p, group]**2), atol=1e-11)
    print('PASS R singlet/triplet lift, spin conservation, R/G poles and weights, empty triplet sectors')

    for nc,nr in [(0,0),(1,0),(0,1),(1,1)]:
        rk = r_kernels(nc,nr)
        gk = kernels(vg,range(2*nc,og),range(og,ng-2*nr))
        dims = [len(rk[i]) for i in [0,2,3,4,5]]
        ip,z = check_adc2(eg,vg,og,2*nc,2*nr,block=1,spatial=v,wk=gk)
        er = arr(np.zeros(n)); zr = arr(np.zeros(n))
        call('r_tda_faddeev_self_energy', True,False,1.,n,nc,o,nr,*dims,e,e,v,*rk,er,zr)
        np.testing.assert_allclose(ip,er[o-1],atol=1e-10)
        np.testing.assert_allclose(z,zr[o-1],atol=1e-9)
        outputs = [arr(np.zeros(1)) for _ in range(4)]
        info = np.zeros(1,dtype=np.int32)
        call('r_tda_faddeev_ip', True,False,1.,n,nc,o,nr,*dims,e,e,v,*rk,*outputs,info)
        assert info[0] == 0 and outputs[3][0] < 1e-8
        np.testing.assert_allclose(outputs[0][0],ip,atol=1e-10)
        np.testing.assert_allclose(outputs[1][0],z,atol=1e-9)
    print('PASS ADC(2) restricted full spectrum and Davidson vs generalized spin reconstruction')


    # Independent ordered-amplitude Faddeev construction: derive V_r on the
    # full ordered tensor space, diagonalize its pair Hamiltonians, and test
    # P M_F = Mbar P for dense, dressed kernels in both branches.
    occ = list(range(o)); vir = list(range(o, n))
    ph = list(it.product(occ, vir)); hh = pairs(occ); ee = pairs(vir)
    wk = []
    for size in [len(ph), len(hh), len(ee)]:
        w = rng.normal(size=(size, size))*.03
        wk.append(arr((w+w.T)/2))
    check_adc2(e,v,o,wk=wk)
    h, minus, plus = build(e, v, o, wk=wk)
    check_action(e, v, o, h, wk)
    omega,z,hw = solve_ip(e, v, o, wk)
    values, vec = np.linalg.eigh(h)
    principal = np.argmax(vec[o-1,:]**2)
    np.testing.assert_allclose(omega, values[principal], atol=1e-10)
    np.testing.assert_allclose(z, sum(vec[:n,principal]**2), atol=1e-9)
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
    hreg,_,_ = build(e,v,o,reg=True,flow=.7)
    wreg = kernels(v,range(o),range(o,n))
    check_action(e,v,o,hreg,wreg,reg=True,flow=.7)
    omega,z,hw = solve_ip(e,v,o,wreg,reg=True,flow=.7)
    values,vec = np.linalg.eigh(hreg)
    principal = np.argmax(vec[o-1,:]**2)
    np.testing.assert_allclose(omega,values[principal],atol=1e-10)

    # A degenerate HF HOMO manifold can split when correlated. Follow the
    # whole manifold and choose the smallest principal IP, not a single
    # arbitrarily selected HF orbital or the nearest-energy satellite.
    edeg = arr([-.7,-.7,.25,.65])
    hdeg,_,_ = build(edeg,v,o)
    wdeg = kernels(v,range(o),range(o,n))
    omega,z,hw = solve_ip(edeg,v,o,wdeg)
    values,vec = np.linalg.eigh(hdeg)
    character = np.sum(vec[:o,:]**2,axis=0)
    roots = np.argsort(character)[-o:]
    np.testing.assert_allclose(omega,np.max(values[roots]),atol=1e-10)
    print('PASS SRG zero-denominator limit and unregularized large-flow limit')
    print('PASS matrix-free products, diagonals and principal-IP Davidson vs dense matrices, R/G, active spaces, dressed kernels, SRG and split HOMO manifold')
    # Force repeated thick restarts in a bounded space, then check that an
    # intentionally exhausted iteration budget returns a failure status.
    nstress, ostress = 8, 4
    stress_factors = rng.normal(size=(6,nstress,nstress))*.05
    stress_factors += stress_factors.transpose(0,2,1)
    vs = arr(np.einsum('xpr,xqs->pqrs', stress_factors, stress_factors))
    es = arr([-1.3,-1.05,-.8,-.55,.25,.4,.65,.9])
    ws = kernels(vs,range(ostress),range(ostress,nstress))
    hs,_,_ = build(es,vs,ostress,wk=ws)
    vals,vecs = np.linalg.eigh(hs)
    target = np.argmax(vecs[ostress-1,:]**2)
    outputs = [arr(np.zeros(1)) for _ in range(4)]
    info = np.zeros(1,dtype=np.int32)
    call('controlled_ip', nstress,ostress,16,6,6,12,256,es,vs,*ws,*outputs,info)
    assert info[0] == 0, (outputs,info)
    np.testing.assert_allclose(outputs[0][0], vals[target], atol=1e-10)
    assert outputs[3][0] < 1e-8
    call('controlled_ip', nstress,ostress,16,6,6,12,1,es,vs,*ws,*outputs,info)
    assert info[0] != 0 and outputs[3][0] > 1e-8
    print('PASS thick restart vs dense root and explicit nonconvergence status')
    def dense_assignment(score):
        # Independent small-target reference by subset dynamic programming.
        nh,ncol=score.shape
        states={0:(0.,[-1]*nh)}
        for col in range(ncol):
            for mask,(value,indices) in list(states.items()):
                for row in range(nh):
                    if mask & (1<<row):
                        continue
                    newmask=mask | (1<<row)
                    candidate=value+score[row,col]
                    if newmask not in states or candidate>states[newmask][0]:
                        selected=indices.copy(); selected[row]=col
                        states[newmask]=(candidate,selected)
        return np.array(states[(1<<nh)-1][1])

    # HF-hole targeted IPs: compare with the distinct dense poles maximizing
    # occupied HF orbital character, rather than the energy-ordered satellites.
    def check_roots(e,v,o,nroots,nc=0,nr=0,block=0,spatial=None,wk=None,adc2=False,reg=False,sources=None):
        n = len(e); no=o-nc; nv=n-nr-o
        wk = kernels(v,range(nc,o),range(o,n-nr)) if wk is None else wk
        h,_,_ = build(e,v,o,nc,nr,block=block,wk=wk,adc2=adc2,reg=reg,flow=.7,sources=sources)
        nphys = n//2 if block else n
        frozen = list(range(nc//2 if block else nc))
        frozen += list(range(nphys-(nr//2 if block else nr),nphys))
        active = [i for i in range(len(h)) if i not in frozen]
        vals,vec = np.linalg.eigh(h[np.ix_(active,active)])
        fermi=(e[o-1]+e[o])/2
        occend = o//2 if block else o
        first=occend-nroots
        ehf=e[::2] if block else e
        while first>(nc//2 if block else nc) and abs(ehf[first-1]-ehf[occend-nroots])<1e-7:
            first-=1
        targets=[active.index(i) for i in range(first,occend)]
        score=vec[targets,:]**2
        score[:,vals>=fermi]=-1.
        matched=dense_assignment(score)
        selected=matched[np.argsort(vals[matched])[::-1]][:nroots]
        outputs=[arr(np.zeros(nroots)) for _ in range(3)]
        info=np.zeros(1,dtype=np.int32)
        eri = v if spatial is None else spatial
        call('g_tda_faddeev_ip_roots',nroots,adc2,reg,.7,block,len(eri),n,nc,o,nr,
             no*nv,no*(no-1)//2,nv*(nv-1)//2,e,e,eri,*wk,*outputs,info, dressed=sources is not None,sources=sources)
        assert info[0]==0 and max(outputs[2])<1e-8, (outputs,info)
        np.testing.assert_allclose(outputs[0],vals[selected],atol=1e-9)
        # Sum over degenerate eigenspaces rather than arbitrary eigenvectors.
        for energy in outputs[0]:
            group=abs(outputs[0]-energy)<1e-8
            exact=abs(vals-energy)<1e-8
            if sum(group)==sum(exact):
                rows=[i for i,j in enumerate(active) if j<nphys]
                np.testing.assert_allclose(sum(outputs[1][group]),np.sum(vec[np.ix_(rows,exact)]**2),atol=1e-8)
        return outputs

    for adc2 in [False,True]:
        for nc,nr in [(0,0),(1,0),(0,1),(1,1)]:
            check_roots(e,v,o,1 if nc else 2,nc,nr,adc2=adc2)
        for nroots in [2,3,4]:
            if nroots<=o:
                check_roots(e,v,o,nroots,adc2=adc2)
            check_roots(eg,vg,og,nroots,adc2=adc2)
            if nroots<=o:
                check_roots(eg,vg,og,nroots,block=1,spatial=v,adc2=adc2)
        rk=r_kernels(0,0); dims=[len(rk[i]) for i in [0,2,3,4,5]]
        expected=check_roots(eg,vg,og,2,block=1,spatial=v,adc2=adc2,reg=True)
        outputs=[arr(np.zeros(2)) for _ in range(3)]; info=np.zeros(1,dtype=np.int32)
        call('r_tda_faddeev_ip_roots',2,adc2,True,.7,n,0,o,0,*dims,e,e,v,*rk,*outputs,info)
        assert info[0]==0
        for got,ref in zip(outputs[:2],expected[:2]):
            np.testing.assert_allclose(got,ref,atol=1e-9)
        # Count 1 must retain the existing principal-HOMO result.
        outputs=[arr(np.zeros(1)) for _ in range(3)]
        wk=kernels(v,range(o),range(o,n))
        principal=solve_ip(e,v,o,wk,adc2=adc2)
        call('g_tda_faddeev_ip_roots',1,adc2,False,1.,0,n,n,0,o,0,o*(n-o),1,1,
             e,e,v,*wk,*outputs,info)
        assert info[0]==0
        np.testing.assert_allclose(outputs[0][0],principal[0],atol=1e-10)
    # A wider dressed problem exercises a bounded multi-root thick restart.
    check_roots(es,vs,ostress,3,wk=ws)
    outputs=[arr(np.zeros(3)) for _ in range(3)]; info=np.zeros(1,dtype=np.int32)
    call('controlled_roots',nstress,ostress,16,6,6,3,16,256,es,vs,*ws,*outputs,info)
    assert info[0]==0 and max(outputs[2])<1e-8, (outputs,info)
    matched=dense_assignment(vecs[ostress-3:ostress,:]**2)
    selected=matched[np.argsort(vals[matched])[::-1]]
    np.testing.assert_allclose(outputs[0],vals[selected],atol=1e-9)
    call('controlled_roots',nstress,ostress,16,6,6,3,16,1,es,vs,*ws,*outputs,info)
    assert info[0]!=0
    # Place two dark satellites at lower IP than both HF holes. HF guesses
    # plus energy ordering would still fail here; character selection must win.
    zero_v=arr(np.zeros_like(v))
    dark_w=[arr(np.zeros((4,4))),arr([[-1.8]]),arr([[0.]])]
    dark=check_roots(e,zero_v,o,2,wk=dark_w)
    np.testing.assert_allclose(dark[0],e[:o][::-1],atol=1e-12)
    np.testing.assert_allclose(dark[1],1.,atol=1e-12)
    # Also test weakly bright satellites, frozen occupied/removed virtuals,
    # and a requested boundary cutting through a degenerate HF manifold.
    check_roots(e,v,o,2,wk=dark_w)
    check_roots(es,vs,ostress,3,nc=1,nr=1)
    check_roots(edeg,v,o,2)
    print('PASS HF-targeted IP roots vs dense orbital-character assignment: R/G, ADC(2)/ADC(3), active spaces, degeneracies, SRG, restart, satellite rejection and failure status')
    print('All algebraic TDA Faddeev checks passed.')


    # New static source approximation: independently solve the linear equations.
    def source(e,o,wk,bph,bpp,nc=0,nr=0,equal=False,reg=False,flow=1.):
        ns,nh,ne=len(wk[0]),len(wk[1]),len(wk[2])
        out=[arr(np.zeros(sh)) for sh in [(ns,ns),(nh,ne),(nh,ne)]]
        call('tda_faddeev_sources',reg,flow,len(e),nc,o,nr,ns,nh,ne,equal,e,*wk,bph,bpp,*out)
        return out

    for nc,nr in [(0,0),(1,0),(0,1),(1,1)]:
        occ=list(range(nc,o)); vir=list(range(o,n-nr)); ph=list(it.product(occ,vir))
        hh=pairs(occ); ee=pairs(vir)
        bare=kernels(v,occ,vir)
        bph=arr([[v[c,a,k,i]-v[c,a,i,k] for i,a in ph] for k,c in ph])
        bpp=arr([[v[a,b,i,j]-v[a,b,j,i] for i,j in hh] for a,b in ee]).reshape((len(ee),len(hh)),order='F')
        wk=[]
        for w in bare:
            delta=arr(rng.normal(size=w.shape)*.003)
            wk.append(arr(w+(delta+delta.T)/2))
        ts=source(e,o,wk,bph,bpp,nc,nr)
        dp=np.diag([e[a]-e[i] for i,a in ph]); dh=np.diag([e[i]+e[j] for i,j in hh])
        de=np.diag([e[a]+e[b] for a,b in ee])
        for got,ref in [(dp@ts[0]+ts[0]@(dp+wk[0]),-bph),
                        ((dh-wk[1])@ts[1]-ts[1]@de,bpp.T),
                        (dh@ts[2]-ts[2]@(de+wk[2]),bpp.T)]:
            np.testing.assert_allclose(got,ref,atol=2e-14)
        h,minus,plus=build(e,v,o,nc,nr,wk=wk,sources=ts)
        old,_,_=build(e,v,o,nc,nr,wk=wk)
        np.testing.assert_array_equal(h[n:,n:],old[n:,n:])
        np.testing.assert_array_equal(h[:n,:n],old[:n,:n])
        check_action(e,v,o,h,wk,nc,nr,sources=ts)
        check_roots(e,v,o,o-nc,nc,nr,wk=wk,sources=ts)
        vals,vec=np.linalg.eigh(h)
        ip,z,_=solve_ip(e,v,o,wk,nc,nr,sources=ts)
        idx=np.argmax(vec[o-1]**2)
        np.testing.assert_allclose(ip,vals[idx],atol=1e-10)
        np.testing.assert_allclose(z,np.sum(vec[:n,idx]**2),atol=1e-9)
        # The first dressed correction differs from ADC(3) at cubic order.
        errors=[]
        for g in [.04,.02,.01]:
            gs=source(e,o,[arr(g*w) for w in bare],arr(g*bph),arr(g*bpp),nc,nr)
            new,_,_=build(e,arr(g*v),o,nc,nr,wk=[arr(g*w) for w in bare],sources=gs)
            adc,_,_=build(e,arr(g*v),o,nc,nr,wk=[arr(g*w) for w in bare])
            errors.append(np.linalg.norm(new[n:,:n]-adc[n:,:n]))
        if errors[-1]>1e-18:
            np.testing.assert_allclose(np.array(errors[:-1])/errors[1:],8,rtol=.03)
        smooth=source(e,o,wk,bph,bpp,nc,nr,reg=True,flow=1e8)
        for got,ref in zip(smooth,ts):np.testing.assert_allclose(got,ref,atol=1e-13)
    zero=source(arr(np.zeros(n)),o,[arr(np.zeros((4,4))),arr(np.zeros((1,1))),arr(np.zeros((1,1)))],
                arr(np.ones((4,4))),arr(np.ones((1,1))),reg=True)
    for tensor in zero:np.testing.assert_array_equal(tensor,0.)
    print('PASS static dressed sources: Sylvester residuals, ADC(3) limit, unchanged charged blocks, dense/Davidson and SRG')

    # Restricted spin reconstruction, including empty triplet sectors.
    for nc,nr in [(0,0),(1,0),(0,1),(1,1)]:
        rk=r_kernels(nc,nr); dims=[len(rk[i]) for i in [0,2,3,4,5]]
        ss=[]
        for trip in [False,True]:
            sp=list(it.combinations(range(nc,o),2)) if trip else list(it.combinations_with_replacement(range(nc,o),2))
            vp=list(it.combinations(range(o,n-nr),2)) if trip else list(it.combinations_with_replacement(range(o,n-nr),2))
            ph=list(it.product(range(nc,o),range(o,n-nr)))
            bph=arr([[(0 if trip else 2*v[a,b,i,j])-v[a,b,j,i] for j,b in ph] for i,a in ph])
            bpp=arr([[(v[a,b,i,j]+(-1 if trip else 1)*v[a,b,j,i])/np.sqrt((1+(i==j))*(1+(a==b)))
                       for i,j in sp] for a,b in vp]).reshape((len(vp),len(sp)),order='F')
            k=[rk[int(trip)],rk[2+int(trip)],rk[4+int(trip)]]
            ss.append(source(e,o,k,bph,bpp,nc,nr,equal=not trip))
        no,nv=2*(o-nc),2*(n-nr-o)
        lifted=[arr(np.zeros(sh)) for sh in [(no*nv,no*nv),(no*(no-1)//2,nv*(nv-1)//2),(no*(no-1)//2,nv*(nv-1)//2)]]
        call('r_tda_faddeev_lift_sources',n,nc,o,nr,*dims,ss[0][0],ss[1][0],ss[0][1],ss[1][1],ss[0][2],ss[1][2],*lifted)
        go=list(range(2*nc,og)); gv=list(range(og,ng-2*nr)); gp=list(it.product(go,gv))
        gh=pairs(go); ge=pairs(gv); gk=kernels(vg,go,gv)
        bph=arr([[vg[c,a,k,i]-vg[c,a,i,k] for i,a in gp] for k,c in gp])
        bpp=arr([[vg[a,b,i,j]-vg[a,b,j,i] for i,j in gh] for a,b in ge]).reshape((len(ge),len(gh)),order='F')
        expected=source(eg,og,gk,bph,bpp,2*nc,2*nr)
        for got,ref in zip(lifted,expected):np.testing.assert_allclose(got,ref,atol=1e-13)
        ha,_,_=build(eg,vg,og,2*nc,2*nr,block=1,wk=gk,sources=lifted)
        check_action(eg,vg,og,ha,gk,2*nc,2*nr,block=1,spatial=v,sources=lifted)
        er=arr(np.zeros(n));zr=er.copy(order='F')
        call('r_tda_faddeev_self_energy',False,False,1.,n,nc,o,nr,*dims,e,e,v,*rk,er,zr,dressed=True,sources=lifted)
        ip,z,_=solve_ip(eg,vg,og,gk,2*nc,2*nr,block=1,spatial=v,sources=lifted)
        np.testing.assert_allclose(ip,er[o-1],atol=1e-10)
        np.testing.assert_allclose(z,zr[o-1],atol=1e-9)
    print('PASS restricted dressed sources vs full spin-orbital calculation, active spaces and empty triplet sectors')

if __name__ == '__main__':
    main()
