#!/usr/bin/env python3
"""R/G molecular integration tests; requires built QuAcK and PySCF.
Build with do_Faddeev_self_energy=.true. and do_ADC_self_energy=.false.
and do_principal_IP=.false. in both R_evParquet.f90 and G_evParquet.f90.
python tests/check_tda_faddeev_molecules.py --workdir /tmp/new-faddeev-molecules
"""
import argparse
import math
import os
import re
from pathlib import Path
import shutil
import subprocess
import sys

from check_gw_parquet import modify


def poles(text):
    marker = 'TDA Faddeev-Parquet dominant poles: orbital, HF (eV), QP (eV), Z'
    assert marker in text, text[-4000:]
    result = []
    for line in text.rsplit(marker, 1)[1].splitlines():
        fields = line.split()
        if len(fields) == 4 and fields[0].isdigit():
            try:
                row = [int(fields[0])] + [float(x) for x in fields[1:]]
            except ValueError:
                break
            assert all(math.isfinite(x) for x in row)
            result.append(row)
        elif result:
            break
    assert result
    assert all(-1e-8 <= row[-1] <= 1+1e-8 for row in result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workdir', type=Path, required=True)
    work = parser.parse_args().workdir.resolve()
    work.mkdir(parents=True, exist_ok=False)
    root = Path(__file__).resolve().parents[1]
    binary = root / 'bin/QuAcK'
    env = dict(os.environ, QUACK_ROOT=str(root),
               OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    base_methods = (root / 'input/methods.default').read_text()
    base_options = (root / 'input/options.default').read_text()

    def inputs(case, kind, iterations=1, reg='F'):
        methods = modify(base_methods, '# RHF ', {0: 'T' if kind == 'R' else 'F',
                                                 2: 'T' if kind == 'G' else 'F'})
        methods = modify(methods, '# evParquet ', {0: 'T'})
        options = modify(base_options, '# Parquet:',
                         {0: 'T', 1: 'T', 2: '1', 3: '1e-8', 4: str(iterations),
                          5: '1e-6', 6: '1', 7: '1', 9: '1000', 10: '1000', 11: reg})
        (case / 'input/methods').write_text(methods)
        (case / 'input/options').write_text(options)

    def run(case, name, command=None, expect_failure=False):
        result = subprocess.run(command or [str(binary), str(case)], cwd=case,
                                env=env, text=True, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, timeout=300)
        text = result.stdout
        (work / f'{name}.log').write_text(text)
        if expect_failure:
            assert 'ERROR STOP' in text
            return text
        assert result.returncode == 0 and 'ERROR STOP' not in text and 'Fortran runtime error' not in text, text[-5000:]
        return text

    water = work / 'water'
    for sub in ['input', 'mol', 'int']:
        (water / sub).mkdir(parents=True, exist_ok=True)
    for source in (root / 'input').iterdir():
        if source.name.endswith('.default'):
            shutil.copy2(source, water / 'input' / source.stem)
        elif source.name in ['hpc_flags', 'geminal', 'dft']:
            shutil.copy2(source, water / 'input' / source.name)
    (water / 'mol/water.xyz').write_text('3\nwater\nO 0 0 0\nH 0.9591 0 0\nH -0.2373 0.9293 0\n')
    inputs(water, 'R')
    text = run(water, 'water-R-bare', [sys.executable, str(root / 'PyDuck.py'), '--working_dir',
                                     str(water), '-x', 'water', '-b', str(root / 'basis/6-31g'), '-nc'])
    refs = {'R': poles(text)}
    original_molecule = (water / 'input/molecule').read_text()
    for kind in ['G', 'R']:
        inputs(water, kind)
        refs[kind] = poles(run(water, f'water-{kind}-bare'))
    for rg, gg in zip(refs['R'], refs['G'][::2]):
        assert abs(rg[2]-gg[2]) < 2e-6, (rg, gg)
        assert abs(rg[3]-gg[3]) < 2e-6, (rg, gg)
    print(f'PASS H2O bare-kernel R/G ADC(3): HOMO {refs["R"][4][2]:.9f} eV', flush=True)

    dressed = {}
    for kind in ['R', 'G']:
        inputs(water, kind, iterations=3, reg='T')
        text = run(water, f'water-{kind}-dressed')
        data = poles(text)
        dressed[kind] = data
        # Fixed iteration snapshots exercise dressed kernels without claiming
        # convergence of the existing parent parquet iteration.
        assert re.search(r'Two-body iteration #\s+3\b', text)
        print(f'PASS H2O {kind} dressed-kernel snapshot (3 pair iterations): HOMO {data[4 if kind=="R" else 9][2]:.9f} eV', flush=True)
    for rg, gg in zip(dressed['R'], dressed['G'][::2]):
        assert abs(rg[2]-gg[2]) < 2e-6, (rg, gg)
        assert abs(rg[3]-gg[3]) < 2e-6, (rg, gg)
    print('PASS R/G agreement with dressed kernels (same fixed iterations, no DIIS)', flush=True)

    # The chemical-potential option changes the pp reference energy only.
    inputs(water, 'R', iterations=3)
    unshifted = poles(run(water, 'water-R-three-iterations'))
    (water / 'use_chem_pot_in_ppRPA').touch()
    shifted = poles(run(water, 'water-R-chemical-potential'))
    (water / 'use_chem_pot_in_ppRPA').unlink()
    assert max(abs(a[2]-b[2]) for a,b in zip(unshifted, shifted)) < 2e-6
    print('PASS restricted pp chemical-potential reference invariance', flush=True)

    molecule = original_molecule.splitlines()
    values = molecule[1].split(); values[3:5] = ['2', '2']; molecule[1] = ' '.join(values)
    (water / 'input/molecule').write_text('\n'.join(molecule)+'\n')
    for kind in ['R','G']:
        inputs(water, kind)
        data = poles(run(water, f'water-{kind}-active-space'))
        assert len(data) == (11 if kind == 'R' else 22)
    print('PASS frozen occupied and removed virtual in both drivers', flush=True)
    (water / 'input/molecule').write_text(original_molecule)
    inputs(water, 'R')
    (water / 'input/options').write_text(modify((water / 'input/options').read_text(), '# Parquet:', {0:'F'}))
    text = run(water, 'reject-non-TDA', expect_failure=True)
    assert 'requires both TDAeh=T and TDApp=T' in text
    print('PASS rejects non-TDA selection', flush=True)

    from pyscf import gto
    from pyscf.gto.basis.parse_nwchem import convert_basis_to_nwchem
    hydrogen = work / 'hydrogen'
    shutil.copytree(water / 'input', hydrogen / 'input')
    (hydrogen / 'mol').mkdir(); (hydrogen / 'int').mkdir()
    (hydrogen / 'mol/hydrogen.xyz').write_text('2\nhydrogen\nH 0 0 0\nH 0 0 0.74\n')
    basis = hydrogen / 'sto-3g.nw'
    basis.write_text(convert_basis_to_nwchem('H', gto.basis.load('sto-3g', 'H')))
    inputs(hydrogen, 'R', iterations=3)
    data = poles(run(hydrogen, 'hydrogen-R', [sys.executable, str(root / 'PyDuck.py'), '--working_dir',
                       str(hydrogen), '-x', 'hydrogen', '-b', str(basis), '-nc']))
    inputs(hydrogen, 'G', iterations=3)
    data_g = poles(run(hydrogen, 'hydrogen-G'))
    assert max(abs(a[2]-b[2]) for a,b in zip(data,data_g[::2])) < 2e-6
    print('PASS H2/STO-3G R/G with empty spatial triplet pair sectors', flush=True)
    print(f'All molecular TDA Faddeev checks passed. Logs: {work}')


if __name__ == '__main__':
    main()
