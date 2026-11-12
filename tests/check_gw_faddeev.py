#!/usr/bin/env python3
"""GW/BSE FTDA molecular checks using a built QuAcK and PySCF.

python tests/check_gw_faddeev.py --workdir /tmp/quack-gw-ftda-checks
"""
import argparse
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys

from check_gw_parquet import modify


def table(text, marker, columns):
    assert marker in text, text[-4000:]
    result = []
    for line in text.rsplit(marker, 1)[1].splitlines()[1:]:
        fields = line.split()
        if len(fields) == columns and fields[0].isdigit():
            row = [float(value.replace('D', 'E')) for value in fields]
            assert all(math.isfinite(value) for value in row)
            result.append(row)
        elif result:
            break
    assert result, marker
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workdir', type=Path, required=True)
    parser.add_argument('--gw-external', action='store_true')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    work = args.workdir.resolve()
    work.mkdir(parents=True, exist_ok=False)
    binary = root / 'bin/QuAcK'
    env = dict(os.environ, QUACK_ROOT=str(root), OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    methods = modify((root / 'input/methods.default').read_text(), '# RHF ', {0: 'T'})
    methods = modify(methods, '# FTDA FRPA', {0: 'T', 1: 'F'})
    options = modify((root / 'input/options.default').read_text(), '# GW:',
                     {3: 'T', 4: '0.0', 5: 'F', 6: 'F'})
    options = modify(options, '# Parquet:', {0: 'T', 1: 'T', 11: 'F'})
    external = 'GW' if args.gw_external else 'HF'
    options = modify(options, '# Faddeev:', {0: '1', 1: '1', 4: 'T' if args.gw_external else 'F'})

    def run(case, name, command, success=True):
        result = subprocess.run(command, cwd=case, env=env, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=180)
        text = result.stdout
        (work / (name + '.log')).write_text(text)
        if success:
            assert result.returncode == 0 and 'ERROR STOP' not in text, text[-4000:]
            assert 'Fortran runtime error' not in text, text[-4000:]
        return text

    from pyscf import gto
    from pyscf.gto.basis.parse_nwchem import convert_basis_to_nwchem
    cases = [('water', '3\nwater\nO 0 0 0\nH 0.9591 0 0\nH -0.2373 0.9293 0\n', ('H', 'O'), 3),
             ('hydrogen', '2\nhydrogen\nH 0 0 0\nH 0 0 0.74\n', ('H',), 1)]
    variants = [('adc2', 2, 'F'), ('adc2-static3', 2, 'T'), ('adc3', 3, 'F'),
                ('adc3-static3', 3, 'T'), ('dressed', 4, 'F')]
    for name, geometry, elements, nroots in cases:
        case = work / name
        for sub in ('input', 'mol', 'int'):
            (case / sub).mkdir(parents=True, exist_ok=True)
        for source in (root / 'input').iterdir():
            if source.name.endswith('.default'):
                shutil.copy2(source, case / 'input' / source.stem)
        basis = case / 'sto-3g.nw'
        basis.write_text('\n'.join(convert_basis_to_nwchem(el, gto.basis.load('sto-3g', el)) for el in elements))
        (case / 'mol/target.xyz').write_text(geometry)
        (case / 'input/options').write_text(options)
        (case / 'input/methods').write_text(modify(methods, '# FTDA FRPA', {0: 'F'}))
        run(case, name + '-integrals', [sys.executable, str(root / 'PyDuck.py'), '--working_dir',
                                      str(case), '-x', 'target', '-b', str(basis), '-nc'])
        (case / 'input/methods').write_text(methods)
        for variant, coupling, static in variants:
            opts = modify(options, '# Faddeev:', {2: str(coupling), 3: static, 5: 'F', 6: str(nroots)})
            (case / 'input/options').write_text(opts)
            dense = run(case, name + '-' + variant + '-dense', [str(binary), str(case)])
            exact = {int(row[0]): row for row in table(dense, 'TDA Faddeev-Parquet dominant poles:', 4)}
            spectrum = table(dense, 'TDA Faddeev-Parquet poles:', 3)
            (case / 'input/options').write_text(modify(opts, '# Faddeev:', {5: 'T'}))
            text = run(case, name + '-' + variant + '-davidson', [str(binary), str(case)])
            assert text.count('* Restricted G0W0 Calculation *') == 1
            assert f'GW-BSE FTDA: {external} external block; frozen GW internal energies' in text
            assert ('FTDA: third-order static self-energy enabled' in text) == (static == 'T')
            assignments = table(text, 'Faddeev HF assignments:', 3)
            poles = table(text, 'Faddeev-Parquet IP roots:', 5)
            assert len(poles) == len(assignments) == nroots
            for assignment, pole in zip(assignments, poles):
                orbital = int(assignment[1])
                assert abs(pole[1] + exact[orbital][2]) < 2e-6, (name, variant, orbital)
                assert abs(assignment[2] - exact[orbital][3]) < 2e-7
                matching = [row for row in spectrum if abs(row[1] + pole[1]) < 2e-6]
                assert len(matching) == 1
                assert abs(pole[3] - matching[0][2]) < 2e-7
                assert pole[4] < 1.0001e-8
            if variant == 'dressed':
                (case / 'use_chem_pot_in_ppRPA').touch()
                shifted = run(case, name + '-shifted-pp', [str(binary), str(case)])
                shifted_poles = table(shifted, 'Faddeev-Parquet IP roots:', 5)
                for a, b in zip(poles, shifted_poles):
                    assert abs(a[1] - b[1]) < 2e-6
                (case / 'use_chem_pot_in_ppRPA').unlink()
            print('PASS', name, variant, 'dense/Davidson agreement')
        nonlinear = modify(opts, '# GW:', {3: 'F'})
        nonlinear = modify(nonlinear, '# Faddeev:', {5: 'T'})
        (case / 'input/options').write_text(nonlinear)
        text = run(case, name + '-gw-root-search', [str(binary), str(case)])
        assert 'Quasiparticle energies obtained by root search' in text
        assert len(table(text, 'Faddeev-Parquet IP roots:', 5)) == nroots
        (case / 'input/methods').write_text(modify(methods, '# evParquet ', {2: 'T'}))
        text = run(case, name + '-reject-mixed-methods', [str(binary), str(case)], success=False)
        assert 'Select either GWParquet or standalone FTDA/FRPA' in text
    print('All GW/BSE FTDA checks passed. Logs:', work)


if __name__ == '__main__':
    main()
