#!/usr/bin/env python3
"""Standalone restricted FTDA/FRPA molecular regression tests."""
import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys

from check_gw_faddeev import table
from check_gw_parquet import modify


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workdir', type=Path, required=True)
    work = parser.parse_args().workdir.resolve()
    work.mkdir(parents=True, exist_ok=False)
    root = Path(__file__).resolve().parents[1]
    env = dict(os.environ, QUACK_ROOT=str(root), OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    binary = root / 'bin/QuAcK'
    methods = modify((root / 'input/methods.default').read_text(), '# RHF ', {0: 'T'})
    options = modify((root / 'input/options.default').read_text(), '# Parquet:', {11: 'F'})

    def run(case, name, command=None, error=None):
        result = subprocess.run(command or [str(binary), str(case)], cwd=case, env=env, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=180)
        text = result.stdout
        (work / (name + '.log')).write_text(text)
        if error:
            assert result.returncode != 0 and error in text, text[-4000:]
        else:
            assert result.returncode == 0 and 'Total wall time for QuAcK' in text, text[-4000:]
        return text

    from pyscf import gto
    from pyscf.gto.basis.parse_nwchem import convert_basis_to_nwchem
    cases = [('water', '3\nwater\nO 0 0 0\nH 0.9591 0 0\nH -0.2373 0.9293 0\n', ('H', 'O'), 3),
             ('hydrogen', '2\nhydrogen\nH 0 0 0\nH 0 0 0.74\n', ('H',), 1)]
    for name, geometry, elements, nroots in cases:
        case = work / name
        for sub in ('input', 'mol', 'int'):
            (case / sub).mkdir(parents=True)
        for source in (root / 'input').glob('*.default'):
            shutil.copy2(source, case / 'input' / source.stem)
        basis = case / 'sto-3g.nw'
        basis.write_text('\n'.join(convert_basis_to_nwchem(el, gto.basis.load('sto-3g', el)) for el in elements))
        (case / 'mol/target.xyz').write_text(geometry)
        (case / 'input/options').write_text(options)
        (case / 'input/methods').write_text(methods)
        run(case, name + '-integrals', [sys.executable, str(root / 'PyDuck.py'), '--working_dir',
                                      str(case), '-x', 'target', '-b', str(basis), '-nc'])
        for index, method in enumerate(('FTDA', 'FRPA')):
            (case / 'input/methods').write_text(modify(methods, '# FTDA FRPA', {index: 'T'}))
            for pairs, coupling in ((0, 2), (0, 3), (0, 4), (0, 5), (2, 3), (2, 5)):
                opts = modify(options, '# Faddeev:', {1: str(pairs), 2: str(coupling), 5: 'F', 6: str(nroots)})
                (case / 'input/options').write_text(opts)
                dense = run(case, f'{name}-{method}-{pairs}-{coupling}-dense')
                if method == 'FRPA' and pairs == 2:
                    assert 'CCD pair interactions with projected RPA metric' in dense
                exact = {orb: energy for orb, _, energy, _ in table(dense, 'dominant poles:', 4)}
                weights = table(dense, 'poles: index, energy', 3)
                (case / 'input/options').write_text(modify(opts, '# Faddeev:', {5: 'T'}))
                iterative = run(case, f'{name}-{method}-{pairs}-{coupling}-davidson')
                assignments = table(iterative, 'Faddeev HF assignments:', 3)
                poles = table(iterative, 'Faddeev-Parquet IP roots:', 5)
                assert len(assignments) == len(poles) == nroots
                for (_, orbital, _), (_, _, pole, z, residual) in zip(assignments, poles):
                    assert abs(exact[orbital] - pole) < 2e-6
                    assert abs(sum(w for _, energy, w in weights if abs(energy-pole) < 2e-6) - z) < 2e-7
                    assert residual < 1.0001e-8
                print(f'PASS {name} {method} pairs {pairs} coupling {coupling} dense/Davidson')
        (case / 'input/methods').write_text(modify(methods, '# evParquet ', {0: 'T'}))
        (case / 'input/options').write_text(modify(options, '# Parquet self-energy:', {0: '2'}))
        run(case, name + '-retired-parquet-selector', error='Faddeev moved to standalone FTDA/FRPA')
        (case / 'input/methods').write_text(modify(methods, '# evParquet ', {2: 'T'}))
        (case / 'input/options').write_text(options + '# GWParquet Faddeev:\nT\n')
        run(case, name + '-retired-gw-selector', error='GW/BSE Faddeev moved to standalone FTDA/FRPA')
    print('All standalone Faddeev molecular checks passed. Logs:', work)


if __name__ == '__main__':
    main()
