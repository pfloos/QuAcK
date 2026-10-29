#!/usr/bin/env python3
"""Molecular smoke/regression tests for a built QuAcK GW-Parquet workflow.

Run with the Python environment used for PyDuck (requires PySCF):
    python tests/check_gw_parquet.py --workdir /tmp/quack-gw-parquet-tests
The selected work directory must be new; logs and integrals are kept there.
"""

import argparse
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys


def modify(text, header, values):
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if line.startswith(header):
            fields = lines[i + 1].split()
            for index, value in values.items():
                fields[index] = value
            lines[i + 1] = '  ' + ' '.join(fields)
            return '\n'.join(lines) + '\n'
    raise AssertionError(f'Missing input record: {header}')


def rows(text):
    marker = 'Native GW-Parquet self-energy and quasiparticle energies'
    assert marker in text, text[-4000:]
    result = []
    for line in text.split(marker, 1)[1].splitlines():
        fields = line.split()
        if len(fields) == 6 and fields[0].isdigit():
            row = [int(fields[0])] + [float(x) for x in fields[1:]]
            assert all(math.isfinite(x) for x in row)
            result.append(row)
    assert result
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workdir', type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    work = args.workdir.resolve()
    work.mkdir(parents=True, exist_ok=False)
    binary = root / 'bin/QuAcK'
    assert binary.is_file(), 'Build QuAcK first'
    env = dict(os.environ, QUACK_ROOT=str(root), OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    methods = modify((root / 'input/methods.default').read_text(), '# RHF ', {0: 'T'})
    methods = modify(methods, '# evParquet ', {2: 'T'})
    # Also enable the ordinary G0W0 switch: the workflow must still run GW once.
    methods = modify(methods, '# G0W0 ', {0: 'T'})
    options = modify((root / 'input/options.default').read_text(), '# GW:', {3: 'T'})
    options = modify(options, '# Parquet:', {0: 'F', 1: 'F', 8: 'T', 9: '1000', 10: '1000', 11: 'T'})

    water = work / 'water'
    for sub in ['input', 'mol', 'int']:
        (water / sub).mkdir(parents=True, exist_ok=True)
    for source in (root / 'input').iterdir():
        if source.name in ['hpc_flags', 'geminal', 'dft']:
            shutil.copy2(source, water / 'input' / source.name)
        elif source.name.endswith('.default'):
            shutil.copy2(source, water / 'input' / source.stem)
    (water / 'input/methods').write_text(methods)
    (water / 'input/options').write_text(options)
    (water / 'mol/water.xyz').write_text('3\nwater\nO 0 0 0\nH 0.9591 0 0\nH -0.2373 0.9293 0\n')

    def run(case, name, command):
        result = subprocess.run(command, cwd=case, env=env, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=180)
        text = result.stdout
        (work / f'{name}.log').write_text(text)
        assert result.returncode == 0, text[-4000:]
        # PyDuck does not propagate a nonzero Fortran exit status.
        assert 'ERROR STOP' not in text and 'Fortran runtime error' not in text, text[-4000:]
        return text

    text = run(water, 'water-full', [sys.executable, str(root / 'PyDuck.py'), '--working_dir',
               str(water), '-x', 'water', '-b', str(root / 'basis/6-31g'), '-nc'])
    reference = rows(text)
    assert len(reference) == 13

    variants = [('water-full', 'F', 'F', 'T', 'T', 'T', '1000'),
                ('water-eh-tda', 'T', 'F', 'T', 'T', 'T', '1000'),
                ('water-pp-tda', 'F', 'T', 'T', 'T', 'T', '1000'),
                ('water-both-tda', 'T', 'T', 'T', 'T', 'T', '1000'),
                ('water-unregularized', 'F', 'F', 'T', 'F', 'T', '1000'),
                ('water-small-flow', 'F', 'F', 'T', 'T', 'T', '0.7'),
                ('water-dyson-roots', 'F', 'F', 'F', 'T', 'T', '1000'),
                ('water-gw-roots', 'F', 'F', 'T', 'T', 'F', '1000')]
    for name, teh, tpp, linear, regularize, gw_linear, flow in variants:
        variant_options = modify(options, '# Parquet:',
                                 {0: teh, 1: tpp, 8: linear, 9: flow, 10: flow, 11: regularize})
        variant_options = modify(variant_options, '# GW:', {3: gw_linear})
        (water / 'input/options').write_text(variant_options)
        text = run(water, name, [str(binary), str(water)])
        data = rows(text)
        assert text.count('* Restricted G0W0 Calculation *') == 1
        assert 'One-body iteration #' not in text and 'Two-body iteration #' not in text
        assert len(data) == 13
        for row, ref in zip(data, reference):
            _, ehf, egw, sigma, z, epa = row
            expected = egw + z * (ehf + sigma - egw) if linear == 'T' else ehf + sigma
            assert abs(epa - expected) < 2e-7, (name, row, expected)
            if gw_linear == 'T':
                assert abs(egw - ref[2]) < 1e-7, 'BSE/PA options fed back into G0W0'
        print(f'PASS {name}: HOMO {data[4][-1]:.9f} eV')

    (water / 'input/options').write_text(options)
    (water / 'use_chem_pot_in_ppRPA').touch()
    shifted = rows(run(water, 'water-pp-chemical-potential', [str(binary), str(water)]))
    (water / 'use_chem_pot_in_ppRPA').unlink()
    assert max(abs(a[-1] - b[-1]) for a, b in zip(shifted, reference)) < 2e-7
    print('PASS absolute pp pole energies with chemical-potential shift')

    (water / 'input/methods').write_text(modify(methods, '# G0W0 ', {0: 'F'}))
    automatic = rows(run(water, 'water-auto-g0w0', [str(binary), str(water)]))
    assert max(abs(a[-1] - b[-1]) for a, b in zip(automatic, reference)) < 2e-7
    (water / 'input/methods').write_text(methods)
    print('PASS automatic G0W0 prerequisite')

    # Frozen core plus a removed virtual: pair dimensions use active counts.
    molecule = (water / 'input/molecule').read_text().splitlines()
    values = molecule[1].split()
    values[3:5] = ['2', '2']
    molecule[1] = ' '.join(values)
    (water / 'input/molecule').write_text('\n'.join(molecule) + '\n')
    (water / 'input/options').write_text(options)
    data = rows(run(water, 'water-active-space', [str(binary), str(water)]))
    assert [r[0] for r in data] == list(range(2, 13))
    print('PASS frozen core and restricted virtual space')

    # Legacy two-field methods record must not consume the next input record.
    legacy = methods.splitlines()
    for i, line in enumerate(legacy):
        if line.startswith('# evParquet '):
            legacy[i] = '# evParquet qsParquet'
            legacy[i + 1] = '  F F'
    (water / 'input/methods').write_text('\n'.join(legacy) + '\n')
    text = run(water, 'legacy-methods', [str(binary), str(water)])
    assert '* Restricted G0W0 Calculation *' in text
    assert 'single-shot GW-Parquet calculation' not in text
    print('PASS legacy two-field methods input')

    # Minimal H2 has empty triplet pair spaces in both branches.
    from pyscf import gto
    from pyscf.gto.basis.parse_nwchem import convert_basis_to_nwchem
    hydrogen = work / 'hydrogen'
    for sub in ['input', 'mol', 'int']:
        (hydrogen / sub).mkdir(parents=True, exist_ok=True)
    for source in (water / 'input').iterdir():
        if source.is_file():
            shutil.copy2(source, hydrogen / 'input' / source.name)
    (hydrogen / 'input/methods').write_text(methods)
    (hydrogen / 'input/options').write_text(options)
    (hydrogen / 'mol/hydrogen.xyz').write_text('2\nhydrogen\nH 0 0 0\nH 0 0 0.74\n')
    basis = hydrogen / 'sto-3g.nw'
    basis.write_text(convert_basis_to_nwchem('H', gto.basis.load('sto-3g', 'H')))
    text = run(hydrogen, 'hydrogen', [sys.executable, str(root / 'PyDuck.py'), '--working_dir',
               str(hydrogen), '-x', 'hydrogen', '-b', str(basis), '-nc'])
    assert len(rows(text)) == 2
    print('PASS empty triplet pair spaces')
    text = run(hydrogen, 'hydrogen-empty-hh-triplet', [sys.executable, str(root / 'PyDuck.py'), '--working_dir',
               str(hydrogen), '-x', 'hydrogen', '-b', str(root / 'basis/6-31g'), '-nc'])
    assert len(rows(text)) == 4
    print('PASS empty hh triplet with nonempty ee triplet space')
    print(f'All GW-Parquet molecular checks passed. Logs: {work}')


if __name__ == '__main__':
    main()
