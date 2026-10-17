"""Small raw-output oracles for ABACUS v3.10/LTS tables."""

import re
from pathlib import Path

import numpy as np
from ase.units import GPa
from ase.stress import full_3x3_to_voigt_6_stress


def assert_band_facts(testcase, bands, running_nscf_path):
    """Compare ASE bands with literal LTS running_nscf.log rows."""
    values = np.asarray(bands.energies)
    testcase.assertEqual(values.ndim, 3)
    nspin, nk, nbands = values.shape
    testcase.assertGreater(nbands, 0)
    testcase.assertEqual(nspin, 1, 'this silicon LTS band smoke is spin-unpolarized')
    lines = Path(running_nscf_path).read_text().splitlines()
    testcase.assertTrue(lines, running_nscf_path)
    headings = []
    for i, line in enumerate(lines):
        match = re.match(r'\s*(\d+)/(\d+)\s+kpoint\s+\(Cartesian\)\s*=.*\(\d+\s+pws\)', line)
        if match:
            headings.append((i, int(match.group(1)), int(match.group(2))))
    testcase.assertEqual(len(headings), nk)
    raw = []
    for i, point, total in headings:
        testcase.assertEqual(total, nk)
        testcase.assertEqual(point, len(raw) + 1)
        energies = []
        for line in lines[i + 1:]:
            match = re.fullmatch(r'\s*(\d+)\s+([-+\d.Ee]+)\s+([-+\d.Ee]+)\s*', line)
            if match:
                testcase.assertEqual(int(match.group(1)), len(energies) + 1)
                energies.append(float(match.group(2)))
                if len(energies) == nbands:
                    break
            elif energies:
                break
        testcase.assertEqual(len(energies), nbands)
        raw.append(energies)
    testcase.assertTrue(np.allclose(values[0], np.asarray(raw), rtol=0, atol=1e-5))


def assert_scf_facts(testcase, calc, directory, natoms, *, forces=False, stress=False):
    """Compare energy, optional EF, forces and stress with literal LTS log."""
    log = Path(directory) / 'OUT.ABACUS' / 'running_scf.log'
    content = log.read_text()
    energies = re.findall(r'!FINAL_ETOT_IS\s+([-+\d.Ee]+)\s+eV', content)
    testcase.assertTrue(energies, f'no final energy in {log}')
    energy = float(calc.results['energy'])
    testcase.assertTrue(np.isfinite(energy))
    testcase.assertAlmostEqual(energy, float(energies[-1]), delta=1e-4)
    fermi = re.findall(r'^\s*E_Fermi\s+[-+\d.Ee]+\s+([-+\d.Ee]+)', content, re.M)
    if fermi:
        testcase.assertIn('fermi_level', calc.results)
        testcase.assertAlmostEqual(float(calc.results['fermi_level']), float(fermi[-1]), delta=1e-4)
    lines = content.splitlines()
    if forces:
        values = np.asarray(calc.results['forces'])
        testcase.assertEqual(values.shape, (natoms, 3))
        testcase.assertTrue(np.all(np.isfinite(values)))
        starts = [i for i, line in enumerate(lines) if re.match(
            r'^\s*TOTAL-FORCE\s*\(eV/Angstrom\)', line, re.I)]
        testcase.assertTrue(starts, f'no LTS force table in {log}')
        rows = []
        for line in lines[starts[-1] + 1:]:
            match = re.fullmatch(r'\s*[A-Z][a-z]?\d+\s+([-+\d.Ee]+)\s+([-+\d.Ee]+)\s+([-+\d.Ee]+)\s*', line)
            if match:
                rows.append([float(value) for value in match.groups()])
                if len(rows) == natoms:
                    break
            elif rows:
                break
        testcase.assertEqual(len(rows), natoms)
        testcase.assertTrue(np.allclose(values, rows, rtol=0, atol=1e-4))
    if stress:
        values = np.asarray(calc.results['stress'])
        testcase.assertEqual(values.shape, (6,))
        testcase.assertTrue(np.all(np.isfinite(values)))
        starts = [i for i, line in enumerate(lines) if re.match(
            r'^\s*TOTAL-STRESS\s*\(KBAR\)', line, re.I)]
        testcase.assertTrue(starts, f'no LTS stress table in {log}')
        rows = []
        for line in lines[starts[-1] + 1:]:
            match = re.fullmatch(r'\s*([-+\d.Ee]+)\s+([-+\d.Ee]+)\s+([-+\d.Ee]+)\s*', line)
            if match:
                rows.append([float(value) for value in match.groups()])
                if len(rows) == 3:
                    break
            elif rows:
                break
        testcase.assertEqual(len(rows), 3)
        tensor = -0.1 * np.asarray(rows, dtype=float) * GPa
        testcase.assertTrue(np.allclose(values, full_3x3_to_voigt_6_stress(tensor),
                                            rtol=0, atol=1e-4))
