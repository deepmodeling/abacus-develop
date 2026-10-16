"""Shared helpers for the ASE interface integration tests."""

import contextlib
import os
import re
import shutil
import tempfile
from pathlib import Path

import numpy as np
from ase.units import GPa
from ase.stress import full_3x3_to_voigt_6_stress


def test_profile_paths(here):
    """Return the engine command and curated assets for a real ASE run."""
    assets = Path(os.environ.get('ASE_ABACUS_ASSETS', here.parent.parent.parent / 'tests' / 'PP_ORB')).resolve()
    if not assets.is_dir():
        raise FileNotFoundError(f'ASE_ABACUS_ASSETS is not a directory: {assets}')
    return os.environ.get('ABACUS_COMMAND', 'mpirun -np 2 abacus'), assets


def assert_scf_facts(testcase, calc, directory, natoms, *, forces=False, stress=False):
    """Compare calculator facts with ABACUS' independently written SCF log."""
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
    if forces:
        values = np.asarray(calc.results['forces'])
        testcase.assertEqual(values.shape, (natoms, 3))
        testcase.assertTrue(np.all(np.isfinite(values)))
        blocks = re.findall(r'#TOTAL-FORCE \(eV/Angstrom\)#(.*?)(?=\n\s*#|\Z)', content, re.S)
        testcase.assertTrue(blocks, f'no raw force table in {log}')
        rows = re.findall(r'^\s*[A-Z][a-z]?\d+\s+([-+\d.Ee]+)\s+([-+\d.Ee]+)\s+([-+\d.Ee]+)\s*$', blocks[-1], re.M)
        testcase.assertEqual(len(rows), natoms)
        testcase.assertTrue(np.allclose(values, np.asarray(rows, dtype=float), rtol=0, atol=1e-4))
    if stress:
        values = np.asarray(calc.results['stress'])
        testcase.assertEqual(values.shape, (6,))
        testcase.assertTrue(np.all(np.isfinite(values)))
        blocks = re.findall(r'#TOTAL-STRESS \(kbar\)#(.*?)(?=\n\s*#|\Z)', content, re.S)
        testcase.assertTrue(blocks, f'no raw stress table in {log}')
        rows = re.findall(r'^\s*([-+\d.Ee]+)\s+([-+\d.Ee]+)\s+([-+\d.Ee]+)\s*$', blocks[-1], re.M)
        testcase.assertGreaterEqual(len(rows), 3)
        tensor = -0.1 * np.asarray(rows[:3], dtype=float) * GPa
        testcase.assertTrue(np.allclose(values, full_3x3_to_voigt_6_stress(tensor), rtol=0, atol=1e-4))


@contextlib.contextmanager
def preserved_tmpdir(prefix):
    """Temporary ABACUS run directory kept on failure for post-mortem CI.

    On success the directory is removed like tempfile.TemporaryDirectory.
    On exception it is copied into ASE_ABACUS_KEEP_DIR (default
    /tmp/abacus_ase_failure) so remote-CI-only failures (e.g. issue #7794)
    can be root-caused from the uploaded artifact.
    """
    tmpdir = tempfile.mkdtemp(prefix=prefix)
    try:
        yield tmpdir
    except Exception:
        keep_root = Path(os.environ.get(
            'ASE_ABACUS_KEEP_DIR', '/tmp/abacus_ase_failure'))
        keep_root.mkdir(parents=True, exist_ok=True)
        keep = keep_root / Path(tmpdir).name
        shutil.copytree(tmpdir, keep, dirs_exist_ok=True)
        print(f'[ase-test] failure, run dir preserved at {keep}')
        raise
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
