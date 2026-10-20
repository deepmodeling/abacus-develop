"""Shared helpers for the ASE interface integration tests."""

import contextlib
import os
import re
import shutil
import tempfile
import warnings
from pathlib import Path

import numpy as np
from ase.units import GPa, fs
from ase.stress import full_3x3_to_voigt_6_stress


def test_profile_paths(here):
    """Return the engine command and curated assets for a real ASE run."""
    assets = Path(os.environ.get('ASE_ABACUS_ASSETS', here.parent.parent.parent / 'tests' / 'PP_ORB')).resolve()
    if not assets.is_dir():
        raise FileNotFoundError(f'ASE_ABACUS_ASSETS is not a directory: {assets}')
    return os.environ.get('ABACUS_COMMAND', 'mpirun -np 2 abacus'), assets


def calculate_socket_frame(calc, atoms):
    """Request all three facts while passing ASE's actual structure changes."""
    changes = calc.check_state(atoms)
    calc.calculate(atoms, properties=['energy', 'forces', 'stress'],
                   system_changes=changes)


def close_socket_checked(calc):
    """Close through ASE's public lifecycle and reject client-exit warnings."""
    with warnings.catch_warnings():
        warnings.filterwarnings('error', message=r'Subprocess exited with status .*',
                                category=UserWarning)
        calc.close()


def assert_band_facts(testcase, bands, eig_occ_path):
    """Check ASE bands against the engine's raw eig_occ rows in eV."""
    if os.environ.get('ASE_ABACUS_FORMAT', 'latest') == 'legacy':
        from legacy_facts import assert_band_facts as assert_legacy_band_facts
        return assert_legacy_band_facts(testcase, bands, eig_occ_path)
    values = np.asarray(bands.energies)
    testcase.assertEqual(values.ndim, 3)
    nspin, nk, nbands = values.shape
    testcase.assertGreater(nbands, 0)
    frames = []
    frame = {}
    key = None
    for line in Path(eig_occ_path).read_text().splitlines():
        if line.strip().endswith('# ionic step'):
            if frame:
                frames.append(frame)
            frame = {}
            key = None
            continue
        header = re.match(r'\s*spin=(\d+)\s+k-point=(\d+)/(\d+)\b', line)
        if header:
            spin, point, total = map(int, header.groups())
            testcase.assertEqual(total, nk)
            key = (spin, point)
            testcase.assertNotIn(key, frame)
            frame[key] = []
            continue
        row = re.fullmatch(r'\s*(\d+)\s+([-+\d.Ee]+)\s+([-+\d.Ee]+)\s*', line)
        if row and key is not None:
            testcase.assertEqual(int(row.group(1)), len(frame[key]) + 1)
            frame[key].append(float(row.group(2)))
    if frame:
        frames.append(frame)
    testcase.assertTrue(frames, f'no raw bands in {eig_occ_path}')
    raw = frames[-1]
    testcase.assertEqual(set(raw), {(spin, point)
                                    for spin in range(1, nspin + 1)
                                    for point in range(1, nk + 1)})
    for spin in range(nspin):
        for point in range(nk):
            row = raw[(spin + 1, point + 1)]
            testcase.assertEqual(len(row), nbands)
            testcase.assertTrue(np.allclose(values[spin, point], row,
                                            rtol=0, atol=1e-5))


def assert_collinear_magmoms(testcase, calc, log_path, natoms):
    """Compare moments in atom order with the last raw ABACUS table (uB)."""
    lines = Path(log_path).read_text().splitlines()
    starts = [index for index, line in enumerate(lines)
              if re.fullmatch(r'\s*Total Magnetism \(uB\)\s*', line, re.I)]
    testcase.assertTrue(starts, f'no per-atom magnetic table in {log_path}')
    rows = []
    for line in lines[starts[-1] + 1:]:
        match = re.fullmatch(r'\s*([A-Z][a-z]?\d+)\s+([-+\d.Ee]+)\s*', line)
        if match:
            rows.append((match.group(1), float(match.group(2))))
        elif rows:
            break
    testcase.assertEqual(len(rows), natoms)
    testcase.assertEqual(len({label for label, _ in rows}), natoms)
    moments = np.asarray(calc.results['magmoms'])
    testcase.assertEqual(moments.shape, (natoms,))
    testcase.assertTrue(np.allclose(moments, [value for _, value in rows],
                                    rtol=0, atol=1e-4))


def assert_native_md_facts(testcase, frames, dump_path, *, expected_frames):
    """Compare ASE frames with the raw MD_dump lattice and Cartesian rows."""
    raw_energies = {}
    ion = None
    for line in Path(dump_path).with_name('running_md.log').read_text().splitlines():
        header = (re.search(r'\bION=\s*(\d+)\b', line) or
                  re.search(r'#ION MOVE#\s*(\d+)\b', line))
        if header:
            ion = int(header.group(1))
        energy = re.match(r'\s*E_KohnSham\s+[-+\d.Ee]+\s+([-+\d.Ee]+)', line)
        if ion is not None and energy:
            raw_energies[ion] = float(energy.group(1))
    testcase.assertEqual(set(raw_energies), set(range(1, expected_frames + 1)))
    sections = re.split(r'(?=^MDSTEP:\s*\d+\s*$)',
                        Path(dump_path).read_text(), flags=re.M)
    sections = [section for section in sections if section.startswith('MDSTEP:')]
    testcase.assertEqual(len(sections), expected_frames)
    testcase.assertEqual(len(frames), expected_frames)
    raw_velocities = []
    for step, (atoms, section) in enumerate(zip(frames, sections), start=1):
        lines = section.splitlines()
        alat = re.search(r'^LATTICE_CONSTANT:\s*([-+\d.Ee]+) Angstrom$',
                         section, re.M)
        testcase.assertIsNotNone(alat)
        lattice = lines.index('LATTICE_VECTORS')
        raw_cell = np.array([[float(value) for value in row.split()]
                             for row in lines[lattice + 1:lattice + 4]])
        rows = [row.split() for row in lines if re.match(r'^\s*\d+\s+[A-Z][a-z]?\s+', row)]
        testcase.assertIn('VELOCITY (Angstrom/fs)', section)
        raw_velocities.append(np.asarray([[float(value) for value in row[8:11]]
                                          for row in rows]))
        testcase.assertEqual(len(rows), len(atoms))
        testcase.assertEqual(atoms.get_chemical_symbols(), [row[1] for row in rows])
        testcase.assertTrue(np.allclose(atoms.cell.array, raw_cell * float(alat.group(1)),
                                        rtol=0, atol=1e-8))
        testcase.assertTrue(np.all(atoms.pbc))
        testcase.assertTrue(np.allclose(atoms.positions,
                                        [[float(value) for value in row[2:5]] for row in rows],
                                        rtol=0, atol=1e-8))
        testcase.assertTrue(np.isfinite(atoms.get_potential_energy()))
        testcase.assertAlmostEqual(atoms.get_potential_energy(), raw_energies[step],
                                   delta=1e-4)
        forces = atoms.get_forces()
        testcase.assertEqual(forces.shape, (len(atoms), 3))
        testcase.assertTrue(np.all(np.isfinite(forces)))
        testcase.assertTrue(np.allclose(forces,
                                        [[float(value) for value in row[5:8]] for row in rows],
                                        rtol=0, atol=1e-4))
    testcase.assertTrue(np.all(np.isfinite(raw_velocities)))
    testcase.assertGreater(np.max(np.abs(raw_velocities[0])), 1e-8,
                           'MD starts with zero velocity')
    cell = frames[0].cell.array
    delta = frames[-1].positions - frames[0].positions
    fractional = np.linalg.solve(cell.T, delta.T).T
    fractional -= np.rint(fractional)
    testcase.assertGreater(np.max(np.linalg.norm(fractional @ cell, axis=1)), 1e-6,
                           'MD trajectory is static')


def assert_stru_velocity_facts(testcase, stru_path, ase_velocities):
    """Check ASE Angstrom/ASE-time velocities converted to Bohr/atomic-time."""
    text = Path(stru_path).read_text()
    tokens = re.findall(r'\bv\s+([-+\d.Ee]+)\s+([-+\d.Ee]+)\s+([-+\d.Ee]+)', text)
    actual = np.asarray(tokens, dtype=float)
    testcase.assertEqual(actual.shape, np.asarray(ase_velocities).shape)
    expected = np.asarray(ase_velocities) * fs * 1.8897270 * 0.02418884326505
    testcase.assertGreater(np.max(np.abs(expected)), 1e-8)
    testcase.assertTrue(np.allclose(actual, expected, rtol=2e-6, atol=1e-12),
                        f'STRU v unit mismatch: actual={actual}, expected={expected}')


def assert_scf_facts(testcase, calc, directory, natoms, *, forces=False, stress=False):
    """Compare calculator facts with ABACUS' independently written SCF log."""
    if os.environ.get('ASE_ABACUS_FORMAT', 'latest') == 'legacy':
        from legacy_facts import assert_scf_facts as assert_legacy_scf_facts
        return assert_legacy_scf_facts(testcase, calc, directory, natoms,
                                       forces=forces, stress=stress)
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
    raw_root = os.environ.get('ASE_ABACUS_RAW_DIR')
    if raw_root:
        Path(raw_root).mkdir(parents=True, exist_ok=True)
    tmpdir = tempfile.mkdtemp(prefix=prefix, dir=raw_root)
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
        if not raw_root:
            shutil.rmtree(tmpdir, ignore_errors=True)
