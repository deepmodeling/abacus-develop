"""Synthetic format and contraction tests; no ABACUS executable is needed."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np


ROOT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("ao_exporter", str(ROOT / "export_ao_potentials.py"))
exporter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(exporter)


def v_bytes(counts, matrix, flag=1, pairs=None):
    pairs = pairs or [(i, j) for i in range(len(counts)) for j in range(i, len(counts))]
    starts = np.cumsum([0] + counts)
    offset = 24 + 4 * len(counts) + 12 * len(pairs)
    header = struct.pack("<6i", -20129433, 1, sum(counts), flag, len(counts), len(pairs))
    header += struct.pack("<%di" % len(counts), *counts)
    records, payload = b"", b""
    all_pairs = [(i, j) for i in range(len(counts)) for j in range(i, len(counts))]
    for i, j in pairs:
        block = np.asarray(matrix[starts[i]:starts[i+1], starts[j]:starts[j+1]],
                           dtype="<c16" if flag else "<f8").tobytes()
        records += struct.pack("<iq", all_pairs.index((i, j)), offset)
        payload += block
        offset += len(block)
    return header + records + payload


def cs_bytes(nao, naux, blocks, pairs=None):
    pairs = pairs or [(i, j) for i in range(len(nao)) for j in range(len(nao))]
    offset = 28 + 36 * len(pairs)
    header = struct.pack("<iiiqq", -10267453, len(nao), 0, len(pairs), len(pairs))
    records, payload = b"", b""
    for i, j in pairs:
        block = np.asarray(blocks[i, j], dtype="<f8").reshape(nao[i], nao[j], naux[i])
        records += struct.pack("<iiiiidq", i+1, j+1, 0, 0, 0, float(np.max(np.abs(block))), offset)
        payload += block.tobytes()
        offset += block.nbytes
    return header + records + payload


class ExportTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=str(ROOT))
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.v = self.root / "v1_coulomb_full_iq_1_rank0.dat"
        self.cs = self.root / "v1_Cs_data_0.txt"
        self.blocks = {(0, 0): [0.5], (0, 1): [2.0], (1, 0): [3.0], (1, 1): [2.0]}
        self.v.write_bytes(v_bytes([1, 1], np.array([[2.0, 0.25], [0.25, 1.0]])))
        self.cs.write_bytes(cs_bytes([1, 1], [1, 1], self.blocks))

    def test_onsite_and_offsite_factor_convention(self):
        actual = exporter.read_cs_v1(self.cs, [1, 1], [1, 1])
        expected = np.array([[[1., 0.], [2., 3.]], [[2., 3.], [0., 4.]]])
        np.testing.assert_array_equal(actual, expected)
        result = exporter.contract(actual, exporter.read_coulomb_v1(self.v)["matrix"])
        expected_potentials = np.array([[[2., 4.75], [4.75, 1.]], [[.25, 3.5], [3.5, 4.]]])
        np.testing.assert_allclose(result, expected_potentials)

    def test_onsite_ao_transpose_not_blind_factor_two(self):
        self.cs.write_bytes(cs_bytes([2], [1], {(0, 0): [1., 2., 3., 4.]}))
        np.testing.assert_array_equal(exporter.read_cs_v1(self.cs, [2], [1])[:, :, 0], [[2., 5.], [5., 8.]])

    def test_complex_coulomb_preserves_lower_conjugation(self):
        matrix = np.array([[2., .2+.3j], [.2-.3j, 1.]])
        self.v.write_bytes(v_bytes([1, 1], matrix))
        np.testing.assert_array_equal(exporter.read_coulomb_v1(self.v)["matrix"], matrix)

    def test_real_coulomb_payload(self):
        matrix = np.array([[2., .25], [.25, 1.]])
        self.v.write_bytes(v_bytes([1, 1], matrix, flag=0))
        np.testing.assert_array_equal(exporter.read_coulomb_v1(self.v)["matrix"], matrix)

    def test_complex_hermitian_contraction_and_rejection(self):
        # Imaginary components cancel only in this deliberately constructed null space.
        antisymmetric = np.array([[0., 1., -1.], [-1., 0., 1.], [1., -1., 0.]])
        matrix = 3*np.eye(3) + .2j*antisymmetric
        coefficients = np.ones((2, 2, 3))
        np.testing.assert_allclose(exporter.contract(coefficients, matrix), 3*np.ones((3, 2, 2)))
        with self.assertRaisesRegex(ValueError, "potential.*Hermitian"):
            exporter.contract(np.array([[[1., 0.]]]), np.array([[2., .2j], [-.2j, 1.]]))

    def test_rejects_nonhermitian_or_nonfinite_values(self):
        for matrix in (np.array([[1., 2.], [0., 1.]]), np.array([[np.nan, 0.], [0., 1.]])):
            self.v.write_bytes(v_bytes([2], matrix))
            with self.assertRaises(ValueError):
                exporter.read_coulomb_v1(self.v)
        with self.assertRaises(ValueError):
            exporter.contract(np.ones((1, 1, 1)), np.zeros((1, 1)))
        self.cs.write_bytes(cs_bytes([1], [1], {(0, 0): [np.inf]}))
        with self.assertRaises(ValueError):
            exporter.read_cs_v1(self.cs, [1], [1])

    def test_rejects_all_truncations_and_trailing_data(self):
        for path, reader in ((self.v, exporter.read_coulomb_v1),
                             (self.cs, lambda p: exporter.read_cs_v1(p, [1, 1], [1, 1]))):
            original = path.read_bytes()
            for length in range(len(original)):
                with self.subTest(file=path.name, length=length):
                    path.write_bytes(original[:length])
                    with self.assertRaises(ValueError):
                        reader(path)
            path.write_bytes(original + b"x")
            with self.assertRaises(ValueError):
                reader(path)

    def test_duplicate_and_missing_records(self):
        self.v.write_bytes(v_bytes([1, 1], np.eye(2), pairs=[(0, 0), (0, 1), (0, 1)]))
        with self.assertRaisesRegex(ValueError, "duplicate"):
            exporter.read_coulomb_v1(self.v)
        self.cs.write_bytes(cs_bytes([1, 1], [1, 1], self.blocks, pairs=[(0, 0), (0, 1), (1, 0), (1, 0)]))
        with self.assertRaisesRegex(ValueError, "duplicate"):
            exporter.read_cs_v1(self.cs, [1, 1], [1, 1])
        self.cs.write_bytes(cs_bytes([1, 1], [1, 1], self.blocks, pairs=[(0, 0)]))
        with self.assertRaises(ValueError):
            exporter.read_cs_v1(self.cs, [1, 1], [1, 1])

    def test_invalid_headers_offsets_translations_and_dimensions(self):
        for path, cases, reader in (
            (self.v, [(0, "i", 0), (4, "i", 2), (12, "i", 3), (16, "i", -1),
                      (24, "i", 0), (32, "i", 8), (36, "q", -1), (36, "q", 69)],
             exporter.read_coulomb_v1),
            (self.cs, [(0, "i", 0), (8, "i", 1), (20, "q", 3), (28, "i", 0),
                       (32, "i", 3), (36, "i", 1), (48, "d", -1.), (48, "d", 42.), (56, "q", 0)],
             lambda p: exporter.read_cs_v1(p, [1, 1], [1, 1]))):
            original = path.read_bytes()
            for offset, fmt, value in cases:
                with self.subTest(file=path.name, offset=offset, value=value):
                    data = bytearray(original)
                    struct.pack_into("<"+fmt, data, offset, value)
                    path.write_bytes(data)
                    with self.assertRaises(ValueError):
                        reader(path)
            path.write_bytes(original)
        with self.assertRaises(ValueError):
            exporter.read_cs_v1(self.cs, [2, 1], [1, 1])

    def run_cli(self, *extra):
        command = [sys.executable, str(ROOT / "export_ao_potentials.py"), "--producer", str(self.root),
                   "--output", str(self.root / "ao.dat"), "--atom-nao", "1", "1"] + list(extra)
        env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
        return subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, env=env)

    def test_cli_output_roundtrip_and_hash_inventory(self):
        (self.root / "INPUT").write_text("INPUT_PARAMETERS\npseudo_dir pp\norbital_dir orbitals\n", encoding="ascii")
        (self.root / "STRU").write_text("ATOMIC_SPECIES\nX 1 X.upf\nNUMERICAL_ORBITAL\nX.orb\nABFS_ORBITAL\norbitals/X.abfs\n", encoding="ascii")
        for directory, name in (("pp", "X.upf"), ("orbitals", "X.orb"), ("orbitals", "X.abfs")):
            (self.root / directory).mkdir(exist_ok=True)
            (self.root / directory / name).write_text("identity fixture\n", encoding="ascii")
        result = self.run_cli("--confirm-serial-isolated-gamma")
        self.assertEqual(result.returncode, 0, result.stderr)
        output = self.root / "ao.dat"
        self.assertEqual(output.read_text().splitlines()[0], "ABACUS_STERNHEIMER_AO_POTENTIALS 1 2 2 Gamma_Hartree")
        values = np.loadtxt(str(output), skiprows=1)
        parsed = (values[:, 0]+1j*values[:, 1]).reshape(2, 2, 2)
        np.testing.assert_allclose(parsed, exporter.contract(exporter.read_cs_v1(self.cs, [1, 1], [1, 1]),
                                                           exporter.read_coulomb_v1(self.v)["matrix"]))
        manifest = json.loads((self.root / "ao.dat.json").read_text())
        self.assertFalse(manifest["external_identity_verified"])
        self.assertEqual(manifest["units"], "Hartree")
        self.assertTrue(any("no unit tag" in item for item in manifest["assumptions"]))
        self.assertEqual(manifest["output"]["sha256"], hashlib.sha256(output.read_bytes()).hexdigest())
        names = {Path(item["path"]).name for item in manifest["inputs"]}
        self.assertTrue({"INPUT", "STRU", "X.upf", "X.orb", "X.abfs", self.v.name, self.cs.name}.issubset(names))
        for item in manifest["inputs"]:
            self.assertEqual(item["sha256"], hashlib.sha256(Path(item["path"]).read_bytes()).hexdigest())
            if Path(item["path"]).name == "X.abfs":
                self.assertIn("STRU ABFS_ORBITAL", item["roles"])
        again = self.run_cli("--confirm-serial-isolated-gamma")
        self.assertNotEqual(again.returncode, 0)

    def test_requires_confirmation_and_rejects_split_or_non_gamma_producer(self):
        self.assertNotEqual(self.run_cli().returncode, 0)
        extra = self.root / "v1_coulomb_full_iq_2_rank0.dat"
        extra.write_bytes(self.v.read_bytes())
        self.assertNotEqual(self.run_cli("--confirm-serial-isolated-gamma").returncode, 0)
        extra.unlink()
        (self.root / "v1_Cs_data_1.txt").write_bytes(self.cs.read_bytes())
        self.assertNotEqual(self.run_cli("--confirm-serial-isolated-gamma").returncode, 0)
        self.assertFalse((self.root / "ao.dat").exists())

    def test_basis_metadata_and_manifest_only_overwrite_refusal(self):
        basis = self.root / "basis_wfc_out"
        basis.write_text("1 2 abacus\n1 1\n1 1\n0\n", encoding="ascii")
        self.assertEqual(exporter.read_basis_counts(basis, [1, 1]), [1, 1])
        basis.write_text("1 2 abacus\n1 2\n1 1\n0\n", encoding="ascii")
        with self.assertRaises(ValueError):
            exporter.read_basis_counts(basis, [1, 1])
        manifest = self.root / "ao.dat.json"
        manifest.write_text("preserve", encoding="ascii")
        self.assertNotEqual(self.run_cli("--confirm-serial-isolated-gamma").returncode, 0)
        self.assertEqual(manifest.read_text(), "preserve")
        self.assertFalse((self.root / "ao.dat").exists())

    def test_basis_file_cli_and_auxiliary_consistency(self):
        for name in ("basis_wfc_out", "basis_aux_out"):
            (self.root / name).write_text("1 2 abacus\n1 1\n1 1\n0\n", encoding="ascii")
        command = [sys.executable, str(ROOT / "export_ao_potentials.py"),
                   "--producer", str(self.root), "--output", str(self.root / "basis.dat"),
                   "--basis-wfc", str(self.root / "basis_wfc_out"),
                   "--basis-aux", str(self.root / "basis_aux_out"), "--atom-types", "1", "1",
                   "--confirm-serial-isolated-gamma"]
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                universal_newlines=True, env=dict(os.environ, PYTHONDONTWRITEBYTECODE="1"))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotEqual(self.run_cli("--atom-naux", "2", "1", "--confirm-serial-isolated-gamma").returncode, 0)
        self.assertNotEqual(self.run_cli("--atom-types", "0", "1", "--confirm-serial-isolated-gamma").returncode, 0)

    def test_present_kpt_must_be_single_gamma(self):
        path = self.root / "KPT"
        for text in ("K_POINTS\n0\nGamma\n1 1 1 0 0 0\n", "K_POINTS\n1\nDirect\n0 0 0 1\n"):
            path.write_text(text, encoding="ascii")
            exporter.validate_gamma_kpt(path)
        for text in ("K_POINTS\n0\nGamma\n2 1 1 0 0 0\n",
                     "K_POINTS\n1\nDirect\n0.1 0 0 1\n", "K_POINTS\n0\nGamma\n1 1 1 nan 0 0\n"):
            path.write_text(text, encoding="ascii")
            self.assertNotEqual(self.run_cli("--confirm-serial-isolated-gamma").returncode, 0)
            self.assertFalse((self.root / "ao.dat").exists())

    def test_heterogeneous_atom_dimensions(self):
        blocks = {(0, 0): np.array([1., 2., 3., 4.]), (0, 1): np.array([5., 6.]),
                  (1, 0): np.array([7., 8., 9., 10.]), (1, 1): np.array([11., 12.])}
        self.cs.write_bytes(cs_bytes([2, 1], [1, 2], blocks))
        result = exporter.read_cs_v1(self.cs, [2, 1], [1, 2])
        self.assertEqual(result.shape, (3, 3, 3))
        np.testing.assert_array_equal(result[:, :, 0], [[2., 5., 5.], [5., 8., 6.], [5., 6., 0.]])
        np.testing.assert_array_equal(result[:, :, 1], [[0., 0., 7.], [0., 0., 9.], [7., 9., 22.]])
        np.testing.assert_array_equal(result[:, :, 2], [[0., 0., 8.], [0., 0., 10.], [8., 10., 24.]])

    def test_refuses_dangling_output_symlink_and_changed_input(self):
        output = self.root / "ao.dat"
        output.symlink_to(self.root / "missing")
        self.assertNotEqual(self.run_cli("--confirm-serial-isolated-gamma").returncode, 0)
        self.assertTrue(output.is_symlink())
        inventory = exporter.InputInventory()
        inventory.add(self.cs, "Cs")
        self.cs.write_bytes(self.cs.read_bytes() + b"changed")
        with self.assertRaisesRegex(ValueError, "changed"):
            inventory.verify_unchanged()

    def test_publication_race_preserves_existing_manifest_and_removes_own_data(self):
        output = self.root / "ao.dat"
        manifest = self.root / "ao.dat.json"
        real_link = os.link

        def competing_link(source, destination):
            if str(destination) == str(manifest):
                manifest.write_text("another writer", encoding="ascii")
            real_link(source, destination)

        with mock.patch.object(exporter.os, "link", side_effect=competing_link):
            with self.assertRaises(FileExistsError):
                exporter.write_outputs(output, np.ones((1, 1, 1), dtype=complex), {}, exporter.InputInventory())
        self.assertFalse(output.exists())
        self.assertEqual(manifest.read_text(), "another writer")
        self.assertEqual(list(self.root.glob(".ao-export-*")), [])


if __name__ == "__main__":
    unittest.main()
