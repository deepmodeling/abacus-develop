#!/usr/bin/env python3
import json
import os
import struct
import subprocess
import sys
import tempfile
import unittest

import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
import analyze_librpa_abf_overlap as analyzer



SCRIPT = os.path.join(ROOT, "analyze_librpa_abf_overlap.py")


def write_basis(directory, type_sizes=(2,), multiplicities=(1,)):
    declared_total = sum(size * multiplicity for size, multiplicity in zip(type_sizes, multiplicities))
    with open(os.path.join(directory, "basis_aux_shrink_out"), "w", encoding="utf-8") as handle:
        handle.write("{:10d}{:10d}    abacus\n".format(len(type_sizes), declared_total))
        for itype, size in enumerate(type_sizes, 1):
            handle.write("{:10d}{:10d}\n".format(itype, size))
        for itype, size in enumerate(type_sizes, 1):
            shells = [0] * size
            handle.write("{:10d}{:10d}\n".format(itype, len(shells)))
            for l_value in shells:
                handle.write("{:10d}\n".format(l_value))


def write_raw(directory, matrix, iq=1, atom_naux=(2,), marker=analyzer.RAW_MARKER):
    matrix = np.asarray(matrix, dtype=np.complex128)
    naux = matrix.shape[0]
    path = os.path.join(directory, "v1_abf_overlap_active_iq_{}.dat".format(iq))
    with open(path, "wb") as handle:
        handle.write(struct.pack("<6i4d", marker, 1, iq, 1, naux, len(atom_naux),
                                 1.0, 0.0, 0.0, 0.0))
        handle.write(struct.pack("<{}i".format(len(atom_naux)), *atom_naux))
        handle.write(np.asarray(matrix, dtype="<c16").tobytes(order="C"))
    return path


def pair_index(i, j, natom):
    return sum(natom - left for left in range(i)) + (j - i)


def write_coulomb(directory, matrix, iq=1, atom_naux=(2,), omit_pair=None, gap_after_first=0):
    matrix = np.asarray(matrix, dtype=np.complex128)
    natom = len(atom_naux)
    naux = sum(atom_naux)
    shifts = np.cumsum([0] + list(atom_naux))
    records = []
    for i in range(natom):
        for j in range(i, natom):
            index = pair_index(i, j, natom)
            if index == omit_pair:
                continue
            payload = matrix[shifts[i]:shifts[i + 1], shifts[j]:shifts[j + 1]]
            records.append((index, payload))
    header_size = 24 + 4 * natom + 12 * len(records)
    offset = header_size
    serial = []
    for record_index, (index, payload) in enumerate(records):
        serial.append((index, offset, payload))
        offset += payload.size * 16
        if record_index == 0:
            offset += gap_after_first
    path = os.path.join(directory, "v1_coulomb_full_iq_{}_rank0.dat".format(iq))
    with open(path, "wb") as handle:
        handle.write(struct.pack("<6i", analyzer.COULOMB_MARKER, iq, naux, 1, natom, len(serial)))
        handle.write(struct.pack("<{}i".format(natom), *atom_naux))
        for index, block_offset, _ in serial:
            handle.write(struct.pack("<iq", index, block_offset))
        for record_index, (_, _, payload) in enumerate(serial):
            handle.write(np.asarray(payload, dtype="<c16").tobytes(order="C"))
            if record_index == 0:
                handle.write(b"\\0" * gap_after_first)
    return path


class AnalyzeLibrpaAbfOverlapTest(unittest.TestCase):
    def fixture(self, S=None, V=None, atom_naux=(2,), type_sizes=(2,), multiplicities=(1,)):
        directory = tempfile.TemporaryDirectory()
        write_basis(directory.name, type_sizes, multiplicities)
        S = np.eye(sum(atom_naux), dtype=np.complex128) if S is None else S
        V = np.eye(sum(atom_naux), dtype=np.complex128) if V is None else V
        write_raw(directory.name, S, atom_naux=atom_naux)
        write_coulomb(directory.name, V, atom_naux=atom_naux)
        return directory

    def test_reports_raw_and_whitened_psd(self):
        S = np.array([[2.0, 0.3j], [-0.3j, 1.0]], dtype=np.complex128)
        V = np.array([[4.0, 0.5], [0.5, 3.0]], dtype=np.complex128)
        with self.fixture(S, V) as directory:
            result = analyzer.analyze_directory(directory)
        point = result["q_points"][0]
        self.assertTrue(point["raw_s"]["psd"])
        self.assertTrue(point["whitened_v"]["psd"])
        self.assertEqual(point["whitened_v"]["rank"], 2)
        self.assertEqual(point["raw_s"]["hermitian_max_residual"], 0.0)
        self.assertGreater(point["whitening"]["cutoff"], 0.0)
        self.assertGreater(point["whitening"]["condition_number"], 1.0)

    def test_one_type_two_atoms_passes_declared_total_check(self):
        atom_naux = (2, 2)
        with self.fixture(atom_naux=atom_naux, type_sizes=(2,), multiplicities=(2,)) as directory:
            result = analyzer.analyze_directory(directory)
        point = result["q_points"][0]
        self.assertEqual(result["naux"], 4)
        self.assertEqual(result["natom"], 2)
        self.assertTrue(point["raw_s"]["psd"])
        self.assertTrue(point["whitened_v"]["psd"])
        self.assertEqual(point["whitened_v"]["rank"], 4)

    def test_two_types_multiple_atoms_passes_declared_total_check(self):
        atom_naux = (1, 1, 2, 2, 2)
        with self.fixture(atom_naux=atom_naux, type_sizes=(1, 2), multiplicities=(2, 3)) as directory:
            result = analyzer.analyze_directory(directory)
        point = result["q_points"][0]
        self.assertEqual(result["naux"], 8)
        self.assertEqual(result["natom"], 5)
        self.assertTrue(point["raw_s"]["psd"])
        self.assertTrue(point["whitened_v"]["psd"])
        self.assertEqual(point["whitened_v"]["rank"], 8)

    def test_reports_non_psd_raw_and_whitened_v(self):
        S = np.array([[1.0, 0.0], [0.0, -0.1]], dtype=np.complex128)
        V = np.array([[-2.0, 0.0], [0.0, 1.0]], dtype=np.complex128)
        with self.fixture(S, V) as directory:
            result = analyzer.analyze_directory(directory)
        point = result["q_points"][0]
        self.assertFalse(point["raw_s"]["psd"])
        self.assertFalse(point["whitened_v"]["psd"])
        self.assertEqual(point["whitened_v"]["rank"], 1)

    def test_rejects_nonhermitian_raw(self):
        S = np.array([[1.0, 1.0], [0.0, 1.0]], dtype=np.complex128)
        with self.fixture(S) as directory:
            with self.assertRaisesRegex(analyzer.AnalysisError, "not Hermitian"):
                analyzer.analyze_directory(directory)

    def test_rejects_missing_coulomb_pair(self):
        with tempfile.TemporaryDirectory() as directory:
            write_basis(directory)
            write_raw(directory, np.eye(2))
            write_coulomb(directory, np.eye(2), omit_pair=0)
            with self.assertRaisesRegex(analyzer.AnalysisError, "missing Coulomb pairs"):
                analyzer.analyze_directory(directory)

    def test_rejects_coulomb_metadata_mismatch(self):
        with tempfile.TemporaryDirectory() as directory:
            write_basis(directory)
            write_raw(directory, np.eye(2))
            write_coulomb(directory, np.eye(1), atom_naux=(1,))
            with self.assertRaisesRegex(analyzer.AnalysisError, "Coulomb/raw metadata mismatch"):
                analyzer.analyze_directory(directory)

    def test_rejects_basis_declared_total_raw_naux_mismatch(self):
        with tempfile.TemporaryDirectory() as directory:
            write_basis(directory, type_sizes=(2,), multiplicities=(2,))
            write_raw(directory, np.eye(2))
            write_coulomb(directory, np.eye(2))
            with self.assertRaisesRegex(analyzer.AnalysisError, "does not match basis declared_total"):
                analyzer.analyze_directory(directory)

    def test_rejects_coulomb_trailing_bytes(self):
        with tempfile.TemporaryDirectory() as directory:
            write_basis(directory)
            write_raw(directory, np.eye(2))
            path = write_coulomb(directory, np.eye(2))
            with open(path, "ab") as handle:
                handle.write(b"unreferenced")
            with self.assertRaisesRegex(analyzer.AnalysisError, "trailing or unreferenced"):
                analyzer.analyze_directory(directory)

    def test_rejects_non_contiguous_coulomb_payload_gap(self):
        with tempfile.TemporaryDirectory() as directory:
            write_basis(directory, type_sizes=(1,), multiplicities=(2,))
            write_raw(directory, np.eye(2), atom_naux=(1, 1))
            path = write_coulomb(directory, np.eye(2), atom_naux=(1, 1), gap_after_first=7)
            with self.assertRaisesRegex(analyzer.AnalysisError, "gap between Coulomb payloads"):
                analyzer.analyze_directory(directory)
            self.assertGreater(os.path.getsize(path), 0)

    def test_rejects_bad_raw_marker_and_truncated_payload(self):
        with tempfile.TemporaryDirectory() as directory:
            write_basis(directory)
            path = write_raw(directory, np.eye(2), marker=0)
            write_coulomb(directory, np.eye(2))
            with self.assertRaisesRegex(analyzer.AnalysisError, "bad raw marker"):
                analyzer.analyze_directory(directory)
            write_raw(directory, np.eye(2))
            with open(path, "r+b") as handle:
                handle.truncate(os.path.getsize(path) - 1)
            with self.assertRaisesRegex(analyzer.AnalysisError, "truncated raw payload"):
                analyzer.analyze_directory(directory)

    def test_cli_json(self):
        with self.fixture() as directory:
            completed = subprocess.run([sys.executable, SCRIPT, directory], check=False,
                                       stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                       universal_newlines=True)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertTrue(json.loads(completed.stdout)["q_points"][0]["raw_s"]["psd"])


if __name__ == "__main__":
    unittest.main()
