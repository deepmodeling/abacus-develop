#!/usr/bin/env python3
"""Unit tests for the ABACUS -> LibRPA producer contract checker."""

import json
import math
import struct
import tempfile
import unittest
from pathlib import Path

from check_librpa_producer import ProducerContractError, check_manifest


def _write_stru(path: Path, symmetry_rows: int = 0, spin_symmetry_rows: int = 0) -> None:
    lines = [
        "1 0 0\n",
        "0 1 0\n",
        "0 0 1\n",
        "6.283185307179586 0 0\n",
        "0 6.283185307179586 0\n",
        "0 0 6.283185307179586\n",
        "1\n",
        "0 0 0 1\n",
    ]
    if symmetry_rows:
        lines.append(f"{symmetry_rows} row\n")
        for _ in range(symmetry_rows):
            lines.append("1 0 0 0 1 0 0 0 1 0 0 0\n")
    if spin_symmetry_rows:
        lines.append("spin_symmetry 1 1\n")
        for _ in range(spin_symmetry_rows):
            lines.append("0 1 0 0 0 1 0 0 0\n")
    path.write_text("".join(lines))


def _write_coulomb(path: Path) -> None:
    # marker, iq, naux, complex flag, natom, nblocks, atom_naux, pair/offset,
    # followed by one complex<double> value.
    payload = struct.pack(
        "<6i i i q 2d",
        -20129433,
        1,
        1,
        1,
        1,
        1,
        1,
        0,
        40,
        1.25,
        -0.5,
    )
    path.write_bytes(payload)


def _write_chi0(path: Path) -> None:
    payload = struct.pack(
        "<6i 2d i i i q 2d",
        -41073291,
        1,
        1,
        1,
        1,
        1,
        0.25,
        0.5,
        1,
        1,
        0,
        60,
        0.75,
        -0.25,
    )
    path.write_bytes(payload)


def _write_lri_coeff(path: Path, value: float = 0.5) -> None:
    # marker, atom/cell dimensions, record table, then one double payload.
    path.write_bytes(
        struct.pack(
            "<3i 2q 5i d q d",
            -10267453,
            1,
            0,
            1,
            1,
            0,
            0,
            0,
            0,
            0,
            1.25,
            64,
            value,
        )
    )


def _write_shrink_sinvs(path: Path, value: float = 0.75) -> None:
    # marker, one record, then one double payload.
    path.write_bytes(
        struct.pack(
            "<2i 7i d q d",
            -30241621,
            1,
            1,
            1,
            0,
            0,
            0,
            0,
            0,
            0.125,
            52,
            value,
        )
    )


def _write_complex_wfc(path: Path, coefficient_count: int = 4) -> None:
    coefficients = " ".join("{:.1f}".format(float(index)) for index in range(coefficient_count))
    path.write_text(
        "1 (index of k points)\n"
        "0.0 0.0 0.0\n"
        "1 (number of bands)\n"
        "2 (number of orbitals)\n"
        "1 (band)\n"
        "-0.5 (Ry)\n"
        "1.0 (Occupations)\n"
        "{}\n".format(coefficients)
    )


def _write_strict_2d_head(path: Path) -> None:
    path.write_text(
        "# ABACUS reader-v1 strict 2D Coulomb head normalization\n"
        "version = 1\n"
        "area_parallel_bohr2 = 10.0\n"
        "multipole_norm_squared = 2.0\n"
        "strict_2d_coulomb_head_coefficient = 3.0\n"
        "strict_2d_sheet_to_raw_scale = 4.0\n"
    )


class ProducerContractTests(unittest.TestCase):
    def test_checks_required_files_and_symmetry_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "OUT.librpa"
            out.mkdir()
            _write_stru(out / "stru_out.txt", symmetry_rows=2)
            (out / "band_out.txt").write_text("1\n1\n2\n2\n0.0\n1 1\n1 1.0 0.0 0.0\n")
            manifest = {
                "output_dir": "OUT.librpa",
                "required": [
                    {"pattern": "stru_out.txt", "kind": "stru", "symmetry_rows": 2},
                    {"pattern": "band_out.txt", "kind": "band"},
                ],
            }
            self.assertEqual(check_manifest(root, manifest), [])

    def test_requires_complete_spin_symmetry_block(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "OUT.librpa"
            out.mkdir()
            _write_stru(out / "stru_out.txt", symmetry_rows=2)
            manifest = {
                "output_dir": "OUT.librpa",
                "required": [
                    {
                        "pattern": "stru_out.txt",
                        "kind": "stru",
                        "symmetry_rows": 2,
                        "spin_symmetry_rows": 2,
                    }
                ],
            }
            with self.assertRaises(ProducerContractError):
                check_manifest(root, manifest)
            _write_stru(out / "stru_out.txt", symmetry_rows=2, spin_symmetry_rows=2)
            self.assertEqual(check_manifest(root, manifest), [])

    def test_missing_required_file_is_an_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "OUT.librpa").mkdir()
            manifest = {
                "output_dir": "OUT.librpa",
                "required": [{"pattern": "stru_out.txt", "kind": "stru"}],
            }
            with self.assertRaises(ProducerContractError):
                check_manifest(root, manifest)

    def test_checks_reader_v1_binary_headers(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "OUT.librpa"
            out.mkdir()
            _write_coulomb(out / "v1_coulomb_full_iq_1_rank0.dat")
            _write_chi0(out / "v1_sternheimer_chi0_iq_1_ifreq_1_rank0.dat")
            manifest = {
                "output_dir": "OUT.librpa",
                "required": [
                    {"pattern": "v1_coulomb_full_iq_*.dat", "kind": "coulomb_v1"},
                    {"pattern": "v1_sternheimer_chi0_*.dat", "kind": "chi0_v1"},
                ],
            }
            self.assertEqual(check_manifest(root, manifest), [])

    def test_checks_strict_2d_head_metadata(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "OUT.librpa"
            out.mkdir()
            _write_strict_2d_head(out / "librpa_2d_coulomb_head.txt")
            manifest = {
                "output_dir": "OUT.librpa",
                "required": [{"pattern": "librpa_2d_coulomb_head.txt", "kind": "strict_2d_head"}],
            }
            self.assertEqual(check_manifest(root, manifest), [])
            (out / "librpa_2d_coulomb_head.txt").write_text("version = 1\n")
            with self.assertRaises(ProducerContractError):
                check_manifest(root, manifest)

    def test_compares_lri_and_shrink_v1_payloads(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "OUT.librpa"
            ref = root / "reference"
            out.mkdir()
            ref.mkdir()
            _write_lri_coeff(out / "v1_Cs_data_1.txt")
            _write_lri_coeff(ref / "v1_Cs_data_1.txt")
            _write_shrink_sinvs(out / "v1_shrink_sinvS_1.txt")
            _write_shrink_sinvs(ref / "v1_shrink_sinvS_1.txt")
            manifest = {
                "output_dir": "OUT.librpa",
                "reference_dir": "reference",
                "required": [
                    {
                        "pattern": "v1_Cs_data_*.txt",
                        "kind": "lri_coeff_v1",
                        "reference": True,
                    },
                    {
                        "pattern": "v1_shrink_sinvS_*.txt",
                        "kind": "shrink_sinvs_v1",
                        "reference": True,
                    },
                ],
            }
            self.assertEqual(check_manifest(root, manifest), [])
            _write_lri_coeff(out / "v1_Cs_data_1.txt", value=0.6)
            with self.assertRaises(ProducerContractError):
                check_manifest(root, manifest)

    def test_gauge_dependent_files_are_presence_only(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "OUT.librpa"
            out.mkdir()
            (out / "KS_eigenvector_0.dat").write_bytes(b"not compared\n")
            (out / "velocity_matrix").write_bytes(b"not compared\n")
            manifest = {
                "output_dir": "OUT.librpa",
                "required": [
                    {"pattern": "KS_eigenvector_*.dat", "kind": "gauge"},
                    {"pattern": "velocity_matrix", "kind": "gauge"},
                ],
            }
            self.assertEqual(check_manifest(root, manifest), [])

    def test_checks_complex_nao_wavefunction_layout(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "OUT.librpa").mkdir()
            out = root / "OUT.autotest"
            out.mkdir()
            wavefunction = out / "wfk1_nao.txt"
            _write_complex_wfc(wavefunction)
            manifest = {
                "output_dir": "OUT.librpa",
                "required": [
                    {
                        "base_dir": ".",
                        "pattern": "OUT.autotest/wfk*_nao.txt",
                        "kind": "wfc_nao",
                    }
                ],
            }
            self.assertEqual(check_manifest(root, manifest), [])
            _write_complex_wfc(wavefunction, coefficient_count=3)
            with self.assertRaises(ProducerContractError):
                check_manifest(root, manifest)

    def test_rejects_non_numeric_wavefunction_coefficients(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "OUT.librpa").mkdir()
            out = root / "OUT.autotest"
            out.mkdir()
            wavefunction = out / "wfk1_nao.txt"
            _write_complex_wfc(wavefunction)
            wavefunction.write_text(wavefunction.read_text().replace("0.0 1.0 2.0 3.0", "0.0 bad 2.0 3.0"))
            manifest = {
                "output_dir": "OUT.librpa",
                "required": [{"base_dir": ".", "pattern": "OUT.autotest/wfk*_nao.txt", "kind": "wfc_nao"}],
            }
            with self.assertRaises(ProducerContractError):
                check_manifest(root, manifest)

    def test_requires_all_declared_wavefunction_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "OUT.librpa").mkdir()
            out = root / "OUT.autotest"
            out.mkdir()
            _write_complex_wfc(out / "wfk1_nao.txt")
            manifest = {
                "output_dir": "OUT.librpa",
                "required": [
                    {"base_dir": ".", "pattern": "OUT.autotest/wfk*_nao.txt", "kind": "wfc_nao", "expected_count": 2}
                ],
            }
            with self.assertRaises(ProducerContractError):
                check_manifest(root, manifest)

    def test_checks_case_relative_vxc_reference(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "OUT.librpa").mkdir()
            output = root / "OUT.autotest"
            reference = root / "reference" / "OUT.autotest"
            output.mkdir()
            reference.mkdir(parents=True)
            (output / "vxc_out.dat").write_text("1.0000000001 2.0\n")
            (reference / "vxc_out.dat").write_text("1.0 2.0\n")
            manifest = {
                "output_dir": "OUT.librpa",
                "reference_dir": "reference",
                "required": [
                    {
                        "base_dir": ".",
                        "pattern": "OUT.autotest/vxc_out.dat",
                        "kind": "text",
                        "reference": True,
                        "abs_tol": 1.0e-8,
                    }
                ],
            }
            self.assertEqual(check_manifest(root, manifest), [])

    def test_numeric_reference_uses_tolerance(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "OUT.librpa"
            ref = root / "reference"
            out.mkdir()
            ref.mkdir()
            (out / "Cs_data_0.txt").write_text("1.0000000001 2.0\n")
            (ref / "Cs_data_0.txt").write_text("1.0 2.0\n")
            manifest = {
                "output_dir": "OUT.librpa",
                "reference_dir": "reference",
                "required": [
                    {
                        "pattern": "Cs_data_*.txt",
                        "kind": "numeric",
                        "reference": True,
                        "abs_tol": 1.0e-8,
                    }
                ],
            }
            self.assertEqual(check_manifest(root, manifest), [])

    def test_nonfinite_numeric_value_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "OUT.librpa"
            out.mkdir()
            (out / "Cs_data_0.txt").write_text("nan 1.0\n")
            manifest = {
                "output_dir": "OUT.librpa",
                "required": [{"pattern": "Cs_data_*.txt", "kind": "numeric"}],
            }
            with self.assertRaises(ProducerContractError):
                check_manifest(root, manifest)


if __name__ == "__main__":
    unittest.main()
