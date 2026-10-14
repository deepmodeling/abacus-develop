"""Reject malformed producer data, including MPI-IO assembly corruption."""

from pathlib import Path
import struct
import tempfile
import unittest

from librpa_wavefunctions import check_ks_nao, check_velocity, read_ks_wfc
from check_librpa_producer import check_manifest, ProducerContractError


class WavefunctionContractTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)
        self.ks = self.root / "KS_wfc_0.dat"
        # Two bands with distinct complex coefficients expose reordering.
        self.data = struct.pack("<6iiq8d", -12345679, 28, 1, 1, 2, 2, 1, 36,
                                1, 0.5, 2, -0.5, 3, 0.25, 4, -0.25)
        self.ks.write_bytes(self.data)
        (self.root / "wfk1_nao.txt").write_text(
            "1 (index of k points)\n0 0 0\n2 (number of bands)\n2 (number of orbitals)\n"
            "1 (band)\n0 (Ry)\n1 (occupation)\n1 0.5 2 -0.5\n"
            "2 (band)\n1 (Ry)\n0 (occupation)\n3 0.25 4 -0.25\n")

    def tearDown(self):
        self.directory.cleanup()

    def test_valid_ks_matches_same_run_text(self):
        self.assertEqual(read_ks_wfc(self.ks)[0], (1, 1, 2, 2))
        check_ks_nao(self.ks, self.root, 1e-7)

    def test_band_header_dimension_order(self):
        band = self.root / "band_out.txt"
        band.write_text("1\n1\n2\n2\n")
        manifest = {"required": [{"pattern": self.ks.name, "kind": "ks_wfc_v1"}]}
        check_manifest(self.root, manifest)
        band.write_text("2\n1\n2\n1\n")
        with self.assertRaisesRegex(ProducerContractError, "dimensions disagree"):
            check_manifest(self.root, manifest)

    def test_spin_blocks_and_global_text_indices(self):
        values = struct.unpack_from("<8d", self.data, 36)
        payload = struct.pack("<16d", *(values + tuple(-v for v in values)))
        header = bytearray(self.data[:36])
        struct.pack_into("<i", header, 12, 2)
        self.ks.write_bytes(header + payload)
        source = (self.root / "wfk1_nao.txt").read_text()
        (self.root / "wfk1s1_nao.txt").write_text(source)
        (self.root / "wfk1s2_nao.txt").write_text(source.replace(
            "1 (index", "2 (index").replace("1 0.5 2 -0.5", "-1 -0.5 -2 0.5").replace(
            "3 0.25 4 -0.25", "-3 -0.25 -4 0.25"))
        check_ks_nao(self.ks, self.root, 1e-7)
        self.ks.write_bytes(header + payload[64:] + payload[:64])
        with self.assertRaisesRegex(ValueError, "disagrees"):
            check_ks_nao(self.ks, self.root, 1e-7)

    def test_truncated_header_and_payload(self):
        for data in (self.data[:20], self.data[:-1], self.data + b"x"):
            self.ks.write_bytes(data)
            with self.assertRaises(ValueError):
                read_ks_wfc(self.ks)

    def test_invalid_marker_dimension_index_and_offset(self):
        for offset, fmt, value in ((0, "i", 0), (12, "i", 3), (20, "i", 0), (24, "i", 2), (28, "q", 40)):
            data = bytearray(self.data)
            struct.pack_into("<" + fmt, data, offset, value)
            self.ks.write_bytes(data)
            with self.assertRaises(ValueError):
                read_ks_wfc(self.ks)

    def test_nonfinite_ks(self):
        data = bytearray(self.data)
        struct.pack_into("<d", data, 36, float("nan"))
        self.ks.write_bytes(data)
        with self.assertRaises(ValueError):
            read_ks_wfc(self.ks)

    def test_finite_but_scrambled_mpi_payload(self):
        data = self.data[:36] + self.data[68:] + self.data[36:68]
        self.ks.write_bytes(data)
        read_ks_wfc(self.ks)
        with self.assertRaisesRegex(ValueError, "disagrees"):
            check_ks_nao(self.ks, self.root, 1e-7)

    def test_velocity_dimensions_indices_payload_and_values(self):
        path = self.root / "velocity_matrix.txt"
        valid = "1\n1\n1\n2\n1 1 1\n1 0\n2 1 1\n2 0\n3 1 1\n3 0\n"
        path.write_text(valid)
        self.assertEqual(check_velocity(path), (1, 1, 1, 2))
        for invalid in (valid[:-4], valid.replace("2 1 1", "1 1 1"), valid.replace("3 0", "nan 0")):
            path.write_text(invalid)
            with self.assertRaises(ValueError):
                check_velocity(path)


if __name__ == "__main__":
    unittest.main()
