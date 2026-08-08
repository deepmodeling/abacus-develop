#!/usr/bin/env python3

import csv
import struct
import tempfile
import unittest
from pathlib import Path

from sternheimer_component_hybrid import (
    ReaderV1,
    build_hybrid_campaign,
    build_hybrid_file,
    ensure_matching_basis_files,
    main,
)


MARKER = -41073291


def write_reader_v1(path, values, *, ifrequency=1, omega=0.5, pair_index=0):
    atom_naux = (2,)
    header_size = 6 * 4 + 2 * 8 + 4 + len(atom_naux) * 4 + 4 + 8
    with path.open("wb") as output:
        output.write(struct.pack("<6i", MARKER, 1, ifrequency, 2, 1, 1))
        output.write(struct.pack("<2d", omega, 0.125))
        output.write(struct.pack("<i", 1))
        output.write(struct.pack("<i", atom_naux[0]))
        output.write(struct.pack("<iq", pair_index, header_size))
        for value in values:
            output.write(struct.pack("<2d", value.real, value.imag))


class SternheimerComponentHybridTest(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)
        self.total40 = self.root / "total40.dat"
        self.component40 = self.root / "component40.dat"
        self.component50 = self.root / "component50.dat"
        self.output = self.root / "hybrid.dat"

    def tearDown(self):
        self.temporary_directory.cleanup()

    def test_builds_hybrid_payload_and_preserves_non_payload_bytes(self):
        total = [1 + 2j, 3 - 1j, -2 + 0.5j, 4 + 0j]
        component40 = [0.5 + 0j, 1 + 1j, -1 + 0j, 0.25 - 0.5j]
        component50 = [0.75 - 0.25j, 0.5 + 2j, -0.5 + 0.5j, 1 + 0j]
        write_reader_v1(self.total40, total)
        write_reader_v1(self.component40, component40)
        write_reader_v1(self.component50, component50)

        result = build_hybrid_file(self.total40, self.component40, self.component50, self.output)

        output = ReaderV1.read(self.output)
        expected = [t - c40 + c50 for t, c40, c50 in zip(total, component40, component50)]
        self.assertEqual(output.values, expected)
        total_bytes = self.total40.read_bytes()
        output_bytes = self.output.read_bytes()
        self.assertEqual(output_bytes[: output.payload_start], total_bytes[: output.payload_start])
        self.assertGreater(result["relative_component_change"], 0.0)

    def test_noop_component_replacement_is_byte_identical(self):
        total = [1 + 0j, 2 + 1j, 3 - 2j, 4 + 0j]
        component = [0.25 + 0j, -0.5 + 0.25j, 0.75 - 1j, 0 + 0j]
        write_reader_v1(self.total40, total)
        write_reader_v1(self.component40, component)
        write_reader_v1(self.component50, component)

        build_hybrid_file(self.total40, self.component40, self.component50, self.output)

        self.assertEqual(self.output.read_bytes(), self.total40.read_bytes())

    def test_rejects_frequency_metadata_mismatch(self):
        values = [1 + 0j] * 4
        write_reader_v1(self.total40, values, ifrequency=1, omega=0.5)
        write_reader_v1(self.component40, values, ifrequency=1, omega=0.5)
        write_reader_v1(self.component50, values, ifrequency=2, omega=0.75)

        with self.assertRaisesRegex(ValueError, "metadata"):
            build_hybrid_file(self.total40, self.component40, self.component50, self.output)

    def test_rejects_block_layout_mismatch(self):
        values = [1 + 0j] * 4
        write_reader_v1(self.total40, values, pair_index=0)
        write_reader_v1(self.component40, values, pair_index=0)
        write_reader_v1(self.component50, values, pair_index=1)

        with self.assertRaisesRegex(ValueError, "layout"):
            build_hybrid_file(self.total40, self.component40, self.component50, self.output)

    def test_rejects_basis_order_mismatch(self):
        basis40 = self.root / "basis40"
        basis50 = self.root / "basis50"
        basis40.write_bytes(b"auxiliary channel 0\nauxiliary channel 1\n")
        basis50.write_bytes(b"auxiliary channel 1\nauxiliary channel 0\n")

        with self.assertRaisesRegex(ValueError, "basis"):
            ensure_matching_basis_files(basis40, basis50)

    def test_builds_frequency_campaign_and_csv_manifest(self):
        run40 = self.root / "run40"
        run50 = self.root / "run50"
        output_root = self.root / "hybrids"
        run40.mkdir()
        run50.mkdir()
        (run40 / "basis_aux_out").write_bytes(b"same ordered auxiliary basis\n")
        (run50 / "basis_aux_out").write_bytes(b"same ordered auxiliary basis\n")
        for ifrequency in (1, 2):
            total = [complex(10 * ifrequency + index, index) for index in range(4)]
            component40 = [complex(ifrequency, index / 4) for index in range(4)]
            component50 = [value + complex(0.5, -0.25) for value in component40]
            write_reader_v1(
                run40 / f"v1_sternheimer_chi0_iq_1_ifreq_{ifrequency}_rank0.dat",
                total,
                ifrequency=ifrequency,
                omega=0.5 * ifrequency,
            )
            for run, values in ((run40, component40), (run50, component50)):
                write_reader_v1(
                    run / f"v1_sternheimer_component_sos_iq_1_ifreq_{ifrequency}_rank0.dat",
                    values,
                    ifrequency=ifrequency,
                    omega=0.5 * ifrequency,
                )

        manifest = build_hybrid_campaign(run40, run50, output_root, components=("sos",))

        self.assertEqual(
            sorted(path.name for path in (output_root / "sos").glob("v1_sternheimer_chi0_*.dat")),
            [
                "v1_sternheimer_chi0_iq_1_ifreq_1_rank0.dat",
                "v1_sternheimer_chi0_iq_1_ifreq_2_rank0.dat",
            ],
        )
        with manifest.open(newline="") as source:
            rows = list(csv.DictReader(source))
        self.assertEqual(len(rows), 2)
        self.assertEqual({row["component"] for row in rows}, {"sos"})
        self.assertEqual({row["ifrequency"] for row in rows}, {"1", "2"})
        self.assertEqual(len({row["basis_sha256"] for row in rows}), 1)

    def test_campaign_command_builds_requested_component(self):
        run40 = self.root / "cli40"
        run50 = self.root / "cli50"
        output_root = self.root / "cli_output"
        run40.mkdir()
        run50.mkdir()
        for run in (run40, run50):
            (run / "basis_aux_out").write_bytes(b"same basis\n")
        write_reader_v1(run40 / "v1_sternheimer_chi0_iq_1_ifreq_1_rank0.dat", [1 + 0j] * 4)
        write_reader_v1(
            run40 / "v1_sternheimer_component_qspace_iq_1_ifreq_1_rank0.dat", [0.25 + 0j] * 4
        )
        write_reader_v1(
            run50 / "v1_sternheimer_component_qspace_iq_1_ifreq_1_rank0.dat", [0.5 + 0j] * 4
        )

        self.assertEqual(
            main(["campaign", str(run40), str(run50), str(output_root), "--components", "qspace"]),
            0,
        )
        self.assertTrue((output_root / "qspace" / "v1_sternheimer_chi0_iq_1_ifreq_1_rank0.dat").is_file())

    def test_campaign_preflight_leaves_no_partial_output_on_key_mismatch(self):
        run40 = self.root / "mismatch40"
        run50 = self.root / "mismatch50"
        output_root = self.root / "mismatch_output"
        run40.mkdir()
        run50.mkdir()
        for run in (run40, run50):
            (run / "basis_aux_out").write_bytes(b"same basis\n")
        write_reader_v1(run40 / "v1_sternheimer_chi0_iq_1_ifreq_1_rank0.dat", [1 + 0j] * 4)
        write_reader_v1(
            run40 / "v1_sternheimer_component_sos_iq_1_ifreq_1_rank0.dat", [0.25 + 0j] * 4
        )

        with self.assertRaisesRegex(ValueError, "keys"):
            build_hybrid_campaign(run40, run50, output_root, components=("sos",))
        self.assertFalse(output_root.exists())


if __name__ == "__main__":
    unittest.main()
