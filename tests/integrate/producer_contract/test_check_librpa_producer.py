#!/usr/bin/env python3
"""Rejection-path tests for the LibRPA producer contract checker.

The checker itself lives in ``tests/integrate/validation_tools`` and keeps the
comparison-only scope of that folder. These negative cases drive its parser
paths with synthetic format fixtures and deliberately contain no physical
model: no lattice, no basis set, and no algorithm-specific data.
"""

import json
import struct
import sys
import tempfile
import unittest
from pathlib import Path

TOOLS = Path(__file__).resolve().parents[1] / "validation_tools"
sys.path.insert(0, str(TOOLS))

from check_librpa_producer import ProducerContractError, check_manifest  # noqa: E402

_LRICOEF_MARKER = -10267453


def _manifest(root, entries):
    manifest = {"output_dir": ".", "reference_dir": "reference", "required": entries}
    (root / "manifest.json").write_text(json.dumps(manifest))
    return manifest


class ProducerRejectionTests(unittest.TestCase):
    def setUp(self):
        self._temporary = tempfile.TemporaryDirectory()
        self.root = Path(self._temporary.name)

    def tearDown(self):
        self._temporary.cleanup()

    def test_missing_required_file_is_rejected(self):
        manifest = _manifest(self.root, [{"pattern": "missing_*.dat", "kind": "text"}])
        with self.assertRaisesRegex(ProducerContractError, "no producer files match"):
            check_manifest(self.root, manifest)

    def test_incomplete_strict_2d_head_metadata_is_rejected(self):
        (self.root / "head.txt").write_text("version = 1\narea_parallel_bohr2 = 1.0\n")
        manifest = _manifest(self.root, [{"pattern": "head.txt", "kind": "strict_2d_head"}])
        with self.assertRaisesRegex(ProducerContractError, "metadata keys"):
            check_manifest(self.root, manifest)

    def test_nonfinite_numeric_data_is_rejected(self):
        (self.root / "values.txt").write_text("1.0 nan\n")
        manifest = _manifest(self.root, [{"pattern": "values.txt", "kind": "numeric"}])
        with self.assertRaisesRegex(ProducerContractError, "non-finite"):
            check_manifest(self.root, manifest)

    def test_lri_payload_mismatch_is_rejected(self):
        (self.root / "reference").mkdir()
        self._write_lri(self.root / "reference" / "Cs_0.dat", payload=1.0)
        self._write_lri(self.root / "Cs_0.dat", payload=3.0)
        manifest = _manifest(
            self.root,
            [{"pattern": "Cs_*.dat", "kind": "lri_coeff_v1", "reference": True, "abs_tol": 1.0e-8}],
        )
        with self.assertRaisesRegex(ProducerContractError, "payload value differs"):
            check_manifest(self.root, manifest)

    @staticmethod
    def _write_lri(path, payload):
        """Write a minimal but format-valid LRI coefficient v1 payload."""
        header = struct.pack("<3i", _LRICOEF_MARKER, 1, 0) + struct.pack("<2q", 1, 1)
        record_size = 36
        payload_offset = 28 + record_size
        record = struct.pack("<5i", 0, 0, 0, 0, 0) + struct.pack("<d", 1.0) + struct.pack("<q", payload_offset)
        path.write_bytes(header + record + struct.pack("<d", payload))


if __name__ == "__main__":
    unittest.main()
