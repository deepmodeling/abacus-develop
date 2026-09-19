"""Unit tests for the LibRPA Sternheimer energy A/B audit."""

import json
from pathlib import Path
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
AUDIT = ROOT / "tools/periodic_sternheimer/audit_rpa_energy_ab.py"


def successful_log(energy: float) -> str:
    return (
        "header\n"
        f"| Total Sternheimer EcRPA: {energy:.16e} 0.0000000000000000e+00\n"
        "libRPA finished successfully\n"
    )


class RpaEnergyAbRoute(unittest.TestCase):
    def run_audit(self, control: str, symmetry: str, tolerance: float = 0.1):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            control_path = root / "control.out"
            symmetry_path = root / "symmetry.out"
            report_path = root / "report.json"
            control_path.write_text(control)
            symmetry_path.write_text(symmetry)
            result = subprocess.run(
                [
                    "python3", str(AUDIT), "--control", str(control_path),
                    "--symmetry", str(symmetry_path), "--output", str(report_path),
                    "--tolerance-mev", str(tolerance),
                ],
                capture_output=True,
                text=True,
            )
            report = json.loads(report_path.read_text()) if report_path.exists() else None
            return result, report

    def test_accepts_finite_equal_energies(self):
        result, report = self.run_audit(successful_log(-1.25), successful_log(-1.25))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(report["energy_gate_passed"])
        self.assertEqual(report["absolute_difference_mev"], 0.0)
        self.assertFalse(report["physical_result"])

    def test_rejects_energy_difference_above_threshold(self):
        result, report = self.run_audit(successful_log(-1.25), successful_log(-1.249), tolerance=0.1)
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(report["energy_gate_passed"])
        self.assertGreater(report["absolute_difference_mev"], 0.1)

    def test_rejects_missing_completion_marker(self):
        result, report = self.run_audit(successful_log(-1.25), "| Total Sternheimer EcRPA: -1 0\n")
        self.assertNotEqual(result.returncode, 0)
        self.assertIsNone(report)


if __name__ == "__main__":
    unittest.main()
