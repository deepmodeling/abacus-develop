"""Regression checks for the plane-wave wavefunction property collector."""

from pathlib import Path
import subprocess
import tempfile
import unittest


TOOLS = Path(__file__).resolve().parent


class WavefunctionMagnitudeTest(unittest.TestCase):
    def run_collector(self, wavefunctions):
        with tempfile.TemporaryDirectory(prefix="abacus-wfc-collector-") as scratch:
            case = Path(scratch)
            (case / "OUT.autotest").mkdir()
            (case / "OUT.autotest/wfk1_pw.txt").write_text(wavefunctions, encoding="utf-8")
            command = (
                'source "$1"; source "$2"; has_wfc_pw=1; nspin=1; '
                'if_separate_k=0; props_result_file=result.out; run_cube_props'
            )
            result = subprocess.run(
                ["bash", "-c", command, "wfc-test", str(TOOLS / "props_common.sh"),
                 str(TOOLS / "props_cube.sh")],
                cwd=case, capture_output=True, text=True, check=False,
            )
            output = (case / "result.out").read_text(encoding="utf-8") if (case / "result.out").exists() else ""
            return result, output

    def test_magnitude_is_invariant_under_global_phase(self):
        wavefunctions = (
            "< Band 1 >\n"
            " 0.6\n"
            " 0.8 -0.6 -0.8\n"
            "< Band 1 >\n"
            "< Band 2 >\n"
            " 1.0 0.0 -1.0 0.0\n"
            "< Band 2 >\n"
        )
        result, output = self.run_collector(wavefunctions)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("Max_wfc_1 1.0000", output)
        self.assertIn("Max_wfc_2 1.0000", output)

    def test_multiple_bands_use_complex_modulus(self):
        wavefunctions = (
            "< Band 1 >\n"
            " 0.6 0.8 -0.6 -0.8\n"
            "< Band 1 >\n"
            "< Band 2 >\n"
            " 0.0 1.0 0.0 -1.0\n"
            "< Band 2 >\n"
            "< Band 3 >\n"
            " 0.8 -0.6 -0.8 0.6\n"
            "< Band 3 >\n"
        )
        result, output = self.run_collector(wavefunctions)
        self.assertEqual(result.returncode, 0, result.stderr)
        for band in (1, 2, 3):
            self.assertIn(f"Max_wfc_{band} 1.0000", output)


if __name__ == "__main__":
    unittest.main()
