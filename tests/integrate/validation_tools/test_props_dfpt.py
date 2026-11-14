"""Regression checks for DFPT phonon frequency collection."""

import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import unittest


TOOLS = Path(__file__).resolve().parent

# minimal STRU for 2-atom Si (the props_dfpt collector reads it for natom)
STRU_2SI = """ATOMIC_SPECIES
Si 28.0855 Si.upf

LATTICE_CONSTANT
1.8897261255

LATTICE_VECTORS
0.00000000 2.71500000 2.71500000
2.71500000 0.00000000 2.71500000
2.71500000 2.71500000 0.00000000

ATOMIC_POSITIONS
Direct

Si
0.0
2
0.00000000 0.00000000 0.00000000 1 1 1
0.25000000 0.25000000 0.25000000 1 1 1
"""

GOOD_FREQ_BLOCK = (
    " DFPT phonon frequencies at q #0 = (0.000000 0.000000 0.000000) (direct) in cm^-1:\n"
    "   mode   0 : -0.000001 cm^-1\n"
    "   mode   1 : 0.000002 cm^-1\n"
    "   mode   2 : 0.000003 cm^-1\n"
    "   mode   3 : 100.000000 cm^-1\n"
    "   mode   4 : 101.000000 cm^-1\n"
    "   mode   5 : 102.000000 cm^-1\n"
)


class DFPTFrequencyCollectorTest(unittest.TestCase):
    def run_collector(self, log_text, stru_text=STRU_2SI):
        with tempfile.TemporaryDirectory(prefix="abacus-dfpt-collector-") as scratch:
            case = Path(scratch)
            (case / "log.txt").write_text(log_text, encoding="utf-8")
            if stru_text is not None:
                (case / "STRU").write_text(stru_text, encoding="utf-8")
            command = (
                'source "$1"; esolver_type=dfpt; '
                'props_result_file=result.out; '
                '> $props_result_file; run_dfpt_props log.txt'
            )
            result = subprocess.run(
                ["bash", "-c", command, "dfpt-collector-test",
                 str(TOOLS / "props_dfpt.sh")],
                cwd=case, capture_output=True, text=True, check=False,
            )
            output = (case / "result.out").read_text(encoding="utf-8") if (case / "result.out").exists() else ""
            return result, output

    # --- successful scenarios ---

    def test_complete_block_passes(self):
        result, output = self.run_collector(GOOD_FREQ_BLOCK)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(output.strip().splitlines()), 6)
        self.assertIn("dfpt_freq_q0_m3 100.000000", output)
        self.assertIn("dfpt_freq_q0_m5 102.000000", output)

    def test_freq_and_loto_blocks_pass(self):
        loto = GOOD_FREQ_BLOCK.replace(
            "phonon frequencies", "LO-TO corrected frequencies")
        result, output = self.run_collector(GOOD_FREQ_BLOCK + loto)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(output.strip().splitlines()), 12)
        self.assertIn("dfpt_loto_q0_m0 -0.000001", output)

    # --- missing mode 0 ---

    def test_missing_mode_zero_fails(self):
        log = GOOD_FREQ_BLOCK.replace("   mode   0 :", "   mode  99 :", 1)
        result, output = self.run_collector(log)
        self.assertNotEqual(result.returncode, 0)

    def test_block_starting_at_mode_one_fails(self):
        lines = GOOD_FREQ_BLOCK.strip().splitlines()
        lines = [ln for ln in lines if "mode   0" not in ln]
        log = "\n".join(lines) + "\n"
        result, output = self.run_collector(log)
        self.assertNotEqual(result.returncode, 0)

    # --- empty blocks ---

    def test_empty_freq_block_fails(self):
        # freq header with no mode rows followed by another block
        log = (" DFPT phonon frequencies at q #0 = test\n"
               + GOOD_FREQ_BLOCK)
        result, output = self.run_collector(log)
        self.assertNotEqual(result.returncode, 0)

    def test_empty_loto_block_fails(self):
        loto_header = (" DFPT LO-TO corrected frequencies at q #0 = test\n")
        result, output = self.run_collector(GOOD_FREQ_BLOCK + loto_header)
        self.assertNotEqual(result.returncode, 0)

    # --- truncated blocks ---

    def test_truncated_freq_block_fails(self):
        # only 3 of 6 modes (STRU says 2 atoms → 6 modes)
        lines = GOOD_FREQ_BLOCK.strip().splitlines()
        log = "\n".join(lines[:4]) + "\n"  # header + modes 0,1,2
        result, output = self.run_collector(log)
        self.assertNotEqual(result.returncode, 0)

    def test_truncated_loto_block_fails(self):
        loto = GOOD_FREQ_BLOCK.replace(
            "phonon frequencies", "LO-TO corrected frequencies")
        lines = loto.strip().splitlines()
        loto_short = "\n".join(lines[:4]) + "\n"  # header + modes 0,1,2
        result, output = self.run_collector(GOOD_FREQ_BLOCK + loto_short)
        self.assertNotEqual(result.returncode, 0)

    # --- nonnumeric token ---

    def test_nonnumeric_token_fails(self):
        log = GOOD_FREQ_BLOCK.replace("100.000000", "abc")
        result, output = self.run_collector(log)
        self.assertNotEqual(result.returncode, 0)

    # --- empty log ---

    def test_empty_log_fails(self):
        result, output = self.run_collector("")
        self.assertNotEqual(result.returncode, 0)

    # --- no STRU → expected_total=0 → falls back to first-block count ---

    def test_no_stru_uses_first_block_count(self):
        result, output = self.run_collector(GOOD_FREQ_BLOCK, stru_text=None)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(output.strip().splitlines()), 6)


if __name__ == "__main__":
    unittest.main()
