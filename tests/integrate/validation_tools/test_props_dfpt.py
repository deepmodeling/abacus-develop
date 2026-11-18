"""Regression checks for DFPT phonon frequency collection."""

import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import unittest


TOOLS = Path(__file__).resolve().parent

# synthetic ELEMENT NATOM table header prepended to every collector log









LOG_HEADER = (
    "                                ELEMENT NATOM\n"
    "                                Si      2\n"
    " ----------------------------------------------------------------\n"
)

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
    def run_collector(self, log_text, header=LOG_HEADER):
        with tempfile.TemporaryDirectory(prefix="abacus-dfpt-collector-") as scratch:
            case = Path(scratch)
            full_log = header + log_text if header else log_text
            (case / "log.txt").write_text(full_log, encoding="utf-8")
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
        # only 3 of 6 modes (log atom table says 2 atoms → 6 modes)
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

    # --- log-derived atom count ---

    def test_multi_species_log_header_passes(self):
        """Si/Ge element table (two species, 3 atoms total)."""
        header = (
            "                                ELEMENT NATOM\n"
            "                                Si      1\n"
            "                                Ge      2\n"
            " ----------------------------------------------------------------\n"
        )
        # 3 atoms → 9 modes; the 6-mode GOOD_FREQ_BLOCK should fail
        result, output = self.run_collector(GOOD_FREQ_BLOCK, header=header)
        self.assertNotEqual(result.returncode, 0)

    def test_no_atom_table_falls_back_to_first_block(self):
        """Log without ELEMENT NATOM table → expected_total=0 → first-block count."""
        result, output = self.run_collector(GOOD_FREQ_BLOCK, header=None)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(output.strip().splitlines()), 6)

    def test_garbage_atom_table_falls_back(self):
        """Garbage after ELEMENT NATOM (no valid numbers) → fallback."""
        header = "                                ELEMENT NATOM\n                                Garbage\n"
        result, output = self.run_collector(GOOD_FREQ_BLOCK, header=header)
        self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == "__main__":
    unittest.main()
