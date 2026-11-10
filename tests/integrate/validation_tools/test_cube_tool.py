#!/usr/bin/env python3
"""Final-frame solvent metrics preserve existing first-frame cube consumers."""
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import cube_tool


def scalar_frame(values, step, voxel_x):
    return (f"Ionic_Step {step}\nscalar field\n1 0 0 0\n"
            f"2 {voxel_x} 0 0\n1 0 1 0\n1 0 0 1\n"
            f"1 1 0 0 0\n{values}\n")


class CubeToolLastFrameTest(unittest.TestCase):
    def setUp(self):
        self.scratch = tempfile.TemporaryDirectory()
        self.addCleanup(self.scratch.cleanup)
        self.path = Path(self.scratch.name) / "sol_cavity.cube"
        self.first = scalar_frame("1 2", 1, 1)
        self.last = scalar_frame("4 5", 2, 2)

    def write_frames(self, text):
        self.path.write_text(text, encoding="utf-8")

    def run_command(self, command):
        tool = Path(cube_tool.__file__).resolve()
        return subprocess.run([sys.executable, str(tool), command, str(self.path)],
                              capture_output=True, text=True, check=False)

    def test_single_frame_keeps_existing_integral(self):
        self.write_frames(self.first)
        self.assertEqual(cube_tool.integrate_cube(self.path), 3.0)
        self.assertEqual(cube_tool.integrate_last_cube(self.path), 3.0)

    def test_last_frame_uses_its_own_values_and_grid(self):
        self.write_frames(self.first + self.last)
        self.assertEqual(cube_tool.integrate_last_cube(self.path), 18.0)
        self.assertEqual(cube_tool.integrate_cube(self.path), 3.0)

    def test_trailing_whitespace_is_accepted(self):
        self.write_frames(self.first + self.last + "\n  \n")
        self.assertEqual(cube_tool.integrate_last_cube(self.path), 18.0)

    def test_cli_selects_last_only_when_requested(self):
        self.write_frames(self.first + self.last)
        old = self.run_command("integrate")
        final = self.run_command("integrate-last")
        self.assertEqual(old.returncode, 0, old.stderr)
        self.assertEqual(final.returncode, 0, final.stderr)
        self.assertEqual(old.stdout.strip(), "3")
        self.assertEqual(final.stdout.strip(), "18")

    def test_truncated_tail_does_not_fall_back_to_first_frame(self):
        self.write_frames(self.first + scalar_frame("4", 2, 2))
        with self.assertRaises(ValueError):
            cube_tool.integrate_last_cube(self.path)
        result = self.run_command("integrate-last")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(result.stdout, "")

    def test_malformed_tail_is_rejected(self):
        for tail in (scalar_frame("bad 5", 2, 2), scalar_frame("nan 5", 2, 2),
                     scalar_frame("4 5 6", 2, 2), "Ionic_Step 2\n"):
            with self.subTest(tail=tail):
                self.write_frames(self.first + tail)
                with self.assertRaises(ValueError):
                    cube_tool.integrate_last_cube(self.path)

    def test_empty_file_is_rejected(self):
        self.write_frames("")
        with self.assertRaises(ValueError):
            cube_tool.integrate_last_cube(self.path)

    def run_solvent_metric(self, frames):
        output_dir = Path(self.scratch.name) / "OUT.autotest"
        output_dir.mkdir()
        (output_dir / "sol_cavity.cube").write_text(frames, encoding="utf-8")
        (output_dir / "sol_eps.cube").write_text(frames, encoding="utf-8")
        tool = Path(cube_tool.__file__).resolve()
        script = tool.with_name("props_basic.sh")
        command = ('source "$1"\nimp_sol=0\nout_sol=1\nassume_isolated=none\n'
                   'CUBE_TOOL="$2"\nprops_result_file=result.out\nrun_basic_props_post_ml')
        return subprocess.run(["bash", "-e", "-c", command, "solvent-metric-test",
                               str(script), str(tool)], cwd=self.scratch.name,
                              capture_output=True, text=True, check=False)

    def test_solvent_metric_invokes_the_final_frame_command(self):
        result = self.run_solvent_metric(self.first + self.last)
        self.assertEqual(result.returncode, 0, result.stderr)
        metric = Path(self.scratch.name) / "result.out"
        self.assertEqual(metric.read_text(), "sol_cavity_integralref 18\n")

    def test_solvent_metric_rejects_an_incomplete_tail(self):
        result = self.run_solvent_metric(self.first + scalar_frame("4", 2, 2))
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse((Path(self.scratch.name) / "result.out").exists())


if __name__ == "__main__":
    unittest.main()
