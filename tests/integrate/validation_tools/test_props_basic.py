"""Regression checks for basic property collection and its shell entry points."""

from importlib.util import find_spec
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


TOOLS = Path(__file__).resolve().parent
INTEGRATE = TOOLS.parent


class ForceCollectionTest(unittest.TestCase):
    def run_collector(self, calculation, atom_count, stru=None, log_force="", mode2=False,
                      stale_stru=None, log_stress="", has_force=1, has_stress=0,
                      input_text="INPUT_PARAMETERS\n"):
        with tempfile.TemporaryDirectory(prefix="abacus-force-collector-") as scratch:
            case = Path(scratch)
            (case / "OUT.autotest").mkdir()
            if stru is not None:
                structure_dir = case / "OUT.autotest/STRUC" if mode2 else case / "OUT.autotest"
                structure_dir.mkdir(exist_ok=True)
                (structure_dir / "STRU_FINAL").write_text(stru, encoding="utf-8")
            if stale_stru is not None:
                (case / "OUT.autotest/STRU_FINAL").write_text(stale_stru, encoding="utf-8")
            (case / "running.log").write_text("ETOT_ 1.0\n" + log_force + log_stress, encoding="utf-8")
            (case / "INPUT").write_text(input_text, encoding="utf-8")
            command = (
                'source "$1"; source "$2"; calculation="$3"; natom="$4"; '
                'is_lr=0; has_force="$5"; has_stress="$6"; esolver_type="$7"; '
                'nspin=1; running_path=running.log; '
                'props_result_file=result.out; run_basic_props'
            )
            result = subprocess.run(
                ["bash", "-c", command, "force-test", str(TOOLS / "props_common.sh"),
                 str(TOOLS / "props_basic.sh"), calculation, str(atom_count),
                 str(has_force), str(has_stress),
                 "ofdft" if "esolver_type ofdft" in input_text else ""],
                cwd=case, capture_output=True, text=True, check=False,
            )
            output = (case / "result.out").read_text(encoding="utf-8") if (case / "result.out").exists() else ""
            return result, output

    def test_md_stress_uses_virial_and_restores_external_pressure(self):
        virial = " MD STEP: 5\n VIRIAL (kbar):\n 1 2 3\n 4 5 6\n 7 8 9\n"
        result, output = self.run_collector(
            "md", 2, stru="# MD step 5\n", log_stress=virial, has_force=0, has_stress=1,
            input_text="INPUT_PARAMETERS\npress1 10\npress2 20\npress3 30\n")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("totalstressref 105.000000", output)

        result, output = self.run_collector(
            "md", 2, stru="# MD step 5\n", log_stress=virial, has_force=0, has_stress=1,
            input_text="INPUT_PARAMETERS\nesolver_type ofdft\npress1 10\npress2 20\npress3 30\n")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("totalstressref 45.000000", output)

        result, output = self.run_collector("md", 2, has_force=0, has_stress=1)
        self.assertNotEqual(result.returncode, 0)
        self.assertNotIn("totalstressref", output)

        earlier = virial.replace("MD STEP: 5", "MD STEP: 4")
        result, output = self.run_collector(
            "md", 2, stru="# MD step 5\n", log_stress=earlier, has_force=0, has_stress=1)
        self.assertNotEqual(result.returncode, 0)
        self.assertNotIn("totalstressref", output)

    def test_non_md_stress_still_uses_total_stress(self):
        table = "#TOTAL-STRESS (kbar)#\nStress_x Stress_y Stress_z\n1 2 3\n4 5 6\n7 8 9\n"
        result, output = self.run_collector("scf", 2, log_stress=table, has_force=0, has_stress=1)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("totalstressref 45.000000", output)

    def test_md_virial_matches_legacy_stress_matrix(self):
        rows = (
            "120.2864521128 3.4986227089 -4.0502591177\n"
            "3.4986227089 120.3390863152 2.3951443025\n"
            "-4.0502591177 2.3951443025 120.3180820327\n"
        )
        md_log = " MD STEP: 6\n VIRIAL (kbar):\n" + rows
        md_result, md_output = self.run_collector(
            "md", 2, stru="# MD step 6\n", log_stress=md_log,
            has_force=0, has_stress=1)
        legacy_log = "#TOTAL-STRESS (kbar)#\n---\nStress_x Stress_y Stress_z\n---\n" + rows
        legacy_result, legacy_output = self.run_collector(
            "scf", 2, log_stress=legacy_log, has_force=0, has_stress=1)
        self.assertEqual(md_result.returncode, 0, md_result.stderr)
        self.assertEqual(legacy_result.returncode, 0, legacy_result.stderr)
        md_stress = next(line for line in md_output.splitlines() if line.startswith("totalstressref "))
        legacy_stress = next(line for line in legacy_output.splitlines() if line.startswith("totalstressref "))
        self.assertEqual(md_stress, legacy_stress)

    def test_final_stru_requires_all_force_rows(self):
        good = (
            "ATOMIC_POSITIONS\nCartesian_angstrom\n\nH\n0.0\n2\n"
            "0 0 0 m 1 1 1 v 0 0 0 f 1 2 3\n"
            "0 0 0 m 1 1 1 v 0 0 0 f 4 5 6\n"
        )
        result, output = self.run_collector("md", 2, good)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("totalforceref 21.000000", output)
        result, output = self.run_collector("md", 2, good, mode2=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("totalforceref 21.000000", output)
        announcement = "Write file: OUT.autotest/STRUC/STRU_FINAL\n"
        stale = good.replace("f 1 2 3", "f 9 9 9")
        result, output = self.run_collector("md", 2, good, log_force=announcement,
                                            mode2=True, stale_stru=stale)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("totalforceref 21.000000", output)

        for stru in (good.split("0 0 0 m")[0], good.rsplit("0 0 0 m", 1)[0],
                     good.replace("f 4 5 6", "")):
            with self.subTest(stru=stru):
                result, output = self.run_collector("md", 2, stru)
                self.assertNotEqual(result.returncode, 0)
                self.assertNotIn("totalforceref", output)

    def test_missing_or_incomplete_legacy_force_fails(self):
        table = "TOTAL-FORCE\nAtoms Force_x Force_y Force_z\nH1 1 2 3\nH2 4 5 6\n"
        result, output = self.run_collector("scf", 2, log_force=table)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("totalforceref 21.000000", output)
        for calculation, count, log in (("md", 2, ""), ("scf", 2, table.split("H2")[0]),
                                        ("scf", 0, table)):
            with self.subTest(calculation=calculation, count=count, log=log):
                result, output = self.run_collector(calculation, count, log_force=log)
                self.assertNotEqual(result.returncode, 0)
                self.assertNotIn("totalforceref", output)

    def test_run_check_propagates_collector_failure(self):
        script = (TOOLS / "run_check.sh").read_text(encoding="utf-8")
        for mode in ("normal", "debug", "reference"):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory(prefix="abacus-run-check-") as scratch:
                root = Path(scratch)
                case = root / "tests/04_FF/01_case"
                tools = root / "tests/integrate/validation_tools"
                case.mkdir(parents=True)
                tools.mkdir(parents=True)
                (case / "run_check.sh").write_text(script, encoding="utf-8")
                failing_collector = '#!/bin/bash\necho partial >"$1"\nexit 1\n'
                (tools / "catch_properties.sh").write_text(failing_collector, encoding="utf-8")
                (tools / "catch_properties.sh").chmod(0o755)
                (case / "validation_tools").mkdir()
                (case / "validation_tools/catch_properties.sh").write_text(failing_collector, encoding="utf-8")
                (case / "validation_tools/catch_properties.sh").chmod(0o755)
                (root / "tests/integrate/general_info").write_text(
                    "EXEC /bin/true\nCHECKACCURACY 6\nNUMBEROFPROCESS 1\n", encoding="utf-8")
                bin_dir = root / "bin"
                bin_dir.mkdir()
                (bin_dir / "mpirun").write_text("#!/bin/bash\nmkdir OUT.autotest\n", encoding="utf-8")
                (bin_dir / "mpirun").chmod(0o755)
                (case / "result.ref").write_text("original\n", encoding="utf-8")
                env = os.environ.copy()
                env["PATH"] = str(bin_dir) + os.pathsep + env["PATH"]
                args = [] if mode == "normal" else ["debug" if mode == "debug" else "ref"]
                result = subprocess.run(["bash", "run_check.sh"] + args, cwd=case, env=env,
                                        capture_output=True, text=True, check=False)
                self.assertNotEqual(result.returncode, 0, result.stdout)
                self.assertEqual((case / "result.ref").read_text(encoding="utf-8"), "original\n")
                self.assertTrue((case / "OUT.autotest").is_dir())

    def test_reference_generation_preserves_existing_reference(self):
        with tempfile.TemporaryDirectory(prefix="abacus-autotest-") as scratch:
            root = Path(scratch) / "tests"
            integrate = root / "integrate"
            integrate.mkdir(parents=True)
            shutil.copy2(INTEGRATE / "Autotest.sh", integrate / "Autotest.sh")
            shutil.copytree(TOOLS, integrate / "validation_tools")
            suite = root / "04_FF"
            suite.mkdir()
            (suite / "CASES_CPU.txt").write_text("01_case\n", encoding="utf-8")
            case = suite / "01_case"
            case.mkdir()
            (case / "INPUT").write_text("calculation md\ncal_force 1\n", encoding="utf-8")
            (case / "result.ref").write_text("original\n", encoding="utf-8")
            fake_abacus = Path(scratch) / "fake_abacus.sh"
            fake_abacus.write_text(
                "#!/bin/bash\nmkdir -p OUT.autotest\n"
                'printf " TOTAL ATOM NUMBER = 2\\n!FINAL_ETOT_IS 1 eV\\n" > OUT.autotest/running_md.log\n',
                encoding="utf-8")
            fake_abacus.chmod(0o755)
            result = subprocess.run(
                ["bash", str(integrate / "Autotest.sh"), "-a", str(fake_abacus),
                 "-n", "0", "-o", "1", "-g"],
                cwd=suite, capture_output=True, text=True, check=False)
            self.assertNotEqual(result.returncode, 0, result.stdout)
            self.assertIn("Reference generation failed", result.stdout)
            self.assertEqual((case / "result.ref").read_text(encoding="utf-8"), "original\n")
            self.assertFalse((case / "result.ref.new").exists())


class PdosCollectionTest(unittest.TestCase):
    @unittest.skipUnless(find_spec("numpy"), "PDOS validation requires NumPy")
    def test_missing_or_corrupt_pdos_fails_through_shell_entry(self):
        fixture = INTEGRATE.parent / "03_NAO_multik/scf_out_dos_spin2"
        with tempfile.TemporaryDirectory() as scratch:
            case = Path(scratch)
            shutil.copy(fixture / "INPUT", case)
            shutil.copy(fixture / "dos.ref", case)
            out = case / "OUT.autotest"
            out.mkdir()
            (out / "running_scf.log").write_text(" TOTAL ATOM NUMBER = 2\nETOT_ -1.0\n", encoding="utf-8")

            for contents, error in ((None, "FileNotFoundError"), (b"invalid NPY", "ValueError")):
                with self.subTest(error=error):
                    if contents is not None:
                        (out / "pdos.npy").write_bytes(contents)
                    result = subprocess.run(["bash", "-e", str(TOOLS / "catch_properties.sh"), "result.out"],
                                            cwd=case, capture_output=True, text=True)
                    self.assertIn(error, result.stderr)
                    self.assertNotEqual(result.returncode, 0, result.stderr)


if __name__ == "__main__":
    unittest.main()
