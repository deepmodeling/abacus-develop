"""Regression checks for MD force collection and its shell entry points."""

import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


TOOLS = Path(__file__).resolve().parent
INTEGRATE = TOOLS.parent


class ForceCollectionTest(unittest.TestCase):
    def run_collector(self, calculation, atom_count, dump=None, log_force=""):
        with tempfile.TemporaryDirectory(prefix="abacus-force-collector-") as scratch:
            case = Path(scratch)
            (case / "OUT.autotest").mkdir()
            if dump is not None:
                (case / "OUT.autotest/MD_dump").write_text(dump, encoding="utf-8")
            (case / "running.log").write_text("ETOT_ 1.0\n" + log_force, encoding="utf-8")
            command = (
                'source "$1"; source "$2"; calculation="$3"; natom="$4"; '
                'is_lr=0; has_force=1; nspin=1; running_path=running.log; '
                'props_result_file=result.out; run_basic_props'
            )
            result = subprocess.run(
                ["bash", "-c", command, "force-test", str(TOOLS / "props_common.sh"),
                 str(TOOLS / "props_basic.sh"), calculation, str(atom_count)],
                cwd=case, capture_output=True, text=True, check=False,
            )
            output = (case / "result.out").read_text(encoding="utf-8") if (case / "result.out").exists() else ""
            return result, output

    def test_last_dump_frame_requires_all_force_rows(self):
        good = (
            "MDSTEP: 0\nINDEX LABEL POSITION (Angstrom) FORCE (eV/Angstrom)\n"
            "0 H 0 0 0 1 2 3\n1 H 0 0 0 4 5 6\n"
        )
        result, output = self.run_collector("md", 2, good)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("totalforceref 21.000000", output)

        for dump in (good + "MDSTEP: 1\n", good.split("1 H")[0],
                     "MDSTEP: 0\nINDEX LABEL POSITION (Angstrom)\n0 H 0 0 0\n"):
            with self.subTest(dump=dump):
                result, output = self.run_collector("md", 2, dump)
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


if __name__ == "__main__":
    unittest.main()
