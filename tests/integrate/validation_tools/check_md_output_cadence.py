"""Check MD structure output and restart cadences independently."""

import argparse
from pathlib import Path
import shutil
import subprocess
import tempfile


def run_case(executable, fixture, out_stru, out_freq_ion, out_app_flag, expected,
             thermo_freq=1, legacy_stru_dir=False):
    with tempfile.TemporaryDirectory(prefix="abacus-md-cadence-") as scratch:
        work = Path(scratch)
        shutil.copyfile(fixture / "STRU", work / "STRU")
        input_text = (fixture / "INPUT").read_text()
        input_text = input_text.replace("out_stru 0", "out_stru " + str(out_stru))
        input_text = input_text.replace("out_freq_ion 2", "out_freq_ion " + str(out_freq_ion))
        input_text = input_text.replace("out_app_flag 0", "out_app_flag " + str(out_app_flag))
        input_text += "\ncal_stress 1\n"
        input_text += "out_thermo_freq {}\n".format(thermo_freq)
        (work / "INPUT").write_text(input_text)
        if legacy_stru_dir:
            (work / "OUT.cadence" / "STRU").mkdir(parents=True)
        result = subprocess.run(
            [executable], cwd=str(work), stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, universal_newlines=True,
            check=False)
        if result.returncode:
            raise AssertionError("ABACUS failed:\n" + result.stdout)

        output = work / "OUT.cadence"
        if (output / "MD_dump").exists():
            raise AssertionError("Removed MD_dump was written")
        running_log = (output / "running_md.log").read_text()
        expected_thermo = 0 if thermo_freq == 0 else len(range(0, 5, thermo_freq))
        if running_log.count("VIRIAL (kbar):") != expected_thermo:
            raise AssertionError("Unexpected virial output cadence")
        if running_log.count("MD TOTAL STRESS (kbar):") != expected_thermo:
            raise AssertionError("Unexpected total stress output cadence")
        if "#TOTAL-STRESS" in running_log or "#TOTAL-FORCE" in running_log:
            raise AssertionError("Solver wrote per-step force or stress outside MD output")
        structure_dir = output / "STRUC" if out_app_flag == 2 else output
        frames = {path.name for path in structure_dir.glob("STRU*") if path.is_file()}
        frames.discard("STRU_restart")
        if frames != expected:
            raise AssertionError("out_stru={}, out_freq_ion={}, out_app_flag={} produced {}, expected {}\n{}".format(
                out_stru, out_freq_ion, out_app_flag, sorted(frames), sorted(expected), result.stdout))
        restart = output / "Restart_md.txt"
        if not restart.is_file() or restart.stat().st_size == 0:
            raise AssertionError("Missing MD restart for out_stru={}".format(out_stru))
        restart_step = restart.read_text().splitlines()[0]
        if restart_step != "3":
            raise AssertionError("out_stru={} restart step is {}, expected 3".format(
                out_stru, restart_step))
        checkpoint = output / "STRU_restart"
        if not checkpoint.is_file():
            raise AssertionError("Missing STRU_restart")
        checkpoint_text = checkpoint.read_text()
        if "# MD step 4" not in checkpoint_text or " v " not in checkpoint_text or " f " in checkpoint_text:
            raise AssertionError("Restart structure lacks matching step or velocities")
        for filename in expected:
            frame_text = (structure_dir / filename).read_text()
            if out_stru == 2:
                if "_atom_site_fract_x" not in frame_text or "_cell_length_a" not in frame_text:
                    raise AssertionError("CIF structure lacks cell or atomic positions: " + filename)
                if " v " in frame_text or " f " in frame_text:
                    raise AssertionError("CIF structure unexpectedly contains STRU atom fields: " + filename)
            elif "Cartesian_angstrom" not in frame_text or " v " not in frame_text or " f " not in frame_text:
                raise AssertionError("Ordinary structure lacks coordinates, velocities or forces: " + filename)

        if legacy_stru_dir:
            legacy_frame = output / "STRU" / "STRU"
            if not legacy_frame.is_file() or "# MD step 5" not in legacy_frame.read_text():
                raise AssertionError("Legacy STRU directory blocked periodic MD output")

        if out_stru == 0:
            restart_input = input_text.replace("md_nstep 5", "md_nstep 7")
            restart_input = restart_input.replace("init_vel 0", "init_vel 1")
            restart_input += "\nmd_restart 1\n"
            (work / "INPUT").write_text(restart_input)
            continuation = subprocess.run(
                [executable], cwd=str(work), stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, universal_newlines=True,
                check=False)
            if continuation.returncode:
                raise AssertionError("MD restart failed:\n" + continuation.stdout)
            if restart.read_text().splitlines()[0] != "6":
                raise AssertionError("MD continuation did not reach restart step 6")
            if "# MD step 7" not in checkpoint.read_text():
                raise AssertionError("MD continuation did not update STRU_restart")


def run_legacy_restart_case(executable, fixture):
    with tempfile.TemporaryDirectory(prefix="abacus-md-legacy-restart-") as scratch:
        work = Path(scratch)
        shutil.copyfile(fixture / "STRU", work / "STRU")
        input_text = (fixture / "INPUT").read_text()
        (work / "INPUT").write_text(input_text)
        initial = subprocess.run([executable], cwd=work, capture_output=True,
                                 text=True, check=False)
        if initial.returncode:
            raise AssertionError("Initial MD run failed:\n" + initial.stdout + initial.stderr)
        output = work / "OUT.cadence"
        legacy_dir = output / "STRU"
        legacy_dir.mkdir()
        (output / "STRU_restart").rename(legacy_dir / "STRU_MD_3")
        restart_input = input_text.replace("md_nstep 5", "md_nstep 7")
        restart_input = restart_input.replace("init_vel 0", "init_vel 1")
        restart_input = restart_input.replace("out_stru 0", "out_stru 1")
        restart_input = restart_input.replace("out_app_flag 0", "out_app_flag 1")
        restart_input += "\nmd_restart 1\n"
        (work / "INPUT").write_text(restart_input)
        continuation = subprocess.run([executable], cwd=work, capture_output=True,
                                      text=True, check=False)
        if continuation.returncode:
            raise AssertionError("Legacy MD restart failed:\n" + continuation.stdout + continuation.stderr)
        if not (legacy_dir / "STRU").is_file() or not (output / "STRU_FINAL").is_file():
            raise AssertionError("Legacy STRU directory blocked current MD output")
        running_log = (output / "running_md.log").read_text()
        if "Read file: " not in running_log or "STRU_MD_3" not in running_log:
            raise AssertionError("Restart source was not recorded in running_md.log")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("executable")
    parser.add_argument("fixture", type=Path)
    args = parser.parse_args()
    run_case(args.executable, args.fixture, 0, 2, 0, set())
    run_case(args.executable, args.fixture, 1, 0, 0, {"STRU_FINAL"})
    run_case(args.executable, args.fixture, 1, 0, 1, {"STRU_FINAL"})
    run_case(args.executable, args.fixture, 1, 0, 2, {"STRU_FINAL"})
    run_case(args.executable, args.fixture, 1, 2, 0, {"STRU", "STRU_FINAL"})
    run_case(args.executable, args.fixture, 1, 2, 1, {"STRU", "STRU_FINAL"})
    run_case(args.executable, args.fixture, 1, 2, 2, {"STRU1", "STRU3", "STRU5", "STRU_FINAL"})
    run_case(args.executable, args.fixture, 1, 2, 1, {"STRU", "STRU_FINAL"}, thermo_freq=3)
    run_case(args.executable, args.fixture, 1, 2, 1, {"STRU", "STRU_FINAL"}, thermo_freq=0)
    run_case(args.executable, args.fixture, 1, 2, 1, {"STRU_FINAL"}, legacy_stru_dir=True)
    run_case(args.executable, args.fixture, 2, 0, 1, {"STRU_FINAL.cif"})
    run_case(args.executable, args.fixture, 2, 2, 0, {"STRU.cif", "STRU_FINAL.cif"})
    run_case(args.executable, args.fixture, 2, 2, 1, {"STRU.cif", "STRU_FINAL.cif"})
    run_case(args.executable, args.fixture, 2, 2, 2,
             {"STRU1.cif", "STRU3.cif", "STRU5.cif", "STRU_FINAL.cif"})
    run_legacy_restart_case(args.executable, args.fixture)


if __name__ == "__main__":
    main()
