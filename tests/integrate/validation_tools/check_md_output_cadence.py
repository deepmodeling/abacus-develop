"""Check MD structure output and restart cadences independently."""

import argparse
from pathlib import Path
import shutil
import subprocess
import tempfile


def run_case(executable, fixture, out_stru, out_freq_ion, out_app_flag, expected):
    with tempfile.TemporaryDirectory(prefix="abacus-md-cadence-") as scratch:
        work = Path(scratch)
        shutil.copyfile(fixture / "STRU", work / "STRU")
        input_text = (fixture / "INPUT").read_text()
        input_text = input_text.replace("out_stru 0", "out_stru " + str(out_stru))
        input_text = input_text.replace("out_freq_ion 2", "out_freq_ion " + str(out_freq_ion))
        input_text = input_text.replace("out_app_flag 0", "out_app_flag " + str(out_app_flag))
        (work / "INPUT").write_text(input_text)
        result = subprocess.run(
            [executable], cwd=str(work), stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, universal_newlines=True, timeout=120,
            check=False)
        if result.returncode:
            raise AssertionError("ABACUS failed:\n" + result.stdout)

        output = work / "OUT.cadence"
        structure_dir = output / "STRUC" if out_app_flag == 2 else output
        frames = {path.name for path in structure_dir.glob("STRU*")}
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
        if "# MD step 4" not in checkpoint_text or " v " not in checkpoint_text:
            raise AssertionError("Restart structure lacks matching step or velocities")
        for filename in expected:
            frame_text = (structure_dir / filename).read_text()
            if "Cartesian_angstrom" not in frame_text or " v " in frame_text:
                raise AssertionError("Ordinary structure has wrong coordinates or velocities: " + filename)

        if out_stru == 0:
            restart_input = input_text.replace("md_nstep 5", "md_nstep 7")
            restart_input = restart_input.replace("init_vel 0", "init_vel 1")
            restart_input += "\nmd_restart 1\n"
            (work / "INPUT").write_text(restart_input)
            continuation = subprocess.run(
                [executable], cwd=str(work), stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, universal_newlines=True, timeout=120,
                check=False)
            if continuation.returncode:
                raise AssertionError("MD restart failed:\n" + continuation.stdout)
            if restart.read_text().splitlines()[0] != "6":
                raise AssertionError("MD continuation did not reach restart step 6")
            if "# MD step 7" not in checkpoint.read_text():
                raise AssertionError("MD continuation did not update STRU_restart")


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


if __name__ == "__main__":
    main()
