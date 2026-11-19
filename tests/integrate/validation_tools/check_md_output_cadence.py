"""Check numbered MD structure frames against unequal output cadences."""

import argparse
from pathlib import Path
import shutil
import subprocess
import tempfile


def run_case(executable, fixture, out_stru, expected):
    with tempfile.TemporaryDirectory(prefix="abacus-md-cadence-") as scratch:
        work = Path(scratch)
        shutil.copyfile(fixture / "STRU", work / "STRU")
        input_text = (fixture / "INPUT").read_text()
        input_text = input_text.replace("out_stru 0", "out_stru " + str(out_stru))
        (work / "INPUT").write_text(input_text)
        result = subprocess.run(
            [executable], cwd=str(work), stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, universal_newlines=True, timeout=120,
            check=False)
        if result.returncode:
            raise AssertionError("ABACUS failed:\n" + result.stdout)

        output = work / "OUT.cadence"
        frames = {path.name for path in output.glob("STRU[0-9]*")}
        if frames != expected:
            raise AssertionError("out_stru={} produced {}, expected {}\n{}".format(
                out_stru, sorted(frames), sorted(expected), result.stdout))
        restart = output / "Restart_md.txt"
        if not restart.is_file() or restart.stat().st_size == 0:
            raise AssertionError("Missing MD restart for out_stru={}".format(out_stru))
        restart_step = restart.read_text().splitlines()[0]
        if restart_step != "3":
            raise AssertionError("out_stru={} restart step is {}, expected 3".format(
                out_stru, restart_step))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("executable")
    parser.add_argument("fixture", type=Path)
    args = parser.parse_args()
    run_case(args.executable, args.fixture, 0, {"STRU1", "STRU4"})
    run_case(args.executable, args.fixture, 1,
             {"STRU1", "STRU3", "STRU4", "STRU5"})


if __name__ == "__main__":
    main()
