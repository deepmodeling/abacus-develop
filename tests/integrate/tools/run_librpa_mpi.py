"""Run the producer regression on two MPI ranks in an isolated case tree."""

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile


def ignore_outputs(directory, names):
    if Path(directory).name == "producer_reference":
        return []
    return [name for name in names if name.startswith("OUT.") or name in ("log.txt", "result.out")]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--abacus", required=True)
    parser.add_argument("--cases", required=True, type=Path)
    args = parser.parse_args()
    source = args.cases.resolve()
    cases = (source / "CASES_LIBRPA_PRODUCER.txt").read_text().split()
    with tempfile.TemporaryDirectory(prefix="abacus-librpa-mpi-") as directory:
        root = Path(directory)
        work = root / "08_RI"
        work.mkdir()
        (root / "PP_ORB").symlink_to(source.parent / "PP_ORB", target_is_directory=True)
        (root / "integrate").symlink_to(source.parent / "integrate", target_is_directory=True)
        shutil.copy2(source / "CASES_LIBRPA_PRODUCER.txt", work)
        for case in cases:
            target = work / case
            shutil.copytree(source / case, target, ignore=ignore_outputs)
            # Binary rank partitions differ from the single-rank reference.
            # Stable text remains compared; KS data use the same-run NAO writer.
            manifest_path = target / "librpa_producer_manifest.json"
            manifest = json.loads(manifest_path.read_text())
            for entry in manifest["required"]:
                if entry["kind"] in ("lri_coeff_v1", "shrink_sinvs_v1", "coulomb_v1"):
                    entry["reference"] = False
            manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
        env = dict(os.environ, OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
        subprocess.run(["bash", "../integrate/Autotest.sh", "-a", args.abacus, "-n", "2", "-j", "1",
                        "-f", "CASES_LIBRPA_PRODUCER.txt"], cwd=work, env=env, check=True, timeout=600)
        for case in cases:
            out = work / case / "OUT.librpa"
            if not (out / "Cs_1.dat").is_file() or not list(out.glob("V_full_*_r1.dat")):
                raise RuntimeError("{} did not produce rank-1 output".format(case))

        # Inject a rank-0-only and a rank-1-only file-open failure. Both must
        # terminate the MPI job promptly, instead of hanging in a collective.
        for name in ("bz_sample.txt", "Cs_1.dat"):
            case = root / ("failure_" + name.replace(".", "_"))
            shutil.copytree(source / cases[0], case, ignore=ignore_outputs)
            # These cases are one directory shallower than the positive runs.
            inp = case / "INPUT"
            inp.write_text(inp.read_text().replace("../../PP_ORB", str(source.parent / "PP_ORB")))
            (case / "OUT.librpa" / name).mkdir(parents=True)
            run = subprocess.run(["mpirun", "-np", "2", args.abacus], cwd=case, env=env,
                                 stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=90)
            if run.returncode == 0 or "RPA producer output failed:" not in run.stdout:
                raise RuntimeError("{} did not report a communicator-wide failure:\n{}".format(name, run.stdout))
            print("MPI output failure check: PASS ({})".format(name))
    print("Two-rank LibRPA producer regression: PASS")


if __name__ == "__main__":
    main()
