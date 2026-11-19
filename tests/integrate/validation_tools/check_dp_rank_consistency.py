"""Compare complete DP outputs for one and four MPI ranks.

Run with the DPA model from issue #8310, for example:
  python3 check_dp_rank_consistency.py --abacus /path/to/abacus \
      --mpi mpirun --model /path/to/frozen_model_op.pth \
      --fixture tests/04_FF/22_DP_rank_consistency/STRU
"""

import argparse
import os
from pathlib import Path
import re
import subprocess
import tempfile


NUMBER = r"[-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?"


def read_observables(log_text):
    energy = re.search(r"#TOTAL ENERGY#\s+(" + NUMBER + r")", log_text)
    force_block = re.search(
        r"#TOTAL-FORCE[^\n]*\n(.*?)#TOTAL-STRESS", log_text, re.S)
    stress_block = re.search(
        r"#TOTAL-STRESS[^\n]*\n(.*?)#TOTAL-PRESSURE", log_text, re.S)
    if not energy or not force_block or not stress_block:
        raise AssertionError("Missing DP energy, forces or stress in running log")
    forces = []
    for line in force_block.group(1).splitlines():
        if re.match(r"\s*Fe\d+\s", line):
            forces.extend(float(value) for value in re.findall(NUMBER, line)[1:])
    stress = []
    for line in stress_block.group(1).splitlines():
        values = re.findall(NUMBER, line)
        if len(values) == 3:
            stress.extend(float(value) for value in values)
    if len(forces) != 6 or len(stress) != 9:
        raise AssertionError("Incomplete DP force or stress matrix")
    return [float(energy.group(1))] + forces + stress


def run_case(args, structure, ranks, case_dir):
    case_dir.mkdir()
    (case_dir / "STRU").write_text(structure)
    (case_dir / "INPUT").write_text(
        "INPUT_PARAMETERS\n"
        "suffix rank_check\n"
        "calculation scf\n"
        "esolver_type dp\n"
        "pot_file " + str(args.model.resolve()) + "\n"
        "cal_force 1\n"
        "cal_stress 1\n"
        "device cpu\n")
    command = [args.mpi, args.mpi_np_flag, str(ranks), args.abacus]
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    result = subprocess.run(
        command, cwd=str(case_dir), env=env, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, universal_newlines=True, timeout=180,
        check=False)
    if result.returncode:
        raise AssertionError("DP run failed on {} ranks:\n{}".format(
            ranks, result.stdout))
    log = case_dir / "OUT.rank_check" / "running_scf.log"
    return read_observables(log.read_text())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--abacus", required=True)
    parser.add_argument("--mpi", required=True)
    parser.add_argument("--mpi-np-flag", default="-np")
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--fixture", type=Path, required=True)
    args = parser.parse_args()
    structure = args.fixture.read_text()
    structures = {
        "original": structure,
        "moved": structure.replace("0.5 0.5 0.5 m", "0.47 0.5 0.5 m"),
        "expanded": structure.replace("2.830954119825", "2.900000000000"),
    }
    with tempfile.TemporaryDirectory(prefix="abacus-dp-ranks-") as scratch:
        baseline = None
        for label, variant in structures.items():
            serial = run_case(args, variant, 1, Path(scratch) / (label + "_np1"))
            parallel = run_case(args, variant, 4, Path(scratch) / (label + "_np4"))
            tolerances = [1e-5] + [1e-4] * 6 + [1e-3] * 9
            for index, (left, right) in enumerate(zip(serial, parallel)):
                if abs(left - right) > tolerances[index]:
                    raise AssertionError("{} observable {} differs: {} vs {}\nnp1={}\nnp4={}".format(
                        label, index, left, right, serial, parallel))
            print("{}: energy, 6 force components and 9 stress components agree; "
                  "energy={}, max_force={}".format(label, serial[0], max(abs(x) for x in serial[1:7])))
            if label == "original":
                baseline = serial
            elif abs(serial[0] - baseline[0]) < 1e-3:
                raise AssertionError("{} geometry update did not change the DP energy".format(label))
            if label == "moved" and max(abs(x) for x in serial[1:7]) < 1e-3:
                raise AssertionError("Coordinate update did not produce a measurable DP force")


if __name__ == "__main__":
    main()
