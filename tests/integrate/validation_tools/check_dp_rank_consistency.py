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


def read_frames(log_text):
    starts = [match.start() for match in re.finditer(r"#TOTAL ENERGY#", log_text)]
    return [read_observables(log_text[start:stop])
            for start, stop in zip(starts, starts[1:] + [len(log_text)])]


def compare_frames(label, serial_frames, parallel_frames, minimum_frames):
    if len(serial_frames) != len(parallel_frames) or len(serial_frames) < minimum_frames:
        raise AssertionError("{} has {} one-rank and {} four-rank evaluations".format(
            label, len(serial_frames), len(parallel_frames)))
    tolerances = [1e-5] + [1e-4] * 6 + [1e-3] * 9
    for step, (serial, parallel) in enumerate(zip(serial_frames, parallel_frames)):
        for index, (left, right) in enumerate(zip(serial, parallel)):
            if abs(left - right) > tolerances[index]:
                raise AssertionError("{} step {} observable {} differs: {} vs {}".format(
                    label, step, index, left, right))


def run_case(args, structure, ranks, case_dir, calculation):
    case_dir.mkdir()
    (case_dir / "STRU").write_text(structure)
    input_text = (
        "INPUT_PARAMETERS\n"
        "suffix rank_check\n"
        "calculation " + calculation + "\n"
        "esolver_type dp\n"
        "pot_file " + str(args.model.resolve()) + "\n"
        "cal_force 1\n"
        "cal_stress 1\n"
        "device cpu\n")
    if calculation != "scf":
        input_text += "relax_nmax 2\nforce_thr_ev 1e-10\nstress_thr 1e-10\n"
    (case_dir / "INPUT").write_text(input_text)
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
    log = case_dir / "OUT.rank_check" / ("running_" + calculation + ".log")
    return read_frames(log.read_text())


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
            serial_frames = run_case(args, variant, 1, Path(scratch) / (label + "_np1"), "scf")
            parallel_frames = run_case(args, variant, 4, Path(scratch) / (label + "_np4"), "scf")
            compare_frames(label, serial_frames, parallel_frames, 1)
            serial = serial_frames[0]
            print("{}: energy, 6 force components and 9 stress components agree; "
                  "energy={}, max_force={}".format(label, serial[0], max(abs(x) for x in serial[1:7])))
            if label == "original":
                baseline = serial
            elif abs(serial[0] - baseline[0]) < 1e-3:
                raise AssertionError("{} geometry update did not change the DP energy".format(label))
            if label == "moved" and max(abs(x) for x in serial[1:7]) < 1e-3:
                raise AssertionError("Coordinate update did not produce a measurable DP force")

        for calculation in ("relax", "cell-relax"):
            variant = structures["moved"]
            serial_frames = run_case(args, variant, 1,
                                     Path(scratch) / (calculation + "_np1"), calculation)
            parallel_frames = run_case(args, variant, 4,
                                       Path(scratch) / (calculation + "_np4"), calculation)
            compare_frames(calculation, serial_frames, parallel_frames, 2)
            print("{}: {} consecutive evaluations agree at one and four ranks".format(
                calculation, len(serial_frames)))


if __name__ == "__main__":
    main()
