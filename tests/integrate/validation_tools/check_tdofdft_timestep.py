#!/usr/bin/env python3
"""Validate matched-time TDOFDFT logs; no calculations or reference writes.

Run the case at dt, dt/2 and dt/4 with the same physical end time, then
provide their running_md.log files. Optionally provide a second MPI layout
at dt. This checks numerical consistency, not a propagator convergence order.
"""

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
import re


def read_log(path, dt):
    if dt <= 0:
        raise ValueError("Timestep must be positive")
    lines = path.read_text().splitlines()
    if any("SCF IS NOT CONVERGED" in line for line in lines):
        raise ValueError("Unconverged initial optimization: {}".format(path))
    frames = []
    for index, line in enumerate(lines):
        match = re.search(r"STEP OF MOLECULAR DYNAMICS:\s*(\d+)", line)
        if match:
            step = int(match.group(1)) - 1
            if step != len(frames):
                raise ValueError("Missing or repeated MD step in {}".format(path))
            frames.append({"step": step})
        if not frames:
            continue
        if "#TOTAL-STRESS (kbar)#" in line:
            values = []
            for row in lines[index + 1:index + 12]:
                tokens = row.split()
                if len(tokens) != 3:
                    continue
                try:
                    numbers = [float(token) for token in tokens]
                except ValueError:
                    continue
                values.extend(numbers)
                if len(values) == 9:
                    break
            if len(values) != 9:
                raise ValueError("Incomplete stress tensor in {}".format(path))
            frames[-1]["stress"] = values
        if "Energy (Ry)" in line and "Potential (Ry)" in line:
            frames[-1]["energy"] = float(lines[index + 1].split()[0])
    if len(frames) < 2:
        raise ValueError("At least two MD outputs required: {}".format(path))
    for frame in frames:
        values = frame["stress"] + [frame["energy"]]
        if not all(math.isfinite(value) for value in values):
            raise ValueError("Non-finite output in {}".format(path))
    return {dt * frame["step"]: frame for frame in frames}


def difference(left, right):
    if max(left) != max(right) or not left.keys() <= right.keys():
        raise ValueError("Trajectories must cover the same end time and aligned outputs")
    component = []
    frobenius = []
    for time, frame in left.items():
        delta = [a - b for a, b in zip(frame["stress"], right[time]["stress"])]
        component.append(max(abs(value) for value in delta))
        frobenius.append(math.sqrt(sum(value * value for value in delta)))
    return {"max_component_kbar": max(component),
            "max_frobenius_kbar": max(frobenius),
            "final_component_kbar": component[-1],
            "final_energy_difference_Ry": abs(left[max(left)]["energy"]
                                              - right[max(right)]["energy"])}


def energy_drift(trajectory):
    initial = trajectory[Fraction(0)]["energy"]
    return max(abs(frame["energy"] - initial) for frame in trajectory.values())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coarse", type=Path, required=True)
    parser.add_argument("--dt", type=Fraction, default=Fraction("0.001"))
    parser.add_argument("--fine", type=Path, required=True)
    parser.add_argument("--finer", type=Path, required=True)
    parser.add_argument("--parallel", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    coarse = read_log(args.coarse, args.dt)
    fine = read_log(args.fine, args.dt / 2)
    finer = read_log(args.finer, args.dt / 4)
    first = difference(coarse, fine)
    second = difference(fine, finer)
    drifts = {name: energy_drift(value) for name, value in
              (("coarse", coarse), ("fine", fine), ("finer", finer))}
    # This case contains two Al atoms. Use a 1 meV/atom trajectory budget;
    # this is separate from Autotest's final-property comparison thresholds.
    drift_mev_per_atom = {name: value * 13.605698 * 1000 / 2
                          for name, value in drifts.items()}
    checks = {
        "same_final_time_0.5_fs": max(coarse) == max(fine) == max(finer) == Fraction("0.5"),
        "energy_drift_below_1_meV_per_atom": max(drift_mev_per_atom.values()) < 1.0,
        "trajectory_stress_difference_below_0.01_kbar": first["max_component_kbar"] < 0.01,
        "stress_difference_decreases_on_refinement":
            second["max_component_kbar"] < first["max_component_kbar"],
        "final_energy_difference_below_1e-5_Ry":
            max(first["final_energy_difference_Ry"], second["final_energy_difference_Ry"]) < 1e-5,
    }
    result = {"end_time_fs": float(max(coarse)), "energy_drift_Ry": drifts,
              "energy_drift_meV_per_atom": drift_mev_per_atom,
              "coarse_vs_fine": first, "fine_vs_finer": second, "checks": checks}
    if args.parallel:
        parallel = read_log(args.parallel, args.dt)
        result["parallel"] = difference(coarse, parallel)
        checks["parallel_stress_difference_below_1e-6_kbar"] = result["parallel"]["max_component_kbar"] < 1e-6
        checks["parallel_energy_drift_below_1_meV_per_atom"] = energy_drift(parallel) * 13.605698 * 1000 / 2 < 1.0
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0 if all(checks.values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
