#!/usr/bin/env python3
"""Run lightweight LCAO relax/MD comparisons for the Verlet neighbor list."""

import argparse
import csv
import hashlib
import math
import os
import re
import shutil
import subprocess
from pathlib import Path


CASES = {
    "relax_si2": ("examples/17_relax/02_lcao_relax_CG_Si2", "INPUT", 5),
    "md_si8": ("examples/18_md/01_lcao_gamma_Si8", "INPUT_0", 10),
}


def replace_or_append(text, key, value):
    pattern = re.compile(rf"(?m)^\s*{re.escape(key)}\s+.*$")
    line = f"{key:<22} {value}"
    return pattern.sub(line, text) if pattern.search(text) else text.rstrip() + "\n" + line + "\n"


def prepare_case(repo, work_root, case_name, skin):
    relative, input_name, steps = CASES[case_name]
    source = repo / relative
    destination = work_root / f"{case_name}_skin_{skin:g}"
    if destination.exists():
        shutil.rmtree(destination)
    shutil.copytree(source, destination)
    input_text = (source / input_name).read_text()
    input_text = replace_or_append(input_text, "pseudo_dir", str(repo / "tests/PP_ORB"))
    input_text = replace_or_append(input_text, "orbital_dir", str(repo / "tests/PP_ORB"))
    input_text = replace_or_append(input_text, "neighbor_skin", skin)
    if case_name == "relax_si2":
        input_text = replace_or_append(input_text, "relax_nmax", steps)
    else:
        input_text = replace_or_append(input_text, "md_nstep", steps)
        input_text = replace_or_append(input_text, "ks_solver", "scalapack_gvx")
    (destination / "INPUT").write_text(input_text)
    return destination


def last_float(pattern, text):
    matches = re.findall(pattern, text)
    return float(matches[-1]) if matches else None


def final_max_force(log_text):
    marker = "#TOTAL-FORCE (eV/Angstrom)#"
    position = log_text.rfind(marker)
    if position < 0:
        return last_float(r"LARGEST GRAD \(eV/Angstrom\)\s+:\s+([-+0-9.eE]+)", log_text)
    maximum = 0.0
    for line in log_text[position:].splitlines()[4:]:
        fields = line.split()
        if len(fields) != 4:
            if maximum > 0.0:
                break
            continue
        try:
            force = [float(value) for value in fields[1:]]
        except ValueError:
            continue
        maximum = max(maximum, math.sqrt(sum(value * value for value in force)))
    return maximum if maximum > 0.0 else None


def run_case(abacus, repo, work_root, case_name, skin, ranks, threads):
    case_dir = prepare_case(repo, work_root, case_name, skin)
    command = ["mpirun", "--mca", "btl", "^tcp", "-np", str(ranks), str(abacus)]
    environment = dict(os.environ)
    environment["OMP_NUM_THREADS"] = str(threads)
    completed = subprocess.run(
        command,
        cwd=case_dir,
        env=environment,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    (case_dir / "validation.stdout").write_text(completed.stdout)

    logs = list(case_dir.glob("OUT.*/running_*.log"))
    log_text = "\n".join(path.read_text(errors="replace") for path in logs)
    decisions = re.findall(r"NEIGHBOR_VERLET .*?decision=(\w+).*?reason=(\w+)", log_text)
    rebuilds = sum(decision == "rebuild" for decision, _ in decisions)
    reuses = sum(decision == "reuse" for decision, _ in decisions)
    final_energy = last_float(r"!FINAL_ETOT_IS\s+([-+0-9.eE]+)", log_text)
    total_time = last_float(r"Total\s+Time\s+:\s+\d+\s+h\s+\d+\s+mins\s+([-+0-9.eE]+)", log_text)
    final_structures = sorted(case_dir.glob("OUT.*/STRU_ION_D"))
    if not final_structures:
        final_structures = sorted(case_dir.glob("OUT.*/STRU.cif"))
    structure_hash = ""
    if final_structures:
        structure_hash = hashlib.sha256(final_structures[-1].read_bytes()).hexdigest()

    return {
        "case": case_name,
        "skin_bohr": skin,
        "ranks": ranks,
        "threads": threads,
        "exit_code": completed.returncode,
        "final_energy_ev": final_energy,
        "final_max_force_ev_per_angstrom": final_max_force(log_text),
        "total_time_s": total_time,
        "rebuild_count": rebuilds,
        "reuse_count": reuses,
        "final_structure_sha256": structure_hash,
        "work_directory": str(case_dir),
    }


def write_report(rows, output):
    output.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0].keys())
    with (output / "neighbor_verlet_validation.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    lines = [
        "# Neighbor Verlet Validation",
        "",
        "| Case | Skin | Ranks | Exit | Rebuild | Reuse | Final energy (eV) | Max force | Time (s) |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['case']} | {row['skin_bohr']} | {row['ranks']} | {row['exit_code']} | "
            f"{row['rebuild_count']} | {row['reuse_count']} | {row['final_energy_ev']} | "
            f"{row['final_max_force_ev_per_angstrom']} | {row['total_time_s']} |"
        )
    (output / "neighbor_verlet_validation.md").write_text("\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--abacus", required=True, type=Path)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--work-dir", type=Path, default=Path("/tmp/abacus-neighbor-verlet"))
    parser.add_argument("--output", type=Path, default=Path("neighbor-verlet-results"))
    parser.add_argument("--ranks", type=int, nargs="+", default=[1, 2])
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    if args.threads <= 0 or any(ranks <= 0 for ranks in args.ranks):
        parser.error("--threads and every --ranks value must be positive")
    if not args.abacus.is_file() or not os.access(args.abacus, os.X_OK):
        parser.error("--abacus must point to an executable file")

    rows = []
    for case_name in CASES:
        for ranks in args.ranks:
            for skin in (0.0, 3.0):
                rows.append(
                    run_case(args.abacus.resolve(),
                             args.repo.resolve(),
                             args.work_dir.resolve(),
                             case_name,
                             skin,
                             ranks,
                             args.threads)
                )
    write_report(rows, args.output.resolve())
    if any(row["exit_code"] != 0 for row in rows):
        raise SystemExit(1)
    for case_name in CASES:
        for ranks in args.ranks:
            baseline = next(row for row in rows
                            if row["case"] == case_name and row["ranks"] == ranks and row["skin_bohr"] == 0.0)
            reused = next(row for row in rows
                          if row["case"] == case_name and row["ranks"] == ranks and row["skin_bohr"] == 3.0)
            required_metrics = (
                baseline["final_energy_ev"],
                reused["final_energy_ev"],
                baseline["final_max_force_ev_per_angstrom"],
                reused["final_max_force_ev_per_angstrom"],
            )
            if any(value is None for value in required_metrics):
                raise SystemExit(f"{case_name} rank {ranks}: validation metrics are incomplete")
            if abs(baseline["final_energy_ev"] - reused["final_energy_ev"]) > 1.0e-8:
                raise SystemExit(f"{case_name} rank {ranks}: final energy mismatch")
            if abs(baseline["final_max_force_ev_per_angstrom"]
                   - reused["final_max_force_ev_per_angstrom"]) > 1.0e-8:
                raise SystemExit(f"{case_name} rank {ranks}: final force mismatch")
            if baseline["final_structure_sha256"] != reused["final_structure_sha256"]:
                raise SystemExit(f"{case_name} rank {ranks}: final structure mismatch")
            if reused["reuse_count"] == 0 or reused["rebuild_count"] >= baseline["rebuild_count"]:
                raise SystemExit(f"{case_name} rank {ranks}: neighbor list was not reused")


if __name__ == "__main__":
    main()
