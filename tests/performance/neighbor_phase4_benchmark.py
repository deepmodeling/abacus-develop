#!/usr/bin/env python3
"""Run reproducible Phase 4 neighbor-search benchmarks and write raw CSV data."""

import argparse
import csv
import hashlib
import json
import os
import platform
import re
import shutil
import statistics
import subprocess
import time
from pathlib import Path


KERNEL_HEADER = [
    "kind",
    "mode",
    "atoms",
    "repeat",
    "storage_us",
    "search_us",
    "page_us",
    "pairs",
    "boxes",
    "pages",
    "page_utilization",
    "page_memory_bytes",
    "peak_rss_kib",
    "pair_hash",
]

MPI_MAX_FIELDS = [
    "exchange_us",
    "search_us",
    "local_atoms",
    "ghost_atoms",
    "neighbor_ranks",
    "nonempty_send_ranks",
    "sent_atoms",
    "received_ghosts",
    "sent_payload",
    "received_payload",
    "pairs",
]

GRID_PATH_HEADER = [
    "kind",
    "implementation",
    "cell",
    "atoms",
    "radius",
    "repeat",
    "pack_us",
    "storage_us",
    "search_us",
    "total_us",
    "pairs",
    "peak_rss_kib",
    "pair_hash",
    "status",
    "reason",
]

VERSION_CASES = {
    "si16": "P100_si16_lcao",
    "si64": "P102_si64_lcao",
    "si128": "P103_si128_lcao",
}


def run(command, cwd=None, env=None, timeout=1200, output=None):
    started = time.monotonic()
    completed = subprocess.run(
        command,
        cwd=cwd,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=timeout,
        check=False,
    )
    elapsed = time.monotonic() - started
    if output:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(completed.stdout)
    if completed.returncode:
        raise RuntimeError(
            f"Command failed ({completed.returncode}): {' '.join(map(str, command))}\n"
            f"{completed.stdout[-4000:]}"
        )
    return completed.stdout, elapsed


def executable(path):
    resolved = path.resolve()
    if not resolved.is_file() or not os.access(resolved, os.X_OK):
        raise ValueError(f"Executable not found: {resolved}")
    return resolved


def write_csv(path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_kernel(args, output):
    binary = executable(args.kernel)
    rows = []
    for atoms in args.kernel_atoms:
        for mode in ("half14", "full27", "half14_ref"):
            text, _ = run(
                [
                    str(binary),
                    "--atoms",
                    str(atoms),
                    "--repeats",
                    str(args.kernel_repeats),
                    "--mode",
                    mode,
                ],
                timeout=args.timeout,
                output=output / "logs" / f"kernel_{mode}_{atoms}.log",
            )
            for line in text.splitlines():
                if not line.startswith("neighbor_kernel,"):
                    continue
                values = line.split(",")
                rows.append(dict(zip(KERNEL_HEADER, values)))
    write_csv(output / "kernel.csv", rows, KERNEL_HEADER)


def run_grid_paths(args, output):
    current_binary = executable(args.current_grid_kernel)
    beta_binary = executable(args.beta_grid_kernel)
    rows = []
    for atoms in args.grid_atoms:
        for cell in args.grid_cells:
            for radius in args.grid_radii:
                current_rows = []
                for mode in ("full27", "half14"):
                    text, _ = run(
                        [
                            str(current_binary),
                            "--atoms",
                            str(atoms),
                            "--repeats",
                            str(args.grid_repeats),
                            "--cell",
                            cell,
                            "--radius",
                            str(radius),
                            "--mode",
                            mode,
                        ],
                        timeout=args.timeout,
                        output=output / "logs" / f"grid_current_{mode}_{cell}_{atoms}_r{radius}.log",
                    )
                    for line in text.splitlines():
                        if line.startswith("neighbor_grid_path,"):
                            row = dict(zip(GRID_PATH_HEADER[:-2], line.split(",")))
                            row.update({"status": "ok", "reason": ""})
                            rows.append(row)
                            current_rows.append(row)

                current_hashes = {row["pair_hash"] for row in current_rows}
                current_pairs = {row["pairs"] for row in current_rows}
                if len(current_hashes) != 1 or len(current_pairs) != 1:
                    raise RuntimeError(
                        f"Current Full27/Half14 mismatch for {cell}, atoms={atoms}, radius={radius}."
                    )

                if atoms > args.beta_grid_max_atoms:
                    rows.append(
                        {
                            "kind": "neighbor_grid_path",
                            "implementation": "beta1_legacy",
                            "cell": cell,
                            "atoms": atoms,
                            "radius": radius,
                            "repeat": "",
                            "pack_us": "",
                            "storage_us": "",
                            "search_us": "",
                            "total_us": "",
                            "pairs": next(iter(current_pairs)),
                            "peak_rss_kib": "",
                            "pair_hash": next(iter(current_hashes)),
                            "status": "skipped",
                            "reason": "legacy O(N^2) exceeds the configured practical size limit",
                        }
                    )
                    continue

                text, _ = run(
                    [
                        str(beta_binary),
                        "--atoms",
                        str(atoms),
                        "--repeats",
                        str(
                            args.beta_grid_large_repeats
                            if atoms > args.beta_grid_full_repeat_max_atoms
                            else args.beta_grid_repeats
                        ),
                        "--cell",
                        cell,
                        "--radius",
                        str(radius),
                    ],
                    timeout=args.timeout,
                    output=output / "logs" / f"grid_beta1_legacy_{cell}_{atoms}_r{radius}.log",
                )
                beta_rows = []
                for line in text.splitlines():
                    if line.startswith("neighbor_grid_path,"):
                        row = dict(zip(GRID_PATH_HEADER[:-2], line.split(",")))
                        row.update({"status": "ok", "reason": ""})
                        rows.append(row)
                        beta_rows.append(row)
                if not beta_rows:
                    raise RuntimeError(
                        f"beta.1 legacy produced no rows for {cell}, atoms={atoms}, radius={radius}."
                    )
                if {row["pair_hash"] for row in beta_rows} != current_hashes:
                    raise RuntimeError(
                        f"beta.1/current pair hash mismatch for {cell}, atoms={atoms}, radius={radius}."
                    )
                if {row["pairs"] for row in beta_rows} != current_pairs:
                    raise RuntimeError(
                        f"beta.1/current pair count mismatch for {cell}, atoms={atoms}, radius={radius}."
                    )
    write_csv(output / "grid_paths.csv", rows, GRID_PATH_HEADER)


def run_mpi(args, output):
    binary = executable(args.mpi_kernel)
    rows = []
    environment = dict(os.environ)
    environment.update(
        {
            "OMP_NUM_THREADS": "1",
            "OMP_PROC_BIND": "true",
            "OMP_PLACES": "cores",
        }
    )
    fields = ["kind", "scaling", "cell", "ranks", "atoms_per_rank", "global_atoms", "repeat", "radius"]
    fields += [f"max_{name}" for name in MPI_MAX_FIELDS]
    fields += [f"sum_{name}" for name in MPI_MAX_FIELDS]
    for scaling in args.mpi_scaling:
        for cell in args.mpi_cells:
            for ranks in args.mpi_ranks:
                command = [
                    "mpirun",
                    "--mca",
                    "btl",
                    "^tcp",
                    "--bind-to",
                    "core",
                    "-np",
                    str(ranks),
                    str(binary),
                    "--atoms-per-rank",
                    str(args.mpi_atoms_per_rank),
                    "--repeats",
                    str(args.mpi_repeats),
                    "--cell",
                    cell,
                ]
                if scaling == "strong":
                    command.extend(
                        [
                            "--shape-x",
                            str(args.mpi_strong_shape[0]),
                            "--shape-y",
                            str(args.mpi_strong_shape[1]),
                            "--shape-z",
                            str(args.mpi_strong_shape[2]),
                        ]
                    )
                text, _ = run(
                    command,
                    env=environment,
                    timeout=args.timeout,
                    output=output / "logs" / f"mpi_{scaling}_{cell}_{ranks}.log",
                )
                for line in text.splitlines():
                    if line.startswith("neighbor_mpi,"):
                        rows.append(dict(zip(fields, line.split(","))))
    write_csv(output / "mpi.csv", rows, fields)


def replace_or_append(text, key, value):
    pattern = re.compile(rf"(?m)^\s*{re.escape(key)}\s+.*$")
    line = f"{key:<22} {value}"
    return pattern.sub(line, text) if pattern.search(text) else text.rstrip() + "\n" + line + "\n"


def prepare_version_case(repo, work_root, version, case_name, ranks):
    source = repo / "tests" / "performance" / VERSION_CASES[case_name]
    destination = work_root / version / case_name
    if destination.exists():
        shutil.rmtree(destination)
    shutil.copytree(source, destination)
    input_path = destination / "INPUT"
    text = input_path.read_text()
    if ranks == 1:
        # The conda ScaLAPACK library is not ABI-compatible with beta.1 when it
        # is linked against MKL. The serial LAPACK solver keeps the comparison
        # inside the same BLAS implementation and avoids that unrelated crash.
        text = replace_or_append(text, "ks_solver", "lapack")
    input_path.write_text(text)
    stru_path = destination / "STRU"
    stru_text = stru_path.read_text()
    stru_text = stru_text.replace("../../PP_ORB/", f"{repo / 'tests' / 'PP_ORB'}/")
    stru_path.write_text(stru_text)
    return destination


def last_float(pattern, text):
    values = re.findall(pattern, text)
    return float(values[-1]) if values else None


def read_abacus_logs(case_dir):
    logs = sorted(case_dir.glob("OUT.*/running_*.log"))
    return "\n".join(path.read_text(errors="replace") for path in logs)


def run_version(args, output):
    configurations = [
        ("beta1", args.beta_repo.resolve(), executable(args.beta_abacus)),
        ("current", args.current_repo.resolve(), executable(args.current_abacus)),
    ]
    rows = []
    fieldnames = [
        "version",
        "case",
        "ranks",
        "repeat",
        "wall_s",
        "peak_rss_kib",
        "final_energy_ev",
        "largest_gradient_ev_per_angstrom",
    ]
    environment = dict(os.environ)
    environment.update(
        {
            "OMP_NUM_THREADS": "1",
            "OMP_PROC_BIND": "true",
            "OMP_PLACES": "cores",
        }
    )
    work_root = output / "version_work"
    for case_name in args.version_cases:
        for ranks in args.version_ranks:
            prepared = {
                version: prepare_version_case(repo, work_root, version, case_name, ranks)
                for version, repo, _ in configurations
            }
            for repeat in range(args.version_repeats + 1):
                for version, repo, abacus in configurations:
                    case_dir = prepared[version]
                    shutil.rmtree(case_dir / "OUT.autotest", ignore_errors=True)
                    command = [
                        "/usr/bin/time",
                        "-v",
                        "mpirun",
                        "--mca",
                        "btl",
                        "^tcp",
                        "--bind-to",
                        "core",
                        "-np",
                        str(ranks),
                        str(abacus),
                    ]
                    text, wall = run(
                        command,
                        cwd=case_dir,
                        env=environment,
                        timeout=args.timeout,
                        output=output / "logs" / f"{version}_{case_name}_r{ranks}_{repeat}.log",
                    )
                    if repeat == 0:
                        continue
                    rss = last_float(r"Maximum resident set size \(kbytes\):\s+([0-9.]+)", text)
                    log_text = read_abacus_logs(case_dir)
                    energy = last_float(r"!FINAL_ETOT_IS\s+([-+0-9.eE]+)", log_text)
                    force = last_float(r"LARGEST GRAD \(eV/Angstrom\)\s+:\s+([-+0-9.eE]+)", log_text)
                    rows.append(
                        {
                            "version": version,
                            "case": case_name,
                            "ranks": ranks,
                            "repeat": repeat - 1,
                            "wall_s": wall,
                            "peak_rss_kib": rss,
                            "final_energy_ev": energy,
                            "largest_gradient_ev_per_angstrom": force,
                        }
                    )
                    # Persist completed repetitions immediately so an unrelated
                    # later timeout does not discard already valid measurements.
                    write_csv(output / "version.csv", rows, fieldnames)


def collect_environment(args, output):
    commands = {
        "uname": ["uname", "-a"],
        "lscpu": ["lscpu"],
        "memory": ["free", "-h"],
        "compiler": ["conda", "run", "-n", "bxcx", "mpicxx", "--version"],
        "mpi": ["conda", "run", "-n", "bxcx", "mpirun", "--version"],
    }
    data = {
        "platform": platform.platform(),
        "current_repo": str(args.current_repo.resolve()),
        "beta_repo": str(args.beta_repo.resolve()),
        "arguments": vars(args),
    }
    data["arguments"] = {key: str(value) if isinstance(value, Path) else value for key, value in data["arguments"].items()}
    for name, command in commands.items():
        try:
            data[name] = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout
        except OSError as error:
            data[name] = str(error)
    for name, repo in (("current_commit", args.current_repo), ("beta_commit", args.beta_repo)):
        data[name] = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=repo, text=True, stdout=subprocess.PIPE, check=True
        ).stdout.strip()
    current_diff = subprocess.run(
        ["git", "diff", "--binary", "HEAD", "--", "source/source_cell/module_neighbor"],
        cwd=args.current_repo,
        text=False,
        stdout=subprocess.PIPE,
        check=True,
    ).stdout
    data["current_worktree_dirty"] = bool(
        subprocess.run(
            ["git", "status", "--short"],
            cwd=args.current_repo,
            text=True,
            stdout=subprocess.PIPE,
            check=True,
        ).stdout.strip()
    )
    data["current_neighbor_diff_sha256"] = hashlib.sha256(current_diff).hexdigest()
    output.mkdir(parents=True, exist_ok=True)
    (output / "environment.json").write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n")


def summarize_csv(path, group_fields, value_fields):
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    groups = {}
    for row in rows:
        key = tuple(row[field] for field in group_fields)
        groups.setdefault(key, []).append(row)
    summary = []
    for key, values in sorted(groups.items()):
        item = dict(zip(group_fields, key))
        for field in value_fields:
            numbers = [
                float(value[field])
                for value in values
                if value.get(field) not in ("", "None", None) and float(value[field]) >= 0.0
            ]
            if numbers:
                item[f"{field}_median"] = statistics.median(numbers)
                item[f"{field}_min"] = min(numbers)
                item[f"{field}_max"] = max(numbers)
        summary.append(item)
    return summary


def write_summary(output):
    summary = {}
    if (output / "grid_paths.csv").exists():
        summary["grid_paths"] = summarize_csv(
            output / "grid_paths.csv",
            ["implementation", "cell", "atoms", "radius", "status"],
            ["pack_us", "storage_us", "search_us", "total_us", "peak_rss_kib"],
        )
    if (output / "kernel.csv").exists():
        summary["kernel"] = summarize_csv(
            output / "kernel.csv", ["mode", "atoms"], ["storage_us", "search_us", "page_us", "peak_rss_kib"]
        )
    if (output / "mpi.csv").exists():
        summary["mpi"] = summarize_csv(
            output / "mpi.csv",
            ["scaling", "cell", "ranks", "atoms_per_rank"],
            ["max_exchange_us", "max_search_us", "sum_ghost_atoms", "sum_sent_payload"],
        )
    if (output / "version.csv").exists():
        summary["version"] = summarize_csv(
            output / "version.csv",
            ["version", "case", "ranks"],
            ["wall_s", "peak_rss_kib", "final_energy_ev", "largest_gradient_ev_per_angstrom"],
        )
    verlet_path = output / "verlet" / "neighbor_verlet_validation.csv"
    if verlet_path.exists():
        summary["verlet"] = summarize_csv(
            verlet_path,
            ["case", "ranks", "skin_bohr"],
            [
                "total_time_s",
                "rebuild_count",
                "reuse_count",
                "final_energy_ev",
                "final_max_force_ev_per_angstrom",
            ],
        )
    (output / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--current-repo", type=Path, required=True)
    parser.add_argument("--beta-repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--kernel", type=Path)
    parser.add_argument("--mpi-kernel", type=Path)
    parser.add_argument("--current-abacus", type=Path)
    parser.add_argument("--beta-abacus", type=Path)
    parser.add_argument("--sections", nargs="+", choices=("grid", "kernel", "mpi", "version"), required=True)
    parser.add_argument("--current-grid-kernel", type=Path)
    parser.add_argument("--beta-grid-kernel", type=Path)
    parser.add_argument("--grid-atoms", type=int, nargs="+", default=[4096, 32768, 131072])
    parser.add_argument("--grid-cells", nargs="+", choices=("orthogonal", "triclinic"), default=["orthogonal", "triclinic"])
    parser.add_argument("--grid-radii", type=float, nargs="+", default=[2.2, 3.4])
    parser.add_argument("--grid-repeats", type=int, default=10)
    parser.add_argument("--beta-grid-repeats", type=int, default=10)
    parser.add_argument("--beta-grid-large-repeats", type=int, default=3)
    parser.add_argument("--beta-grid-full-repeat-max-atoms", type=int, default=4096)
    parser.add_argument("--beta-grid-max-atoms", type=int, default=32768)
    parser.add_argument("--kernel-atoms", type=int, nargs="+", default=[4096, 32768, 131072])
    parser.add_argument("--kernel-repeats", type=int, default=10)
    parser.add_argument("--mpi-atoms-per-rank", type=int, default=4096)
    parser.add_argument("--mpi-ranks", type=int, nargs="+", default=[1, 2, 4, 8, 16])
    parser.add_argument("--mpi-repeats", type=int, default=5)
    parser.add_argument("--mpi-scaling", nargs="+", choices=("strong", "weak"), default=["strong", "weak"])
    parser.add_argument("--mpi-cells", nargs="+", choices=("orthogonal", "triclinic"), default=["orthogonal", "triclinic"])
    parser.add_argument("--mpi-strong-shape", type=int, nargs=3, default=[64, 32, 32])
    parser.add_argument("--version-cases", nargs="+", choices=tuple(VERSION_CASES), default=["si16", "si64", "si128"])
    parser.add_argument("--version-ranks", type=int, nargs="+", default=[1, 4])
    parser.add_argument("--version-repeats", type=int, default=3)
    parser.add_argument("--timeout", type=int, default=1200)
    args = parser.parse_args()

    if "kernel" in args.sections and not args.kernel:
        parser.error("--kernel is required for the kernel section")
    if "grid" in args.sections and (not args.current_grid_kernel or not args.beta_grid_kernel):
        parser.error("--current-grid-kernel and --beta-grid-kernel are required for the grid section")
    if "mpi" in args.sections and not args.mpi_kernel:
        parser.error("--mpi-kernel is required for the mpi section")
    if "version" in args.sections and (not args.current_abacus or not args.beta_abacus):
        parser.error("--current-abacus and --beta-abacus are required for the version section")
    if any(
        value <= 0
        for value in args.grid_atoms + args.kernel_atoms + args.mpi_ranks + args.version_ranks
    ):
        parser.error("Atom and rank values must be positive")

    output = args.output.resolve()
    collect_environment(args, output)
    if "grid" in args.sections:
        run_grid_paths(args, output)
    if "kernel" in args.sections:
        run_kernel(args, output)
    if "mpi" in args.sections:
        run_mpi(args, output)
    if "version" in args.sections:
        run_version(args, output)
    write_summary(output)


if __name__ == "__main__":
    main()
