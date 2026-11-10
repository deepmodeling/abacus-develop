"""Run an isolated Fe DeltaSpin benchmark and retain numerical diagnostics."""

import argparse
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?"


def digest(path: Path) -> str:
    """Return the SHA256 of a file."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def last_number(text: str, pattern: str) -> float | None:
    values = re.findall(pattern, text)
    return float(values[-1]) if values else None


def vectors(text: str, heading: str, count: int) -> list[list[float]]:
    """Read the final three-component table, with or without atom labels."""
    if heading not in text or count <= 0:
        return []
    pattern = re.compile(rf"\s*(?:(?:[A-Za-z]+\d*|ATOM\s+\d+)\s+)?({NUMBER})\s+({NUMBER})\s+({NUMBER})\s*")
    rows = []
    for line in text.rsplit(heading, 1)[1].splitlines():
        match = pattern.fullmatch(line)
        if match:
            rows.append([float(value) for value in match.groups()])
        elif rows:
            break
        if len(rows) == count:
            break
    return rows


def summarize(directory: Path) -> dict:
    """Extract runtime, convergence, and physical outputs independently."""
    stdout = (directory / "run.log").read_text()
    log_path = directory / "OUT.ABACUS/running_scf.log"
    log = log_path.read_text() if log_path.exists() else ""
    iterations = re.findall(r"^\s*DS(\d+)\s+", stdout, re.MULTILINE)
    energy = last_number(log, rf"!FINAL_ETOT_IS\s+({NUMBER})")
    atom_counts = re.findall(r"TOTAL ATOM NUMBER\s*=\s*(\d+)", log)
    atom_count = int(atom_counts[-1]) if atom_counts else 0
    fields = {
        "atom_count": atom_count,
        "energy_eV": energy,
        "forces_eV_per_A": vectors(log, "TOTAL-FORCE (eV/Angstrom)", atom_count),
        "stress_kbar": vectors(log, "TOTAL-STRESS (KBAR)", 3),
        "moments_uB": vectors(log, "Total Magnetism (uB)", atom_count),
        "magnetic_forces_eV_per_uB": vectors(log, "Magnetic force (eV/uB)", atom_count),
        "target_moments_uB": vectors(stdout, "target spin (uB):", atom_count),
        "native_scf_converged": "charge density convergence is achieved" in log,
        "scf_steps": int(iterations[-1]) if iterations else 0,
        "density_error": last_number(log, rf"Density error is\s+({NUMBER})"),
        "spin_rms_uB": last_number(stdout, rf"Current RMS:\s*({NUMBER})"),
        "inner_loop_limit_hits": stdout.count("Reach maximum number of steps ("),
    }
    expected = {"forces_eV_per_A": atom_count, "stress_kbar": 3, "moments_uB": atom_count,
                "magnetic_forces_eV_per_uB": atom_count, "target_moments_uB": atom_count}
    finite = atom_count > 0 and energy is not None and math.isfinite(energy)
    for key, count in expected.items():
        finite = finite and len(fields[key]) == count
        finite = finite and all(math.isfinite(value) for row in fields[key] for value in row)
    fields["finite_outputs"] = finite
    fields["runtime_complete"] = bool(re.search(r"FINISH\s+Time", stdout)) and finite
    fields["timers"] = []
    if "TIME STATISTICS" in log:
        timer_pattern = re.compile(rf"^\s*(\w+)\s+(.+?)\s+({NUMBER})\s+(\d+)\s+({NUMBER})\s+({NUMBER})\s*$")
        for line in log.rsplit("TIME STATISTICS", 1)[1].splitlines():
            match = timer_pattern.fullmatch(line)
            if match:
                fields["timers"].append({"class": match[1], "name": match[2],
                                         "seconds": float(match[3]), "calls": int(match[4]),
                                         "percent": float(match[6])})
    return fields


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--case", choices=sorted(p.name for p in (ROOT / "inputs").iterdir()), required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--threads", type=int, default=12)
    parser.add_argument("--ranks", type=int, default=1)
    parser.add_argument("--kpar", type=int)
    parser.add_argument("--steps", type=int)
    parser.add_argument("--profile", action="store_true")
    parser.add_argument("--cpu-profile", type=Path, help="perf binary for sampling the main application thread")
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", args.label):
        parser.error("label must contain only letters, digits, underscores, periods, or hyphens")
    if min(args.threads, args.ranks) < 1 or (args.steps is not None and args.steps < 1):
        parser.error("threads, ranks, and steps must be positive")
    if args.cpu_profile and (args.profile or args.ranks != 1):
        parser.error("CPU sampling requires one rank and cannot be combined with --profile")
    binary = args.binary.resolve(strict=True)
    cpu_profiler = args.cpu_profile.resolve(strict=True) if args.cpu_profile else None
    directory = ROOT / "results" / args.label
    directory.mkdir(exist_ok=False)
    for name in ("INPUT", "STRU", "KPT"):
        shutil.copy2(ROOT / "inputs" / args.case / name, directory / name)
    for asset in (ROOT / "assets").iterdir():
        (directory / asset.name).symlink_to(asset)
    input_path = directory / "INPUT"
    content = input_path.read_text()
    if args.steps is not None:
        content, count = re.subn(r"(?m)^scf_nmax\s+.*$", f"scf_nmax {args.steps}", content)
        assert count == 1
    if args.kpar is not None:
        if args.kpar < 1 or args.ranks % args.kpar:
            parser.error("kpar must be positive and divide ranks")
        content += f"\nkpar {args.kpar}\n"
    input_path.write_text(content)
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS=str(args.threads), OPENBLAS_NUM_THREADS="1", CUDA_VISIBLE_DEVICES="0",
               OMPI_ALLOW_RUN_AS_ROOT="1", OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1")
    command = ["mpirun", "--bind-to", "none", "-np", str(args.ranks), str(binary)]
    if cpu_profiler:
        command = command[:-1] + [str(cpu_profiler), "record", "-e", "cpu-clock:u", "-F", "99",
                                 "--no-inherit", "--per-thread", "--call-graph", "dwarf,2048", "-m", "8",
                                 "-o", "cpu.perf.data", "--", str(binary)]
    if args.profile:
        command = ["nsys", "profile", "--trace=cuda,nvtx,osrt", "--sample=none",
                   "--cpuctxsw=none", "--output=profile"] + command
    metadata = {
        "case": args.case, "label": args.label, "binary": str(binary), "binary_sha256": digest(binary),
        "started_at": datetime.now(timezone.utc).isoformat(), "command": command,
        "threads": args.threads, "ranks": args.ranks, "kpar": args.kpar or 1,
        "omp_wait_policy": env.get("OMP_WAIT_POLICY", "runtime_default"),
        "gomp_spincount": env.get("GOMP_SPINCOUNT", "runtime_default"),
        "step_limit_override": args.steps, "profiled": args.profile or cpu_profiler is not None,
        "profile_tool": "perf" if cpu_profiler else "nsys" if args.profile else None,
        "cpu_sample_scope": "main_thread" if cpu_profiler else None,
        "input_sha256": {name: digest(directory / name) for name in ("INPUT", "STRU", "KPT")},
        "asset_sha256": {asset.name: digest(asset) for asset in (ROOT / "assets").iterdir()},
    }
    (directory / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    with (directory / "gpu.csv").open("w") as gpu_log:
        monitor = subprocess.Popen(
            ["nvidia-smi", "--id=0", "--query-gpu=timestamp,memory.used,utilization.gpu",
             "--format=csv,noheader,nounits", "--loop-ms=100"],
            stdout=gpu_log, stderr=subprocess.STDOUT)
        try:
            start = time.monotonic()
            with (directory / "run.log").open("w") as stream:
                completed = subprocess.run(command, cwd=directory, env=env, stdout=stream, stderr=subprocess.STDOUT)
            elapsed = time.monotonic() - start
        finally:
            monitor.terminate()
            monitor.wait(timeout=5)
    memory_samples = [float(line.split(",")[1]) for line in (directory / "gpu.csv").read_text().splitlines()
                      if len(line.split(",")) == 3]
    result = {**metadata, "elapsed_seconds": elapsed,
              "sampled_peak_gpu_memory_MiB": max(memory_samples) if memory_samples else None,
              "gpu_memory_sample_interval_ms": 100,
              "return_code": completed.returncode, **summarize(directory)}
    (directory / "result.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({key: result[key] for key in ("label", "elapsed_seconds", "return_code",
                                                 "runtime_complete", "native_scf_converged", "scf_steps")}), flush=True)
    if completed.returncode or not result["runtime_complete"]:
        raise SystemExit(completed.returncode or 1)


if __name__ == "__main__":
    main()
