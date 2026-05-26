#!/usr/bin/env python3
"""Run the 05_pexsi.md section 5.1 solver-scaling benchmark.

The benchmark compares the conventional LCAO generalized eigensolver
(`scalapack_gvx`) with `pexsi` on the five systems listed in section 5.1:
Si64 Gamma, Al128 4x4x4, Cu256 2x2x2, TiO2-192 2x2x2, and H2O Gamma.

Each case is run for np=1,2,4,8,16 by default. Results are written after every
job so interrupted runs can be resumed.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import json
import math
import os
import re
import shutil
import signal
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


RY_TO_EV = 13.605693009
ENERGY_TOL_RY = 1.0e-5
FORCE_TOL_EV_A = 1.0e-4
STRESS_TOL_GPA = 1.0e-3

DEFAULT_NP = (1, 2, 4, 8, 16)
DEFAULT_SOLVERS = ("scalapack_gvx", "pexsi")

FINAL_ETOT_RE = re.compile(r"!\s*FINAL_ETOT_IS\s+([-+0-9.eE]+)\s+eV")
TOTAL_TIME_RE = re.compile(r"Total\s+Time\s*:\s*([-+0-9.eE]+)")
TIMER_RE = re.compile(r"^\s*([A-Za-z0-9_./ -]+?)\s+([-+0-9.eE]+)\s+(?:s|sec|seconds)?\s*$")


@dataclass(frozen=True)
class CaseSpec:
    name: str
    natom: int
    kmesh: tuple[int, int, int]
    generator: str
    scf_nmax: int = 100
    ecutwfc: int = 100
    mixing_beta: float = 0.3
    nbands: int | None = None
    smearing_sigma: float = 0.015

    @property
    def nks(self) -> int:
        return self.kmesh[0] * self.kmesh[1] * self.kmesh[2]


CASES: dict[str, CaseSpec] = {
    "si64_gamma": CaseSpec(
        name="si64_gamma",
        natom=64,
        kmesh=(1, 1, 1),
        generator="copy_si64",
        scf_nmax=100,
        ecutwfc=60,
        mixing_beta=0.4,
        nbands=200,
    ),
    "al128_444": CaseSpec(
        name="al128_444",
        natom=128,
        kmesh=(4, 4, 4),
        generator="gen_al128",
        mixing_beta=0.3,
    ),
    "cu256_222": CaseSpec(
        name="cu256_222",
        natom=256,
        kmesh=(2, 2, 2),
        generator="gen_cu256",
        mixing_beta=0.4,
    ),
    "tio2_192_222": CaseSpec(
        name="tio2_192_222",
        natom=192,
        kmesh=(2, 2, 2),
        generator="gen_tio2_192",
        mixing_beta=0.3,
    ),
    "h2o_gamma": CaseSpec(
        name="h2o_gamma",
        natom=3,
        kmesh=(1, 1, 1),
        generator="copy_h2o",
        scf_nmax=50,
        ecutwfc=100,
        mixing_beta=0.4,
        nbands=6,
        smearing_sigma=0.02,
    ),
}


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def timestamp() -> str:
    return _dt.datetime.now().strftime("%Y%m%d_%H%M%S")


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def format_vec(vec: Iterable[float]) -> str:
    return " ".join(f"{x:.10f}" for x in vec)


def direct_supercell_positions(
    base_positions: Iterable[tuple[float, float, float]],
    reps: tuple[int, int, int],
) -> list[tuple[float, float, float]]:
    rx, ry, rz = reps
    positions: list[tuple[float, float, float]] = []
    for ix in range(rx):
        for iy in range(ry):
            for iz in range(rz):
                for x, y, z in base_positions:
                    positions.append(((x + ix) / rx, (y + iy) / ry, (z + iz) / rz))
    return positions


def write_kpt(case_dir: Path, kmesh: tuple[int, int, int]) -> None:
    write_text(
        case_dir / "KPT",
        "K_POINTS\n"
        "0\n"
        "Gamma\n"
        f"{kmesh[0]} {kmesh[1]} {kmesh[2]} 0 0 0\n",
    )


def copy_case_files(src: Path, dst: Path, kmesh: tuple[int, int, int]) -> None:
    dst.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src / "STRU", dst / "STRU")
    write_kpt(dst, kmesh)


def write_al128(dst: Path) -> None:
    reps = (4, 4, 2)
    base = (
        (0.0, 0.0, 0.0),
        (0.5, 0.5, 0.0),
        (0.5, 0.0, 0.5),
        (0.0, 0.5, 0.5),
    )
    positions = direct_supercell_positions(base, reps)
    assert len(positions) == 128
    lines = [
        "ATOMIC_SPECIES",
        "Al 26.982 Al.PD04.PBE.UPF",
        "",
        "NUMERICAL_ORBITAL",
        "Al_gga_10au_100Ry_3s3p2d.orb",
        "",
        "LATTICE_CONSTANT",
        "1.8897261254578281",
        "",
        "LATTICE_VECTORS",
        format_vec((4.04495 * reps[0], 0.0, 0.0)),
        format_vec((0.0, 4.04495 * reps[1], 0.0)),
        format_vec((0.0, 0.0, 4.04495 * reps[2])),
        "",
        "ATOMIC_POSITIONS",
        "Direct",
        "",
        "Al",
        "0.0",
        str(len(positions)),
    ]
    lines.extend(f"{x:.10f} {y:.10f} {z:.10f} 1 1 1" for x, y, z in positions)
    write_text(dst / "STRU", "\n".join(lines) + "\n")
    write_kpt(dst, (4, 4, 4))


def write_cu256(dst: Path) -> None:
    reps = (8, 8, 4)
    positions = direct_supercell_positions(((0.0, 0.0, 0.0),), reps)
    assert len(positions) == 256
    base_vectors = (
        (0.5, 0.5, 0.0),
        (0.5, 0.0, 0.5),
        (0.0, 0.5, 0.5),
    )
    lines = [
        "ATOMIC_SPECIES",
        "Cu 1.000 Cu.LDA.UPF",
        "",
        "NUMERICAL_ORBITAL",
        "Cu_lda_7.0au_100Ry_2s2p2d",
        "",
        "LATTICE_CONSTANT",
        "6.91640",
        "",
        "LATTICE_VECTORS",
    ]
    for vec, rep in zip(base_vectors, reps):
        lines.append(format_vec(tuple(rep * x for x in vec)))
    lines.extend(["", "ATOMIC_POSITIONS", "Direct", "", "Cu", "0.0", str(len(positions))])
    lines.extend(f"{x:.10f} {y:.10f} {z:.10f} 1 1 1" for x, y, z in positions)
    write_text(dst / "STRU", "\n".join(lines) + "\n")
    write_kpt(dst, (2, 2, 2))


def write_tio2_192(dst: Path) -> None:
    reps = (4, 4, 2)
    u = 0.305
    ti_base = ((0.0, 0.0, 0.0), (0.5, 0.5, 0.5))
    o_base = (
        (u, u, 0.0),
        (1.0 - u, 1.0 - u, 0.0),
        (0.5 + u, 0.5 - u, 0.5),
        (0.5 - u, 0.5 + u, 0.5),
    )
    ti_positions = direct_supercell_positions(ti_base, reps)
    o_positions = direct_supercell_positions(o_base, reps)
    assert len(ti_positions) + len(o_positions) == 192
    lines = [
        "ATOMIC_SPECIES",
        "Ti 47.867 Ti_ONCV_PBE-1.0.upf",
        "O 15.9994 O_ONCV_PBE-1.0.upf",
        "",
        "NUMERICAL_ORBITAL",
        "Ti_gga_8au_100Ry_4s2p2d1f.orb",
        "O_gga_8au_100Ry_2s2p1d.orb",
        "",
        "LATTICE_CONSTANT",
        "1.889716",
        "",
        "LATTICE_VECTORS",
        format_vec((4.594 * reps[0], 0.0, 0.0)),
        format_vec((0.0, 4.594 * reps[1], 0.0)),
        format_vec((0.0, 0.0, 2.959 * reps[2])),
        "",
        "ATOMIC_POSITIONS",
        "Direct",
        "",
        "Ti",
        "0.0",
        str(len(ti_positions)),
    ]
    lines.extend(f"{x:.10f} {y:.10f} {z:.10f} 1 1 1" for x, y, z in ti_positions)
    lines.extend(["", "O", "0.0", str(len(o_positions))])
    lines.extend(f"{x:.10f} {y:.10f} {z:.10f} 1 1 1" for x, y, z in o_positions)
    write_text(dst / "STRU", "\n".join(lines) + "\n")
    write_kpt(dst, (2, 2, 2))


def prepare_structure(repo: Path, case: CaseSpec, case_dir: Path) -> None:
    if case.generator == "copy_si64":
        copy_case_files(repo / "examples/33_pexsi/02_scf_Si64", case_dir, case.kmesh)
    elif case.generator == "copy_h2o":
        copy_case_files(repo / "stage3_pexsi_tests/01_h2o_pexsi", case_dir, case.kmesh)
    elif case.generator == "gen_al128":
        write_al128(case_dir)
    elif case.generator == "gen_cu256":
        write_cu256(case_dir)
    elif case.generator == "gen_tio2_192":
        write_tio2_192(case_dir)
    else:
        raise ValueError(f"unknown generator: {case.generator}")


def input_text(
    repo: Path,
    case: CaseSpec,
    solver: str,
    np: int,
    suffix: str,
    pexsi_npole: int,
    pexsi_elec_thr: float,
    pexsi_mu_thr: float,
) -> str:
    gamma_only = 1 if case.kmesh == (1, 1, 1) else 0
    kpar = min(np, case.nks)
    lines = [
        "INPUT_PARAMETERS",
        f"suffix              {suffix}",
        f"pseudo_dir          {repo / 'tests/PP_ORB'}",
        f"orbital_dir         {repo / 'tests/PP_ORB'}",
        "calculation         scf",
        "basis_type          lcao",
        f"gamma_only          {gamma_only}",
        "symmetry            0",
        f"ecutwfc             {case.ecutwfc}",
        "lcao_dr             1e-3",
        "scf_thr             1e-6",
        f"scf_nmax            {case.scf_nmax}",
        "cal_force           1",
        "cal_stress          1",
        "smearing_method     fd",
        f"smearing_sigma      {case.smearing_sigma}",
        "mixing_type         broyden",
        f"mixing_beta         {case.mixing_beta}",
    ]
    if case.nbands is not None:
        lines.append(f"nbands              {case.nbands}")
    if kpar > 1:
        lines.append(f"kpar                {kpar}")
    lines.append(f"ks_solver           {solver}")
    if solver == "pexsi":
        lines.extend(
            [
                f"pexsi_npole         {pexsi_npole}",
                "pexsi_nproc_pole    1",
                f"pexsi_temp          {case.smearing_sigma}",
                "pexsi_gap           0",
                "pexsi_delta_e       20",
                "pexsi_mu_lower      -10",
                "pexsi_mu_upper      10",
                "pexsi_mu            0",
                f"pexsi_mu_thr        {pexsi_mu_thr}",
                "pexsi_mu_expand     0.3",
                "pexsi_mu_guard      0.2",
                f"pexsi_elec_thr      {pexsi_elec_thr}",
                "pexsi_zero_thr      1e-10",
            ]
        )
    return "\n".join(lines) + "\n"


def parse_numeric_tail(line: str, n: int) -> list[float] | None:
    vals: list[float] = []
    for token in line.replace(",", " ").split():
        try:
            vals.append(float(token))
        except ValueError:
            continue
    if len(vals) < n:
        return None
    return vals[-n:]


def parse_block_vectors(lines: list[str], marker: str, width: int, max_rows: int | None) -> list[list[float]]:
    vectors: list[list[float]] = []
    active = False
    for line in lines:
        if marker in line:
            active = True
            vectors = []
            continue
        if not active:
            continue
        if not line.strip():
            if vectors:
                active = False
            continue
        values = parse_numeric_tail(line, width)
        if values is None:
            continue
        vectors.append(values)
        if max_rows is not None and len(vectors) >= max_rows:
            active = False
    return vectors


def parse_log(case_dir: Path, suffix: str) -> dict[str, object]:
    log_path = case_dir / f"OUT.{suffix}" / "running_scf.log"
    result: dict[str, object] = {
        "log": str(log_path),
        "converged": False,
        "energy_ev": None,
        "total_time_s": None,
        "max_force_abs": None,
        "max_stress_abs": None,
        "pexsidft_s": None,
        "hsolver_s": None,
        "transmat2ccs_s": None,
    }
    if not log_path.exists():
        return result
    text = log_path.read_text(errors="replace")
    lines = text.splitlines()
    result["converged"] = "#SCF IS CONVERGED" in text or "SCF IS CONVERGED" in text
    energies = FINAL_ETOT_RE.findall(text)
    if energies:
        result["energy_ev"] = float(energies[-1])
    total_times = TOTAL_TIME_RE.findall(text)
    if total_times:
        result["total_time_s"] = float(total_times[-1])

    forces = parse_block_vectors(lines, "TOTAL-FORCE", 3, None)
    if forces:
        result["max_force_abs"] = max(abs(v) for row in forces for v in row)
    stress = parse_block_vectors(lines, "TOTAL-STRESS", 3, 3)
    if stress:
        result["max_stress_abs"] = max(abs(v) for row in stress for v in row)

    for line in lines:
        if "PEXSIDFT" in line:
            nums = parse_numeric_tail(line, 1)
            if nums:
                result["pexsidft_s"] = nums[-1]
        if "HSolverLCAO" in line and "solve" in line:
            nums = parse_numeric_tail(line, 1)
            if nums:
                result["hsolver_s"] = nums[-1]
        if "TransMAT2CCS" in line:
            nums = parse_numeric_tail(line, 1)
            if nums:
                result["transmat2ccs_s"] = nums[-1]
    return result


def run_command(cmd: list[str], cwd: Path, env: dict[str, str], timeout_s: int) -> tuple[int, bool, float]:
    stdout = cwd / "abacus_stdout.log"
    stderr = cwd / "abacus_stderr.log"
    start = time.monotonic()
    with stdout.open("w", encoding="utf-8") as out, stderr.open("w", encoding="utf-8") as err:
        proc = subprocess.Popen(cmd, cwd=cwd, env=env, stdout=out, stderr=err, start_new_session=True)
        try:
            returncode = proc.wait(timeout=timeout_s)
            timed_out = False
        except subprocess.TimeoutExpired:
            timed_out = True
            os.killpg(proc.pid, signal.SIGTERM)
            try:
                proc.wait(timeout=30)
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid, signal.SIGKILL)
                proc.wait()
            returncode = 124
    elapsed = time.monotonic() - start
    return returncode, timed_out, elapsed


def result_path(run_root: Path, case: str, solver: str, np: int) -> Path:
    return run_root / case / solver / f"np{np}" / "result.json"


def load_result(path: Path) -> dict[str, object] | None:
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def write_summaries(run_root: Path, cases: list[str], nps: list[int], solvers: list[str]) -> None:
    rows: list[dict[str, object]] = []
    by_key: dict[tuple[str, str, int], dict[str, object]] = {}
    for case in cases:
        for solver in solvers:
            for np in nps:
                data = load_result(result_path(run_root, case, solver, np))
                if data is None:
                    continue
                by_key[(case, solver, np)] = data
                rows.append(data)

    for data in rows:
        case = str(data["case"])
        solver = str(data["solver"])
        np = int(data["np"])
        data["delta_e_ev"] = None
        data["delta_e_ry"] = None
        data["energy_pass"] = None
        data["delta_force_max"] = None
        data["force_pass"] = None
        data["delta_stress_max"] = None
        data["stress_pass"] = None
        if solver != "pexsi":
            continue
        ref = by_key.get((case, "scalapack_gvx", np))
        if not ref:
            continue
        energy = data.get("energy_ev")
        ref_energy = ref.get("energy_ev")
        if isinstance(energy, (int, float)) and isinstance(ref_energy, (int, float)):
            delta_ev = float(energy) - float(ref_energy)
            data["delta_e_ev"] = delta_ev
            data["delta_e_ry"] = delta_ev / RY_TO_EV
            data["energy_pass"] = abs(delta_ev / RY_TO_EV) <= ENERGY_TOL_RY
        force = data.get("max_force_abs")
        ref_force = ref.get("max_force_abs")
        if isinstance(force, (int, float)) and isinstance(ref_force, (int, float)):
            delta_force = abs(float(force) - float(ref_force))
            data["delta_force_max"] = delta_force
            data["force_pass"] = delta_force <= FORCE_TOL_EV_A
        stress = data.get("max_stress_abs")
        ref_stress = ref.get("max_stress_abs")
        if isinstance(stress, (int, float)) and isinstance(ref_stress, (int, float)):
            delta_stress = abs(float(stress) - float(ref_stress))
            data["delta_stress_max"] = delta_stress
            data["stress_pass"] = delta_stress <= STRESS_TOL_GPA

    summary_json = run_root / "summary.json"
    write_text(summary_json, json.dumps(rows, indent=2, sort_keys=True) + "\n")

    headers = [
        "case",
        "solver",
        "np",
        "status",
        "returncode",
        "wall_s",
        "energy_ev",
        "delta_e_ev",
        "delta_e_ry",
        "energy_pass",
        "max_force_abs",
        "delta_force_max",
        "force_pass",
        "max_stress_abs",
        "delta_stress_max",
        "stress_pass",
        "total_time_s",
        "hsolver_s",
        "pexsidft_s",
        "transmat2ccs_s",
        "case_dir",
    ]
    lines = ["\t".join(headers)]
    for data in sorted(rows, key=lambda x: (str(x["case"]), int(x["np"]), str(x["solver"]))):
        lines.append("\t".join(str(data.get(h, "")) for h in headers))
    write_text(run_root / "summary.tsv", "\n".join(lines) + "\n")


def parse_csv_arg(value: str, allowed: Iterable[str] | None = None) -> list[str]:
    items = [item.strip() for item in value.split(",") if item.strip()]
    if allowed is not None:
        allowed_set = set(allowed)
        bad = [item for item in items if item not in allowed_set]
        if bad:
            raise argparse.ArgumentTypeError(f"unknown value(s): {', '.join(bad)}")
    return items


def main(argv: list[str] | None = None) -> int:
    repo = repo_root()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, default=None)
    parser.add_argument("--abacus-bin", type=Path, default=repo / "build-pexsi-conda/abacus_basic_para")
    parser.add_argument("--mpirun", default="mpirun")
    parser.add_argument("--timeout", type=int, default=3600)
    parser.add_argument("--pexsi-npole", type=int, default=120)
    parser.add_argument("--pexsi-elec-thr", type=float, default=0.001)
    parser.add_argument("--pexsi-mu-thr", type=float, default=0.05)
    parser.add_argument("--np-list", default=",".join(str(x) for x in DEFAULT_NP))
    parser.add_argument("--solvers", default=",".join(DEFAULT_SOLVERS))
    parser.add_argument("--cases", default=",".join(CASES.keys()))
    parser.add_argument("--force", action="store_true", help="rerun jobs with an existing result.json")
    parser.add_argument("--prepare-only", action="store_true", help="only generate input directories")
    args = parser.parse_args(argv)

    cases = parse_csv_arg(args.cases, CASES.keys())
    solvers = parse_csv_arg(args.solvers, DEFAULT_SOLVERS)
    nps = [int(x) for x in parse_csv_arg(args.np_list)]
    run_root = args.run_root or Path(f"/tmp/abacus_pexsi_51_scaling_{timestamp()}")
    run_root.mkdir(parents=True, exist_ok=True)

    meta = {
        "created_at": _dt.datetime.now().isoformat(timespec="seconds"),
        "repo": str(repo),
        "abacus_bin": str(args.abacus_bin),
        "timeout_s": args.timeout,
        "cases": cases,
        "solvers": solvers,
        "np_list": nps,
        "pexsi_npole": args.pexsi_npole,
        "pexsi_elec_thr": args.pexsi_elec_thr,
        "pexsi_mu_thr": args.pexsi_mu_thr,
        "energy_tolerance_ry": ENERGY_TOL_RY,
        "force_tolerance_ev_per_angstrom": FORCE_TOL_EV_A,
        "stress_tolerance_gpa": STRESS_TOL_GPA,
    }
    write_text(run_root / "metadata.json", json.dumps(meta, indent=2, sort_keys=True) + "\n")

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = env.get("OMP_NUM_THREADS", "1")
    env["OPENBLAS_NUM_THREADS"] = env.get("OPENBLAS_NUM_THREADS", "1")
    env["MKL_NUM_THREADS"] = env.get("MKL_NUM_THREADS", "1")
    env["ABACUS_PEXSI_TRACE"] = env.get("ABACUS_PEXSI_TRACE", "stage")

    conda_bin = repo / ".conda/pexsi-env/bin"
    conda_lib = repo / ".conda/pexsi-env/lib"
    if conda_bin.exists():
        env["PATH"] = f"{conda_bin}:{env.get('PATH', '')}"
    if conda_lib.exists():
        env["LD_LIBRARY_PATH"] = f"{conda_lib}:{env.get('LD_LIBRARY_PATH', '')}"

    for case_name in cases:
        case = CASES[case_name]
        for np in nps:
            for solver in solvers:
                job_dir = run_root / case.name / solver / f"np{np}"
                result_file = job_dir / "result.json"
                if result_file.exists() and not args.force:
                    print(f"skip existing {case.name} {solver} np={np}")
                    continue

                if job_dir.exists() and args.force:
                    shutil.rmtree(job_dir)
                job_dir.mkdir(parents=True, exist_ok=True)
                prepare_structure(repo, case, job_dir)
                suffix = f"p51_{case.name}_{solver}_np{np}"
                write_text(
                    job_dir / "INPUT",
                    input_text(
                        repo,
                        case,
                        solver,
                        np,
                        suffix,
                        args.pexsi_npole,
                        args.pexsi_elec_thr,
                        args.pexsi_mu_thr,
                    ),
                )
                if args.prepare_only:
                    continue

                cmd = [args.mpirun, "-np", str(np), str(args.abacus_bin)]
                print(f"run case={case.name} solver={solver} np={np} timeout={args.timeout}s", flush=True)
                returncode, timed_out, wall_s = run_command(cmd, job_dir, env, args.timeout)
                parsed = parse_log(job_dir, suffix)
                status = "timeout" if timed_out else ("converged" if parsed["converged"] else "failed")
                data: dict[str, object] = {
                    "case": case.name,
                    "natom": case.natom,
                    "kmesh": "x".join(str(x) for x in case.kmesh),
                    "solver": solver,
                    "np": np,
                    "status": status,
                    "returncode": returncode,
                    "timed_out": timed_out,
                    "wall_s": wall_s,
                    "case_dir": str(job_dir),
                }
                data.update(parsed)
                write_text(result_file, json.dumps(data, indent=2, sort_keys=True) + "\n")
                write_summaries(run_root, cases, nps, solvers)

    write_summaries(run_root, cases, nps, solvers)
    print(f"summary: {run_root / 'summary.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
