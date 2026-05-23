#!/usr/bin/env python3
"""Validate ABACUS LCAO 10-case outputs against the official energy baseline."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


REFERENCE_ENERGIES_EV = {
    "001": ("001_4GaAs", -7836.15655820),
    "002": ("002_C2H6O", -665.55500011),
    "003": ("003_4MoS2", -9699.00659511),
    "004": ("004_12Pt111", -39600.74431945),
    "005": ("005_3BaTiO3", -10749.40029987),
    "006": ("006_16Na", -18524.76009654),
    "007": ("007_27Fe", -86945.53067515),
    "008": ("008_32H2O", -14950.45084341),
    "009": ("009_Battery108", -124070.65411079),
    "010": ("010_216Si", -23123.69990200),
}

FINAL_ETOT_RE = re.compile(r"!\s*FINAL_ETOT_IS\s+([-+0-9.eE]+)\s+eV")


def find_case_dir(run_root: Path, case_id: str) -> Path | None:
    matches = sorted(path for path in run_root.iterdir() if path.is_dir() and path.name.startswith(f"{case_id}_"))
    return matches[0] if matches else None


def find_running_log(case_dir: Path) -> Path | None:
    matches = sorted(case_dir.glob("OUT.*/running_scf.log"))
    return matches[0] if matches else None


def parse_log(log_path: Path) -> tuple[bool, float | None]:
    text = log_path.read_text(errors="replace")
    converged = "#SCF IS CONVERGED#" in text
    matches = FINAL_ETOT_RE.findall(text)
    final_energy = float(matches[-1]) if matches else None
    return converged, final_energy


def validate(run_root: Path, tolerance_ev: float) -> int:
    if not run_root.is_dir():
        print(f"ERROR: run root does not exist: {run_root}", file=sys.stderr)
        return 2

    failures: list[str] = []
    print(f"{'case':<18} {'status':<10} {'energy(eV)':>18} {'ref(eV)':>18} {'|dE|(eV)':>12}")
    for case_id, (case_name, reference_energy) in REFERENCE_ENERGIES_EV.items():
        case_dir = find_case_dir(run_root, case_id)
        if case_dir is None:
            failures.append(f"{case_name}: missing case directory")
            print(f"{case_name:<18} {'MISSING':<10}")
            continue

        log_path = find_running_log(case_dir)
        if log_path is None:
            failures.append(f"{case_name}: missing OUT.*/running_scf.log")
            print(f"{case_name:<18} {'NO_LOG':<10}")
            continue

        converged, final_energy = parse_log(log_path)
        if final_energy is None:
            failures.append(f"{case_name}: missing FINAL_ETOT_IS")
            print(f"{case_name:<18} {'NO_ETOT':<10}")
            continue

        delta = abs(final_energy - reference_energy)
        status = "PASS" if converged and delta <= tolerance_ev else "FAIL"
        if status != "PASS":
            reason = "not converged" if not converged else f"|dE|={delta:.3e} > {tolerance_ev:.3e}"
            failures.append(f"{case_name}: {reason}")
        print(f"{case_name:<18} {status:<10} {final_energy:18.10f} {reference_energy:18.10f} {delta:12.3e}")

    if failures:
        print("\nFailures:", file=sys.stderr)
        for failure in failures:
            print(f"  - {failure}", file=sys.stderr)
        return 1
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_root", type=Path, help="directory containing the 001_* ... 010_* LCAO case outputs")
    parser.add_argument(
        "--tolerance-ev",
        type=float,
        default=5.0e-7,
        help="maximum allowed absolute total-energy error in eV",
    )
    args = parser.parse_args()
    return validate(args.run_root, args.tolerance_ev)


if __name__ == "__main__":
    raise SystemExit(main())
