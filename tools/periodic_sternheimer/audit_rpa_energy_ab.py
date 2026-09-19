#!/usr/bin/env python3
"""Audit equality of two LibRPA Sternheimer correlation-energy runs."""

import argparse
import json
import math
from pathlib import Path
import re


HARTREE_TO_MEV = 27.211386245988 * 1000.0
ENERGY_PATTERN = re.compile(
    r"^\| Total Sternheimer EcRPA:\s*([+\-0-9.eE]+)\s+([+\-0-9.eE]+)\s*$",
    re.MULTILINE,
)


def read_energy(path: Path) -> tuple[float, float]:
    text = path.read_text(encoding="utf-8", errors="replace")
    if "libRPA finished successfully" not in text:
        raise ValueError(f"LibRPA did not finish successfully: {path}")
    matches = ENERGY_PATTERN.findall(text)
    if len(matches) != 1:
        raise ValueError(f"expected one Total Sternheimer EcRPA record: {path}")
    real, imaginary = (float(value) for value in matches[0])
    if not math.isfinite(real) or not math.isfinite(imaginary):
        raise ValueError(f"non-finite Sternheimer EcRPA: {path}")
    return real, imaginary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--control", type=Path, required=True)
    parser.add_argument("--symmetry", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tolerance-mev", type=float, default=0.1)
    args = parser.parse_args()
    if not math.isfinite(args.tolerance_mev) or args.tolerance_mev < 0.0:
        raise ValueError("tolerance must be finite and non-negative")

    control_real, control_imaginary = read_energy(args.control)
    symmetry_real, symmetry_imaginary = read_energy(args.symmetry)
    difference_mev = abs(symmetry_real - control_real) * HARTREE_TO_MEV
    report = {
        "control_log": str(args.control),
        "symmetry_log": str(args.symmetry),
        "control_ecrpa_hartree": control_real,
        "symmetry_ecrpa_hartree": symmetry_real,
        "control_ecrpa_imaginary_hartree": control_imaginary,
        "symmetry_ecrpa_imaginary_hartree": symmetry_imaginary,
        "absolute_difference_mev": difference_mev,
        "tolerance_mev": args.tolerance_mev,
        "energy_gate_passed": difference_mev <= args.tolerance_mev,
        "physical_result": False,
    }
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="ascii")
    if not report["energy_gate_passed"]:
        raise SystemExit("RPA energy A/B gate failed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
