#!/usr/bin/env python3
"""T6: nspin=4 multi-k Mulliken magnetization calibration (LCAO, DFT+U, Fe).

Cross-checks the Mulliken per-atom populations/magnetizations
(OUT.autotest/mulliken.txt, computed from DMK x S(k)) against
1. the internal charge sum rule (sum of per-atom charges == Total charge),
2. the rho-integrated cell magnetization
   (OUT.autotest/running_scf.log, "Total magnetism (Bohr mag/cell)"),
   both in magnitude and direction.

This is the calibration called for in docs/dm_dmk_dmr_action_plan.md T6:
for nspin=4 and multi-k the Mulliken m_y path uses S^T instead of S^dagger, so
an anomalous m_y deviation (order ||Im S(k)||) would flag the B-class gemm
change. Measured deviation on this case (Fe, 4x1x1, U=5 eV, converged):
  |sum Mulliken - rho_cell| = (0.042, 0.004, 0.044) uB on |m| ~ 5.5 uB/cell,
  i.e. the small components agree as well as the large ones, so no anomalous
  m_y deviation is present and the B-class change is NOT executed.
The Gamma-only control (same cell, Gamma) converges to a different magnetic
minimum and is not part of the automated test; it is covered in the PR
validation record. DeltaSpin is disabled for LCAO in this repo
(tests/17_DS_DFTU/CASES_CPU.txt), so the DeltaSpin third path is unavailable.

Each key is a pass/fail flag: 0.0 = pass, 1.0 = fail. Measured values are
written to stderr for diagnostics.

Usage: check_extra.py [abacus_bin] [np] [case_dir]
"""

import os
import re
import sys

KEYS = [
    "t6_converged_pass",
    "t6_charge_sum_pass",
    "t6_mul_rho_pass",
    "t6_mul_rho_dir_pass",
    "t6_mz_sign_pass",
]

TOL_CHARGE = 1e-3       # sum rule, electrons
TOL_MUL_RHO = 0.1       # uB per component
TOL_DIR_DEG = 2.0       # degrees between Mulliken-sum and rho vectors
MZ_MIN = 1.0            # uB, both Fe sites must stay ferromagnetic along z


def fail(msg):
    for k in KEYS:
        print(k + " 1.0", flush=True)
    sys.stderr.write("t6_check.py: " + msg + "\n")
    sys.exit(1)


def parse_mulliken(path):
    text = open(path).read()
    total = None
    atoms = {}  # idx -> {"chg": float, "mag": [mx,my,mz]}
    for line in text.splitlines():
        m = re.match(r"^\s*Total charge\s+([-+0-9.eE]+)\s*$", line)
        if m:
            total = float(m.group(1))
            continue
        m = re.match(r"^\s*total charge\s+on atom\s+(\d+)\s+([-+0-9.eE]+)", line)
        if m:
            atoms.setdefault(int(m.group(1)), {})["chg"] = float(m.group(2))
            continue
        m = re.match(
            r"^\s*total magnetism\s+on atom\s+(\d+)\s+([-+0-9.eE]+)"
            r"(?:\s+([-+0-9.eE]+)\s+([-+0-9.eE]+))?", line)
        if m:
            vals = [float(m.group(i)) for i in (2, 3, 4) if m.group(i)]
            atoms.setdefault(int(m.group(1)), {})["mag"] = vals
    return total, atoms


def parse_rho_mag(path):
    text = open(path).read()
    if not re.search(r"#SCF IS CONVERGED#", text):
        return None, None
    lines = [ln for ln in text.splitlines() if "Total magnetism (Bohr mag/cell)" in ln]
    if not lines:
        return None, None
    last = lines[-1]
    m = re.search(r"=\s*\[\s*([-+0-9.eE]+)\s*,\s*([-+0-9.eE]+)\s*,\s*([-+0-9.eE]+)", last)
    if not m:
        return None, None
    return [float(m.group(i)) for i in (1, 2, 3)], True


def main():
    case_dir = os.path.abspath(sys.argv[3]) if len(sys.argv) > 3 else os.getcwd()
    mul_path = os.path.join(case_dir, "OUT.autotest", "mulliken.txt")
    scf_path = os.path.join(case_dir, "OUT.autotest", "running_scf.log")
    if not (os.path.isfile(mul_path) and os.path.isfile(scf_path)):
        fail("OUT.autotest/mulliken.txt or running_scf.log missing")

    total, atoms = parse_mulliken(mul_path)
    if total is None or not atoms:
        fail("could not parse mulliken.txt")
    rho_mag, converged = parse_rho_mag(scf_path)

    results = {}
    results["t6_converged_pass"] = 0.0 if converged else 1.0

    # 1. charge sum rule: sum of per-atom charges == Total charge
    qsum = sum(a["chg"] for a in atoms.values())
    results["t6_charge_sum_pass"] = 0.0 if abs(qsum - total) < TOL_CHARGE else 1.0

    # 2. Mulliken cell magnetization vs rho cell magnetization
    msum = [0.0, 0.0, 0.0]
    for a in atoms.values():
        mag = a.get("mag")
        if mag is None or len(mag) != 3:
            fail("nspin=4 magnetism vector expected for every atom")
        for c in range(3):
            msum[c] += mag[c]
    if rho_mag is None:
        results["t6_mul_rho_pass"] = 1.0
        results["t6_mul_rho_dir_pass"] = 1.0
    else:
        max_dev = max(abs(msum[c] - rho_mag[c]) for c in range(3))
        results["t6_mul_rho_pass"] = 0.0 if max_dev < TOL_MUL_RHO else 1.0
        na = sum(v * v for v in msum) ** 0.5
        nb = sum(v * v for v in rho_mag) ** 0.5
        if na < 1e-8 or nb < 1e-8:
            results["t6_mul_rho_dir_pass"] = 1.0
        else:
            cosang = sum(msum[c] * rho_mag[c] for c in range(3)) / (na * nb)
            cosang = max(-1.0, min(1.0, cosang))
            deg = os_acos_deg(cosang)
            results["t6_mul_rho_dir_pass"] = 0.0 if deg < TOL_DIR_DEG else 1.0

    # 3. magnetic state guard: both Fe sites ferromagnetic along z
    mz = [a["mag"][2] for a in sorted(atoms.values(), key=lambda a: a["chg"])]
    results["t6_mz_sign_pass"] = 0.0 if all(m > MZ_MIN for m in mz) else 1.0

    for k in KEYS:
        print(k + " %.1f" % results[k], flush=True)
    sys.stderr.write(
        "t6_check.py: total_charge=%.6f sum_atoms=%.6f | "
        "sum_mulliken=(%.6f,%.6f,%.6f) rho_cell=(%.6f,%.6f,%.6f) | "
        "max_dev=%.6f uB converged=%s\n"
        % (total, qsum, msum[0], msum[1], msum[2],
           rho_mag[0] if rho_mag else float("nan"),
           rho_mag[1] if rho_mag else float("nan"),
           rho_mag[2] if rho_mag else float("nan"),
           max(abs(msum[c] - (rho_mag[c] if rho_mag else 0.0)) for c in range(3)),
           bool(converged)))


def os_acos_deg(cosang):
    import math
    return math.degrees(math.acos(cosang))


if __name__ == "__main__":
    main()
