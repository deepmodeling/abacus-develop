#!/usr/bin/env python3
"""T7: collinear (nspin=2) AFM spin-channel check (LCAO, DFT+U, FeO).

The Mulliken output for nspin=2 reports per-atom charge and magnetization for
each spin channel separately (mulliken.txt "Total charge of spin1/spin2" plus
the per-atom "sum lmz" / "total magnetism on atom i"). AFM order must show up
in the channels individually: the two Fe sites must carry opposite-sign
magnetizations of comparable magnitude, the O sites must stay (almost)
non-magnetic, and the spin-channel charges must add up to the total charge.
Checking the channels separately (rather than only the total charge/force)
guards against an accidental cancellation of the two channel errors.

Each key is a pass/fail flag: 0.0 = pass, 1.0 = fail. Measured values are
written to stderr for diagnostics.

Usage: check_extra.py [abacus_bin] [np] [case_dir]
"""

import os
import re
import sys

KEYS = [
    "t7_converged_pass",
    "t7_spin_sum_pass",
    "t7_charge_sum_pass",
    "t7_fe_opposite_sign_pass",
    "t7_fe_mag_ratio_pass",
    "t7_o_small_pass",
    "t7_mul_rho_pass",
]

TOL_CHARGE = 1e-3       # electrons
TOL_MUL_RHO = 0.05      # uB: sum(Mulliken per-atom m) vs rho cell magnetization
FE_MIN = 1.0            # uB, Fe sites must be non-vanishing
RATIO_LO, RATIO_HI = 0.3, 3.0
O_MAX = 0.5             # uB, O sites stay (almost) non-magnetic


def fail(msg):
    for k in KEYS:
        print(k + " 1.0", flush=True)
    sys.stderr.write("t7_check.py: " + msg + "\n")
    sys.exit(1)


def parse_mulliken(path):
    text = open(path).read()
    total = None
    spins = {}
    atoms = {}  # idx -> {"chg": float, "mag": float}
    for line in text.splitlines():
        m = re.match(r"^\s*Total charge\s+([-+0-9.eE]+)\s*$", line)
        if m:
            total = float(m.group(1))
            continue
        m = re.match(r"^\s*Total charge of spin(\d+)\s+([-+0-9.eE]+)", line)
        if m:
            spins[int(m.group(1))] = float(m.group(2))
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
    return total, spins, atoms


def parse_rho_mag(path):
    text = open(path).read()
    if not re.search(r"#SCF IS CONVERGED#", text):
        return None, False
    lines = [ln for ln in text.splitlines() if "Total magnetism (Bohr mag/cell)" in ln]
    if not lines:
        return None, True
    # nspin=2: scalar cell magnetization "Total magnetism (Bohr mag/cell) = v"
    m = re.search(r"=\s*([-+0-9.eE]+)\s*$", lines[-1])
    if m:
        return float(m.group(1)), True
    return None, True


def main():
    case_dir = os.path.abspath(sys.argv[3]) if len(sys.argv) > 3 else os.getcwd()
    mul_path = os.path.join(case_dir, "OUT.autotest", "mulliken.txt")
    scf_path = os.path.join(case_dir, "OUT.autotest", "running_scf.log")
    if not (os.path.isfile(mul_path) and os.path.isfile(scf_path)):
        fail("OUT.autotest/mulliken.txt or running_scf.log missing")

    total, spins, atoms = parse_mulliken(mul_path)
    if total is None or not atoms:
        fail("could not parse mulliken.txt")
    rho_mag, converged = parse_rho_mag(scf_path)

    results = {}
    results["t7_converged_pass"] = 0.0 if converged else 1.0

    # 1. spin-channel sum rule: spin1 + spin2 == total charge
    if len(spins) == 2:
        results["t7_spin_sum_pass"] = 0.0 if abs(spins[1] + spins[2] - total) < TOL_CHARGE else 1.0
    else:
        results["t7_spin_sum_pass"] = 1.0

    # 2. per-atom charge sum rule
    qsum = sum(a["chg"] for a in atoms.values())
    results["t7_charge_sum_pass"] = 0.0 if abs(qsum - total) < TOL_CHARGE else 1.0

    # 3-6. channel-wise AFM physics (index order follows STRU: Fe, Fe, O, O)
    mags = []
    for idx in sorted(atoms):
        m = atoms[idx].get("mag")
        if m is None or len(m) != 1:
            fail("nspin=2 scalar magnetism expected for every atom")
        mags.append(m[0])
    if len(mags) != 4:
        fail("expected 4 atoms (2 Fe, 2 O)")
    m_fe1, m_fe2, m_o1, m_o2 = mags

    results["t7_fe_opposite_sign_pass"] = (
        0.0 if m_fe1 > FE_MIN and m_fe2 < -FE_MIN else 1.0)
    ratio = abs(m_fe1) / abs(m_fe2) if abs(m_fe2) > 1e-8 else 1e6
    results["t7_fe_mag_ratio_pass"] = (
        0.0 if RATIO_LO <= ratio <= RATIO_HI else 1.0)
    results["t7_o_small_pass"] = (
        0.0 if abs(m_o1) < O_MAX and abs(m_o2) < O_MAX else 1.0)

    # 7. Mulliken cell sum vs rho cell magnetization
    if rho_mag is None:
        results["t7_mul_rho_pass"] = 1.0
    else:
        msum = sum(mags)
        results["t7_mul_rho_pass"] = 0.0 if abs(msum - rho_mag) < TOL_MUL_RHO else 1.0

    for k in KEYS:
        print(k + " %.1f" % results[k], flush=True)
    sys.stderr.write(
        "t7_check.py: total=%.6f spin1=%.6f spin2=%.6f sum_atoms=%.6f | "
        "m_Fe1=%.6f m_Fe2=%.6f m_O1=%.6f m_O2=%.6f | "
        "sum_mulliken=%.6f rho_cell=%.6f converged=%s\n"
        % (total, spins.get(1, float("nan")), spins.get(2, float("nan")), qsum,
           m_fe1, m_fe2, m_o1, m_o2, sum(mags),
           rho_mag if rho_mag is not None else float("nan"), bool(converged)))


if __name__ == "__main__":
    main()
