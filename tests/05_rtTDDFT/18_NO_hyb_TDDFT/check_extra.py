#!/usr/bin/env python3
"""T5: rt-TDDFT current R-pairing / phase checks (multi-k, hybrid gauge Si).

Appended to result.out / result.ref by the Autotest.sh check_extra.py hook.
Each key is a pass/fail flag (0.0 = pass, 1.0 = fail); the measured values are
written to stderr for diagnostics. All sub-runs use a 4x1x1 Gamma-centred mesh,
which contains the non-time-reversal-invariant pair k=(1/4,0,0) / -k=(3/4,0,0).

Checks
------
1. zero_field_current_max: max |j(t)| with td_vext=0. Si is time-reversal
   symmetric, so the total current must vanish to ~1e-8 (the +/-k pairs cancel
   at the level of the total). A nonzero value indicates a broken R pairing /
   phase convention in the DMR -> current contraction.
2. current_path_max_dev: max |j_full - sum_k j_k| over steps/directions with a
   finite field (td_vext=1). The full path contracts the all-k DMR
   (cal_DMR()); the per-k path contracts each cal_DMR(ik) and sums the
   k-weighted contributions. The two constructions must agree to ~1e-10.
3. current_pmk_cancel_max: max |j(k) + j(-k)| over the +/-k pair and steps.
   Each per-k current is not physical alone (single-k DMR is not
   Hermitian-paired), but time reversal gives j(-k) = -j(k); the residual
   must be ~1e-8.

NOTE on sensitivity: the current is a closed Frobenius trace with a full R sum,
so it is protected against the e^{-ikR}/e^{+ikR} DMR phase choice (flipping the
sign leaves all keys at roundoff). These checks guard against gross
pairing/phase breakage (e.g. a mismatched
phase_hybrid in the hybrid gauge). The DMR sign sentinel is
DMTest.T1_fourier_round_trip.

Usage: check_extra.py [abacus_bin] [np] [case_dir]
"""

import os
import re
import shutil
import subprocess
import sys
import tempfile

# ABACUS integration-test convention: OMP_NUM_THREADS=1 keeps the SCF and the
# two DMR current paths bit-deterministic (multi-threaded summation order makes
# the full-k and per-k paths disagree at ~1e-8 instead of roundoff).
os.environ.setdefault("OMP_NUM_THREADS", "1")

KPT_TEXT = "K_POINTS\n0\nGamma\n4 1 1  0 0 0\n"
TOL_PATH = 1e-10
TOL_ZERO = 1e-8
TOL_PMK = 1e-8


def fail(msg):
    print("current_path_max_dev 1.0", flush=True)
    print("zero_field_current_max 1.0", flush=True)
    print("current_pmk_cancel_max 1.0", flush=True)
    sys.stderr.write("check_extra.py: " + msg + "\n")
    sys.exit(1)


def get_input_value(text, key, default=None):
    for line in text.splitlines():
        line = line.strip()
        if line.startswith("#"):
            continue
        parts = line.split()
        if parts and parts[0] == key:
            return parts[1]
    return default


def set_input_value(text, key, value):
    if re.search(rf"^\s*{key}\s", text, flags=re.M):
        return re.sub(rf"^\s*{key}\s+\S+", f"{key}            {value}", text, flags=re.M)
    return text + f"\n{key}            {value}\n"


def read_current_rows(path):
    rows = []
    for line in open(path):
        parts = line.split()
        if len(parts) == 4:
            rows.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return rows


def read_kpoints(running_log):
    """Parse the K-point list (direct coords) from KPT.info in the log dir."""
    ks = []
    in_list = False
    for line in open(running_log):
        if "K-POINTS DIRECT COORDINATES" in line:
            in_list = True
            continue
        if in_list:
            parts = line.split()
            if len(parts) == 5 and parts[0].isdigit():
                ks.append([float(parts[1]), float(parts[2]), float(parts[3])])
            elif ks and len(parts) != 5:
                break
    return ks


def run_case(case_dir, workdir, abacus, np_, mods, kpt_text=KPT_TEXT):
    """Copy the case into workdir, apply INPUT modifications + KPT, run ABACUS,
    and return (current_tot rows, {ik: rows} for per-k currents, k vectors)."""
    os.makedirs(workdir, exist_ok=True)
    inp = open(os.path.join(case_dir, "INPUT")).read()
    inp = inp.replace("../../PP_ORB", os.path.abspath(os.path.join(case_dir, "../../PP_ORB")))
    for key, value in mods.items():
        inp = set_input_value(inp, key, value)
    with open(os.path.join(workdir, "INPUT"), "w") as f:
        f.write(inp)
    for name in ("STRU",):
        shutil.copy(os.path.join(case_dir, name), os.path.join(workdir, name))
    with open(os.path.join(workdir, "KPT"), "w") as f:
        f.write(kpt_text)
    with open(os.path.join(workdir, "run.log"), "w") as f:
        ret = subprocess.run(["mpirun", "-np", str(np_), abacus], cwd=workdir,
                             stdout=f, stderr=subprocess.STDOUT)
    if ret.returncode != 0:
        fail("ABACUS run failed in " + workdir)
    out = os.path.join(workdir, "OUT.autotest")
    tot = os.path.join(out, "current_tot.txt")
    if not os.path.isfile(tot):
        fail("current_tot.txt missing in " + workdir)
    tot_rows = read_current_rows(tot)
    if not tot_rows:
        fail("current_tot.txt empty in " + workdir)
    perk = {}
    for name in sorted(os.listdir(out)):
        m = re.fullmatch(r"current_s(\d+)k(\d+)\.txt", name)
        if m:
            ik = int(m.group(2))
            perk[ik] = read_current_rows(os.path.join(out, name))
    ks = read_kpoints(os.path.join(out, "KPT.info"))
    if not ks:
        fail("could not parse K-point list from KPT.info in " + workdir)
    return tot_rows, perk, ks


def max_diff(a_rows, b_rows):
    if len(a_rows) != len(b_rows):
        fail("step count mismatch between current paths")
    return max(abs(x - y)
               for ra, rb in zip(a_rows, b_rows)
               for x, y in zip(ra, rb))


def main():
    abacus = sys.argv[1] if len(sys.argv) > 1 else "abacus"
    np_ = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    case_dir = os.path.abspath(sys.argv[3]) if len(sys.argv) > 3 else os.getcwd()

    with tempfile.TemporaryDirectory(prefix="t5_current_") as tmp:
        # finite field: nonzero current; compare full-k DMR path vs per-k DMR path
        j_full, _, _ = run_case(case_dir, os.path.join(tmp, "full"), abacus, np_,
                                {"out_current": "1", "out_current_k": "0", "td_vext": "1"})
        j_eachk, perk, ks = run_case(case_dir, os.path.join(tmp, "eachk"), abacus, np_,
                                     {"out_current": "1", "out_current_k": "1", "td_vext": "1"})
        path_dev = max_diff(j_full, j_eachk)

        # sum of per-k files must reproduce the per-k-path total
        sum_dev = 0.0
        if perk:
            nstep = len(next(iter(perk.values())))
            acc = [[0.0, 0.0, 0.0] for _ in range(nstep)]
            for rows in perk.values():
                for step, row in enumerate(rows):
                    if step < nstep:
                        for d in range(3):
                            acc[step][d] += row[d]
            for step, row in enumerate(acc):
                if step < len(j_eachk):
                    sum_dev = max(sum_dev, max(abs(row[d] - j_eachk[step][d]) for d in range(3)))
        if sum_dev > TOL_PATH:
            fail("per-k file sum deviates from current_tot (per-k path): %.3e" % sum_dev)

        # +/-k cancellation on the per-k currents
        pmk_dev = 0.0
        nk = len(ks)
        for i in range(nk):
            for j in range(i + 1, nk):
                if sum(abs(ks[i][d] + ks[j][d]) for d in range(3)) < 1e-10:
                    ri = perk.get(i + 1)
                    rj = perk.get(j + 1)
                    if ri is None or rj is None:
                        fail("missing per-k current file for k pair %d/%d" % (i + 1, j + 1))
                    pmk_dev = max(pmk_dev, max(abs(a + b)
                                               for ra, rb in zip(ri, rj)
                                               for a, b in zip(ra, rb)))

        # zero external field: time-reversal symmetric ground state
        j_zero, _, _ = run_case(case_dir, os.path.join(tmp, "zero"), abacus, np_,
                                {"out_current": "1", "out_current_k": "0", "td_vext": "0"})
        zero_dev = max(abs(v) for row in j_zero for v in row)

    sys.stderr.write("check_extra.py: path_dev=%.3e sum_dev=%.3e pmk_dev=%.3e zero_dev=%.3e\n"
                     % (path_dev, sum_dev, pmk_dev, zero_dev))
    print("current_path_max_dev %.1f" % (1.0 if path_dev > TOL_PATH else 0.0), flush=True)
    print("zero_field_current_max %.1f" % (1.0 if zero_dev > TOL_ZERO else 0.0), flush=True)
    print("current_pmk_cancel_max %.1f" % (1.0 if pmk_dev > TOL_PMK else 0.0), flush=True)


if __name__ == "__main__":
    main()
