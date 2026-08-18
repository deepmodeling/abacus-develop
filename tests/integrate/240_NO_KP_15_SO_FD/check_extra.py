#!/usr/bin/env python3
"""T4: SOC nonlocal force finite-difference validation.

Compares the analytic LCAO force against a central finite difference of the
total energy for every atom and Cartesian direction:

    F_fd = -(E(+delta) - E(-delta)) / (2*delta),   delta = 1e-3 Bohr

and reports two pass/fail keys for result.out (0.0 = pass, 1.0 = fail):

    fd_force_pass  max |F_analytic - F_fd| over all atoms/directions < 1e-4 eV/Bohr
    acoustic_sum_pass  max |sum_i F_analytic(i,d)| over directions d < 1e-4 eV/Bohr

The acoustic sum (translational invariance) is the first quantity that breaks
if the R pairing of the nonlocal force is wrong; the per-component comparison
prevents accidental cancellations between force terms. The measured max_dev is
written to stderr for diagnostics (reference ~5.8e-5 eV/Bohr on this case,
deterministic at scf_thr=1e-8); comparing raw FD metrics run-to-run at 1e-6
would be machine-fragile, so the harness comparison uses the pass/fail flags.

NOTE on sensitivity: this non-magnetic GaAs SOC case exercises the nonlocal
SOC force path, but its force is (numerically) identical under the
e^{-ikR}/e^{+ikR} DMR phase choice (closed-trace protection, see
docs/dm_dmk_dmr_action_plan.md). The sentinel for the DMR Fourier sign itself
is the unit test DMTest.T1_fourier_round_trip
(source/source_estate/module_dm/test/test_cal_dm_R.cpp), which turns red when
the sign is flipped. This case guards FD-consistency of the analytic force,
including the SOC nonlocal contribution.

Usage: check_extra.py [abacus_bin] [np] [case_dir]
Keys are printed to stdout so the integrate harness (Autotest.sh) appends them
to result.out / result.ref.
"""

import math
import os
import re
import shutil
import subprocess
import sys
import tempfile

BOHR_TO_ANG = 0.529177210903


def fail(msg):
    print("fd_force_pass 1.0", flush=True)
    print("acoustic_sum_pass 1.0", flush=True)
    sys.stderr.write("fd_force.py: " + msg + "\n")
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


def parse_stru(text):
    """Return (lattice_matrix, atom_list) from an ABACUS STRU file.

    lattice_matrix: 3x3 rows in Bohr.
    atom_list: list of dicts {label, idx, frac (3 floats), extra (rest of line)}.
    """
    lines = text.splitlines()
    i = 0
    n = len(lines)

    def next_data():
        nonlocal i
        while i < n:
            ln = lines[i].strip()
            i += 1
            if ln and not ln.startswith("#"):
                return ln
        return None

    # ATOMIC_SPECIES block
    while True:
        ln = next_data()
        if ln is None:
            raise RuntimeError("ATOMIC_SPECIES not found")
        if ln.startswith("ATOMIC_SPECIES"):
            break
    lat_const = 1.0
    while True:
        ln = next_data()
        if ln is None:
            raise RuntimeError("LATTICE_CONSTANT not found")
        if ln.startswith("LATTICE_CONSTANT"):
            lat_const = float(next_data())
            break

    # skip NUMERICAL_ORBITAL if present
    while True:
        ln = next_data()
        if ln is None:
            raise RuntimeError("LATTICE_VECTORS not found")
        if ln.startswith("LATTICE_VECTORS"):
            break
    latvec = []
    for _ in range(3):
        latvec.append([float(x) for x in next_data().split()[:3]])
    lat = [[lat_const * v for v in row] for row in latvec]

    while True:
        ln = next_data()
        if ln is None:
            raise RuntimeError("ATOMIC_POSITIONS not found")
        if ln.startswith("ATOMIC_POSITIONS"):
            break
    coord_type = next_data()

    atoms = []
    while True:
        ln = next_data()
        if ln is None:
            break
        if ln.startswith("#") or ln == "":
            continue
        label = ln.strip()
        if next_data() is None:
            break  # magnetism line
        nat_line = next_data()
        if nat_line is None:
            break
        nat = int(nat_line)
        for ia in range(nat):
            pos = next_data()
            if pos is None:
                raise RuntimeError("unexpected end of STRU atom block")
            fields = pos.split()
            frac = [float(fields[0]), float(fields[1]), float(fields[2])]
            atoms.append({"label": label, "idx": ia, "frac": frac, "extra": fields[3:]})
    return lat, atoms


def rewrite_stru(text, lat, atoms, delta_vec_bohr, target_idx):
    """Return STRU text with atom positions shifted by delta (Cartesian, Bohr)."""
    inv = invert3(lat)
    # fraction-coordinate shift dfrac = delta . A^{-1} (row-vector convention:
    # Cartesian delta = dfrac . A). A is not symmetric for non-orthogonal
    # lattices, so the index order below matters.
    dfrac = [0.0, 0.0, 0.0]
    for a in range(3):
        for b in range(3):
            dfrac[a] += delta_vec_bohr[b] * inv[b][a]

    lines = text.splitlines()
    out = []
    i = 0
    n = len(lines)
    atom_iter = iter(atoms)
    shifted = {}
    for idx, atom in enumerate(atoms):
        if idx == target_idx:
            shifted[id(atom)] = [atom["frac"][a] + dfrac[a] for a in range(3)]
        else:
            shifted[id(atom)] = list(atom["frac"])

    # copy everything up to and including the ATOMIC_POSITIONS header line
    while i < n:
        out.append(lines[i])
        if lines[i].strip().startswith("ATOMIC_POSITIONS"):
            i += 1
            if i < n:
                out.append(lines[i])  # Direct/Cartesian
                i += 1
            break
        i += 1

    # element blocks: label, magnetism, nat, nat coordinate lines
    while i < n:
        ln = lines[i]
        stripped = ln.strip()
        if stripped == "" or stripped.startswith("#"):
            out.append(ln)
            i += 1
            continue
        label = stripped
        out.append(ln)
        i += 1
        # magnetism line
        while i < n and (lines[i].strip() == "" or lines[i].strip().startswith("#")):
            out.append(lines[i])
            i += 1
        if i >= n:
            break
        out.append(lines[i])
        i += 1
        # atom count line
        while i < n and (lines[i].strip() == "" or lines[i].strip().startswith("#")):
            out.append(lines[i])
            i += 1
        if i >= n:
            break
        out.append(lines[i])
        nat = int(lines[i].strip())
        i += 1
        for ia in range(nat):
            while i < n and (lines[i].strip() == "" or lines[i].strip().startswith("#")):
                out.append(lines[i])
                i += 1
            if i >= n:
                break
            pos = lines[i].split()
            target = None
            for atom in atoms:
                if atom["label"] == label and atom["idx"] == ia:
                    target = atom
                    break
            if target is None:
                raise RuntimeError("atom not found: %s %d" % (label, ia))
            frac = shifted[id(target)]
            fields = ["%14.9f" % v for v in frac]
            fields += pos[3:]
            out.append(" ".join(fields))
            i += 1
    return "\n".join(out)


def invert3(m):
    a, b, c = m[0][0], m[0][1], m[0][2]
    d, e, f = m[1][0], m[1][1], m[1][2]
    g, h, k = m[2][0], m[2][1], m[2][2]
    det = a * (e * k - f * h) - b * (d * k - f * g) + c * (d * h - e * g)
    if abs(det) < 1e-30:
        raise RuntimeError("singular lattice matrix")
    inv = [
        [(e * k - f * h) / det, (c * h - b * k) / det, (b * f - c * e) / det],
        [(f * g - d * k) / det, (a * k - c * g) / det, (c * d - a * f) / det],
        [(d * h - e * g) / det, (b * g - a * h) / det, (a * e - b * d) / det],
    ]
    return inv


def run_scf(workdir, abacus, np, tag):
    """Run ABACUS SCF in workdir, return final total energy in eV."""
    log = os.path.join(workdir, "run_%s.log" % tag)
    with open(log, "w") as f:
        subprocess.run(["mpirun", "-np", str(np), abacus],
                       cwd=workdir, stdout=f, stderr=subprocess.STDOUT)
    # find the OUT.<suffix> directory
    inp = open(os.path.join(workdir, "INPUT")).read()
    suffix = get_input_value(inp, "suffix", "autotest")
    outdir = os.path.join(workdir, "OUT." + suffix)
    scf = os.path.join(outdir, "running_scf.log")
    if not os.path.isfile(scf):
        raise RuntimeError("missing " + scf)
    text = open(scf).read()
    m = re.search(r"!FINAL_ETOT_IS\s+([-0-9.eE+]+)", text)
    if not m:
        raise RuntimeError("FINAL_ETOT_IS not found in " + scf)
    return float(m.group(1))


def read_analytic_force(case_dir, suffix):
    scf = os.path.join(case_dir, "OUT." + suffix, "running_scf.log")
    if not os.path.isfile(scf):
        return None
    text = open(scf).read()
    m = re.search(r"#TOTAL-FORCE.*?\n(.*?)(?:\n#|$)", text, re.S)
    if not m:
        return None
    block = m.group(1)
    forces = []
    for line in block.splitlines():
        parts = line.split()
        if len(parts) == 4 and re.match(r"^[A-Za-z]+\d+$", parts[0]):
            forces.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return forces


def main():
    abacus = sys.argv[1] if len(sys.argv) > 1 else "abacus"
    np_ = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    case_dir = os.path.abspath(sys.argv[3]) if len(sys.argv) > 3 else os.getcwd()
    delta = 1.0e-3  # Bohr

    inp_path = os.path.join(case_dir, "INPUT")
    stru_path = os.path.join(case_dir, "STRU")
    if not (os.path.isfile(inp_path) and os.path.isfile(stru_path)):
        fail("INPUT/STRU missing in " + case_dir)

    inp_text = open(inp_path).read()
    suffix = get_input_value(inp_text, "suffix", "autotest")
    pseudo_dir = get_input_value(inp_text, "pseudo_dir", "../../PP_ORB")
    orbital_dir = get_input_value(inp_text, "orbital_dir", "../../PP_ORB")
    pseudo_abs = os.path.abspath(os.path.join(case_dir, pseudo_dir))
    orbital_abs = os.path.abspath(os.path.join(case_dir, orbital_dir))

    stru_text = open(stru_path).read()
    lat, atoms = parse_stru(stru_text)
    natom = len(atoms)
    if natom == 0:
        fail("no atoms in STRU")

    # analytic force from the base run (already produced by the harness) or run
    # a fresh base SCF if missing.
    base_forces = read_analytic_force(case_dir, suffix)
    if base_forces is None:
        shutil.copy(inp_path, case_dir + ".INPUT.bak")
        try:
            shutil.copy(stru_path, case_dir + ".STRU.bak")
            run_scf(case_dir, abacus, np_, "base")
            base_forces = read_analytic_force(case_dir, suffix)
        finally:
            shutil.move(case_dir + ".INPUT.bak", inp_path)
            shutil.move(case_dir + ".STRU.bak", stru_path)
    if base_forces is None or len(base_forces) != natom:
        fail("could not read analytic forces (%d atoms)" % natom)

    fd_force = [[0.0, 0.0, 0.0] for _ in range(natom)]
    with tempfile.TemporaryDirectory(prefix="t4_fd_") as tmp:
        for ia in range(natom):
            for d in range(3):
                shift = [0.0, 0.0, 0.0]
                shift[d] = delta
                energies = []
                for sgn in (1.0, -1.0):
                    wd = os.path.join(tmp, "a%d_d%d_%s" % (ia, d, "+" if sgn > 0 else "-"))
                    os.makedirs(wd)
                    inp_mod = re.sub(r"^\s*pseudo_dir.*$", "pseudo_dir " + pseudo_abs,
                                     inp_text, flags=re.M)
                    inp_mod = re.sub(r"^\s*orbital_dir.*$", "orbital_dir " + orbital_abs,
                                     inp_mod, flags=re.M)
                    stru_mod = rewrite_stru(stru_text, lat, atoms,
                                            [sgn * v for v in shift], ia)
                    with open(os.path.join(wd, "INPUT"), "w") as f:
                        f.write(inp_mod)
                    with open(os.path.join(wd, "STRU"), "w") as f:
                        f.write(stru_mod)
                    shutil.copy(os.path.join(case_dir, "KPT"), os.path.join(wd, "KPT"))
                    energies.append(run_scf(wd, abacus, np_, "fd"))
                # F = -dE/dx: FD force = -(E(+d)-E(-d)) / (2*d)
                fd_force[ia][d] = (energies[1] - energies[0]) / (2.0 * delta)

    # convert analytic force eV/Angstrom -> eV/Bohr (1 Bohr = 0.529177 Angstrom,
    # so dE/dx[Bohr] = dE/dx[Ang] * 0.529177)
    max_dev = 0.0
    acoustic = [0.0, 0.0, 0.0]
    for ia in range(natom):
        for d in range(3):
            fa = base_forces[ia][d] * BOHR_TO_ANG
            dev = abs(fa - fd_force[ia][d])
            max_dev = max(max_dev, dev)
            acoustic[d] += fa
    acoustic_max = max(abs(v) for v in acoustic)

    print("fd_force_pass %.1f" % (1.0 if max_dev > 1e-4 else 0.0), flush=True)
    print("acoustic_sum_pass %.1f" % (1.0 if acoustic_max > 1e-4 else 0.0), flush=True)
    sys.stderr.write(
        "fd_force.py: max_dev=%.3e eV/Bohr acoustic=%.3e eV/Bohr (threshold 1e-4)\n"
        % (max_dev, acoustic_max))


if __name__ == "__main__":
    main()
