"""Compare uneven EXX/ACE k-point pools against a single-pool calculation."""
import argparse
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile


def read_states(path):
    states = {}
    spin = None
    ik = None
    for line in path.read_text().splitlines():
        header = re.search(r"spin=(\d+) k-point=(\d+)/(\d+)", line)
        if header:
            spin = int(header[1])
            ik = int(header[2])
            continue
        row = re.fullmatch(r"\s*(\d+)\s+(\S+)\s+(\S+)\s*", line)
        if spin is not None and row:
            states[(spin, ik, int(row[1]))] = (float(row[2]), float(row[3]))
    assert states, path
    return states


def compare(reference, actual):
    energy_ref, states_ref = reference
    energy, states = actual
    assert states_ref.keys() == states.keys(), "Missing or duplicate band rows"
    energy_delta = abs(energy - energy_ref)
    occupation_delta = max(abs(states[key][1] - states_ref[key][1]) for key in states)
    eigenvalue_delta = max(abs(states[key][0] - states_ref[key][0]) for key in states)
    assert energy_delta < 1e-6, energy_delta
    assert occupation_delta < 1e-6, occupation_delta
    assert eigenvalue_delta < 1e-4, eigenvalue_delta
    return energy_delta, occupation_delta, eigenvalue_delta


def run_case(binary, launcher, source, root, functional, nspin, nk, ranks, kpar, smearing, low_symmetry=False):
    prefix = "triclinic-" if low_symmetry else ""
    case = root / f"{prefix}{functional}-{smearing}-spin{nspin}-nk{nk}-np{ranks}-kpar{kpar}"
    case.mkdir()
    metallic = smearing == "gaussian"
    fixture = "001_PW_UPF100_Al" if metallic else "097_PW_PBE0"
    shutil.copyfile(source / "01_PW" / fixture / "STRU", case / "STRU")
    lattice = "latname fcc\n" if metallic else ""
    if low_symmetry:
        (case / "STRU").write_text("""ATOMIC_SPECIES
Al 26.98 Al.pbe-rrkj.UPF
LATTICE_CONSTANT
7.5
LATTICE_VECTORS
-0.50 0.01 0.52
0.02 0.47 0.49
-0.48 0.51 0.03
ATOMIC_POSITIONS
Direct
Al
0
1
0 0 0 0 0 0
""")
        lattice = ""
    cutoff = 30 if metallic else 10
    nbands = 4 if metallic else 3
    magnetization = "nupdown 1\n" if metallic and nspin == 2 else ""
    transverse_mesh = 3 if metallic else 1
    gamma_extrapolation = "true" if metallic else "false"
    (case / "KPT").write_text(f"K_POINTS\n0\nGamma\n{nk} {transverse_mesh} {transverse_mesh} 0 0 0\n")
    (case / "INPUT").write_text(f"""INPUT_PARAMETERS
suffix test
basis_type pw
calculation scf
pseudo_dir {source / 'PP_ORB'}
nbands {nbands}
init_wfc random
pw_seed 1
ecutwfc {cutoff}
{lattice}
symmetry -1
smearing_method {smearing}
smearing_sigma 0.02
scf_thr 1e-8
scf_nmax 200
mixing_beta 0.2
exx_hybrid_step 200
dft_functional {functional}
nspin {nspin}
{magnetization}
kpar {kpar}
exx_separate_loop true
exx_gamma_extrapolation {gamma_extrapolation}
exx_thr_type density
cal_force 0
cal_stress 0
""")
    env = dict(os.environ, OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    command = [launcher, "-n", str(ranks), str(binary)]
    result = subprocess.run(command, cwd=case, env=env, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=300)
    (case / "stdout.log").write_text(result.stdout)
    if result.returncode:
        raise AssertionError(f"{case.name}: exit {result.returncode}\n{result.stdout}")
    log = (case / "OUT.test/running_scf.log").read_text()
    assert "#SCF IS CONVERGED#" in log, case.name
    assert "!!EXX IS NOT CONVERGED!!" not in log, case.name
    assert "construct_ace" in result.stdout, f"ACE was not exercised: {case.name}"
    if metallic:
        entropy = float(re.findall(r"E_entropy\(-TS\)\s+(\S+)", log)[-1])
        assert abs(entropy) > 1e-6, f"No partial occupations exercised: {case.name}"
    energy = float(re.findall(r"!FINAL_ETOT_IS\s+(\S+)", log)[-1])
    assert math.isfinite(energy), case.name
    states = read_states(case / "OUT.test/eig_occ.txt")
    nk_total = nk * transverse_mesh * transverse_mesh
    assert len(states) == nk_total * nspin * nbands, case.name
    electron_number = sum(state[1] for state in states.values())
    expected_electrons = 3 if metallic else 2
    assert abs(electron_number - expected_electrons) < 1e-6, case.name
    if metallic:
        band_weights = [value[1] for (spin, ik, band), value in states.items()
                        if spin == 1 and band == 2]
        assert max(band_weights) - min(band_weights) > 1e-3, case.name
    return energy, states


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--launcher", default="mpirun")
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path)
    parser.add_argument("--low-symmetry-only", action="store_true")
    args = parser.parse_args()
    binary = args.binary.resolve()
    source = args.source.resolve()
    with tempfile.TemporaryDirectory(prefix="abacus-exx-kpar-") as directory:
        root = args.work_dir.resolve() if args.work_dir else Path(directory)
        root.mkdir(parents=True, exist_ok=True)
        smearings = () if args.low_symmetry_only else ("fixed", "gaussian")
        for smearing in smearings:
            for functional in ("pbe0", "hse"):
                for nspin in (1, 2):
                    # Local FFT/batched and distributed FFT paths. The H2
                    # mesh tests one/two idle pools; the Al 3x3x3 mesh tests
                    # metallic partial occupations and a divisible control.
                    layouts = ((3, 2, 2), (3, 4, 2), (4, 6, 3), (5, 6, 3))
                    if smearing == "gaussian":
                        layouts = ((3, 2, 2), (3, 4, 2), (3, 6, 3))
                    for nk, ranks, kpar in layouts:
                        reference = run_case(binary, args.launcher, source, root,
                                             functional, nspin, nk, ranks, 1, smearing)
                        actual = run_case(binary, args.launcher, source, root,
                                          functional, nspin, nk, ranks, kpar, smearing)
                        delta, occupation_delta, eigenvalue_delta = compare(reference, actual)
                        print(f"PASS {functional} {smearing} nspin={nspin} mesh={nk}x{3 if smearing == 'gaussian' else 1}x{3 if smearing == 'gaussian' else 1} ranks={ranks} "
                              f"kpar={kpar}: |delta E|={delta:.3e} eV, "
                              f"max |delta wg|={occupation_delta:.3e}", flush=True)

        for functional in ("pbe0", "hse"):
            for nspin in (1, 2):
                reference = run_case(binary, args.launcher, source, root,
                                     functional, nspin, 3, 4, 1, "gaussian", True)
                actual = run_case(binary, args.launcher, source, root,
                                  functional, nspin, 3, 4, 2, "gaussian", True)
                delta, occupation_delta, eigenvalue_delta = compare(reference, actual)
                print(f"PASS triclinic {functional} nspin={nspin} nk=27 ranks=4 kpar=2: "
                      f"|delta E|={delta:.3e} eV, max |delta wg|={occupation_delta:.3e}", flush=True)



if __name__ == "__main__":
    main()
