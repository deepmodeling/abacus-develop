#!/usr/bin/env python3
"""Exercise a frozen HSE ensemble on independent PW target k points."""
import argparse
import hashlib
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

REPO = Path(__file__).resolve().parents[3]
COMMON = """INPUT_PARAMETERS
suffix hybrid
basis_type pw
nbands 3
ecutwfc 10
scf_thr 1e-9
scf_nmax 100
pw_diag_thr 1e-10
pw_seed 1
smearing_method fixed
dft_functional HSE
symmetry -1
exxace false
exx_gamma_extrapolation false
out_band 1
cal_force 0
cal_stress 0
"""
MESH = "K_POINTS\n0\nGamma\n2 1 1 0 0 0\n"
PATH = "K_POINTS\n3\nDirect\n0 0 0 1\n0.125 0 0 1\n0.25 0 0 1\n"


def create_case(root, name, extra, points=MESH):
    case = root / name
    case.mkdir()
    shutil.copy(REPO / "tests/01_PW/097_PW_PBE0/STRU", case / "STRU")
    text = COMMON + f"pseudo_dir {REPO / 'tests/PP_ORB'}\n" + extra
    (case / "INPUT").write_text(text)
    (case / "KPT").write_text(points)
    return case


def run(executable, case, launcher, expected_error=None):
    env = dict(os.environ, OMP_NUM_THREADS="1")
    with (case / "run.log").open("w") as log:
        result = subprocess.run(launcher + [str(executable)], cwd=case, env=env,
                                stdout=log, stderr=subprocess.STDOUT, timeout=180)
    text = (case / "run.log").read_text()
    if expected_error:
        assert expected_error in text, text[-4000:]
    else:
        assert result.returncode == 0 and "TOTAL  Time" in text, text[-4000:]


def bands(case, spin=0):
    filename = "band.txt" if spin == 0 else f"bands{spin}.txt"
    lines = (case / "OUT.hybrid" / filename).read_text().splitlines()
    return [[float(value) for value in line.split()[2:]] for line in lines if line.strip()]


def max_difference(reference, target):
    assert len(reference) == len(target)
    return max(abs(a - b) for row, other in zip(reference, target) for a, b in zip(row, other))


def digest(directory):
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in directory.iterdir() if p.is_file()}


def verify(executable, root, mpi_ranks):
    scf = create_case(root, "scf", "calculation scf\nout_chg 1\nout_wfc_pw 2\n")
    run(executable, scf, [])
    source = scf / "OUT.hybrid"
    assert (source / "EXX_SOURCE").exists()
    original = digest(source)
    extra = f"calculation nscf\nread_file_dir {source}\nout_chg 0\nout_wfc_pw 0\n"
    # More target bands than source bands, initialized independently.
    same = create_case(root, "same", extra)
    (same / "INPUT").write_text((same / "INPUT").read_text().replace("nbands 3", "nbands 5"))
    run(executable, same, [])
    delta = max_difference(bands(scf), bands(same))
    assert delta < 2e-4, f"SCF/NSCF eigenvalue mismatch: {delta} eV"
    print(f"same mesh, source 3 bands / target 5 bands: max difference {delta:.3g} eV")

    path = create_case(root, "path", extra, PATH)
    run(executable, path, [])
    values = bands(path)
    assert len(values) == 3 and all(math.isfinite(x) for row in values for x in row)
    gamma_delta = max_difference([bands(same)[0]], [values[0]])
    assert gamma_delta < 2e-4
    assert abs(values[1][0] - values[0][0]) > 0.1
    print(f"independent 3-point target / 2-point source: shared Gamma difference {gamma_delta:.3g} eV")

    if mpi_ranks > 1:
        parallel = create_case(root, "parallel", extra)
        (parallel / "INPUT").write_text((parallel / "INPUT").read_text().replace("nbands 3", "nbands 5"))
        run(executable, parallel, ["mpirun", "-np", str(mpi_ranks)])
        delta = max_difference(bands(same), bands(parallel))
        assert delta < 2e-6, delta
        print(f"serial SCF restart into {mpi_ranks} MPI ranks: max difference {delta:.3g} eV")

    spin_scf = create_case(root, "spin_scf", "calculation scf\nout_chg 1\nout_wfc_pw 2\nnspin 2\n")
    run(executable, spin_scf, [])
    spin_source = spin_scf / "OUT.hybrid"
    spin = create_case(root, "spin_nscf", f"calculation nscf\nread_file_dir {spin_source}\nnspin 2\n")
    run(executable, spin, [])
    for channel in (1, 2):
        delta = max_difference(bands(spin_scf, channel), bands(spin, channel))
        assert delta < 2e-4, delta
        print(f"spin channel {channel}: max difference {delta:.3g} eV")
    assert digest(source) == original, "NSCF modified the SCF restart ensemble"

    broken_source = root / "broken_source"
    shutil.copytree(source, broken_source)
    (broken_source / "EXX_SOURCE").unlink()
    broken = create_case(root, "missing", extra.replace(str(source), str(broken_source)))
    run(executable, broken, [], "EXX source checkpoint")
    (broken_source / "EXX_SOURCE").write_text((source / "EXX_SOURCE").read_text()[:-6])
    truncated = create_case(root, "truncated", extra.replace(str(source), str(broken_source)))
    run(executable, truncated, [], "Invalid EXX source occupation")
    unsupported = create_case(root, "unsupported_ace", extra)
    (unsupported / "INPUT").write_text((unsupported / "INPUT").read_text().replace("exxace false", "exxace true"))
    run(executable, unsupported, [], "Hybrid NSCF currently requires")
    print("missing/truncated checkpoint and unsupported ACE rejected")

    # A failed SCF must invalidate an older companion before replacing orbitals.
    failed = create_case(root, "failed_scf", "calculation scf\nout_chg 1\nout_wfc_pw 2\nexx_hybrid_step 1\n")
    shutil.copytree(source, failed / "OUT.hybrid")
    text = (failed / "INPUT").read_text().replace("scf_nmax 100", "scf_nmax 1")
    (failed / "INPUT").write_text(text)
    run(executable, failed, [])
    assert not (failed / "OUT.hybrid/EXX_SOURCE").exists(), "Failed SCF left a stale checkpoint"
    print("failed SCF invalidates the previous source checkpoint")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("--workdir", type=Path)
    parser.add_argument("--mpi-ranks", type=int, default=2)
    args = parser.parse_args()
    executable = args.executable.resolve()
    if args.workdir:
        args.workdir.mkdir(parents=True, exist_ok=True)
        verify(executable, args.workdir.resolve(), args.mpi_ranks)
    else:
        with tempfile.TemporaryDirectory(prefix="abacus-hybrid-nscf-") as directory:
            verify(executable, Path(directory), args.mpi_ranks)
    print("PASS: screened hybrid NSCF regression")


if __name__ == "__main__":
    main()
