#!/usr/bin/env python3
"""Exercise a frozen HSE ensemble on independent PW target k points."""
import argparse
import hashlib
import math
import os
from pathlib import Path
import shutil
import subprocess
import struct
import tempfile

REPO = Path(__file__).resolve().parents[3]
COMMON = """INPUT_PARAMETERS
suffix hybrid
basis_type pw
device cpu
precision double
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
exx_gamma_extra false
out_band 1
cal_force 0
cal_stress 0
"""
MESH = "K_POINTS\n0\nGamma\n2 1 1 0 0 0\n"
PATH = "K_POINTS\n3\nDirect\n0 0 0 1\n0.125 0 0 1\n0.25 0 0 1\n"


def create_case(root, name, extra, points=MESH):
    case = root / name
    case.mkdir()
    shutil.copy(REPO / "tests/01_PW/scf_pbe0_spin1/STRU", case / "STRU")
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
    assert (source / "eig_occ.txt").exists()
    assert not (source / "EXX_SOURCE").exists()
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
    (broken_source / "eig_occ.txt").unlink()
    broken = create_case(root, "missing", extra.replace(str(source), str(broken_source)))
    run(executable, broken, [], "EXX source eig_occ.txt")
    saved_occupations = (source / "eig_occ.txt").read_text()
    (broken_source / "eig_occ.txt").write_text(saved_occupations[:len(saved_occupations) // 2])
    truncated = create_case(root, "truncated", extra.replace(str(source), str(broken_source)))
    run(executable, truncated, [], "EXX source")
    (broken_source / "eig_occ.txt").write_bytes((source / "eig_occ.txt").read_bytes())
    wavefunction = (source / "wfk1_pw.dat").read_bytes()
    (broken_source / "wfk1_pw.dat").write_bytes(wavefunction[:90])
    broken_wave = create_case(root, "truncated_wave", extra.replace(str(source), str(broken_source)))
    run(executable, broken_wave, [], "EXX source wavefunction")
    invalid_miller = bytearray(wavefunction)
    struct.pack_into("i", invalid_miller, 164, 2**30)
    (broken_source / "wfk1_pw.dat").write_bytes(invalid_miller)
    bad_grid = create_case(root, "bad_source_grid", extra.replace(str(source), str(broken_source)))
    run(executable, bad_grid, [], "EXX source Miller index is incompatible with FFT grid")
    unsupported = create_case(root, "unsupported_ace", extra)
    (unsupported / "INPUT").write_text((unsupported / "INPUT").read_text().replace("exxace false", "exxace true"))
    run(executable, unsupported, [], "Hybrid NSCF currently requires")
    print("missing/truncated source, invalid Miller mapping and unsupported ACE rejected")

    equivalent = create_case(root, "equivalent_settings", extra + "exx_erfc_alpha 2.5e-1\nexx_erfc_omega 1.1e-1\n")
    text = (equivalent / "INPUT").read_text().replace("dft_functional HSE", "dft_functional hse")
    (equivalent / "INPUT").write_text(text)
    run(executable, equivalent, [])
    assert "EXX source configuration differs" not in (equivalent / "OUT.hybrid/warning.log").read_text()
    assert max_difference(bands(same), bands(equivalent)) < 2e-4
    print("equivalent numeric notation and functional case do not trigger mismatch warnings")

    changed_exchange = create_case(root, "changed_exchange", extra + "exx_erfc_alpha 0.3\n")
    run(executable, changed_exchange, [])
    assert "EXX source configuration differs for exx_erfc_alpha" in (changed_exchange / "OUT.hybrid/warning.log").read_text()
    print("changed exchange fraction is accepted with a warning")

    legacy_source = root / "legacy_source"
    shutil.copytree(source, legacy_source)
    (legacy_source / "INPUT.info").unlink()
    legacy = create_case(root, "legacy", extra.replace(str(source), str(legacy_source)))
    run(executable, legacy, [])
    assert "EXX source configuration is incomplete" in (legacy / "OUT.hybrid/warning.log").read_text()
    assert max_difference(bands(same), bands(legacy)) < 2e-4
    print("legacy source without INPUT.info is accepted with a warning")


def verify_correction(executable, root, gpu=False):
    """Check the shared finite-limit scheme and restart convention on each device."""
    device = "gpu" if gpu else "cpu"
    scf = create_case(root, f"limits_{device}_scf",
                      "calculation scf\nout_chg 1\nout_wfc_pw 2\n")
    # Exercise the defaults, including automatic disabling of gamma extrapolation.
    text = (scf / "INPUT").read_text().replace("exx_gamma_extra false\n", "")
    text = text.replace("device cpu", f"device {device}")
    (scf / "INPUT").write_text(text)
    run(executable, scf, [])
    source = scf / "OUT.hybrid"
    original = digest(source)
    extra = f"calculation nscf\nread_file_dir {source}\nexx_singularity_correction limits\n"
    near_gamma = "K_POINTS\n3\nDirect\n0 0 0 1\n0.001 0 0 1\n-0.001 0 0 1\n"
    target = create_case(root, f"limits_{device}_continuity", extra, near_gamma)
    text = (target / "INPUT").read_text().replace("device cpu", f"device {device}")
    (target / "INPUT").write_text(text)
    run(executable, target, [])
    values = bands(target)
    delta = max(max_difference([values[0]], [row]) for row in values[1:])
    assert delta < 5e-3, f"Finite-limit discontinuity near Gamma: {delta} eV"
    assert max_difference([bands(scf)[0]], [values[0]]) < 2e-4
    assert digest(source) == original
    print(f"{device} default/explicit limits, Gamma +/-1e-3: max difference {delta:.3g} eV")

    mismatch = create_case(root, f"{device}_correction_mismatch",
                           extra.replace("correction limits", "correction gygi"))
    text = (mismatch / "INPUT").read_text().replace("device cpu", f"device {device}")
    (mismatch / "INPUT").write_text(text)
    run(executable, mismatch, [])
    assert "EXX source configuration differs for exx_singularity_correction" in (mismatch / "OUT.hybrid/warning.log").read_text()
    print(f"{device} correction mismatch is accepted with a warning")
    bad = create_case(root, f"{device}_bad_correction", extra.replace("correction limits", "correction spencer"))
    run(executable, bad, [], "PW exx_singularity_correction must be limits or gygi")
    old_name = create_case(root, f"{device}_obsolete_correction",
                           extra.replace("correction limits", "correction auxiliary"))
    run(executable, old_name, [], "PW exx_singularity_correction must be limits or gygi")
    gamma = create_case(root, f"{device}_limits_gamma", extra)
    text = (gamma / "INPUT").read_text().replace("exx_gamma_extra false", "exx_gamma_extra true")
    (gamma / "INPUT").write_text(text)
    run(executable, gamma, [], "PW limits requires screened exchange")
    unscreened = create_case(root, f"{device}_limits_fock", extra)
    text = (unscreened / "INPUT").read_text().replace("dft_functional HSE", "dft_functional PBE0")
    (unscreened / "INPUT").write_text(text)
    run(executable, unscreened, [], "PW limits requires screened exchange")

    gygi = create_case(root, f"gygi_{device}_scf",
                            "calculation scf\nout_chg 1\nout_wfc_pw 2\nexx_singularity_correction gygi\n")
    text = (gygi / "INPUT").read_text().replace("device cpu", f"device {device}")
    (gygi / "INPUT").write_text(text)
    run(executable, gygi, [])
    aux_source = gygi / "OUT.hybrid"
    aux_target = create_case(root, f"gygi_{device}_nscf",
                             f"calculation nscf\nread_file_dir {aux_source}\nexx_singularity_correction gygi\n")
    text = (aux_target / "INPUT").read_text().replace("device cpu", f"device {device}")
    (aux_target / "INPUT").write_text(text)
    run(executable, aux_target, [])
    assert max_difference(bands(gygi), bands(aux_target)) < 2e-4
    # Cross the zero-transfer threshold with the same source: the historical
    # gygi scheme must exhibit the jump this test is intended to detect.
    aux_near = create_case(root, f"gygi_{device}_continuity",
                           f"calculation nscf\nread_file_dir {aux_source}\nexx_singularity_correction gygi\n",
                           near_gamma)
    text = (aux_near / "INPUT").read_text().replace("device cpu", f"device {device}")
    (aux_near / "INPUT").write_text(text)
    run(executable, aux_near, [])
    aux_values = bands(aux_near)
    aux_jump = max(max_difference([aux_values[0]], [row]) for row in aux_values[1:])
    assert aux_jump > 5e-3, "Auxiliary control failed to expose the zero-transfer jump"
    shift = max_difference(bands(scf), bands(gygi))
    assert shift > 1e-5, "limits and gygi unexpectedly produced identical SCF bands"
    print(f"{device} gygi control jump near Gamma: {aux_jump:.3g} eV")
    print(f"{device} gygi SCF/NSCF matched; limits/gygi SCF difference {shift:.3g} eV")
    print("unsupported scheme, gamma and unscreened limits rejected")


def verify_gpu(executable, root):
    """Compare devices using the same frozen ensemble, then reverse the restart."""
    cpu_scf = create_case(root, "cpu_scf", "calculation scf\nout_chg 1\nout_wfc_pw 2\n")
    run(executable, cpu_scf, [])
    source = cpu_scf / "OUT.hybrid"
    original = digest(source)
    extra = f"calculation nscf\nread_file_dir {source}\n"
    for name, points in (("same", MESH), ("path", PATH)):
        reference = create_case(root, f"cpu_{name}", extra, points)
        target = create_case(root, f"gpu_{name}", extra, points)
        for case in (reference, target):
            text = (case / "INPUT").read_text().replace("nbands 3", "nbands 5")
            if case == target:
                text = text.replace("device cpu", "device gpu")
            (case / "INPUT").write_text(text)
            run(executable, case, [])
        delta = max_difference(bands(reference), bands(target))
        assert delta < 2e-6, f"CPU/GPU {name} mismatch: {delta} eV"
        print(f"CPU source, CPU/GPU {name} targets: max difference {delta:.3g} eV")
    single = create_case(root, "gpu_single", extra, PATH)
    text = (single / "INPUT").read_text().replace("device cpu", "device gpu")
    text = text.replace("precision double", "precision single")
    text = text.replace("pw_diag_thr 1e-10", "pw_diag_thr 1e-6")
    (single / "INPUT").write_text(text)
    run(executable, single, [])
    reference = [row[:3] for row in bands(root / "cpu_path")]
    delta = max_difference(reference, bands(single))
    assert delta < 2e-4, f"Single-precision GPU restart mismatch: {delta} eV"
    print(f"CPU double source, GPU single 3-band path: max difference {delta:.3g} eV")
    assert digest(source) == original, "GPU NSCF modified CPU source files"

    gpu_scf = create_case(root, "gpu_scf", "calculation scf\nout_chg 1\nout_wfc_pw 2\nnspin 2\n")
    text = (gpu_scf / "INPUT").read_text().replace("device cpu", "device gpu")
    (gpu_scf / "INPUT").write_text(text)
    run(executable, gpu_scf, [])
    source = gpu_scf / "OUT.hybrid"
    assert (source / "eig_occ.txt").exists()
    assert not (source / "EXX_SOURCE").exists()
    original = digest(source)
    targets = []
    for device in ("cpu", "gpu"):
        target = create_case(root, f"{device}_spin", f"calculation nscf\nread_file_dir {source}\nnspin 2\n", PATH)
        text = (target / "INPUT").read_text().replace("device cpu", f"device {device}")
        (target / "INPUT").write_text(text)
        run(executable, target, [])
        targets.append(target)
    for channel in (1, 2):
        delta = max_difference(bands(targets[0], channel), bands(targets[1], channel))
        assert delta < 2e-6, delta
        gamma_delta = max_difference([bands(gpu_scf, channel)[0]], [bands(targets[1], channel)[0]])
        assert gamma_delta < 2e-4, gamma_delta
        print(f"GPU source, spin {channel}, independent CPU/GPU targets: {delta:.3g} eV")
    assert digest(source) == original, "NSCF modified GPU source files"
    unsupported = create_case(root, "gpu_unsupported_ace", extra)
    text = (unsupported / "INPUT").read_text().replace("device cpu", "device gpu")
    text = text.replace("exxace false", "exxace true")
    (unsupported / "INPUT").write_text(text)
    run(executable, unsupported, [], "Hybrid NSCF currently requires")
    print("CUDA NSCF rejects ACE")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("--workdir", type=Path)
    parser.add_argument("--gpu", action="store_true", help="verify CUDA using one rank and double precision")
    parser.add_argument("--mpi-ranks", type=int, default=2)
    args = parser.parse_args()
    executable = args.executable.resolve()
    verify_case = verify_gpu if args.gpu else lambda exe, root: verify(exe, root, args.mpi_ranks)
    if args.workdir:
        args.workdir.mkdir(parents=True, exist_ok=True)
        verify_case(executable, args.workdir.resolve())
        verify_correction(executable, args.workdir.resolve(), args.gpu)
    else:
        with tempfile.TemporaryDirectory(prefix="abacus-hybrid-nscf-") as directory:
            verify_case(executable, Path(directory))
            verify_correction(executable, Path(directory), args.gpu)
    print("PASS: screened hybrid NSCF regression")


if __name__ == "__main__":
    main()
