"""Run the spinor LCAO interface and check its formatted Pauli matrices."""

import re
import shutil
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np


def read_coefficients(path):
    with path.open("rb") as stream:
        _, _, _, _, bands, orbitals = struct.unpack("=i3d2i", stream.read(36))
        records = np.fromfile(stream, dtype=[("band", "=i4"), ("e", "=f8"),
                                             ("occ", "=f8"), ("c", "=c16", (orbitals,))])
    np.testing.assert_array_equal(records["band"], np.arange(1, bands + 1))
    return records["c"].T


def read_overlap(path):
    with path.open("rb") as stream:
        size, = struct.unpack("=i", stream.read(4))
        overlap = np.zeros((size, size), dtype=complex)
        overlap[np.triu_indices(size)] = np.fromfile(stream, dtype="=c16")
    return overlap + np.triu(overlap, 1).conj().T


def check_spin(out):
    pauli = [np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1, -1])]
    with (out / "oxygen.spn").open() as stream:
        next(stream)
        bands, kpoints = map(int, next(stream).split())
        assert (bands, kpoints) == (13, 2)
        data = np.loadtxt(stream)
    assert data.shape == (3 * kpoints * bands * (bands + 1) // 2, 2)
    actual = data[:, 0] + 1j * data[:, 1]
    expected = []
    overlaps = []
    for ik in range(1, kpoints + 1):
        c = read_coefficients(next(out.rglob(f"wfk{ik}s4_nao.dat")))
        c = np.delete(c, 2, axis=1)  # exclude_bands = 3
        s = read_overlap(next(out.rglob(f"sk{ik}_nao.dat")))[::2, ::2]
        overlaps.append(s)
        matrices = [c.conj().T @ np.kron(s, sigma) @ c for sigma in pauli]
        for m in range(bands):
            for n in range(m + 1):
                expected.extend(mat[n, m] for mat in matrices)
    assert np.max(abs(overlaps[0] - overlaps[1])) > 1e-4
    np.testing.assert_allclose(actual, expected, atol=1e-10, rtol=0, equal_nan=False)
    print(f"Verified {actual.size} spin matrix elements, max error {max(abs(actual - expected)):.3e}")


def check_spn_io_failure(work, abacus, logs, error):
    result = subprocess.run(abacus, cwd=work, capture_output=True, text=True, timeout=60)
    assert result.returncode != 0, "ABACUS accepted an unwritable .spn file"
    diagnostic = f"{error} OUT.ABACUS/oxygen.spn"
    for log in logs:
        assert diagnostic in log.read_text(), f"{log.name} missing '{diagnostic}'\n{result.stdout}{result.stderr}"
    print(f"Verified '{diagnostic}' on {len(logs)} MPI rank(s)")


def main(work, wannier, postw90, abacus):
    tests = Path(__file__).resolve().parents[3] / "tests"
    case = tests / "03_NAO_multik" / "wannier_spn"
    shutil.copy(case / "STRU", work)
    shutil.copy(case / "KPT", work)
    shutil.copy(case / "oxygen.win", work)
    common = (case / "INPUT").read_text() + f"\npseudo_dir {tests}/PP_ORB\norbital_dir {tests}/PP_ORB\n"
    subprocess.run([wannier, "-pp", "oxygen"], cwd=work, check=True)
    scf = common + "calculation scf\nout_chg 1\n"
    (work / "INPUT").write_text(scf + "out_wannier_spn 1\n")
    result = subprocess.run(abacus, cwd=work, capture_output=True, text=True)
    assert result.returncode != 0 and "Invalid INPUT parameter 'out_wannier_spn'" in result.stdout + result.stderr
    (work / "INPUT").write_text(scf)
    subprocess.run(abacus, cwd=work, check=True)
    nscf = common + "calculation nscf\ninit_chg file\ntowannier90 1\nwannier_method 2\nnnkpfile oxygen.nnkp\n"
    out = work / "OUT.ABACUS"
    assert not (out / "oxygen.spn").exists()
    (work / "INPUT").write_text(nscf + "out_wannier_spn 0\nout_wannier_mmn 0\nout_wannier_amn 0\n")
    subprocess.run(abacus, cwd=work, check=True)
    assert (out / "oxygen.eig").is_file() and not (out / "oxygen.spn").exists()
    (work / "INPUT").write_text(nscf + "out_wfc_lcao 2\nout_hsk 2\n")
    subprocess.run(abacus, cwd=work, check=True)
    running_log = (out / "running_nscf.log").read_text()
    match = re.search(r"^ Wrote Wannier90 spin matrices to (.+)$", running_log, re.MULTILINE)
    assert match, "Missing .spn path in running_nscf.log"
    logged_spn = Path(match.group(1))
    assert logged_spn.is_absolute() and logged_spn.samefile(out / "oxygen.spn"), logged_spn
    check_spin(out)
    for ext in ["eig", "mmn", "amn", "spn"]:
        shutil.copy(out / f"oxygen.{ext}", work)
    subprocess.run([wannier, "oxygen"], cwd=work, check=True)
    subprocess.run([postw90, "oxygen"], cwd=work, check=True)
    output = (work / "oxygen.wpout").read_text()
    assert "Reading spin matrices from oxygen.spn" in output and "All done" in output
    moment = [float(x) for x in re.findall(r"[xyz] component:\s+(\S+)", output)]
    assert len(moment) == 3 and np.isfinite(moment).all()
    print(f"postw90 spin moment (Bohr magnetons/cell): {moment}")

    # Native ABACUS can export spin matrices independently of the Python workflow.
    companions = [out / f"oxygen.{ext}" for ext in ("mmn", "amn", "eig")]
    spn = out / "oxygen.spn"
    for path in [spn, *companions]:
        path.unlink()
    (work / "INPUT").write_text(nscf + "out_wannier_spn 1\nout_wannier_mmn 0\n"
                                "out_wannier_amn 0\nout_wannier_eig 0\nout_alllog 1\n")
    subprocess.run(abacus, cwd=work, check=True)
    assert spn.stat().st_size > 0 and not any(path.exists() for path in companions)
    logs = sorted(out.glob("running_nscf_*.log"))
    assert logs, "out_alllog did not produce per-rank NSCF logs"
    spn.unlink()
    spn.mkdir()
    check_spn_io_failure(work, abacus, logs, "Cannot open")
    spn.rmdir()
    if Path("/dev/full").is_char_device():
        spn.symlink_to("/dev/full")
        check_spn_io_failure(work, abacus, logs, "Error writing")


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as scratch:
        main(Path(scratch), sys.argv[1], sys.argv[2], sys.argv[3:])
