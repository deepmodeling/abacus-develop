"""Check first-step PW DFT+U output with and without a loaded matrix."""

import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def replace_setting(text, old, new):
    if old not in text:
        raise AssertionError(f"Missing fixture setting: {old}")
    return text.replace(old, new, 1)


def first_matrix_element(text):
    lines = text.splitlines()
    for index, line in enumerate(lines[:-1]):
        if line.strip().startswith("spin 1 nelec"):
            return float(lines[index + 1].split()[0])
    raise AssertionError("First spin matrix is missing")


def run_case(executable, fixture, work, init_occ_mat):
    for name in ("KPT", "STRU"):
        shutil.copy2(fixture / name, work / name)
    if init_occ_mat:
        shutil.copy2(fixture / "occ_mat.txt", work / "occ_mat.txt")

    input_text = (fixture / "INPUT").read_text()
    input_text = replace_setting(input_text, "calculation     relax", "calculation     scf")
    input_text = replace_setting(input_text, "scf_nmax        70", "scf_nmax        1")
    input_text = replace_setting(input_text, "out_freq_ion    2", "out_freq_ion    1")
    if init_occ_mat != 2:
        input_text = replace_setting(input_text, "init_occ_mat    2", f"init_occ_mat    {init_occ_mat}")
    data_dir = fixture.parent.parent / "PP_ORB"
    input_text = replace_setting(input_text, "pseudo_dir      ../../PP_ORB/", f"pseudo_dir      {data_dir}/")
    input_text = replace_setting(input_text, "orbital_dir     ../../PP_ORB/", f"orbital_dir     {data_dir}/")
    (work / "INPUT").write_text(input_text)

    output_dir = work / "OUT.autotest"
    output_dir.mkdir()
    numbered = output_dir / "occ_matg1.txt"
    numbered.write_text("# Electronic step 99\nSTALE SECTION\n")

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    result = subprocess.run(
        [str(executable)], cwd=work, env=env, capture_output=True, text=True, timeout=180
    )
    if result.returncode != 0:
        raise AssertionError(
            f"ABACUS failed with init_occ_mat={init_occ_mat}: {result.returncode}\n"
            + (result.stdout + result.stderr)[-4000:]
        )

    assert numbered.exists(), "Numbered occupation output is missing"
    snapshot = numbered.read_text()
    assert "STALE SECTION" not in snapshot
    assert snapshot.count("# Electronic step 1") == 1
    latest = output_dir / "occ_mat.txt"
    if init_occ_mat:
        assert latest.exists(), "Loaded matrix was not written at the first PW step"
        latest_text = latest.read_text()
        assert "# Electronic step 1" in latest_text
        assert "Fe Atom 1 L 2" in latest_text
        assert "Fe Atom 1 L 2" in snapshot
        assert "N/A" not in latest_text + snapshot
        source_value = first_matrix_element((fixture / "occ_mat.txt").read_text())
        assert abs(first_matrix_element(latest_text) - source_value) < 1e-6
        assert abs(first_matrix_element(snapshot) - source_value) < 1e-6
    else:
        assert not latest.exists(), "Fresh PW run wrote an unavailable first-step matrix"
        assert "N/A" in snapshot
        assert "Fe Atom 1 L 2" not in snapshot


def main():
    executable = Path(sys.argv[1]).resolve()
    fixture = Path(sys.argv[2]).resolve()
    with tempfile.TemporaryDirectory(prefix="abacus_dftu_pw_first_step_") as scratch:
        root = Path(scratch)
        for init_occ_mat in (2, 0):
            work = root / f"init_occ_mat_{init_occ_mat}"
            work.mkdir()
            run_case(executable, fixture, work, init_occ_mat)
    print("PW DFT+U first-step output checks passed for init_occ_mat=0 and 2")


if __name__ == "__main__":
    main()
