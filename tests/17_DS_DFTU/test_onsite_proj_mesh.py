"""Check onsite projection with unequal orbital meshes and an excessive radius."""

import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def run_case(binary, reference, radius, root):
    case = root / str(radius)
    case.mkdir()
    for name in ("INPUT", "STRU", "KPT"):
        shutil.copy2(reference / name, case / name)

    orbital_dir = reference.parent.parent / "PP_ORB"
    inp = (case / "INPUT").read_text()
    inp = inp.replace("../../PP_ORB", str(orbital_dir))
    inp = inp.replace("scf_nmax    100", "scf_nmax    1")
    inp = inp.replace("onsite_radius   3.0", "onsite_radius   " + str(radius))
    inp = "\n".join(line for line in inp.splitlines()
                    if not line.lstrip().startswith(("dft_plus_u", "orbital_corr", "hubbard_u"))) + "\n"
    (case / "INPUT").write_text(inp)

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    return subprocess.run([str(binary)], cwd=case, env=env,
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          universal_newlines=True)


def main():
    binary = Path(sys.argv[1]).resolve()
    reference = Path(sys.argv[2]).resolve()
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        valid = run_case(binary, reference, 3.0, root)
        if valid.returncode != 0 or "INIT SCF" not in valid.stdout:
            raise AssertionError("Valid onsite mesh failed:\n" + valid.stdout[-3000:])

        invalid = run_case(binary, reference, 6.2, root)
        expected = "onsite_radius must lie inside the orbital radial grid: Fe_gga_6au_100Ry_4s2p2d1f.orb"
        if invalid.returncode == 0 or expected not in invalid.stdout:
            raise AssertionError("Invalid onsite mesh was not rejected:\n" + invalid.stdout[-3000:])


if __name__ == "__main__":
    main()
