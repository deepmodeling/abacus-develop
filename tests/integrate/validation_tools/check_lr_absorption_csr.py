#!/usr/bin/env python3
"""Compare HSE LR spectra from in-memory and exported optical matrices."""
import argparse
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

from check_lr_absorption import compare_spectra, read_spectrum


def run_case(directory, fixture, launch, parameters):
    directory.mkdir()
    for name in ("STRU", "KPT"):
        shutil.copyfile(fixture / name, directory / name)
    lines = [f"{key} {value}\n" for key, value in parameters.items()]
    (directory / "INPUT").write_text("INPUT_PARAMETERS\n" + "".join(lines))
    env = dict(os.environ, OMP_NUM_THREADS="1")
    result = subprocess.run(launch, cwd=directory, env=env, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True, check=False)
    if result.returncode != 0:
        tail = "\n".join(result.stdout.splitlines()[-40:])
        raise ValueError(f"{directory.name} failed with exit code {result.returncode}\n{tail}")
    return directory / "OUT.autotest"


def compare_workflows(args, root):
    fixture = args.case.resolve()
    orbitals = fixture.parent.parent / "PP_ORB"
    parameters = dict(line.split(maxsplit=1) for line in (fixture / "INPUT").read_text().splitlines()
                      if len(line.split(maxsplit=1)) == 2)
    parameters.update(pseudo_dir=orbitals, orbital_dir=orbitals)
    launch = [args.mpiexec, args.mpi_np_flag, "4", *args.mpi_preflag,
              str(args.abacus.resolve()), *args.mpi_postflag]

    producer = dict(parameters, out_hsr="1 16", out_chg="1 16",
                    out_wfc_lcao=1, out_app_flag=0)
    ks_out = run_case(root / "in-memory", fixture, launch, producer)
    consumer = dict(parameters, calculation="nscf", esolver_type="lr",
                    init_chg="file", init_wfc="file", read_file_dir=ks_out)
    csr_out = run_case(root / "standalone", fixture, launch, consumer)
    compare_spectra(ks_out / "absorption.dat", csr_out / "absorption.dat")

    direct_inputs = dict(consumer, abs_gauge="velocity")
    direct_out = run_case(root / "direct", fixture, launch, direct_inputs)
    commutator = read_spectrum(csr_out / "absorption.dat")
    direct = read_spectrum(direct_out / "absorption.dat")
    if any(full[:2] != bare[:2] for full, bare in zip(commutator, direct)):
        raise ValueError("The commutator and direct velocity spectra use different grids.")
    peak = max(row[2] for row in commutator)
    difference = max(abs(full[2] - bare[2]) for full, bare in zip(commutator, direct))
    if peak <= 0 or difference <= 0.05 * peak:
        raise ValueError("The HSE commutator spectrum must differ from direct velocity by over 5% of peak height.")
    print(f"In-memory and CSR spectra agree. Direct velocity differs by {difference / peak:.1%} of peak height.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--abacus", type=Path, required=True)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--mpiexec", required=True)
    parser.add_argument("--mpi-np-flag", required=True)
    parser.add_argument("--mpi-preflag", action="append", default=[])
    parser.add_argument("--mpi-postflag", action="append", default=[])
    args = parser.parse_args()
    try:
        with tempfile.TemporaryDirectory(prefix="abacus-lr-absorption-") as directory:
            compare_workflows(args, Path(directory))
    except (OSError, ValueError) as error:
        parser.exit(1, f"{error}\n")


if __name__ == "__main__":
    main()
