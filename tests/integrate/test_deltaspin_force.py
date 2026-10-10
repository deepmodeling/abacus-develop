#!/usr/bin/env python3
"""Check the LCAO DeltaSpin force component against a frozen-state derivative reference."""

import argparse
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import tempfile


NUMBER = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[Ee][-+]?\d+)?"
# Independent frozen-D/lambda energy derivative; derivation and evidence: PR #8117.
REFERENCE_FORCE = [[-0.003514237992, 0.003514237992, 0.0],
                   [0.003514237992, -0.003514237992, 0.0]]
FORCE_TOLERANCE = 1e-6


def read_result(directory):
    log = (directory / "OUT.test/running_scf.log").read_text()
    stdout = (directory / "stdout.log").read_text()
    headers = next(line.split() for line in stdout.splitlines()
                   if line.lstrip().startswith("ITER "))
    rows = [line.split() for line in stdout.splitlines()
            if re.match(r"\s*[A-Z]+\d+\s", line)]
    last = dict(zip(headers, rows[-1]))
    drho = float(last["DRHO"])
    rms = float(last["RMS"])
    delta_energy = float(re.findall(r"DeltaE_womix\s*=\s*(" + NUMBER + ")", log)[-1])
    energy = float(re.findall(r"!FINAL_ETOT_IS\s+(" + NUMBER + ")", log)[-1])
    force_text = log.rsplit("#DeltaSpin  FORCE#", 1)[1]
    forces = []
    for line in force_text.splitlines():
        fields = line.split()
        if len(fields) == 4 and re.fullmatch(r"Fe[12]", fields[0]):
            forces.append([float(value) for value in fields[1:]])
        if len(forces) == 2:
            break
    finite = all(math.isfinite(value) for value in [energy, drho, rms, delta_energy]
                 + [value for row in forces for value in row])
    if not (finite and len(forces) == 2 and "#SCF IS CONVERGED#" in log
            and drho < 1e-10 and rms < 1.01e-10 and abs(delta_energy) < 1e-9):
        raise RuntimeError("Unconverged or invalid result in " + str(directory))
    return {"energy_eV": energy, "force_eV_per_A": forces, "drho": drho,
            "moment_rms_uB": rms, "delta_energy_eV": delta_energy,
            "final_scf_row": last}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--abacus", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--device", choices=("cpu", "gpu"), default="cpu")
    parser.add_argument("--launcher", default="", help="e.g. 'mpirun -np 1'")
    args = parser.parse_args()
    executable = args.abacus.resolve()
    fixture = Path(__file__).resolve().parents[1] / "17_DS_DFTU/65_LCAO_DS_S4_SO_CF"
    assets = Path(__file__).resolve().parents[1] / "PP_ORB"
    args.work_dir.mkdir(parents=True, exist_ok=True)
    root = Path(tempfile.mkdtemp(prefix="deltaspin-force-", dir=str(args.work_dir.resolve())))
    print("Results: " + str(root), flush=True)
    solver = "cusolver" if args.device == "gpu" else "scalapack_gvx"
    input_text = (fixture / "INPUT").read_text()
    input_text += "device {}\nks_solver {}\npseudo_dir {}\norbital_dir {}\n".format(
        args.device, solver, assets, assets)
    environment = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    results = {}
    name = "base"
    directory = root / name
    directory.mkdir()
    (directory / "INPUT").write_text(input_text)
    structure = (fixture / "STRU").read_text()
    (directory / "STRU").write_text(structure)
    shutil.copyfile(str(fixture / "KPT"), str(directory / "KPT"))
    with (directory / "stdout.log").open("w") as out, (directory / "stderr.log").open("w") as err:
        subprocess.run(shlex.split(args.launcher) + [str(executable)], cwd=str(directory),
                       env=environment, stdout=out, stderr=err, check=True)
    results[name] = read_result(directory)
    print(name + ": converged", flush=True)
    force = results["base"]["force_eV_per_A"]
    error = max(abs(force[i][j] - REFERENCE_FORCE[i][j])
                for i in range(2) for j in range(3))
    summary = {"cases": results, "reference_eV_per_A": REFERENCE_FORCE,
               "max_error_eV_per_A": error, "tolerance_eV_per_A": FORCE_TOLERANCE,
               "passed": error < FORCE_TOLERANCE}
    (root / "result.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    if not summary["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
