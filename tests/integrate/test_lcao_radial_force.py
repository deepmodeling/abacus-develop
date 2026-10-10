#!/usr/bin/env python3
"""Compare an unconstrained LCAO Fe2 force with a converged energy derivative."""

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
STEP_ANGSTROM = 0.0005
FORCE_TOLERANCE = 1e-4
# Match ModuleBase::BOHR_TO_A and the fixture's LATTICE_CONSTANT.
COORDINATE_TO_ANGSTROM = 0.5291770 * 1.8897261254578281


def read_result(directory):
    log = (directory / "OUT.test/running_scf.log").read_text()
    stdout = (directory / "stdout.log").read_text()
    headers = next(line.split() for line in stdout.splitlines()
                   if line.lstrip().startswith("ITER "))
    rows = [line.split() for line in stdout.splitlines()
            if re.match(r"\s*[A-Z]+\d+\s", line)]
    last = dict(zip(headers, rows[-1]))
    drho = float(last["DRHO"])
    delta_energy = float(re.findall(r"DeltaE_womix\s*=\s*(" + NUMBER + ")", log)[-1])
    energy = float(re.findall(r"!FINAL_ETOT_IS\s+(" + NUMBER + ")", log)[-1])
    force_text = log.rsplit("#TOTAL-FORCE (eV/Angstrom)#", 1)[1]
    forces = []
    for line in force_text.splitlines():
        fields = line.split()
        if len(fields) == 4 and re.fullmatch(r"Fe[12]", fields[0]):
            forces.append([float(value) for value in fields[1:]])
        if len(forces) == 2:
            break
    net = re.findall(r"Net force vector.*?\[([^\]]+)\]", log)[-1]
    net_force = [float(value) for value in re.findall(NUMBER, net)]
    finite = all(math.isfinite(value) for value in [energy, drho, delta_energy]
                 + net_force + [value for row in forces for value in row])
    if not (finite and len(net_force) == 3 and len(forces) == 2 and "#SCF IS CONVERGED#" in log
            and drho < 1e-10 and abs(delta_energy) < 1e-9):
        raise RuntimeError("Unconverged or invalid result in " + str(directory))
    return {"energy_eV": energy, "force_eV_per_A": forces, "drho": drho,
            "net_force_eV_per_A": net_force, "delta_energy_eV": delta_energy,
            "final_scf_row": last}


def displaced_structure(text, displacement):
    lines = text.splitlines()
    atom = 0
    for index, line in enumerate(lines):
        fields = line.split()
        if "magmom" not in fields:
            continue
        # Displace Fe2 only; keep the atom at the cell origin fixed.
        direction = 0.0 if atom == 0 else 1.0
        fields[0] = "{:.16f}".format(float(fields[0])
                                    + direction * displacement / COORDINATE_TO_ANGSTROM)
        lines[index] = " ".join(fields)
        atom += 1
    if atom != 2:
        raise RuntimeError("Expected two Fe atoms in fixture")
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cutoff", type=int, choices=(100, 200), required=True)
    parser.add_argument("--abacus", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--device", choices=("cpu", "gpu"), default="cpu")
    parser.add_argument("--launcher", default="", help="e.g. 'mpirun -np 1'")
    args = parser.parse_args()
    executable = args.abacus.resolve()
    fixture = Path(__file__).resolve().parent / "fixtures/lcao_radial_force"
    assets = Path(__file__).resolve().parents[1] / "PP_ORB"
    args.work_dir.mkdir(parents=True, exist_ok=True)
    root = Path(tempfile.mkdtemp(prefix="lcao-radial-force-", dir=str(args.work_dir.resolve())))
    print("Results: " + str(root), flush=True)
    solver = "cusolver" if args.device == "gpu" else "scalapack_gvx"
    input_text = (fixture / "INPUT").read_text()
    input_text = input_text.replace("ecutwfc 100", "ecutwfc " + str(args.cutoff))
    input_text += "device {}\nks_solver {}\npseudo_dir {}\norbital_dir {}\n".format(
        args.device, solver, assets, assets)
    environment = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    results = {}
    for name, displacement in (("base", 0.0), ("minus", -STEP_ANGSTROM),
                               ("plus", STEP_ANGSTROM)):
        directory = root / name
        directory.mkdir()
        (directory / "INPUT").write_text(input_text)
        structure = displaced_structure((fixture / "STRU").read_text(), displacement)
        (directory / "STRU").write_text(structure)
        shutil.copyfile(str(fixture / "KPT"), str(directory / "KPT"))
        with (directory / "stdout.log").open("w") as out, (directory / "stderr.log").open("w") as err:
            subprocess.run(shlex.split(args.launcher) + [str(executable)], cwd=str(directory),
                           env=environment, stdout=out, stderr=err, check=True)
        results[name] = read_result(directory)
        print(name + ": converged", flush=True)
    force = results["base"]["force_eV_per_A"]
    # Undo the uniform net-force subtraction in ABACUS printed forces.
    analytical = force[1][0] + results["base"]["net_force_eV_per_A"][0] / 2.0
    numerical = -(results["plus"]["energy_eV"] - results["minus"]["energy_eV"]) / (2 * STEP_ANGSTROM)
    error = numerical - analytical
    summary = {"cases": results, "step_A": STEP_ANGSTROM,
               "analytical_eV_per_A": analytical, "finite_difference_eV_per_A": numerical,
               "error_eV_per_A": error, "tolerance_eV_per_A": FORCE_TOLERANCE,
               "passed": abs(error) < FORCE_TOLERANCE}
    summary["cutoff_Ry"] = args.cutoff
    summary["sc_mag_switch"] = 0
    (root / "result.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0 if summary["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
