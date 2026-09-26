#!/usr/bin/env python3
"""Small rVV10 SCF and supported-configuration regressions; Python stdlib only."""

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile


# !FINAL_ETOT_IS uses this ABACUS conversion, not a newer CODATA value.
EV_PER_RY = 13.605698
ENERGY_ATOL_RY = 2.0e-6
EFFECT_ATOL_RY = 2.0e-7
SUCCESS = "INPUT parameters have been successfully checked!"
FUNCTIONALS = {"rvv10": "GGA_X_RPW86+GGA_C_PBE", "base": "GGA_X_RPW86+GGA_C_PBE", "pbe": "PBE"}
NONLOCAL = {"rvv10": "rvv10", "base": "none", "pbe": "none"}
# Converged ABACUS references, independently cross-checked with QE. See README.
REFERENCES_RY = {
    "he-rvv10": -5.6120307435324417,
    "he-base": -5.6290681676412015,
    "he-pbe": -5.5813045054136214,
    "si-rvv10": -15.705104918677593,
    "si-base": -15.734085857260649,
}
EFFECTS_RY = {"he": 0.0170374241087598, "si": 0.0289809385830559}
PSEUDOS = {
    "he": ("He_ONCV_PBE-1.0.upf", "f0af35831c74f9b6a53dcca409ffc3ca9ad8e1e7144b7e627321ff12817710e9"),
    "si": ("Si_dojo_nsoc.upf", "39822757f53f36e3bf3bfb779356152a8d3f21199c7db9dd5a931e5d18c45282"),
    "na": ("Na.pbe-spn-rrkjus_psl.0.2.UPF", "1ccbb6048a1c0bcd537e14cf92a1f54d09c98dcec8ba0bbc4ecd9de480572912"),
}

UNSUPPORTED_VARIANTS = (
    ("noncollinear", {"nspin": "4"}),
    ("force", {"cal_force": "1"}),
    ("stress", {"cal_stress": "1"}),
    ("relax", {"calculation": "relax"}),
    ("dispersion", {"vdw_method": "d3_bj"}),
    ("basis", {"basis_type": "lcao", "ks_solver": "lapack"}),
    ("gpu", {"device": "gpu"}),
)


def parse_energy(text):
    if "#SCF IS CONVERGED#" not in text:
        raise ValueError("SCF was not converged")
    rows = re.findall(r"(?m)^\s*!FINAL_ETOT_IS\s+([^\n]+)", text)
    if len(rows) != 1:
        raise ValueError("Expected exactly one final energy")
    words = rows[0].split()
    if len(words) != 2 or words[1] != "eV":
        raise ValueError("Final energy must explicitly carry eV units")
    energy = float(words[0]) / EV_PER_RY
    if not math.isfinite(energy):
        raise ValueError("Final energy is not finite")
    return energy


def check_input_result(text, returncode, rejection=None):
    if any(marker in text for marker in ("ERROR: AddressSanitizer", "ERROR: LeakSanitizer",
                                         "AddressSanitizer:DEADLYSIGNAL", "runtime error:")):
        raise ValueError("Sanitizer failure during INPUT or runtime validation")
    if rejection is None:
        if returncode != 0 or SUCCESS not in text:
            raise ValueError("Valid INPUT did not finish its check successfully")
    elif returncode not in (0, 1) or rejection not in text or SUCCESS in text:
        # WARNING_QUIT exit codes alone cannot distinguish intentional rejection
        # from a crash, an unrelated failure, or a silently accepted INPUT.
        raise ValueError("Expected controlled rejection: " + rejection)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def check_rejected_run(text, returncode, reason):
    check_input_result(text, returncode, reason)
    if any(marker in text for marker in
           ("Start SCF iteration", "#ELEC ITER#", "#SCF IS CONVERGED#", "!FINAL_ETOT_IS")):
        raise ValueError("Rejected run nevertheless started SCF")


def prepare_case(root, name, system, pseudo_dir, changes):
    case = root / name
    case.mkdir()
    fixture = Path(__file__).resolve().parent / "fixtures" / system
    values = {}
    for line in (fixture / "INPUT").read_text().splitlines()[1:]:
        if line.strip():
            key, value = line.split(None, 1)
            values[key] = value
    values.update(changes)
    values["pseudo_dir"] = str(pseudo_dir)
    (case / "INPUT").write_text(
        "INPUT_PARAMETERS\n" + "".join(key + " " + value + "\n" for key, value in values.items())
    )
    for filename in ("STRU", "KPT"):
        shutil.copyfile(fixture / filename, case / filename)
    return case


def execute(executable, case, flags, logfile, env, timeout):
    command = executable if isinstance(executable, list) else [str(executable)]
    with (case / logfile).open("w") as output:
        result = subprocess.run(
            command + flags, cwd=str(case), env=env,
            stdout=output, stderr=subprocess.STDOUT, timeout=timeout, check=False,
        )
    return result.returncode, (case / logfile).read_text(errors="replace")


def run(args, root, summary):
    executable = Path(args.executable).resolve(strict=True)
    if args.launch is not None and str(executable) not in args.launch:
        raise ValueError("The launcher command must include the resolved --executable path")
    command = args.launch or [str(executable)]
    pseudo_dir = Path(args.pseudo_dir).resolve(strict=True)
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    summary["executable"] = str(executable)
    summary["executable_sha256"] = sha256(executable)
    summary["command"] = command
    code, version = execute(command, root, ["--version"], "version.log", env, 30)
    if code != 0 or not version.strip():
        raise ValueError("Executable identity check failed")
    summary["version"] = version.strip()
    systems = {"all": ("he", "si", "na"), "scf": ("he", "si"), "mpi-scf": ("he",),
               "uspp": ("he", "na")}.get(args.mode, ("he",))
    for system in systems:
        filename, expected = PSEUDOS[system]
        if sha256(pseudo_dir / filename) != expected:
            raise ValueError("Pseudopotential differs from reference: " + filename)

    if args.mode in ("no-libxc", "mpi-spin-reject"):
        name = "reject-" + args.mode
        reason = ("xc_nonlocal=rvv10 requires a LIBXC-enabled build" if args.mode == "no-libxc"
                  else "xc_nonlocal=rvv10 supports CPU double PW nspin=1/2 SCF")
        values = {"nspin": "4"} if args.mode == "mpi-spin-reject" else {}
        case = prepare_case(root, name, "he", pseudo_dir, values)
        code, text = execute(command, case, ["--check-input"], "check.log", env, 30)
        check_rejected_run(text, code, reason)
        summary["checks"].append({"case": name, "passed": True})
        print("PASS " + name, flush=True)
        return

    changes = ({"dft_functional": "PBE", "xc_nonlocal": "none"}
               if args.mode == "mpi-pbe" else {})
    case = prepare_case(root, "valid-input", "he", pseudo_dir, changes)
    code, text = execute(command, case, ["--check-input"], "check.log", env, 30)
    check_input_result(text, code)
    summary["checks"].append({"case": "valid-input", "passed": True})

    if args.mode in ("all", "reject"):
        for name, changes in UNSUPPORTED_VARIANTS:
            case = prepare_case(root, "reject-" + name, "he", pseudo_dir, changes)
            code, text = execute(command, case, ["--check-input"], "check.log", env, 30)
            check_input_result(text, code, "xc_nonlocal=rvv10 supports CPU double PW nspin=1/2 SCF")
            summary["checks"].append({"case": "reject-" + name, "passed": True})
            print("PASS unsupported INPUT rejected: " + name, flush=True)

        # The nonlocal component is not coupled to the RPW86+PBE parser path.
        # This is an INPUT-level check only; a production parameter set still
        # needs its own numerical regression before being advertised.
        case = prepare_case(root, "accept-pbe-rvv10", "he", pseudo_dir,
                            {"dft_functional": "PBE", "xc_nonlocal": "rvv10",
                             "rvv10_b": "6.3", "rvv10_c": "0.0093"})
        code, text = execute(command, case, ["--check-input"], "check.log", env, 30)
        check_input_result(text, code)
        summary["checks"].append({"case": "accept-pbe-rvv10", "passed": True})
        print("PASS independent PBE+rVV10 INPUT", flush=True)

        # A Libxc functional that already contains a nonlocal term must not be
        # combined with the independent rVV10 component.
        for functional in ("GGA_XC_BEEF_VDW", "GGA_XC_VV10", "MGGA_X_SCAN+MGGA_C_SCAN_VV10"):
            name = "reject-" + functional
            case = prepare_case(root, name, "he", pseudo_dir,
                                {"dft_functional": functional, "xc_nonlocal": "rvv10"})
            code, text = execute(command, case, [], "run.log", env, args.timeout)
            for log in (case / "OUT.rvv10_regression").glob("*.log"):
                text += "\n" + log.read_text(errors="replace")
            check_rejected_run(text, code, "xc_nonlocal=rvv10 cannot be combined with a functional that already contains vdW or VV10")
            summary["checks"].append({"case": name, "passed": True})
            print("PASS unsupported nonlocal functional rejected: " + functional, flush=True)

    if args.mode in ("all", "uspp"):
        # USPP is identified while reading the actual UPF, after --check-input.
        case = prepare_case(root, "reject-uspp", "he", pseudo_dir,
                            {"nbands": "6", "scf_nmax": "1"})
        shutil.copyfile(Path(__file__).resolve().parent / "fixtures" / "na" / "STRU", case / "STRU")
        code, text = execute(command, case, [], "run.log", env, args.timeout)
        for log in (case / "OUT.rvv10_regression").glob("*.log"):
            text += "\n" + log.read_text(errors="replace")
        check_rejected_run(text, code, "rVV10 currently supports norm-conserving pseudopotentials only")
        summary["checks"].append({"case": "reject-uspp", "passed": True})
        print("PASS actual USPP rejected before SCF", flush=True)

    if args.mode in ("spin-scf", "mpi-spin-scf"):
        # A zero-magnetization nspin=2 run should reproduce the nspin=1
        # reference for the same unpolarized He system.  Running both cases
        # in one invocation avoids freezing a second absolute reference and
        # checks the production SCF adapter rather than only the component
        # unit test.
        energies = {}
        for spin in (1, 2):
            name = "he-rvv10-spin{}".format(spin)
            case = prepare_case(root, name, "he", pseudo_dir,
                                {"nspin": str(spin), "nupdown": "0.0"})
            code, unused = execute(command, case, [], "run.log", env, args.timeout)
            if code != 0:
                raise ValueError(name + ": SCF exited with code " + str(code))
            log = case / "OUT.rvv10_regression" / "running_scf.log"
            actual = parse_energy(log.read_text(errors="replace"))
            energies[spin] = actual
            summary["checks"].append({"case": name, "energy_Ry": actual})
            print("PASS {}: E = {:.13f} Ry".format(name, actual), flush=True)
        delta = energies[2] - energies[1]
        summary["checks"].append({"case": "nspin2-vs-nspin1", "delta_Ry": delta,
                                  "tolerance_Ry": ENERGY_ATOL_RY})
        if abs(delta) > ENERGY_ATOL_RY:
            raise ValueError("nspin=2 and nspin=1 energies differ by {:.9g} Ry (limit {:.3g} Ry)".format(
                delta, ENERGY_ATOL_RY))
        print("PASS nspin2-vs-nspin1: delta = {:.3g} Ry".format(delta), flush=True)
        return

    if args.mode in ("all", "scf", "mpi-scf", "mpi-pbe"):
        energies = {}
        if args.mode == "mpi-pbe":
            references = {"he-pbe": REFERENCES_RY["he-pbe"]}
        elif args.mode == "mpi-scf":
            references = {"he-rvv10": REFERENCES_RY["he-rvv10"]}
        else:
            references = REFERENCES_RY
        for name, expected in references.items():
            system, method = name.split("-")
            case = prepare_case(root, name, system, pseudo_dir,
                                {"dft_functional": FUNCTIONALS[method],
                                 "xc_nonlocal": NONLOCAL[method]})
            code, unused = execute(command, case, [], "run.log", env, args.timeout)
            if code != 0:
                raise ValueError(name + ": SCF exited with code " + str(code))
            log = case / "OUT.rvv10_regression" / "running_scf.log"
            actual = parse_energy(log.read_text(errors="replace"))
            delta = actual - expected
            summary["checks"].append({"case": name, "energy_Ry": actual,
                                      "reference_Ry": expected, "delta_Ry": delta,
                                      "tolerance_Ry": ENERGY_ATOL_RY})
            if abs(delta) > ENERGY_ATOL_RY:
                raise ValueError("{}: energy differs by {:.9g} Ry (limit {:.3g} Ry)".format(
                    name, delta, ENERGY_ATOL_RY))
            energies[name] = actual
            print("PASS {}: E = {:.13f} Ry, delta = {:.3g} Ry".format(name, actual, delta), flush=True)
        for system, expected in ({} if args.mode in ("mpi-pbe", "mpi-scf") else EFFECTS_RY).items():
            actual = energies[system + "-rvv10"] - energies[system + "-base"]
            summary["checks"].append({"case": system + "-nonlocal-scf-effect", "effect_Ry": actual,
                                      "reference_Ry": expected, "tolerance_Ry": EFFECT_ATOL_RY})
            if abs(actual - expected) > EFFECT_ATOL_RY:
                raise ValueError(system + ": rVV10-minus-base SCF energy change is incorrect")
            print("PASS {}: rVV10-minus-base = {:.13f} Ry".format(system, actual), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", required=True)
    parser.add_argument("--pseudo-dir", required=True)
    parser.add_argument("--artifacts-dir", help="Parent for a fresh, always-preserved result directory")
    parser.add_argument(
        "--mode",
        choices=("all", "scf", "spin-scf", "reject", "uspp", "no-libxc", "mpi-scf",
                 "mpi-spin-scf", "mpi-spin-reject", "mpi-pbe"),
        default="all",
    )
    parser.add_argument("--timeout", type=int, default=300, help="Maximum seconds per SCF (default: 300)")
    parser.add_argument("--launch", nargs=argparse.REMAINDER,
                        help="Complete launcher command with resolved --executable path; must be the last option")
    args = parser.parse_args()
    if args.timeout <= 0:
        parser.error("--timeout must be positive")
    if args.artifacts_dir:
        Path(args.artifacts_dir).mkdir(parents=True, exist_ok=True)
    root = Path(tempfile.mkdtemp(prefix="rvv10-regression-", dir=args.artifacts_dir)).resolve()
    print("Artifacts: " + str(root), flush=True)
    summary = {"mode": args.mode, "passed": False, "checks": [], "energy_unit": "Ry"}
    try:
        run(args, root, summary)
        summary["passed"] = True
    except (OSError, ValueError, subprocess.TimeoutExpired) as error:
        summary["error"] = str(error)
        print("FAIL: " + str(error), file=sys.stderr, flush=True)
    finally:
        (root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    if summary["passed"]:
        print("PASS rVV10 regression (" + args.mode + ")", flush=True)
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())
