"""Check that a SNAP error on a nonzero MPI rank remains visible."""

import argparse
import os
from pathlib import Path
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--abacus", type=Path, required=True)
    parser.add_argument("--mpiexec", type=Path, required=True)
    parser.add_argument("--numproc-flag", default="-n")
    parser.add_argument("--models", type=Path, required=True)
    args = parser.parse_args()

    coefficient = (args.models / "snap_102.coeff").resolve()
    with tempfile.TemporaryDirectory(prefix="abacus-snap-rank-") as scratch:
        case = Path(scratch)
        (case / "INPUT").write_text(
            "INPUT_PARAMETERS\n"
            "calculation md\n"
            "suffix rank_error\n"
            "esolver_type snap\n"
            f"pot_file {coefficient}\n"
            "snap_param_file snap.param\n"
            "cal_force 1\n"
            "md_nstep 1\n"
            "md_type nve\n"
            "init_vel 0\n",
            encoding="utf-8",
        )
        (case / "snap.param").write_text(
            "rcutfac 4.9\n"
            "twojmax 8\n"
            "rfac0 0.99363\n"
            "rmin0 0.2\n"
            "bzeroflag 0\n",
            encoding="utf-8",
        )
        # Rank 0 owns the distant atom; rank 1 owns the close pair.
        (case / "STRU").write_text(
            "ATOMIC_SPECIES\nSi 28.0855 auto auto\n\n"
            "LATTICE_CONSTANT\n1.0\n\n"
            "LATTICE_VECTORS\n20.0 0.0 0.0\n0.0 20.0 0.0\n0.0 0.0 20.0\n\n"
            "ATOMIC_POSITIONS\nCartesian\n\nSi\n0\n3\n"
            "2.0 0.0 0.0 m 1 1 1\n"
            "14.0 0.0 0.0 m 1 1 1\n"
            "14.2 0.0 0.0 m 1 1 1\n",
            encoding="utf-8",
        )
        environment = os.environ.copy()
        environment["OMP_NUM_THREADS"] = "1"
        environment["OMPI_ALLOW_RUN_AS_ROOT"] = "1"
        environment["OMPI_ALLOW_RUN_AS_ROOT_CONFIRM"] = "1"
        command = [str(args.mpiexec), args.numproc_flag, "2", str(args.abacus)]
        result = subprocess.run(
            command, cwd=str(case), env=environment, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=90,
        )
        expected = "ESolver_SNAP rank 1: SNAP evaluation failed: " \
                   "SNAP neighbor distance must be greater than rmin0."
        coefficient_messages = result.stdout.count(" Read SNAP coefficient file: ")
        parameter_messages = result.stdout.count(" Read SNAP parameter file: ")
        if (result.returncode == 0 or expected not in result.stderr
                or "SNAP evaluation failed" in result.stdout
                or coefficient_messages != 1 or parameter_messages != 1):
            raise AssertionError(
                "Expected one rank-0 file confirmation, a rank-1 SNAP error, "
                f"and nonzero exit; got {result.returncode}.\n"
                f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
            )


if __name__ == "__main__":
    main()
