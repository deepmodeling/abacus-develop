#!/usr/bin/env python3
"""Run the native DFTB3 C2N case and compare its bands with DFTB+."""

import argparse
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


EXPECTED_KPOINTS = 121
EXPECTED_BANDS = 144
MAX_DIFFERENCE_EV = 7.5e-4
RMS_DIFFERENCE_EV = 3.5e-4


def read_dftbplus_bands(path):
    bands = []
    current_kpoint = 0
    last_band = 0
    for line_number, line in enumerate(path.read_text().splitlines(), start=1):
        fields = line.split()
        if not fields:
            continue
        if fields[0] == "KPT":
            current_kpoint += 1
            last_band = 0
            continue
        if len(fields) != 3 or current_kpoint == 0:
            continue
        try:
            band_index = int(fields[0])
            energy_ev = float(fields[1])
            float(fields[2])  # Occupation column; parse it to reject malformed rows.
        except ValueError:
            continue
        if band_index != last_band + 1 or not math.isfinite(energy_ev):
            raise ValueError(f"Invalid DFTB+ band row at {path}:{line_number}")
        last_band = band_index
        bands.append((current_kpoint, band_index, energy_ev))
    return bands


def read_abacus_bands(path):
    bands = []
    for line_number, line in enumerate(path.read_text().splitlines(), start=1):
        fields = line.split()
        if not fields or fields[0].startswith("#"):
            continue
        if len(fields) != 9:
            raise ValueError(f"Unexpected ABACUS band row at {path}:{line_number}")
        try:
            kpoint = int(fields[0])
            band_index = int(fields[6])
            energy_ev = float(fields[7])
        except ValueError as error:
            raise ValueError(f"Invalid ABACUS band row at {path}:{line_number}") from error
        if not math.isfinite(energy_ev):
            raise ValueError(f"Non-finite ABACUS eigenvalue at {path}:{line_number}")
        bands.append((kpoint, band_index, energy_ev))
    return bands


def reverse_stru_species_order(path):
    lines = path.read_text().splitlines()
    species_header = lines.index("ATOMIC_SPECIES")
    positions_header = lines.index("ATOMIC_POSITIONS")

    species_end = species_header + 1
    while species_end < positions_header and lines[species_end].strip():
        species_end += 1
    species_records = lines[species_header + 1 : species_end]
    species_names = [record.split()[0] for record in species_records]
    if len(species_names) < 2 or len(set(species_names)) != len(species_names):
        raise ValueError("Cannot reverse an invalid ATOMIC_SPECIES table")

    cursor = positions_header + 1
    while cursor < len(lines) and not lines[cursor].strip():
        cursor += 1
    if cursor >= len(lines):
        raise ValueError("STRU has no ATOMIC_POSITIONS mode")
    mode = lines[cursor]
    cursor += 1

    atom_blocks = []
    while cursor < len(lines):
        while cursor < len(lines) and not lines[cursor].strip():
            cursor += 1
        if cursor == len(lines):
            break
        block_start = cursor
        species = lines[cursor].split()[0]
        if species not in species_names:
            raise ValueError(f"Unexpected species block in STRU: {species}")
        if cursor + 2 >= len(lines):
            raise ValueError(f"Incomplete STRU block for {species}")
        try:
            atom_count = int(lines[cursor + 2].split()[0])
        except (IndexError, ValueError) as error:
            raise ValueError(f"Invalid atom count for STRU species {species}") from error
        cursor += 3 + atom_count
        if cursor > len(lines):
            raise ValueError(f"Incomplete atom coordinates for STRU species {species}")
        atom_blocks.append((species, lines[block_start:cursor]))

    block_names = [species for species, _ in atom_blocks]
    if block_names != species_names:
        raise ValueError(
            "This regression expects ATOMIC_POSITIONS blocks to follow ATOMIC_SPECIES order"
        )

    reordered_species = list(reversed(species_records))
    updated_lines = lines[: species_header + 1] + reordered_species + lines[species_end:]
    positions_header = updated_lines.index("ATOMIC_POSITIONS")
    reordered_blocks = list(reversed(atom_blocks))
    rebuilt_positions = updated_lines[: positions_header + 1] + [mode, ""]
    for index, (_, block) in enumerate(reordered_blocks):
        if index:
            rebuilt_positions.append("")
        rebuilt_positions.extend(block)

    return "\n".join(rebuilt_positions) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--abacus", required=True, type=Path)
    parser.add_argument("--case-dir", required=True, type=Path)
    parser.add_argument("--mpi-exec")
    parser.add_argument("--mpi-num-procs-flag")
    parser.add_argument("--mpi-num-procs", type=int, default=1)
    parser.add_argument("--mpi-preflag", action="append", default=[])
    parser.add_argument("--mpi-postflag", action="append", default=[])
    parser.add_argument("--reverse-species-order", action="store_true")
    args = parser.parse_args()

    if args.mpi_num_procs < 1:
        raise ValueError("--mpi-num-procs must be positive")
    if args.mpi_num_procs > 1 and not args.mpi_exec:
        raise ValueError("Multiple MPI ranks require --mpi-exec")

    executable = args.abacus.resolve()
    case_dir = args.case_dir.resolve()
    reference = read_dftbplus_bands(case_dir / "reference" / "band.out")
    expected_rows = EXPECTED_KPOINTS * EXPECTED_BANDS
    if len(reference) != expected_rows:
        raise ValueError(f"DFTB+ reference has {len(reference)} rows; expected {expected_rows}")

    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = "1"
    environment["OMP_THREAD_LIMIT"] = "1"
    with tempfile.TemporaryDirectory(prefix="abacus-dftb-c2n-") as temporary_directory:
        working_directory = Path(temporary_directory)
        for filename in ("INPUT", "STRU", "KPT", "dftb_native.in", "dftb_band_path.in"):
            shutil.copy2(case_dir / filename, working_directory / filename)
        shutil.copytree(case_dir / "parameters", working_directory / "parameters")
        if args.reverse_species_order:
            stru_file = working_directory / "STRU"
            stru_file.write_text(reverse_stru_species_order(stru_file))

        command = [str(executable)]
        if args.mpi_exec:
            if not args.mpi_num_procs_flag:
                raise ValueError("--mpi-num-procs-flag is required with --mpi-exec")
            command = (
                [args.mpi_exec]
                + args.mpi_preflag
                + [args.mpi_num_procs_flag, str(args.mpi_num_procs), str(executable)]
                + args.mpi_postflag
            )
        print(
            f"Running C2N regression with {args.mpi_num_procs if args.mpi_exec else 1} MPI rank(s); "
            f"reversed species order={args.reverse_species_order}"
        )
        try:
            result = subprocess.run(
                command,
                cwd=working_directory,
                env=environment,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                timeout=900,
                check=False,
            )
        except subprocess.TimeoutExpired as error:
            output = error.stdout or ""
            raise RuntimeError(f"ABACUS timed out after 900 s. Output tail:\n{output[-5000:]}") from error
        if result.returncode != 0:
            raise RuntimeError(
                f"ABACUS exited with status {result.returncode}. Output tail:\n{result.stdout[-5000:]}"
            )

        output_directory = working_directory / "OUT.dftb_native"
        dftb_log = output_directory / "dftb.log"
        band_file = output_directory / "band.txt"
        if not dftb_log.is_file() or not band_file.is_file():
            raise RuntimeError("ABACUS did not write OUT.dftb_native/dftb.log and band.txt")
        if "# SCC converged: true" not in dftb_log.read_text():
            raise RuntimeError("The C2N DFTB3 SCC calculation did not report convergence")

        calculated = read_abacus_bands(band_file)
        if len(calculated) != expected_rows:
            raise ValueError(f"ABACUS output has {len(calculated)} rows; expected {expected_rows}")
        differences = []
        for expected, actual in zip(reference, calculated):
            if expected[:2] != actual[:2]:
                raise ValueError(f"Band ordering mismatch: DFTB+ {expected[:2]}, ABACUS {actual[:2]}")
            differences.append(actual[2] - expected[2])

    maximum_difference = max(abs(value) for value in differences)
    mean_absolute_difference = sum(abs(value) for value in differences) / len(differences)
    rms_difference = math.sqrt(sum(value * value for value in differences) / len(differences))
    print(
        f"Compared {len(differences)} eigenvalues: max |ΔE|={maximum_difference:.6e} eV, "
        f"mean |ΔE|={mean_absolute_difference:.6e} eV, RMS ΔE={rms_difference:.6e} eV"
    )
    if maximum_difference > MAX_DIFFERENCE_EV or rms_difference > RMS_DIFFERENCE_EV:
        raise AssertionError(
            "ABACUS and DFTB+ band differences exceed the regression limits "
            f"({MAX_DIFFERENCE_EV:.2e} eV max, {RMS_DIFFERENCE_EV:.2e} eV RMS)"
        )


if __name__ == "__main__":
    try:
        main()
    except Exception as error:  # Keep one clear failure message in CTest output.
        print(f"DFTB3 C2N regression failed: {error}", file=sys.stderr)
        sys.exit(1)
