#!/usr/bin/env python3
"""Utilities for cube files used by ABACUS integration tests.

This file intentionally keeps the cube parser and all cube-specific regression
checks together. Use one of the following subcommands:

``integrate CUBE``
    Integrate one scalar cube file over real space. This is suitable for charge
    density, partial charge, wavefunction modulus, potential, and other scalar
    cube outputs. The command prints one number: ``sum(values) * voxel_volume``.

``compare-wfc CAL_RE CAL_IM REF_RE REF_IM [CAL_RE CAL_IM REF_RE REF_IM ...]``
    Compare one complex wavefunction state after removing its arbitrary global
    U(1) phase. Pass one group of four files for a scalar wavefunction or for an
    independently aligned spin channel. Pass two groups for an nspin=4 spinor;
    all groups then share the same fitted phase. The command prints the maximum
    pointwise complex error after alignment. It does not perform unitary
    alignment within a degenerate subspace.

``fingerprint-wfc RE IM [RE IM ...]``
    Generate a compact, deterministic fingerprint of one complex wavefunction
    state. One real/imaginary pair represents a scalar component; two pairs in
    upper/lower order represent an nspin=4 spinor. The fingerprint contains
    grid metadata, the RMS amplitude, and phase-invariant products of eight
    fixed complex projections. It is insensitive to one common U(1) phase but
    remains sensitive to spatial changes and relative spinor phases. This is a
    regression summary, not a proof of pointwise equality and not an invariant
    of unitary rotations within a degenerate band subspace.

    For flattened values ``psi[i]``, let ``E = sum(abs(psi[i])**2)`` and
    ``q[j] = sum(w[j,i] * psi[i]) / sqrt(E)``. Eight deterministic SplitMix64
    sketches map each ``(j, i)`` to one of ``{1, -1, 1j, -1j}``. The reported
    powers ``abs(q[j])**2`` and cyclic products
    ``q[j] * conj(q[(j + 1) % 8])`` are unchanged by one common phase applied
    to every component. The constants and flattening order are part of the
    reference format and must not be changed without regenerating references.

``check-spinor DIRECTORY [--tolerance VALUE]``
    Check pointwise Pauli identities among separate-k PW nspin=4 outputs in one
    ABACUS output directory. This requires wavefunction modulus, spinor real and
    imaginary parts, and all four partial-charge components. It also checks that
    their grids agree with each other and, when available, with the dense grid
    reported in the running log. Successful checks produce no output.

Examples:

    python3 cube_tool.py integrate OUT.autotest/charge.cube
    python3 cube_tool.py compare-wfc cal_re.cube cal_im.cube ref_re.cube ref_im.cube
    python3 cube_tool.py fingerprint-wfc state_re.cube state_im.cube
    python3 cube_tool.py check-spinor OUT.autotest --tolerance 1e-8

All commands use only the Python standard library. Invalid or missing inputs are
reported on standard error and return a nonzero exit status.
"""

import argparse
import math
import re
import sys
from pathlib import Path


FINGERPRINT_SKETCHES = 8
UINT64_MASK = (1 << 64) - 1
SPLITMIX_GOLDEN = 0x9E3779B97F4A7C15
SPLITMIX_MIX1 = 0xBF58476D1CE4E5B9
SPLITMIX_MIX2 = 0x94D049BB133111EB
FINGERPRINT_WEIGHTS = (1.0 + 0.0j, -1.0 + 0.0j, 0.0 + 1.0j, 0.0 - 1.0j)


def read_cube(path):
    """Return the grid shape, voxel volume, and scalar data from a cube file."""
    input_path = Path(path)
    with input_path.open("r", encoding="utf-8") as stream:
        stream.readline()
        stream.readline()

        atom_line = stream.readline().split()
        if len(atom_line) < 4:
            raise ValueError(f"invalid atom header in {input_path}")
        natom = int(atom_line[0])

        shape = []
        vectors = []
        for _ in range(3):
            grid_line = stream.readline().split()
            if len(grid_line) < 4:
                raise ValueError(f"invalid grid header in {input_path}")
            shape.append(int(grid_line[0]))
            vectors.append([float(value) for value in grid_line[1:4]])

        v1, v2, v3 = vectors
        val0 = v2[1] * v3[2] - v2[2] * v3[1]
        val1 = v2[0] * v3[2] - v2[2] * v3[0]
        val2 = v2[0] * v3[1] - v2[1] * v3[0]
        volume = abs(v1[0] * val0 - v1[1] * val1 + v1[2] * val2)

        # Preserve the existing ABACUS cube reader behavior for negative atom counts.
        for _ in range(max(natom, 0)):
            stream.readline()

        count = shape[0] * shape[1] * shape[2]
        values = []
        for line in stream:
            values.extend(float(value) for value in line.split())
            if len(values) >= count:
                break

    if len(values) < count:
        raise ValueError(f"expected {count} grid values in {input_path}, found {len(values)}")
    return tuple(shape), volume, values[:count]


def integrate_cube(path):
    _, volume, values = read_cube(path)
    return sum(values) * volume


def splitmix64(value):
    value = (value + SPLITMIX_GOLDEN) & UINT64_MASK
    value = ((value ^ (value >> 30)) * SPLITMIX_MIX1) & UINT64_MASK
    value = ((value ^ (value >> 27)) * SPLITMIX_MIX2) & UINT64_MASK
    return value ^ (value >> 31)


def read_wfc_components(paths):
    if len(paths) == 0 or len(paths) % 2 != 0:
        raise ValueError("expected one or more RE IM cube path pairs")

    common_shape = None
    common_volume = None
    components = []
    for index in range(0, len(paths), 2):
        real_path = paths[index]
        imag_path = paths[index + 1]
        real_shape, real_volume, real_values = read_cube(real_path)
        imag_shape, imag_volume, imag_values = read_cube(imag_path)
        if imag_shape != real_shape:
            raise ValueError(
                f"grid shape mismatch in {imag_path}: expected {real_shape}, found {imag_shape}"
            )
        if imag_volume != real_volume:
            raise ValueError(
                f"voxel volume mismatch in {imag_path}: expected {real_volume}, found {imag_volume}"
            )
        if common_shape is None:
            common_shape = real_shape
            common_volume = real_volume
        elif real_shape != common_shape:
            raise ValueError(
                f"component grid shape mismatch in {real_path}: "
                f"expected {common_shape}, found {real_shape}"
            )
        elif real_volume != common_volume:
            raise ValueError(
                f"component voxel volume mismatch in {real_path}: "
                f"expected {common_volume}, found {real_volume}"
            )
        components.append(
            [complex(real, imag) for real, imag in zip(real_values, imag_values)]
        )
    return common_shape, common_volume, components


def wfc_fingerprint(paths):
    shape, volume, components = read_wfc_components(paths)
    point_count = sum(len(component) for component in components)
    energy = sum(
        value.real * value.real + value.imag * value.imag
        for component in components
        for value in component
    )
    if not math.isfinite(energy) or energy == 0.0:
        raise ValueError("wavefunction has zero or non-finite norm")

    projections = [0.0j] * FINGERPRINT_SKETCHES
    flat_index = 0
    for component in components:
        for value in component:
            for sketch in range(FINGERPRINT_SKETCHES):
                weight_index = splitmix64(flat_index + sketch * SPLITMIX_GOLDEN) & 3
                projections[sketch] += FINGERPRINT_WEIGHTS[weight_index] * value
            flat_index += 1

    norm = math.sqrt(energy)
    normalized = [projection / norm for projection in projections]
    metrics = [
        ("components", len(components)),
        ("nx", shape[0]),
        ("ny", shape[1]),
        ("nz", shape[2]),
        ("voxel", volume),
        ("rms", math.sqrt(energy / point_count)),
    ]
    metrics.extend(
        (f"power_{sketch}", abs(normalized[sketch]) ** 2)
        for sketch in range(FINGERPRINT_SKETCHES)
    )
    for sketch in range(FINGERPRINT_SKETCHES):
        next_sketch = (sketch + 1) % FINGERPRINT_SKETCHES
        cross = normalized[sketch] * normalized[next_sketch].conjugate()
        metrics.append((f"cross_{sketch}_re", cross.real))
        metrics.append((f"cross_{sketch}_im", cross.imag))

    if any(not math.isfinite(float(value)) for _, value in metrics):
        raise ValueError("wavefunction fingerprint contains a non-finite value")
    return metrics


def format_metric(value):
    if isinstance(value, int):
        return str(value)
    return f"{value:.12e}"


def read_wfc_component(calculated_re, calculated_im, reference_re, reference_im):
    paths = (calculated_re, calculated_im, reference_re, reference_im)
    cubes = [read_cube(path) for path in paths]
    shape = cubes[0][0]
    for path, cube in zip(paths[1:], cubes[1:]):
        if cube[0] != shape:
            raise ValueError(f"grid shape mismatch in {path}: expected {shape}, found {cube[0]}")

    calculated = [complex(real, imag) for real, imag in zip(cubes[0][2], cubes[1][2])]
    reference = [complex(real, imag) for real, imag in zip(cubes[2][2], cubes[3][2])]
    return shape, calculated, reference


def phase_aligned_error(file_groups):
    calculated_components = []
    reference_components = []
    common_shape = None
    for group in file_groups:
        shape, calculated, reference = read_wfc_component(*group)
        if common_shape is None:
            common_shape = shape
        elif shape != common_shape:
            raise ValueError(
                f"spinor component grid shape mismatch in {group[0]}: "
                f"expected {common_shape}, found {shape}"
            )
        calculated_components.append(calculated)
        reference_components.append(reference)

    calculated_norm = sum(
        abs(value) ** 2 for component in calculated_components for value in component
    )
    reference_norm = sum(
        abs(value) ** 2 for component in reference_components for value in component
    )
    if not math.isfinite(calculated_norm) or calculated_norm == 0.0:
        raise ValueError("calculated wavefunction has zero or non-finite norm")
    if not math.isfinite(reference_norm) or reference_norm == 0.0:
        raise ValueError("reference wavefunction has zero or non-finite norm")

    overlap = sum(
        reference.conjugate() * calculated
        for calculated_component, reference_component in zip(
            calculated_components, reference_components
        )
        for calculated, reference in zip(calculated_component, reference_component)
    )
    if (
        not math.isfinite(overlap.real)
        or not math.isfinite(overlap.imag)
        or abs(overlap) == 0.0
    ):
        raise ValueError("wavefunction overlap is zero or non-finite; the global phase is undefined")

    phase = overlap.conjugate() / abs(overlap)
    return max(
        abs(phase * calculated - reference)
        for calculated_component, reference_component in zip(
            calculated_components, reference_components
        )
        for calculated, reference in zip(calculated_component, reference_component)
    )


def parse_wfc_groups(paths):
    if len(paths) % 4 != 0:
        raise ValueError("expected four cube paths per component: CAL_RE CAL_IM REF_RE REF_IM")
    return [paths[index : index + 4] for index in range(0, len(paths), 4)]


def read_dense_shape(directory):
    pattern = re.compile(
        r"fft grid for dense charge/potential\s*=\s*\[\s*([0-9]+),\s*([0-9]+),\s*([0-9]+)\s*\]"
    )
    for log_path in sorted(directory.glob("running_*.log")):
        with log_path.open("r", encoding="utf-8") as stream:
            for line in stream:
                match = pattern.search(line)
                if match:
                    return tuple(int(value) for value in match.groups())
    return None


def max_abs_difference(actual, expected):
    return max(abs(lhs - rhs) for lhs, rhs in zip(actual, expected))


def check_cube_metadata(path):
    with path.open("r", encoding="utf-8") as stream:
        stream.readline()
        metadata = stream.readline()
    if "Fermi energy" in metadata:
        raise ValueError(f"unexpected Fermi-energy metadata in {path}")


def check_spinor_state(directory, band, kpoint, tolerance, expected_shape):
    prefix = f"wfi{band}"
    pchg_prefix = f"pchgi{band}"
    paths = {
        "norm": directory / f"{prefix}s1k{kpoint}.cube",
        "up_re": directory / f"{prefix}s1k{kpoint}re.cube",
        "up_im": directory / f"{prefix}s1k{kpoint}im.cube",
        "down_re": directory / f"{prefix}s2k{kpoint}re.cube",
        "down_im": directory / f"{prefix}s2k{kpoint}im.cube",
        "rho0": directory / f"{pchg_prefix}s1k{kpoint}.cube",
        "mx": directory / f"{pchg_prefix}s2k{kpoint}.cube",
        "my": directory / f"{pchg_prefix}s3k{kpoint}.cube",
        "mz": directory / f"{pchg_prefix}s4k{kpoint}.cube",
    }
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise ValueError("missing spinor output files: " + ", ".join(missing))

    shape = None
    data = {}
    for name, path in paths.items():
        check_cube_metadata(path)
        current_shape, _, data[name] = read_cube(path)
        if shape is None:
            shape = current_shape
        elif current_shape != shape:
            raise ValueError(f"grid shape mismatch in {path}")
    if expected_shape is not None and shape != expected_shape:
        raise ValueError(f"output grid {shape} does not match dense charge grid {expected_shape}")

    rho0_expected = []
    mx_expected = []
    my_expected = []
    mz_expected = []
    norm_squared = []
    for up_re, up_im, down_re, down_im, norm in zip(
        data["up_re"], data["up_im"], data["down_re"], data["down_im"], data["norm"]
    ):
        up_norm = up_re * up_re + up_im * up_im
        down_norm = down_re * down_re + down_im * down_im
        rho0_expected.append(up_norm + down_norm)
        mx_expected.append(2.0 * (up_re * down_re + up_im * down_im))
        my_expected.append(2.0 * (up_re * down_im - down_re * up_im))
        mz_expected.append(up_norm - down_norm)
        norm_squared.append(norm * norm)

    errors = {
        "norm": max_abs_difference(norm_squared, rho0_expected),
        "rho0": max_abs_difference(data["rho0"], rho0_expected),
        "mx": max_abs_difference(data["mx"], mx_expected),
        "my": max_abs_difference(data["my"], my_expected),
        "mz": max_abs_difference(data["mz"], mz_expected),
    }
    failed = {
        name: error
        for name, error in errors.items()
        if not math.isfinite(error) or error > tolerance
    }
    if failed:
        details = ", ".join(f"{name}={error:.6e}" for name, error in failed.items())
        raise ValueError(f"band {band}, k-point {kpoint}: pointwise identity failure ({details})")


def check_spinor_directory(directory, tolerance):
    pattern = re.compile(r"^pchgi([0-9]+)s1k([0-9]+)[.]cube$")
    states = []
    for path in directory.iterdir():
        match = pattern.match(path.name)
        if match:
            states.append((int(match.group(1)), int(match.group(2))))

    if not states:
        raise ValueError(f"no separate-k PW spinor partial-density files found in {directory}")

    expected_shape = read_dense_shape(directory)
    for band, kpoint in sorted(states):
        check_spinor_state(directory, band, kpoint, tolerance, expected_shape)


def run_integrate(args):
    print(f"{integrate_cube(args.cube):.10g}")


def run_compare_wfc(args):
    print(f"{phase_aligned_error(parse_wfc_groups(args.cube_paths)):.10g}")


def run_fingerprint_wfc(args):
    for name, value in wfc_fingerprint(args.cube_paths):
        print(f"{name} {format_metric(value)}")


def run_check_spinor(args):
    check_spinor_directory(args.directory, args.tolerance)


def build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    subparsers = parser.add_subparsers(dest="command")

    integrate_parser = subparsers.add_parser(
        "integrate", help="integrate one scalar cube file"
    )
    integrate_parser.add_argument("cube", type=Path)
    integrate_parser.set_defaults(handler=run_integrate)

    compare_parser = subparsers.add_parser(
        "compare-wfc", help="compare one complex state after global phase alignment"
    )
    compare_parser.add_argument(
        "cube_paths",
        nargs="+",
        help="CAL_RE CAL_IM REF_RE REF_IM, repeated for spinor components",
    )
    compare_parser.set_defaults(handler=run_compare_wfc)

    fingerprint_parser = subparsers.add_parser(
        "fingerprint-wfc",
        help="generate a deterministic global-phase-invariant wavefunction fingerprint",
    )
    fingerprint_parser.add_argument(
        "cube_paths",
        nargs="+",
        help="RE IM, repeated in component order for spinors",
    )
    fingerprint_parser.set_defaults(handler=run_fingerprint_wfc)

    spinor_parser = subparsers.add_parser(
        "check-spinor", help="check pointwise PW nspin=4 Pauli identities"
    )
    spinor_parser.add_argument("directory", type=Path)
    spinor_parser.add_argument("--tolerance", type=float, default=1.0e-8)
    spinor_parser.set_defaults(handler=run_check_spinor)
    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()
    if not hasattr(args, "handler"):
        parser.error("a subcommand is required")
    try:
        args.handler(args)
    except (OSError, ValueError) as error:
        print(error, file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
