#!/usr/bin/env python3

import argparse
import csv
import math
import re
from pathlib import Path


VERSION_LINE = "ABACUS_STERNHEIMER_GRID_DIAGNOSTICS 1"
OPERATOR_PATTERN = re.compile(r"^STERNHEIMER_DELTA_GRID_MATRICES_spin_(\d+)\.dat$")
PERTURBATION_PATTERN = re.compile(r"^STERNHEIMER_DELTA_PERTURBATION_spin_(\d+)\.dat$")
OPERATOR_QUANTITIES = (
    "overlap",
    "kinetic",
    "local_potential",
    "nonlocal",
    "hamiltonian",
    "occupied_virtual_overlap",
)


def _parse_grid_file(path):
    path = Path(path)
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if not lines or lines[0] != VERSION_LINE:
        raise ValueError("invalid grid diagnostic version in {}".format(path))

    metadata = {}
    sections = {}
    current = None
    for line in lines[1:]:
        fields = line.split()
        key = fields[0]
        if key == "matrix":
            if len(fields) != 2:
                raise ValueError("invalid matrix header in {}".format(path))
            current = fields[1]
            if current in sections:
                raise ValueError("duplicate section {} in {}".format(current, path))
            sections[current] = {}
            continue
        if key == "tensor":
            if fields != ["tensor", "perturbation", "occupied", "virtual", "auxiliary"]:
                raise ValueError("invalid tensor header in {}".format(path))
            current = "perturbation"
            sections[current] = {}
            continue
        if current is None:
            if key == "grid" and len(fields) == 4:
                metadata[key] = tuple(int(value) for value in fields[1:])
            elif key in ("spin", "occupied", "virtuals", "auxiliaries") and len(fields) == 2:
                metadata[key] = int(fields[1])
            elif key == "volume_element" and len(fields) == 2:
                metadata[key] = float(fields[1])
            elif key in ("kind", "energy_unit", "potential_unit", "storage") and len(fields) == 2:
                metadata[key] = fields[1]
            else:
                raise ValueError("invalid metadata line in {}: {}".format(path, line))
            continue

        coordinate_size = 3 if current == "perturbation" else 2
        if len(fields) != coordinate_size + 2:
            raise ValueError("invalid data row in {}: {}".format(path, line))
        coordinate = tuple(int(value) for value in fields[:coordinate_size])
        if coordinate in sections[current]:
            raise ValueError("duplicate coordinate in {}: {}".format(path, coordinate))
        sections[current][coordinate] = complex(float(fields[-2]), float(fields[-1]))

    required_metadata = ("grid", "spin", "occupied", "virtuals", "auxiliaries", "volume_element", "kind")
    if any(key not in metadata for key in required_metadata):
        raise ValueError("incomplete metadata in {}".format(path))
    expected_sections = ("perturbation",) if metadata["kind"] == "delta_perturbation_tensor" else OPERATOR_QUANTITIES
    if set(sections) != set(expected_sections):
        raise ValueError("unexpected sections in {}".format(path))
    for section, values in sections.items():
        if section == "perturbation":
            expected = metadata["occupied"] * metadata["virtuals"] * metadata["auxiliaries"]
        elif section == "occupied_virtual_overlap":
            expected = metadata["occupied"] * metadata["virtuals"]
        else:
            expected = metadata["virtuals"] * metadata["virtuals"]
        if len(values) != expected:
            raise ValueError("invalid {} size in {}".format(section, path))
    return metadata, sections


def _indexed_files(directory, pattern):
    indexed = {}
    for path in Path(directory).iterdir():
        match = pattern.match(path.name)
        if match is None:
            continue
        spin = int(match.group(1))
        if spin in indexed:
            raise ValueError("duplicate spin {} in {}".format(spin, directory))
        indexed[spin] = path
    return indexed


def _validate_pair(metadata40, metadata50, values40, values50, path40, path50):
    for key in ("spin", "occupied", "virtuals", "auxiliaries", "kind"):
        if metadata40[key] != metadata50[key]:
            raise ValueError("{} mismatch: {} vs {}".format(key, path40, path50))
    if set(values40) != set(values50):
        raise ValueError("coordinate mismatch: {} vs {}".format(path40, path50))


def _invariant_profile(quantity, values):
    if quantity == "perturbation":
        squared_norms = {}
        for (occupied, virtual, auxiliary), value in values.items():
            key = (occupied, auxiliary)
            squared_norms[key] = squared_norms.get(key, 0.0) + abs(value) ** 2
        return "occupied_auxiliary_norm", [
            math.sqrt(squared_norms[key]) for key in sorted(squared_norms)
        ]
    if quantity == "occupied_virtual_overlap":
        squared_norms = {}
        for (occupied, virtual), value in values.items():
            squared_norms[occupied] = squared_norms.get(occupied, 0.0) + abs(value) ** 2
        return "occupied_norm", [math.sqrt(squared_norms[key]) for key in sorted(squared_norms)]
    diagonal = [value for (row, column), value in sorted(values.items()) if row == column]
    return "diagonal", diagonal


def _comparison_row(kind, spin, quantity, metadata40, metadata50, values40, values50):
    ordered_coordinates = sorted(values40)
    vector40 = [values40[index] for index in ordered_coordinates]
    vector50 = [values50[index] for index in ordered_coordinates]
    differences = [value50 - value40 for value40, value50 in zip(vector40, vector50)]
    norm40 = math.sqrt(sum(abs(value) ** 2 for value in vector40))
    norm50 = math.sqrt(sum(abs(value) ** 2 for value in vector50))
    difference_norm = math.sqrt(sum(abs(value) ** 2 for value in differences))
    relative_difference = difference_norm / norm40 if norm40 > 0.0 else difference_norm
    relative_norm_change = abs(norm50 - norm40) / norm40 if norm40 > 0.0 else norm50
    profile_kind40, profile40 = _invariant_profile(quantity, values40)
    profile_kind50, profile50 = _invariant_profile(quantity, values50)
    if profile_kind40 != profile_kind50 or len(profile40) != len(profile50):
        raise ValueError("diagnostic profile mismatch for {}".format(quantity))
    profile_difference_norm = math.sqrt(
        sum(abs(value50 - value40) ** 2 for value40, value50 in zip(profile40, profile50))
    )
    profile_norm40 = math.sqrt(sum(abs(value) ** 2 for value in profile40))
    profile_relative_difference = (
        profile_difference_norm / profile_norm40 if profile_norm40 > 0.0 else profile_difference_norm
    )
    return {
        "kind": kind,
        "spin": spin,
        "quantity": quantity,
        "nvalues": len(vector40),
        "norm40": "{:.17e}".format(norm40),
        "norm50": "{:.17e}".format(norm50),
        "difference_norm": "{:.17e}".format(difference_norm),
        "relative_difference": "{:.17e}".format(relative_difference),
        "relative_norm_change": "{:.17e}".format(relative_norm_change),
        "profile_kind": profile_kind40,
        "profile_relative_difference": "{:.17e}".format(profile_relative_difference),
        "max_abs_difference": "{:.17e}".format(max((abs(value) for value in differences), default=0.0)),
        "grid40": "x".join(str(value) for value in metadata40["grid"]),
        "grid50": "x".join(str(value) for value in metadata50["grid"]),
        "volume_element40": "{:.17e}".format(metadata40["volume_element"]),
        "volume_element50": "{:.17e}".format(metadata50["volume_element"]),
    }


def compare_campaign(run40, run50, output_path):
    run40 = Path(run40)
    run50 = Path(run50)
    output_path = Path(output_path)
    if output_path.exists():
        raise FileExistsError(output_path)

    rows = []
    for kind, pattern in (("operator", OPERATOR_PATTERN), ("tensor", PERTURBATION_PATTERN)):
        files40 = _indexed_files(run40, pattern)
        files50 = _indexed_files(run50, pattern)
        if not files40 or set(files40) != set(files50):
            raise ValueError("{} spin files mismatch".format(kind))
        for spin in sorted(files40):
            metadata40, sections40 = _parse_grid_file(files40[spin])
            metadata50, sections50 = _parse_grid_file(files50[spin])
            if set(sections40) != set(sections50):
                raise ValueError("section mismatch: {} vs {}".format(files40[spin], files50[spin]))
            for quantity in sections40:
                _validate_pair(
                    metadata40, metadata50, sections40[quantity], sections50[quantity], files40[spin], files50[spin]
                )
                rows.append(
                    _comparison_row(
                        kind, spin, quantity, metadata40, metadata50, sections40[quantity], sections50[quantity]
                    )
                )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0])
    with output_path.open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def main(argv=None):
    parser = argparse.ArgumentParser(description="Compare ABACUS Delta-ST grid diagnostic text files.")
    parser.add_argument("run40", type=Path)
    parser.add_argument("run50", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args(argv)
    print(compare_campaign(args.run40, args.run50, args.output))


if __name__ == "__main__":
    main()
