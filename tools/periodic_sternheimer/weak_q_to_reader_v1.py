#!/usr/bin/env python3
"""Package WEAK_Q potential-potential response matrices for LibRPA.

WEAK_Q writes M = V chi0 V.  LibRPA's dedicated Sternheimer RPA task forms
Pi = V^(-1/2) M V^(-1/2), so this module deliberately applies no Coulomb
transformation.  Only the numerical Hermitian projection required by the
reader-v1 upper atom-block format is applied and audited.
"""

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import struct

import numpy as np


CHI0_V1_MARKER = -41073291
COULOMB_V1_MARKER = -20129433
COMPLEX_FLAG = 1
HARTREE_TO_EV = 27.211386245988


class ReaderV1:
    def __init__(self, iq, ifrequency, omega, weight, atom_naux, matrix):
        self.iq = iq
        self.ifrequency = ifrequency
        self.omega = omega
        self.weight = weight
        self.atom_naux = tuple(atom_naux)
        self.matrix = matrix


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _pairs(natoms):
    return tuple((i, j) for i in range(natoms) for j in range(i, natoms))


def _offsets(atom_naux):
    result = [0]
    for count in atom_naux:
        _require(type(count) is int and count > 0, "atom_naux entries must be positive integers")
        result.append(result[-1] + count)
    return result


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _matrix_summary(matrix):
    norm = float(np.linalg.norm(matrix))
    anti = matrix - matrix.conj().T
    return {
        "frobenius_norm": norm,
        "relative_hermiticity_defect": float(np.linalg.norm(anti)) / norm if norm else 0.0,
        "max_abs": float(np.max(np.abs(matrix))) if matrix.size else 0.0,
        "max_abs_antihermitian": float(np.max(np.abs(anti))) if matrix.size else 0.0,
        "max_abs_imag_diagonal": float(np.max(np.abs(np.diag(matrix).imag))) if matrix.size else 0.0,
    }


def _read_header_and_matrix(path, expected_marker):
    path = Path(path)
    data = path.read_bytes()
    position = 0

    def unpack(fmt, label):
        nonlocal position
        size = struct.calcsize("<" + fmt)
        _require(position + size <= len(data), f"truncated {label}: {path}")
        values = struct.unpack_from("<" + fmt, data, position)
        position += size
        return values[0] if len(values) == 1 else values

    if expected_marker == CHI0_V1_MARKER:
        marker, iq, ifrequency, naux, value_flag, natoms = unpack("6i", "response header")
        omega, weight = unpack("2d", "response frequency metadata")
        nblocks = unpack("i", "response block count")
    else:
        marker, iq, naux, value_flag, natoms, nblocks = unpack("6i", "Coulomb header")
        ifrequency, omega, weight = 0, 0.0, 1.0
    _require(marker == expected_marker, f"invalid reader-v1 marker: {path}")
    _require(iq > 0 and naux > 0 and natoms > 0, f"invalid reader-v1 dimensions: {path}")
    _require(value_flag in (0, COMPLEX_FLAG), f"invalid reader-v1 value flag: {path}")
    if expected_marker == CHI0_V1_MARKER:
        _require(ifrequency > 0 and value_flag == COMPLEX_FLAG,
                 f"invalid response metadata: {path}")
        _require(math.isfinite(omega) and math.isfinite(weight) and weight > 0.0,
                 f"invalid response frequency metadata: {path}")

    atom_naux = unpack(f"{natoms}i", "atom auxiliary dimensions")
    if natoms == 1:
        atom_naux = (atom_naux,)
    atom_naux = tuple(int(value) for value in atom_naux)
    offsets = _offsets(atom_naux)
    _require(offsets[-1] == naux, f"atom_naux dimension mismatch: {path}")
    pairs = _pairs(natoms)
    _require(0 <= nblocks <= len(pairs), f"invalid reader-v1 block count: {path}")

    records = []
    seen = set()
    for _ in range(nblocks):
        pair_index = unpack("i", "pair index")
        byte_offset = unpack("q", "block offset")
        _require(pair_index not in seen and 0 <= pair_index < len(pairs),
                 f"invalid or duplicate pair index: {path}")
        seen.add(pair_index)
        records.append((pair_index, byte_offset))

    matrix = np.zeros((naux, naux), dtype=np.complex128)
    ranges = []
    value_dtype = np.dtype("<c16") if value_flag == COMPLEX_FLAG else np.dtype("<f8")
    for pair_index, byte_offset in records:
        iatom, jatom = pairs[pair_index]
        ni, nj = atom_naux[iatom], atom_naux[jatom]
        byte_count = ni * nj * value_dtype.itemsize
        _require(position <= byte_offset <= len(data) - byte_count,
                 f"invalid reader-v1 payload offset: {path}")
        ranges.append((byte_offset, byte_offset + byte_count))
        block = np.frombuffer(data, dtype=value_dtype, count=ni * nj, offset=byte_offset)
        block = block.astype(np.complex128, copy=False).reshape((ni, nj))
        i0, j0 = offsets[iatom], offsets[jatom]
        matrix[i0:i0 + ni, j0:j0 + nj] = block
        if iatom != jatom:
            matrix[j0:j0 + nj, i0:i0 + ni] = block.conj().T
    ranges.sort()
    _require(all(right <= next_left for (_, right), (next_left, _) in zip(ranges, ranges[1:])),
             f"overlapping reader-v1 payloads: {path}")
    _require(np.isfinite(matrix).all(), f"reader-v1 matrix contains non-finite values: {path}")
    return ReaderV1(iq, ifrequency, omega, weight, atom_naux, matrix)


def read_reader_v1(path):
    return _read_header_and_matrix(path, CHI0_V1_MARKER)


def read_coulomb_v1(path):
    return _read_header_and_matrix(path, COULOMB_V1_MARKER)


def read_full_coulomb_text(path):
    path = Path(path)
    rows = None
    columns = None
    matrix = None
    seen = None
    record_count = 0
    complete = False

    with path.open(encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, start=1):
            fields = line.split()
            if not fields:
                continue
            if fields[0] == "full_matrix_rows":
                _require(len(fields) == 2 and rows is None,
                         f"invalid row metadata at line {line_number}")
                rows = int(fields[1])
            elif fields[0] == "full_matrix_columns":
                _require(len(fields) == 2 and columns is None,
                         f"invalid column metadata at line {line_number}")
                columns = int(fields[1])
            elif fields[0] == "coulomb_integral":
                _require(rows is not None and columns is not None,
                         f"matrix dimensions must precede records at line {line_number}")
                _require(rows > 0 and rows == columns,
                         "matrix dimensions must be positive and square")
                if matrix is None:
                    matrix = np.empty((rows, columns), dtype=np.complex128)
                    seen = np.zeros((rows, columns), dtype=np.bool_)
                _require(len(fields) == 5, f"invalid matrix record at line {line_number}")
                row, column = int(fields[1]) - 1, int(fields[2]) - 1
                _require(0 <= row < rows and 0 <= column < columns,
                         f"matrix record is out of bounds at line {line_number}")
                _require(not seen[row, column], f"duplicate matrix record at line {line_number}")
                value = complex(float(fields[3]), float(fields[4]))
                _require(math.isfinite(value.real) and math.isfinite(value.imag),
                         f"non-finite matrix record at line {line_number}")
                matrix[row, column] = value
                seen[row, column] = True
                record_count += 1
            elif fields[0] == "full_matrix_complete":
                _require(fields == ["full_matrix_complete", "yes"],
                         f"invalid completion marker at line {line_number}")
                complete = True

    _require(rows is not None and columns is not None and rows > 0 and rows == columns,
             "matrix dimensions must be positive and square")
    _require(complete, "matrix completion marker is missing")
    _require(matrix is not None and record_count == rows * columns and bool(np.all(seen)),
             "matrix coverage is incomplete")
    return matrix


def package_matrix(matrix, output, *, iq, ifrequency, omega, weight, atom_naux,
                   hermitian_mode="project"):
    output = Path(output)
    if output.exists():
        raise FileExistsError(output)
    matrix = np.asarray(matrix, dtype=np.complex128)
    atom_naux = tuple(atom_naux)
    offsets = _offsets(atom_naux)
    _require(iq > 0, "iq must be positive")
    _require(ifrequency > 0, "ifrequency must be positive")
    _require(math.isfinite(omega) and omega >= 0.0, "omega must be finite and non-negative")
    _require(math.isfinite(weight) and weight > 0.0, "weight must be finite and positive")
    _require(matrix.ndim == 2 and matrix.shape[0] == matrix.shape[1], "matrix must be square")
    _require(matrix.shape == (offsets[-1], offsets[-1]), "matrix dimension does not match atom_naux")
    _require(np.isfinite(matrix).all(), "matrix values must be finite")

    _require(hermitian_mode in ("project", "sum"), "invalid Hermitian completion mode")
    hermitian = (0.5 * (matrix + matrix.conj().T)
                 if hermitian_mode == "project" else matrix + matrix.conj().T)
    pairs = _pairs(len(atom_naux))
    header_size = 7 * 4 + 2 * 8 + 4 * len(atom_naux) + 12 * len(pairs)
    records = []
    byte_offset = header_size
    for pair_index, (iatom, jatom) in enumerate(pairs):
        records.append((pair_index, byte_offset))
        byte_offset += atom_naux[iatom] * atom_naux[jatom] * 16

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.tmp.{os.getpid()}")
    try:
        with temporary.open("xb") as stream:
            stream.write(struct.pack("<6i2di", CHI0_V1_MARKER, iq, ifrequency, offsets[-1],
                                     COMPLEX_FLAG, len(atom_naux), omega, weight, len(pairs)))
            stream.write(struct.pack(f"<{len(atom_naux)}i", *atom_naux))
            for pair_index, payload_offset in records:
                stream.write(struct.pack("<iq", pair_index, payload_offset))
            for iatom, jatom in pairs:
                i0, i1 = offsets[iatom], offsets[iatom + 1]
                j0, j1 = offsets[jatom], offsets[jatom + 1]
                block = np.ascontiguousarray(hermitian[i0:i1, j0:j1], dtype="<c16")
                stream.write(block.tobytes(order="C"))
            stream.flush()
            os.fsync(stream.fileno())
        os.link(temporary, output)
        temporary.unlink()
    except Exception:
        temporary.unlink(missing_ok=True)
        raise

    restored = read_reader_v1(output)
    denominator = float(np.linalg.norm(hermitian))
    round_trip = float(np.linalg.norm(restored.matrix - hermitian)) / denominator if denominator else 0.0
    return {
        "path": str(output),
        "sha256": _sha256(output),
        "matrix_kind": "potential_potential_response",
        "mathematical_object": "M = V chi0 V",
        "coulomb_transform": "none",
        "projection": "(M + M^H) / 2" if hermitian_mode == "project" else "M + M^H",
        "raw": _matrix_summary(matrix),
        "packaged": _matrix_summary(restored.matrix),
        "round_trip_relative_frobenius": round_trip,
    }


def _fixed_q_routes(routes, *, iq, full_kpoint_count):
    _require(iq > 0 and full_kpoint_count > 0, "partial routes require positive q and full-k dimensions")
    by_member = {}
    representatives = set()
    normalized = []
    for raw in routes:
        _require(isinstance(raw, dict), "fixed-q route must be a mapping")
        try:
            representative = int(raw["representative_ik_full"])
            member = int(raw["member_ik_full"])
            spatial_isym = int(raw["spatial_isym"])
            time_reversal = bool(raw["time_reversal"])
            fold = tuple(int(value) for value in raw["fold_G"])
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError("invalid fixed-q route") from error
        _require(0 <= representative < full_kpoint_count and 0 <= member < full_kpoint_count,
                 "fixed-q route k index is outside the full grid")
        _require(spatial_isym >= 0 and len(fold) == 3, "invalid fixed-q route symmetry metadata")
        _require(member not in by_member, "duplicate fixed-q route member")
        item = {
            "representative_ik_full": representative,
            "member_ik_full": member,
            "spatial_isym": spatial_isym,
            "time_reversal": time_reversal,
            "fold_G": fold,
        }
        by_member[member] = item
        representatives.add(representative)
        normalized.append(item)
    _require(set(by_member) == set(range(full_kpoint_count)),
             "fixed-q routes do not cover the full k grid")
    for representative in representatives:
        _require(representative in by_member
                 and by_member[representative]["representative_ik_full"] == representative,
                 "fixed-q route has no representative identity member")
    return tuple(sorted(normalized, key=lambda item: item["member_ik_full"])), tuple(sorted(representatives))


def _write_full_kpoints(path, full_kpoints):
    with Path(path).open("x", encoding="utf-8") as stream:
        stream.write("# ik_full kx ky kz\n")
        for ik, point in enumerate(full_kpoints):
            stream.write(f"{ik} {point[0]:.17g} {point[1]:.17g} {point[2]:.17g}\n")


def _write_fixed_q_routes(path, *, iq, routes):
    with Path(path).open("x", encoding="utf-8") as stream:
        stream.write("version 1\n")
        stream.write("# iq representative_ik member_ik spatial_isym time_reversal fold_Gx fold_Gy fold_Gz\n")
        for route in routes:
            fold = route["fold_G"]
            stream.write(
                f"{iq} {route['representative_ik_full']} {route['member_ik_full']} "
                f"{route['spatial_isym']} {int(route['time_reversal'])} {fold[0]} {fold[1]} {fold[2]}\n")


def package_partial_k_response(matrices, output_dir, *, iq, frequencies, atom_naux,
                               full_kpoints, fixed_q_routes, qpoint, qweight):
    """Write a fixed-q representative-k response bundle for LibRPA restoration.

    `matrices[(ik_full, ifrequency)]` must contain only full occupied-subspace
    contributions from fixed-q representatives.  The caller is responsible for
    combining weak-unit band and column shards before invoking this function.
    """
    output_dir = Path(output_dir)
    _require(not output_dir.exists(), f"output directory already exists: {output_dir}")
    _require(isinstance(matrices, dict) and matrices, "partial response matrices are empty")
    _require(isinstance(frequencies, dict) and frequencies, "partial response frequencies are empty")
    _require(len(qpoint) == 3 and all(math.isfinite(float(value)) for value in qpoint),
             "qpoint must contain three finite coordinates")
    _require(math.isfinite(qweight) and qweight > 0.0, "qweight must be finite and positive")
    full_kpoints = tuple(tuple(float(value) for value in point) for point in full_kpoints)
    _require(full_kpoints and all(len(point) == 3 and all(math.isfinite(value) for value in point)
                                  for point in full_kpoints),
             "full k-point manifest is invalid")
    routes, representatives = _fixed_q_routes(
        fixed_q_routes, iq=iq, full_kpoint_count=len(full_kpoints))
    normalized_frequencies = {}
    for ifrequency, values in frequencies.items():
        ifrequency = int(ifrequency)
        _require(ifrequency > 0 and len(values) == 2, "invalid partial response frequency metadata")
        omega, weight = float(values[0]), float(values[1])
        _require(math.isfinite(omega) and omega >= 0 and math.isfinite(weight) and weight > 0,
                 "non-finite partial response frequency metadata")
        normalized_frequencies[ifrequency] = (omega, weight)
    expected = {(representative, ifrequency)
                for representative in representatives for ifrequency in normalized_frequencies}
    supplied = {(int(key[0]), int(key[1])) for key in matrices}
    _require(supplied == expected,
             "partial response matrices must cover every frequency of each representative and no nonrepresentative k point")

    output_dir.mkdir(parents=True, exist_ok=False)
    try:
        records = []
        reports = []
        for ik_full, ifrequency in sorted(expected):
            matrix = matrices[(ik_full, ifrequency)]
            omega, weight = normalized_frequencies[ifrequency]
            filename = f"v1_sternheimer_chi0_iq_{iq}_ik_{ik_full}_ifreq_{ifrequency}.dat"
            report = package_matrix(matrix, output_dir / filename, iq=iq,
                                    ifrequency=ifrequency, omega=omega, weight=weight,
                                    atom_naux=atom_naux, hermitian_mode="sum")
            report.update({"ik_full": ik_full, "ifrequency": ifrequency, "file": filename})
            reports.append(report)
            records.append((ik_full, ifrequency, filename))

        manifest = output_dir / f"v1_sternheimer_partial_manifest_iq_{iq}.dat"
        with manifest.open("x", encoding="utf-8") as stream:
            stream.write("# iq ik_full ifreq response_file\n")
            for ik_full, ifrequency, filename in records:
                stream.write(f"{iq} {ik_full} {ifrequency} {filename}\n")
        _write_full_kpoints(output_dir / "v1_sternheimer_full_kpoints.dat", full_kpoints)
        _write_fixed_q_routes(output_dir / f"v1_sternheimer_symmetry_routes_iq_{iq}.dat",
                              iq=iq, routes=routes)
        qpoint_path = output_dir / f"v1_sternheimer_qpoint_iq_{iq}.dat"
        with qpoint_path.open("x", encoding="utf-8") as stream:
            stream.write(f"{iq} {qpoint[0]:.17g} {qpoint[1]:.17g} {qpoint[2]:.17g} {qweight:.17g}\n")
    except Exception:
        for path in output_dir.glob("*"):
            path.unlink()
        output_dir.rmdir()
        raise

    return {
        "partial_response_complete": True,
        "physical_result": False,
        "iq": iq,
        "representative_kpoints": list(representatives),
        "full_kpoint_count": len(full_kpoints),
        "frequency_count": len(normalized_frequencies),
        "fixed_q_route_count": len(routes),
        "partial_manifest": manifest.name,
        "symmetry_routes": f"v1_sternheimer_symmetry_routes_iq_{iq}.dat",
        "full_kpoints_manifest": "v1_sternheimer_full_kpoints.dat",
        "qpoint_fragment": qpoint_path.name,
        "outputs": reports,
    }


def _read_audit_records(path):
    values = {}
    lists = {}
    with Path(path).open(encoding="utf-8") as stream:
        for line in stream:
            fields = line.split()
            if not fields or fields[0].startswith("#"):
                continue
            key = fields[0]
            if key in {"column_file", "fixed_q_route"}:
                lists.setdefault(key, []).append(fields[1:])
            else:
                _require(len(fields) >= 2, f"invalid audit record: {path}")
                values[key] = fields[1:]
    return values, lists


def _audit_value(values, key, path):
    _require(key in values and len(values[key]) == 1, f"missing audit field {key}: {path}")
    return values[key][0]


def _read_weak_columns(path, channels):
    values, _ = _read_audit_records(path)
    ifrequency = int(_audit_value(values, "ifrequency", path))
    omega = float(_audit_value(values, "omega_Ha", path))
    weight = float(_audit_value(values, "weight_Ha", path))
    input_hash = _audit_value(values, "input_manifest_sha256", path)
    matrices = {}
    current_column = None
    seen = set()
    with Path(path).open(encoding="utf-8") as stream:
        for line in stream:
            fields = line.split()
            if not fields or fields[0].startswith("#"):
                continue
            if fields[0] == "owned_column":
                _require(len(fields) == 2, f"invalid owned column record: {path}")
                current_column = int(fields[1]) - 1
                _require(0 <= current_column < channels and current_column not in matrices,
                         f"invalid or duplicate owned column: {path}")
                matrices[current_column] = np.zeros((channels, channels), dtype=np.complex128)
            elif fields[0] == "response":
                _require(current_column is not None and len(fields) == 5,
                         f"response appears outside an owned column block: {path}")
                row, column = int(fields[1]) - 1, int(fields[2]) - 1
                _require(column == current_column and 0 <= row < channels,
                         f"invalid response column or row: {path}")
                key = (current_column, row)
                _require(key not in seen, f"duplicate weak response element: {path}")
                value = complex(float(fields[3]), float(fields[4]))
                _require(math.isfinite(value.real) and math.isfinite(value.imag),
                         f"non-finite weak response element: {path}")
                matrices[current_column][row, current_column] = value
                seen.add(key)
            elif fields[0] == "columns_complete":
                _require(fields == ["columns_complete", "yes"],
                         f"invalid weak column completion marker: {path}")
    for column in matrices:
        _require(sum(1 for item in seen if item[0] == column) == channels,
                 f"weak column coverage is incomplete: {path}")
    return ifrequency, omega, weight, input_hash, matrices


def _read_weak_input_manifest(path):
    values, _ = _read_audit_records(path)
    full_kpoint_count = int(_audit_value(values, "full_kpoints", path))
    iq = int(_audit_value(values, "iq", path))
    qpoint = tuple(float(value) for value in values["qpoint"])
    _require(len(qpoint) == 3 and all(math.isfinite(value) for value in qpoint),
             f"invalid input-manifest q point: {path}")
    kpoints = {}
    channel_atoms = {}
    with Path(path).open(encoding="utf-8") as stream:
        for line in stream:
            fields = line.split()
            if not fields or fields[0].startswith("#"):
                continue
            if fields[0] == "kpoint":
                _require(len(fields) >= 5, f"invalid input-manifest k point: {path}")
                ik = int(fields[1]) - 1
                _require(0 <= ik < full_kpoint_count and ik not in kpoints,
                         f"invalid or duplicate input-manifest k point: {path}")
                point = tuple(float(value) for value in fields[2:5])
                _require(all(math.isfinite(value) for value in point),
                         f"non-finite input-manifest k point: {path}")
                kpoints[ik] = point
            elif fields[0] == "auxiliary_channel":
                _require(len(fields) >= 4 and fields[2] == "atom_index",
                         f"invalid auxiliary-channel manifest record: {path}")
                channel = int(fields[1]) - 1
                atom = int(fields[3])
                _require(channel >= 0 and atom >= 0 and channel not in channel_atoms,
                         f"invalid or duplicate auxiliary channel: {path}")
                channel_atoms[channel] = atom
    _require(set(kpoints) == set(range(full_kpoint_count)),
             f"input-manifest k-point coverage is incomplete: {path}")
    _require(channel_atoms and set(channel_atoms) == set(range(len(channel_atoms))),
             f"input-manifest auxiliary-channel coverage is incomplete: {path}")
    ordered_atoms = [channel_atoms[index] for index in range(len(channel_atoms))]
    _require(ordered_atoms == sorted(ordered_atoms),
             "auxiliary channels are not ordered by atom in the input manifest")
    atom_naux = []
    for atom in sorted(set(ordered_atoms)):
        _require(ordered_atoms.count(atom) > 0, "empty atom auxiliary block in input manifest")
        atom_naux.append(ordered_atoms.count(atom))
    return {
        "iq": iq,
        "qpoint": qpoint,
        "full_kpoints": tuple(kpoints[index] for index in range(full_kpoint_count)),
        "atom_naux": tuple(atom_naux),
        "channels": len(channel_atoms),
    }


def merge_weak_q_unit_columns(unit_audits, output_dir, *, qweight):
    """Merge complete occupied-band WEAK_Q units into representative-k v1 files."""
    _require(unit_audits, "at least one weak-q unit audit is required")
    bundles = []
    for audit_path in map(Path, unit_audits):
        values, lists = _read_audit_records(audit_path)
        _require(_audit_value(values, "unit_complete", audit_path) == "yes",
                 f"weak unit is incomplete: {audit_path}")
        _require(_audit_value(values, "all_source_bands", audit_path) == "yes",
                 f"weak unit does not contain the complete occupied subspace: {audit_path}")
        _require(values.get("symmetry_skipped_nonrepresentative", ["no"]) != ["yes"],
                 f"nonrepresentative weak unit cannot be merged: {audit_path}")
        input_path = audit_path.parent / _audit_value(values, "input_manifest", audit_path)
        input_manifest = _read_weak_input_manifest(input_path)
        source_k = int(_audit_value(values, "source_k", audit_path)) - 1
        iq = int(_audit_value(values, "iq", audit_path))
        routes = []
        for fields in lists.get("fixed_q_route", []):
            _require(len(fields) == 8, f"invalid fixed-q route in audit: {audit_path}")
            routes.append({
                "representative_ik_full": int(fields[1]),
                "member_ik_full": int(fields[2]),
                "spatial_isym": int(fields[3]),
                "time_reversal": bool(int(fields[4])),
                "fold_G": tuple(int(value) for value in fields[5:8]),
            })
        _require(routes, f"weak unit has no fixed-q symmetry routes: {audit_path}")
        column_files = lists.get("column_file", [])
        _require(column_files, f"weak unit has no column files: {audit_path}")
        columns = []
        for fields in column_files:
            _require(len(fields) == 2, f"invalid weak column-file record: {audit_path}")
            column_path = audit_path.parent / fields[1]
            columns.append(_read_weak_columns(column_path, input_manifest["channels"]))
        bundles.append((audit_path, values, input_manifest, source_k, iq, routes, columns))

    first = bundles[0]
    iq = first[4]
    input_manifest = first[2]
    _require(all(bundle[4] == iq and bundle[2] == input_manifest for bundle in bundles),
             "weak units do not share one physical q/input manifest")
    frequencies = {}
    matrices = {}
    matrix_columns = {}
    all_routes = {}
    for audit_path, values, _, source_k, _, routes, columns in bundles:
        _require(source_k >= 0, f"invalid weak source k: {audit_path}")
        for route in routes:
            key = route["member_ik_full"]
            _require(key not in all_routes or all_routes[key] == route,
                     f"conflicting fixed-q route: {audit_path}")
            all_routes[key] = route
        for ifrequency, omega, weight, input_hash, column_map in columns:
            expected_hash = _audit_value(values, "input_manifest_sha256", audit_path)
            _require(input_hash == expected_hash, f"column/input manifest mismatch: {audit_path}")
            _require(ifrequency not in frequencies or frequencies[ifrequency] == (omega, weight),
                     "conflicting weak frequency metadata")
            frequencies[ifrequency] = (omega, weight)
            key = (source_k, ifrequency)
            matrix = matrices.setdefault(
                key, np.zeros((input_manifest["channels"], input_manifest["channels"]), dtype=np.complex128))
            seen_columns = matrix_columns.setdefault(key, set())
            for column, values_by_row in column_map.items():
                _require(column not in seen_columns, "duplicate weak response column bundle")
                matrix[:, column] = values_by_row[:, column]
                seen_columns.add(column)
    for key, seen_columns in matrix_columns.items():
        _require(seen_columns == set(range(input_manifest["channels"])),
                 f"weak response column coverage is incomplete for source/frequency {key}")
    routes = tuple(all_routes.values())
    _require(routes, "weak merged route set is empty")
    report = package_partial_k_response(
        matrices,
        output_dir,
        iq=iq,
        frequencies=frequencies,
        atom_naux=input_manifest["atom_naux"],
        full_kpoints=input_manifest["full_kpoints"],
        fixed_q_routes=routes,
        qpoint=input_manifest["qpoint"],
        qweight=qweight,
    )
    report["unit_audits"] = [str(Path(path).resolve()) for path in unit_audits]
    report["unit_count"] = len(unit_audits)
    (Path(output_dir) / "PARTIAL_K_PACKAGING_COMPLETE.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return report


def package_coulomb_matrix(matrix, output, *, iq, atom_naux):
    output = Path(output)
    if output.exists():
        raise FileExistsError(output)
    matrix = np.asarray(matrix, dtype=np.complex128)
    atom_naux = tuple(atom_naux)
    offsets = _offsets(atom_naux)
    _require(iq > 0, "iq must be positive")
    _require(matrix.ndim == 2 and matrix.shape[0] == matrix.shape[1], "matrix must be square")
    _require(matrix.shape == (offsets[-1], offsets[-1]), "matrix dimension does not match atom_naux")
    _require(np.isfinite(matrix).all(), "matrix values must be finite")

    hermitian = 0.5 * (matrix + matrix.conj().T)
    pairs = _pairs(len(atom_naux))
    header_size = 6 * 4 + 4 * len(atom_naux) + 12 * len(pairs)
    records = []
    byte_offset = header_size
    for pair_index, (iatom, jatom) in enumerate(pairs):
        records.append((pair_index, byte_offset))
        byte_offset += atom_naux[iatom] * atom_naux[jatom] * 16

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.tmp.{os.getpid()}")
    try:
        with temporary.open("xb") as stream:
            stream.write(struct.pack("<6i", COULOMB_V1_MARKER, iq, offsets[-1], COMPLEX_FLAG,
                                     len(atom_naux), len(pairs)))
            stream.write(struct.pack(f"<{len(atom_naux)}i", *atom_naux))
            for pair_index, payload_offset in records:
                stream.write(struct.pack("<iq", pair_index, payload_offset))
            for iatom, jatom in pairs:
                i0, i1 = offsets[iatom], offsets[iatom + 1]
                j0, j1 = offsets[jatom], offsets[jatom + 1]
                block = np.ascontiguousarray(hermitian[i0:i1, j0:j1], dtype="<c16")
                stream.write(block.tobytes(order="C"))
            stream.flush()
            os.fsync(stream.fileno())
        os.link(temporary, output)
        temporary.unlink()
    except Exception:
        temporary.unlink(missing_ok=True)
        raise

    restored = read_coulomb_v1(output)
    denominator = float(np.linalg.norm(hermitian))
    round_trip = float(np.linalg.norm(restored.matrix - hermitian)) / denominator if denominator else 0.0
    return {
        "path": str(output),
        "sha256": _sha256(output),
        "matrix_kind": "finite_part_coulomb",
        "projection": "(V + V^H) / 2",
        "raw": _matrix_summary(matrix),
        "packaged": _matrix_summary(restored.matrix),
        "round_trip_relative_frobenius": round_trip,
    }


def sternheimer_inverse_sqrt(coulomb, threshold):
    coulomb = np.asarray(coulomb, dtype=np.complex128)
    _require(coulomb.ndim == 2 and coulomb.shape[0] == coulomb.shape[1],
             "Coulomb matrix must be square")
    coulomb = 0.5 * (coulomb + coulomb.conj().T)
    eigenvalues, eigenvectors = np.linalg.eigh(coulomb)
    active = eigenvalues > threshold
    _require(np.any(active), "Coulomb active subspace is empty")
    vectors = eigenvectors[:, active]
    return (vectors / np.sqrt(eigenvalues[active])) @ vectors.conj().T


def sternheimer_pi_from_inverse_sqrt(inverse_sqrt, response):
    inverse_sqrt = np.asarray(inverse_sqrt, dtype=np.complex128)
    response = np.asarray(response, dtype=np.complex128)
    _require(inverse_sqrt.shape == response.shape and response.ndim == 2
             and response.shape[0] == response.shape[1], "Coulomb/response dimension mismatch")
    response = 0.5 * (response + response.conj().T)
    pi = inverse_sqrt @ response @ inverse_sqrt
    return 0.5 * (pi + pi.conj().T)


def trace_log_from_pi(pi):
    sign, logabs = np.linalg.slogdet(np.eye(pi.shape[0], dtype=np.complex128) - pi)
    _require(abs(sign) > 0.0 and math.isfinite(logabs), "non-finite RPA log determinant")
    return complex(np.log(sign) + logabs + np.trace(pi))


def sternheimer_pi(coulomb, response, threshold):
    return sternheimer_pi_from_inverse_sqrt(sternheimer_inverse_sqrt(coulomb, threshold), response)


def sternheimer_trace_log(coulomb, response, threshold):
    return trace_log_from_pi(sternheimer_pi(coulomb, response, threshold))


def package_summary(summary_path, coulomb_path, output_dir, threshold=1.0e-10,
                    energy_gate_ev=1.0e-6):
    summary_path = Path(summary_path).resolve()
    coulomb_path = Path(coulomb_path).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=False)
    summary = json.loads(summary_path.read_text())
    _require(summary.get("merge_complete") is True and summary.get("physical_result") is False,
             "input summary is not a complete diagnostic merge")
    coulomb = read_coulomb_v1(coulomb_path)
    _require(coulomb.ifrequency == 0, "expected a Coulomb reader-v1 file")
    inverse_sqrt = sternheimer_inverse_sqrt(coulomb.matrix, threshold)

    reports = []
    energy_change_hartree = 0.0
    try:
        for item in summary.get("outputs", []):
            matrix_path = summary_path.parent / item["file"]
            _require(_sha256(matrix_path) == item["sha256"], f"matrix SHA256 mismatch: {matrix_path}")
            matrix = np.load(matrix_path, allow_pickle=False)
            output = output_dir / (f"v1_sternheimer_chi0_iq_{coulomb.iq}_ifreq_"
                                   f"{item['ifrequency']}_rank0.dat")
            report = package_matrix(
                matrix,
                output,
                iq=coulomb.iq,
                ifrequency=int(item["ifrequency"]),
                omega=float(item["omega_Ha"]),
                weight=float(item["weight_Ha"]),
                atom_naux=coulomb.atom_naux,
            )
            restored = read_reader_v1(output)
            raw_integrand = trace_log_from_pi(sternheimer_pi_from_inverse_sqrt(inverse_sqrt, matrix))
            packaged_integrand = trace_log_from_pi(
                sternheimer_pi_from_inverse_sqrt(inverse_sqrt, restored.matrix))
            change = (packaged_integrand - raw_integrand) * restored.weight / (2.0 * math.pi)
            energy_change_hartree += change.real
            report.update({
                "ifrequency": restored.ifrequency,
                "omega_Ha": restored.omega,
                "weight_Ha": restored.weight,
                "raw_integrand": [raw_integrand.real, raw_integrand.imag],
                "packaged_integrand": [packaged_integrand.real, packaged_integrand.imag],
                "weighted_energy_change_eV": change.real * HARTREE_TO_EV,
            })
            reports.append(report)
        _require(len(reports) == len(summary["outputs"]) and len(reports) > 0,
                 "response frequency coverage is empty or incomplete")
        total_change_ev = energy_change_hartree * HARTREE_TO_EV
        _require(abs(total_change_ev) <= energy_gate_ev,
                 f"Hermitian packaging energy change {total_change_ev:.16e} eV exceeds gate")
        final = {
            "schema_version": 1,
            "packaging_complete": True,
            "physical_result": False,
            "matrix_kind": "potential_potential_response",
            "mathematical_object": "M = V chi0 V",
            "librpa_transform": "Pi = V^(-1/2) M V^(-1/2)",
            "source_summary": str(summary_path),
            "source_summary_sha256": _sha256(summary_path),
            "coulomb_file": str(coulomb_path),
            "coulomb_sha256": _sha256(coulomb_path),
            "implementation_sha256": _sha256(Path(__file__).resolve()),
            "iq": coulomb.iq,
            "atom_naux": list(coulomb.atom_naux),
            "frequency_count": len(reports),
            "energy_projection_change_eV": total_change_ev,
            "energy_projection_gate_eV": energy_gate_ev,
            "outputs": reports,
        }
        report_path = output_dir / "PACKAGING_COMPLETE.json"
        report_path.write_text(json.dumps(final, indent=2, allow_nan=False) + "\n")
        return final
    except Exception as exc:
        failure = {"packaging_complete": False, "error": str(exc)}
        (output_dir / "FAILURE.json").write_text(json.dumps(failure, indent=2) + "\n")
        raise


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary")
    parser.add_argument("--coulomb")
    parser.add_argument("--unit-audit", action="append", default=[])
    parser.add_argument("--output", required=True)
    parser.add_argument("--sqrt-coulomb-threshold", type=float, default=1.0e-10)
    parser.add_argument("--energy-gate-ev", type=float, default=1.0e-6)
    parser.add_argument("--qweight", type=float)
    args = parser.parse_args(argv)
    if args.unit_audit:
        _require(not args.summary and not args.coulomb,
                 "--unit-audit mode cannot be combined with --summary or --coulomb")
        _require(args.qweight is not None and math.isfinite(args.qweight) and args.qweight > 0,
                 "--unit-audit mode requires a positive finite --qweight")
        merge_weak_q_unit_columns(args.unit_audit, args.output, qweight=args.qweight)
    else:
        _require(args.summary and args.coulomb and args.qweight is None,
                 "summary mode requires --summary and --coulomb and forbids --qweight")
        package_summary(args.summary, args.coulomb, args.output,
                        threshold=args.sqrt_coulomb_threshold, energy_gate_ev=args.energy_gate_ev)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
