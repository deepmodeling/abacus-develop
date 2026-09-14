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


def package_matrix(matrix, output, *, iq, ifrequency, omega, weight, atom_naux):
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

    hermitian = 0.5 * (matrix + matrix.conj().T)
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
        "projection": "(M + M^H) / 2",
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
    parser.add_argument("--summary", required=True)
    parser.add_argument("--coulomb", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--sqrt-coulomb-threshold", type=float, default=1.0e-10)
    parser.add_argument("--energy-gate-ev", type=float, default=1.0e-6)
    args = parser.parse_args(argv)
    package_summary(args.summary, args.coulomb, args.output,
                    threshold=args.sqrt_coulomb_threshold, energy_gate_ev=args.energy_gate_ev)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
