#!/usr/bin/env python3
"""Validate the ABACUS producer-side file contract consumed by LibRPA.

The checker deliberately separates stable producer data from gauge-dependent
wavefunction data.  It validates names, headers, dimensions, finite values,
and optional tolerant references; it does not require KS eigenvector or
velocity-matrix values to be byte-identical.
"""

import argparse
import glob
import json
import math
import re
import struct
from pathlib import Path


class ProducerContractError(ValueError):
    """Raised when a producer output violates the declared contract."""


_COULOMB_MARKER = -20129433
_CHI0_MARKER = -41073291
_LRICOEF_MARKER = -10267453
_SHRINK_SINVS_MARKER = -30241621


def _matches(root, pattern):
    return sorted(Path(path) for path in glob.glob(str(root / pattern), recursive=True))


def _number(token):
    try:
        value = float(token)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(value):
        raise ProducerContractError("non-finite numeric value: {}".format(token))
    return value


def _required_number(path, token, context):
    value = _number(token)
    if value is None:
        raise ProducerContractError("{} has a non-numeric {}: {}".format(path, context, token))
    return value


def _numeric_tokens(path):
    values = []
    for line_number, line in enumerate(path.read_text().splitlines(), 1):
        for token in line.split():
            value = _number(token)
            if value is not None:
                values.append((line_number, token, value))
    if not values:
        raise ProducerContractError("{} contains no numeric values".format(path))
    return values


def _compare_text(actual, reference, tolerance):
    actual_lines = actual.read_text().splitlines()
    reference_lines = reference.read_text().splitlines()
    if len(actual_lines) != len(reference_lines):
        raise ProducerContractError(
            "{} has {} lines; reference has {}".format(actual, len(actual_lines), len(reference_lines))
        )
    for line_number, (aline, rline) in enumerate(zip(actual_lines, reference_lines), 1):
        atokens = aline.split()
        rtokens = rline.split()
        if len(atokens) != len(rtokens):
            raise ProducerContractError("{} line {} token count differs".format(actual, line_number))
        for column, (atoken, rtoken) in enumerate(zip(atokens, rtokens), 1):
            avalue = _number(atoken)
            rvalue = _number(rtoken)
            if avalue is None or rvalue is None:
                if atoken != rtoken:
                    raise ProducerContractError("{} line {} column {} differs".format(actual, line_number, column))
            elif abs(avalue - rvalue) > tolerance:
                raise ProducerContractError(
                    "{} line {} column {} differs by {} (tol {})".format(
                        actual, line_number, column, abs(avalue - rvalue), tolerance
                    )
                )


def _check_stru(path, minimum_symmetry_rows, minimum_spin_symmetry_rows=None):
    lines = path.read_text().splitlines()
    if len(lines) < 8:
        raise ProducerContractError("{} is too short for stru_out".format(path))
    for line in lines[:6]:
        if len(line.split()) != 3:
            raise ProducerContractError("{} has an invalid lattice block".format(path))
        for token in line.split():
            _number(token)
    try:
        natom = int(lines[6].split()[0])
    except (IndexError, ValueError):
        raise ProducerContractError("{} has an invalid atom count".format(path))
    if natom <= 0 or len(lines) < 7 + natom:
        raise ProducerContractError("{} has an invalid atom block".format(path))
    for line in lines[7 : 7 + natom]:
        if len(line.split()) != 4:
            raise ProducerContractError("{} has an invalid atom row".format(path))
        for token in line.split()[:3]:
            _number(token)
        try:
            int(line.split()[3])
        except ValueError:
            raise ProducerContractError("{} has an invalid atom type".format(path))
    row_match = None
    row_index = None
    for index, line in enumerate(lines[7 + natom :], 7 + natom):
        match = re.match(r"^\s*(\d+)\s+row\s*$", line)
        if match:
            row_match = int(match.group(1))
            row_index = index
            break
    if minimum_symmetry_rows is not None:
        if row_match is None or row_match < int(minimum_symmetry_rows):
            raise ProducerContractError("{} has too few symmetry rows".format(path))
    if row_match is not None:
        if len(lines) < row_index + 1 + row_match:
            raise ProducerContractError("{} symmetry block is truncated".format(path))
        for line in lines[row_index + 1 : row_index + 1 + row_match]:
            tokens = line.split()
            if len(tokens) != 12:
                raise ProducerContractError("{} has an invalid symmetry row".format(path))
            for token in tokens[:9]:
                try:
                    int(token)
                except ValueError:
                    raise ProducerContractError("{} has a non-integer symmetry rotation".format(path))
            for token in tokens[9:]:
                _number(token)
    if minimum_spin_symmetry_rows is not None:
        if row_match is None:
            raise ProducerContractError("{} has no spatial symmetry block for spin symmetry".format(path))
        spin_header_index = row_index + 1 + row_match
        if spin_header_index >= len(lines):
            raise ProducerContractError("{} has no spin_symmetry block".format(path))
        header = re.match(r"^\s*spin_symmetry\s+([01])\s+([01])\s*$", lines[spin_header_index])
        if header is None or header.group(2) != "1":
            raise ProducerContractError("{} has an invalid spin_symmetry header".format(path))
        spin_rows = [line for line in lines[spin_header_index + 1 :] if line.strip()]
        if len(spin_rows) < int(minimum_spin_symmetry_rows) or len(spin_rows) != row_match:
            raise ProducerContractError("{} has an incomplete spin_symmetry block".format(path))
        for line in spin_rows:
            tokens = line.split()
            if len(tokens) != 9:
                raise ProducerContractError("{} has an invalid spin symmetry row".format(path))
            try:
                antiunitary = int(tokens[0])
            except ValueError:
                raise ProducerContractError("{} has an invalid spin antiunitary flag".format(path))
            if antiunitary not in (0, 1):
                raise ProducerContractError("{} has an invalid spin antiunitary flag".format(path))
            for token in tokens[1:]:
                _number(token)


def _check_wfc_nao(path):
    """Validate the text LCAO wavefunction format consumed by GW preprocessing."""
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if len(lines) < 6:
        raise ProducerContractError("{} is too short for an LCAO wavefunction".format(path))

    complex_layout = "(index of k points)" in lines[0]
    header_index = 0
    coefficient_width = 1
    if complex_layout:
        try:
            if int(lines[0].split()[0]) <= 0:
                raise ValueError
        except (IndexError, ValueError):
            raise ProducerContractError("{} has an invalid k-point index".format(path))
        coordinates = lines[1].split()
        if len(coordinates) != 3:
            raise ProducerContractError("{} has an invalid k-point coordinate".format(path))
        for coordinate in coordinates:
            _number(coordinate)
        header_index = 2
        coefficient_width = 2

    try:
        nbands = int(lines[header_index].split()[0])
        nlocal = int(lines[header_index + 1].split()[0])
    except (IndexError, ValueError):
        raise ProducerContractError("{} has an invalid wavefunction dimension header".format(path))
    if nbands <= 0 or nlocal <= 0:
        raise ProducerContractError("{} has non-positive wavefunction dimensions".format(path))

    first_band_index = header_index + 2
    band_indices = [
        index for index, line in enumerate(lines[first_band_index:], first_band_index) if line.endswith("(band)")
    ]
    if len(band_indices) != nbands:
        raise ProducerContractError("{} has {} bands; expected {}".format(path, len(band_indices), nbands))

    for iband, band_index in enumerate(band_indices, 1):
        try:
            band_number = int(lines[band_index].split()[0])
        except (IndexError, ValueError):
            raise ProducerContractError("{} has an invalid band index".format(path))
        if band_number != iband:
            raise ProducerContractError("{} has a non-contiguous band index".format(path))
        if band_index + 2 >= len(lines):
            raise ProducerContractError("{} has an incomplete band record".format(path))
        try:
            energy = lines[band_index + 1].split()[0]
            occupation = lines[band_index + 2].split()[0]
        except IndexError:
            raise ProducerContractError("{} has an incomplete band record".format(path))
        _required_number(path, energy, "band energy")
        _required_number(path, occupation, "band occupation")
        next_band_index = band_indices[iband] if iband < nbands else len(lines)
        coefficients = []
        for line in lines[band_index + 3 : next_band_index]:
            for token in line.split():
                coefficients.append(_required_number(path, token, "wavefunction coefficient"))
        expected_coefficients = coefficient_width * nlocal
        if len(coefficients) != expected_coefficients:
            raise ProducerContractError(
                "{} band {} has {} coefficients; expected {}".format(
                    path, iband, len(coefficients), expected_coefficients
                )
            )


def _upper_pair(pair_index, natom):
    current = 0
    for iatom in range(natom):
        for jatom in range(iatom, natom):
            if current == pair_index:
                return iatom, jatom
            current += 1
    raise ProducerContractError("invalid upper-triangular atom pair index {}".format(pair_index))


def _check_v1_binary(path, kind):
    data = path.read_bytes()
    if kind == "coulomb_v1":
        if len(data) < 24:
            raise ProducerContractError("{} is too short for a Coulomb v1 header".format(path))
        marker, iq, naux, value_flag, natom, nblocks = struct.unpack_from("<6i", data, 0)
        header_end = 24
        if marker != _COULOMB_MARKER or iq <= 0 or naux <= 0 or natom <= 0 or nblocks < 0:
            raise ProducerContractError("{} has an invalid Coulomb v1 header".format(path))
    else:
        if len(data) < 44:
            raise ProducerContractError("{} is too short for a chi0 v1 header".format(path))
        marker, iq, ifrequency, naux, value_flag, natom = struct.unpack_from("<6i", data, 0)
        omega, weight = struct.unpack_from("<2d", data, 24)
        nblocks = struct.unpack_from("<i", data, 40)[0]
        header_end = 44
        if marker != _CHI0_MARKER or iq <= 0 or ifrequency <= 0 or naux <= 0 or natom <= 0 or nblocks < 0:
            raise ProducerContractError("{} has an invalid chi0 v1 header".format(path))
        if not math.isfinite(omega) or not math.isfinite(weight):
            raise ProducerContractError("{} has non-finite chi0 metadata".format(path))
    if value_flag != 1:
        raise ProducerContractError("{} is not a complex v1 file".format(path))
    atom_bytes = 4 * natom
    block_bytes = 12 * nblocks
    table_end = header_end + atom_bytes + block_bytes
    if len(data) < table_end:
        raise ProducerContractError("{} has a truncated v1 table".format(path))
    atom_naux = struct.unpack_from("<{}i".format(natom), data, header_end)
    if any(value <= 0 for value in atom_naux):
        raise ProducerContractError("{} has invalid atom auxiliary dimensions".format(path))
    table_offset = header_end + atom_bytes
    blocks = []
    seen = set()
    for iblock in range(nblocks):
        pair_index, offset = struct.unpack_from("<iq", data, table_offset + 12 * iblock)
        if pair_index in seen:
            raise ProducerContractError("{} has duplicate v1 atom-pair blocks".format(path))
        seen.add(pair_index)
        iatom, jatom = _upper_pair(pair_index, natom)
        payload_end = offset + 16 * atom_naux[iatom] * atom_naux[jatom]
        if offset < table_end or payload_end > len(data):
            raise ProducerContractError("{} has an invalid v1 payload offset".format(path))
        blocks.append((offset, payload_end))
    if blocks and any(right > next_left for (_, right), (next_left, _) in zip(sorted(blocks), sorted(blocks)[1:])):
        raise ProducerContractError("{} has overlapping v1 payload blocks".format(path))


def _v1_layout(path, kind):
    """Read the stable metadata/table/payload regions of a v1 binary file."""
    data = path.read_bytes()
    if kind == "coulomb_v1":
        if len(data) < 24:
            raise ProducerContractError("{} is too short for a Coulomb v1 header".format(path))
        header = struct.unpack_from("<6i", data, 0)
        marker, iq, naux, value_flag, natom, nblocks = header
        floats = ()
        header_end = 24
    else:
        if len(data) < 44:
            raise ProducerContractError("{} is too short for a chi0 v1 header".format(path))
        header = struct.unpack_from("<6i", data, 0)
        marker, iq, ifrequency, naux, value_flag, natom = header
        floats = struct.unpack_from("<2d", data, 24)
        nblocks = struct.unpack_from("<i", data, 40)[0]
        header_end = 44
    atom_naux = struct.unpack_from("<{}i".format(natom), data, header_end)
    table_offset = header_end + 4 * natom
    table = [struct.unpack_from("<iq", data, table_offset + 12 * index) for index in range(nblocks)]
    payload_start = table_offset + 12 * nblocks
    return data, header, floats, atom_naux, table, payload_start


def _compare_v1_binary(actual, reference, kind, tolerance):
    adata, aheader, afloat, aatom_naux, atable, apayload_start = _v1_layout(actual, kind)
    rdata, rheader, rfloat, ratom_naux, rtable, rpayload_start = _v1_layout(reference, kind)
    if aheader != rheader or aatom_naux != ratom_naux or atable != rtable:
        raise ProducerContractError("{} v1 metadata/table differs from {}".format(actual, reference))
    if len(adata) != len(rdata) or len(afloat) != len(rfloat):
        raise ProducerContractError("{} v1 payload size differs from {}".format(actual, reference))
    for index, (avalue, rvalue) in enumerate(zip(afloat, rfloat)):
        if not math.isfinite(avalue) or not math.isfinite(rvalue) or abs(avalue - rvalue) > tolerance:
            raise ProducerContractError("{} v1 header value {} differs from reference".format(actual, index))
    if (len(adata) - apayload_start) % 8 != 0 or (len(rdata) - rpayload_start) % 8 != 0:
        raise ProducerContractError("{} v1 payload is not a sequence of doubles".format(actual))
    if apayload_start != rpayload_start:
        raise ProducerContractError("{} v1 payload layout differs from {}".format(actual, reference))
    for index in range(apayload_start, len(adata), 8):
        avalue = struct.unpack_from("<d", adata, index)[0]
        rvalue = struct.unpack_from("<d", rdata, index)[0]
        if not math.isfinite(avalue) or not math.isfinite(rvalue) or abs(avalue - rvalue) > tolerance:
            raise ProducerContractError("{} v1 payload value differs at byte {} (tol {})".format(actual, index, tolerance))


def _lri_layout(path, kind):
    data = path.read_bytes()
    if kind == "lri_coeff_v1":
        if len(data) < 28:
            raise ProducerContractError("{} is too short for an LRI coefficient v1 header".format(path))
        marker, natom, ncell = struct.unpack_from("<3i", data, 0)
        nrecords, nrecords_max = struct.unpack_from("<2q", data, 12)
        if marker != _LRICOEF_MARKER or natom <= 0 or ncell < 0 or nrecords < 0 or nrecords > nrecords_max:
            raise ProducerContractError("{} has an invalid LRI coefficient v1 header".format(path))
        record_size = 36
        table_start = 28
        table_end = table_start + record_size * nrecords_max
        records = []
        for index in range(nrecords_max):
            offset = table_start + record_size * index
            ia1, ia2, r0, r1, r2 = struct.unpack_from("<5i", data, offset)
            max_abs = struct.unpack_from("<d", data, offset + 20)[0]
            payload_offset = struct.unpack_from("<q", data, offset + 28)[0]
            records.append((ia1, ia2, r0, r1, r2, max_abs, payload_offset))
        return data, (marker, natom, ncell, nrecords, nrecords_max), records, table_end
    if len(data) < 8:
        raise ProducerContractError("{} is too short for a shrink_sinvS v1 header".format(path))
    marker, nrecords = struct.unpack_from("<2i", data, 0)
    if marker != _SHRINK_SINVS_MARKER or nrecords < 0:
        raise ProducerContractError("{} has an invalid shrink_sinvS v1 header".format(path))
    record_size = 44
    table_start = 8
    table_end = table_start + record_size * nrecords
    records = []
    for index in range(nrecords):
        offset = table_start + record_size * index
        ints = struct.unpack_from("<7i", data, offset)
        qweight = struct.unpack_from("<d", data, offset + 28)[0]
        payload_offset = struct.unpack_from("<q", data, offset + 36)[0]
        records.append((ints, qweight, payload_offset))
    return data, (marker, nrecords), records, table_end


def _compare_lri_binary(actual, reference, kind, tolerance):
    adata, aheader, arecords, atable_end = _lri_layout(actual, kind)
    rdata, rheader, rrecords, rtable_end = _lri_layout(reference, kind)
    if aheader != rheader or len(adata) != len(rdata) or atable_end != rtable_end:
        raise ProducerContractError("{} v1 metadata/table differs from {}".format(actual, reference))
    if kind == "lri_coeff_v1":
        if len(arecords) != len(rrecords):
            raise ProducerContractError("{} v1 record count differs from {}".format(actual, reference))
        for index, (arecord, rrecord) in enumerate(zip(arecords, rrecords)):
            if arecord[:5] != rrecord[:5] or arecord[6] != rrecord[6]:
                raise ProducerContractError("{} v1 record {} differs from {}".format(actual, index, reference))
            if abs(arecord[5] - rrecord[5]) > tolerance:
                raise ProducerContractError("{} v1 record {} differs from reference".format(actual, index))
    else:
        for index, (arecord, rrecord) in enumerate(zip(arecords, rrecords)):
            if arecord[0] != rrecord[0] or arecord[2] != rrecord[2]:
                raise ProducerContractError("{} v1 record {} differs from {}".format(actual, index, reference))
            if abs(arecord[1] - rrecord[1]) > tolerance:
                raise ProducerContractError("{} v1 record {} differs from reference".format(actual, index))
    start = atable_end
    if start != rtable_end:
        raise ProducerContractError("{} v1 payload layout differs from {}".format(actual, reference))
    for index in range(start, len(adata), 8):
        avalue = struct.unpack_from("<d", adata, index)[0]
        rvalue = struct.unpack_from("<d", rdata, index)[0]
        if not math.isfinite(avalue) or not math.isfinite(rvalue) or abs(avalue - rvalue) > tolerance:
            raise ProducerContractError("{} v1 payload value differs at byte {} (tol {})".format(actual, index, tolerance))


def _check_file(path, entry):
    kind = entry.get("kind", "numeric")
    if not path.is_file() or path.stat().st_size == 0:
        raise ProducerContractError("required producer file is missing or empty: {}".format(path))
    if kind == "gauge" or kind == "presence":
        return
    if kind == "stru":
        _check_stru(path, entry.get("symmetry_rows"), entry.get("spin_symmetry_rows"))
    elif kind == "wfc_nao":
        _check_wfc_nao(path)
    elif kind == "band":
        lines = path.read_text().splitlines()
        if len(lines) < 5:
            raise ProducerContractError("{} is too short for band_out".format(path))
        for index in range(4):
            try:
                if int(lines[index].split()[0]) <= 0:
                    raise ValueError
            except (IndexError, ValueError):
                raise ProducerContractError("{} has an invalid band_out header".format(path))
        _number(lines[4].split()[0])
        for line in lines[5:]:
            for token in line.split():
                _number(token)
    elif kind in ("coulomb_v1", "chi0_v1"):
        _check_v1_binary(path, kind)
    elif kind in ("lri_coeff_v1", "shrink_sinvs_v1"):
        _lri_layout(path, kind)
    elif kind in ("numeric", "text"):
        _numeric_tokens(path)
    else:
        raise ProducerContractError("unknown producer contract kind: {}".format(kind))


def check_manifest(root, manifest):
    """Return diagnostics (empty on success) or raise on contract failure."""
    root = Path(root)
    output_root = root / manifest.get("output_dir", ".")
    if not output_root.is_dir():
        raise ProducerContractError("producer output directory is missing: {}".format(output_root))
    reference_root = root / manifest["reference_dir"] if manifest.get("reference_dir") else None
    diagnostics = []
    for entry in manifest.get("required", []):
        entry_root = root / entry.get("base_dir", manifest.get("output_dir", "."))
        matches = _matches(entry_root, entry["pattern"])
        if not matches:
            raise ProducerContractError("no producer files match {}".format(entry["pattern"]))
        expected_count = entry.get("expected_count")
        if expected_count is not None and len(matches) != int(expected_count):
            raise ProducerContractError(
                "{} matches {} files; expected {}".format(entry["pattern"], len(matches), expected_count)
            )
        for path in matches:
            _check_file(path, entry)
            if entry.get("reference"):
                if reference_root is None:
                    raise ProducerContractError("reference requested without reference_dir")
                reference = reference_root / path.relative_to(entry_root)
                if not reference.is_file():
                    raise ProducerContractError("missing producer reference: {}".format(reference))
                tolerance = float(entry.get("abs_tol", 1.0e-8))
                if entry.get("kind", "numeric") in ("coulomb_v1", "chi0_v1"):
                    _compare_v1_binary(path, reference, entry["kind"], tolerance)
                elif entry.get("kind", "numeric") in ("lri_coeff_v1", "shrink_sinvs_v1"):
                    _compare_lri_binary(path, reference, entry["kind"], tolerance)
                else:
                    _compare_text(path, reference, tolerance)
    for entry in manifest.get("optional", []):
        entry_root = root / entry.get("base_dir", manifest.get("output_dir", "."))
        for path in _matches(entry_root, entry["pattern"]):
            _check_file(path, entry)
    return diagnostics


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    args = parser.parse_args(argv)
    manifest = json.loads(args.manifest.read_text())
    try:
        check_manifest(args.root, manifest)
    except ProducerContractError as error:
        parser.error(str(error))
    print("LibRPA producer contract: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
