#!/usr/bin/env python3
"""Export isolated-Gamma AO Coulomb potentials; Python 3.6+ and NumPy only."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import shlex
import struct
import sys
import tempfile

import numpy as np


COULOMB_MARKER = -20129433
CS_MARKER = -10267453
HERMITIAN_TOLERANCE = 1.e-10


def positive_counts(values, label):
    if not values or any(not isinstance(n, (int, np.integer)) or n <= 0 for n in values):
        raise ValueError(label + " must contain positive integers")
    return [int(n) for n in values]


def unpack(data, fmt, offset, label):
    if offset < 0 or offset + struct.calcsize(fmt) > len(data):
        raise ValueError(label + ": truncated header or record")
    return struct.unpack_from(fmt, data, offset)


def check_hermitian(values, label):
    if not np.all(np.isfinite(values)):
        raise ValueError(label + ": non-finite values")
    with np.errstate(over="ignore", invalid="ignore"):
        norm = float(np.linalg.norm(values))
        difference = float(np.linalg.norm(values - values.conj().swapaxes(-1, -2)))
    if not math.isfinite(norm) or norm == 0.0 or not math.isfinite(difference):
        raise ValueError(label + ": zero or non-finite Frobenius norm")
    if difference > HERMITIAN_TOLERANCE * norm:
        raise ValueError(label + " is not Hermitian")
    return difference / norm


def read_coulomb_v1(path):
    """Read complete serial upper atom-pair blocks; never repair diagonal blocks."""
    data = Path(path).read_bytes()
    marker, iq, naux, flag, natom, nblocks = unpack(data, "<6i", 0, "Coulomb")
    if marker != COULOMB_MARKER or iq != 1 or flag not in (0, 1):
        raise ValueError("Coulomb requires v1 marker, iq=1 and real/complex payload flag")
    if natom <= 0 or naux <= 0 or nblocks != natom * (natom + 1) // 2:
        raise ValueError("Coulomb requires complete serial atom-pair coverage and positive dimensions")
    header_end = 24 + 4 * natom + 12 * nblocks
    if header_end > len(data):
        raise ValueError("Coulomb: truncated metadata")
    counts = positive_counts(list(unpack(data, "<%di" % natom, 24, "Coulomb")), "Coulomb atom_naux")
    if sum(counts) != naux:
        raise ValueError("Coulomb total auxiliary dimension differs from atom counts")
    pairs = [(i, j) for i in range(natom) for j in range(i, natom)]
    records, seen, expected_offset = [], set(), header_end
    dtype = np.dtype("<c16" if flag else "<f8")
    for record in range(nblocks):
        pair, offset = unpack(data, "<iq", 24 + 4*natom + 12*record, "Coulomb")
        if pair < 0 or pair >= nblocks:
            raise ValueError("Coulomb atom-pair index out of range")
        if pair in seen:
            raise ValueError("Coulomb duplicate atom-pair record")
        seen.add(pair)
        i, j = pairs[pair]
        size = counts[i] * counts[j]
        if offset != expected_offset or offset + size*dtype.itemsize > len(data):
            raise ValueError("Coulomb offset or block size mismatch")
        records.append((i, j, offset, size))
        expected_offset += size*dtype.itemsize
    if expected_offset != len(data):
        raise ValueError("Coulomb trailing data or payload length mismatch")
    starts = np.cumsum([0] + counts).tolist()
    matrix = np.zeros((naux, naux), dtype=np.complex128)
    for i, j, offset, size in records:
        block = np.frombuffer(data, dtype=dtype, count=size, offset=offset).reshape(counts[i], counts[j])
        matrix[starts[i]:starts[i+1], starts[j]:starts[j+1]] = block
        if i != j:
            matrix[starts[j]:starts[j+1], starts[i]:starts[i+1]] = block.conj().T
    error = check_hermitian(matrix, "Coulomb V")
    return {"matrix": matrix, "atom_naux": counts, "natom": natom, "naux": naux,
            "iq": iq, "value_flag": flag, "nblocks": nblocks, "hermitian_relative_error": error}


def read_cs_v1(path, atom_nao, atom_naux):
    """Reconstruct C[p,q,nu] with both AO orientations, including onsite records."""
    atom_nao = positive_counts(atom_nao, "atom_nao")
    atom_naux = positive_counts(atom_naux, "atom_naux")
    data = Path(path).read_bytes()
    marker, natom, ncell, nrecords, nblocks = unpack(data, "<iiiqq", 0, "Cs")
    if marker != CS_MARKER or natom != len(atom_nao) or natom != len(atom_naux) or ncell != 0:
        raise ValueError("Cs marker, atom layout or ncell=0 contract differs")
    if nrecords != natom*natom or nblocks != nrecords:
        raise ValueError("Cs requires complete serial directed atom-pair coverage")
    header_end = 28 + 36*nrecords
    if header_end > len(data):
        raise ValueError("Cs: truncated metadata")
    records, seen, expected_offset = [], set(), header_end
    for record in range(nrecords):
        i, j, r0, r1, r2, max_abs, offset = unpack(data, "<iiiiidq", 28 + 36*record, "Cs")
        if not 1 <= i <= natom or not 1 <= j <= natom:
            raise ValueError("Cs atom index out of range")
        if (r0, r1, r2) != (0, 0, 0):
            raise ValueError("Cs molecular export requires zero translations")
        i, j = i-1, j-1
        if (i, j) in seen:
            raise ValueError("Cs duplicate directed atom-pair record")
        seen.add((i, j))
        shape = (atom_nao[i], atom_nao[j], atom_naux[i])
        size = shape[0]*shape[1]*shape[2]
        if offset != expected_offset or offset + 8*size > len(data):
            raise ValueError("Cs offset or block size mismatch")
        block = np.frombuffer(data, dtype="<f8", count=size, offset=offset).reshape(shape)
        if not math.isfinite(max_abs) or max_abs < 0.0 or not np.all(np.isfinite(block)):
            raise ValueError("Cs contains non-finite values or invalid max_abs")
        actual_max = float(np.max(np.abs(block)))
        if abs(actual_max - max_abs) > 1.e-12 * max(1.0, actual_max):
            raise ValueError("Cs record max_abs differs from payload")
        records.append((i, j, block))
        expected_offset += 8*size
    if expected_offset != len(data):
        raise ValueError("Cs trailing data or payload length mismatch")
    ao = np.cumsum([0] + atom_nao).tolist()
    aux = np.cumsum([0] + atom_naux).tolist()
    coefficients = np.zeros((ao[-1], ao[-1], aux[-1]))
    for i, j, block in records:
        coefficients[ao[i]:ao[i+1], ao[j]:ao[j+1], aux[i]:aux[i+1]] += block
        coefficients[ao[j]:ao[j+1], ao[i]:ao[i+1], aux[i]:aux[i+1]] += block.transpose(1, 0, 2)
    if not np.all(np.isfinite(coefficients)):
        raise ValueError("Cs reconstruction produced non-finite values")
    return coefficients


def contract(coefficients, coulomb):
    if coefficients.ndim != 3 or coefficients.shape[0] != coefficients.shape[1]:
        raise ValueError("AO coefficient dimensions differ")
    naux = coefficients.shape[2]
    if coulomb.shape != (naux, naux) or not np.all(np.isfinite(coefficients)):
        raise ValueError("Auxiliary dimensions differ or Cs is non-finite")
    check_hermitian(coulomb, "Coulomb V")
    potentials = np.einsum("pqv,vm->mpq", coefficients, coulomb, optimize=True)
    check_hermitian(potentials, "AO potentials")
    return potentials


def fields(line):
    tokens = shlex.split(line, comments=True)
    for index, token in enumerate(tokens):
        if token.startswith("//"):
            return tokens[:index]
    return tokens


def read_basis_counts(path, atom_types):
    """Read ABACUS split basis_wfc_out/basis_aux_out, checking all shell counts."""
    atom_types = positive_counts(atom_types, "atom_types (one-based)")
    tokens = [word for line in Path(path).read_text(encoding="ascii").splitlines() for word in fields(line)]
    if len(tokens) < 3 or tokens[2] != "abacus":
        raise ValueError("Basis metadata requires 'ntype total abacus' header")
    ntypes, total = int(tokens[0]), int(tokens[1])
    if ntypes <= 0 or total <= 0 or 3 + 4*ntypes > len(tokens):
        raise ValueError("Basis metadata dimensions or length invalid")
    position, sizes, shells_seen = 3, {}, set()
    for _ in range(ntypes):
        type_id, size = int(tokens[position]), int(tokens[position+1])
        position += 2
        if not 1 <= type_id <= ntypes or type_id in sizes or size <= 0:
            raise ValueError("Basis type size record invalid or duplicate")
        sizes[type_id] = size
    for _ in range(ntypes):
        if position + 2 > len(tokens):
            raise ValueError("Basis metadata truncated shell header")
        type_id, nshell = int(tokens[position]), int(tokens[position+1])
        position += 2
        if type_id not in sizes or type_id in shells_seen or nshell <= 0 or position+nshell > len(tokens):
            raise ValueError("Basis shell record invalid, duplicate or truncated")
        shells_seen.add(type_id)
        angular = [int(value) for value in tokens[position:position+nshell]]
        position += nshell
        if any(l < 0 for l in angular) or sum(2*l+1 for l in angular) != sizes[type_id]:
            raise ValueError("Basis angular shell layout differs from type size")
    if position != len(tokens) or any(t not in sizes for t in atom_types):
        raise ValueError("Basis trailing data or unknown atom type")
    counts = [sizes[t] for t in atom_types]
    if sum(counts) != total:
        raise ValueError("Basis total differs from requested atom-type order")
    return counts


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            digest.update(block)
    return digest.hexdigest()


class InputInventory:
    def __init__(self):
        self.files = {}
        self.missing = set()

    def add(self, path, role, required=True):
        path = Path(path).resolve()
        if not path.is_file():
            if required:
                raise ValueError("Missing input file: " + str(path))
            self.missing.add(str(path))
            return
        key = str(path)
        if key not in self.files:
            self.files[key] = {"path": key, "bytes": path.stat().st_size, "sha256": sha256(path), "roles": []}
        if role not in self.files[key]["roles"]:
            self.files[key]["roles"].append(role)

    def verify_unchanged(self):
        for entry in self.files.values():
            if sha256(entry["path"]) != entry["sha256"]:
                raise ValueError("Input changed during export: " + entry["path"])


def validate_gamma_kpt(path):
    rows = [fields(line) for line in Path(path).read_text(encoding="ascii").splitlines()]
    rows = [row for row in rows if row]
    if len(rows) != 4 or rows[0] not in (["K_POINTS"], ["KPOINTS"]) or len(rows[1]) != 1 or len(rows[2]) != 1:
        raise ValueError("Producer KPT must explicitly describe a single Gamma point")
    count, mode = int(rows[1][0]), rows[2][0].lower()
    if count == 0 and mode in ("gamma", "mp") and len(rows[3]) == 6:
        mesh = [int(value) for value in rows[3][:3]]
        shifts = [float(value) for value in rows[3][3:]]
        if mesh == [1, 1, 1] and shifts == [0.0, 0.0, 0.0]:
            return
    if count == 1 and mode in ("direct", "cartesian") and len(rows[3]) == 4:
        coordinates = [float(value) for value in rows[3]]
        if coordinates[:3] == [0.0, 0.0, 0.0] and math.isfinite(coordinates[3]) and coordinates[3] > 0.0:
            return
    raise ValueError("Producer KPT is not an unshifted single Gamma point")


def collect_provenance(producer, inventory):
    input_path = producer / "INPUT"
    inventory.add(input_path, "producer INPUT", required=False)
    options = {}
    if input_path.is_file():
        for line in input_path.read_text(encoding="utf-8").splitlines():
            words = fields(line)
            if len(words) >= 2:
                options[words[0].lower()] = words[1]
    stru = producer / options.get("stru_file", "STRU")
    inventory.add(stru, "producer STRU", required=False)
    kpt = producer / options.get("kpoint_file", "KPT")
    inventory.add(kpt, "producer KPT", required=False)
    if kpt.is_file():
        validate_gamma_kpt(kpt)
    directories = {"ATOMIC_SPECIES": producer / options.get("pseudo_dir", "."),
                   "NUMERICAL_ORBITAL": producer / options.get("orbital_dir", "."),
                   "ABFS_ORBITAL": producer}
    if stru.is_file():
        section = None
        for line in stru.read_text(encoding="utf-8").splitlines():
            words = fields(line)
            if not words:
                continue
            if len(words) == 1 and re.match(r"^[A-Z_]+$", words[0]):
                section = words[0]
                continue
            if section in directories:
                index = 2 if section == "ATOMIC_SPECIES" else 0
                if len(words) <= index:
                    raise ValueError("Incomplete producer STRU provenance section: " + section)
                inventory.add(directories[section] / words[index], "STRU " + section, required=False)
    suffixes = {".upf", ".vwr", ".orb", ".abfs", ".psp", ".psp8", ".gth"}
    for path in producer.rglob("*"):
        if path.is_file() and (path.suffix.lower() in suffixes or path.name.startswith("basis_")
                               or path.name.startswith("INPUT_")):
            inventory.add(path, "producer basis or input identity")
    return kpt.is_file()


def check_serial_directory(producer):
    for pattern, allowed in (("v1_coulomb_full_iq_*_rank*.dat", "v1_coulomb_full_iq_1_rank0.dat"),
                              ("v1_Cs_data_*.txt", "v1_Cs_data_0.txt")):
        for path in producer.glob(pattern):
            if path.name != allowed:
                raise ValueError("Unexpected q/rank file in serial molecular producer: " + str(path))


def write_outputs(output, potentials, metadata, inventory):
    manifest = Path(str(output) + ".json")
    for path in (output, manifest):
        if os.path.lexists(str(path)):
            raise FileExistsError("Refusing to overwrite output: " + str(path))
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary, published = [], []
    try:
        fd, name = tempfile.mkstemp(prefix=".ao-export-", dir=str(output.parent))
        temporary.append(Path(name))
        with os.fdopen(fd, "w", encoding="ascii") as stream:
            naux, nao, _ = potentials.shape
            stream.write("ABACUS_STERNHEIMER_AO_POTENTIALS 1 {} {} Gamma_Hartree\n".format(nao, naux))
            # Stream one auxiliary matrix at a time; avoid another full tensor-sized copy.
            for matrix in potentials:
                values = matrix.reshape(-1)
                np.savetxt(stream, np.column_stack((values.real, values.imag)), fmt="%.17e")
        metadata["output"] = {"path": str(output), "sha256": sha256(temporary[0]),
                              "bytes": temporary[0].stat().st_size, "shape": list(potentials.shape)}
        fd, name = tempfile.mkstemp(prefix=".ao-export-", dir=str(output.parent))
        temporary.append(Path(name))
        with os.fdopen(fd, "w", encoding="ascii") as stream:
            json.dump(metadata, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
        inventory.verify_unchanged()
        # Hard links publish complete files without overwriting an existing name, even in a race.
        for source, destination in zip(temporary, (output, manifest)):
            os.link(str(source), str(destination))
            published.append((source, destination))
    except BaseException:
        for source, destination in published:
            if destination.exists() and os.path.samefile(str(source), str(destination)):
                destination.unlink()
        raise
    finally:
        for path in temporary:
            path.unlink()
    return manifest


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--producer", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--atom-nao", nargs="+", type=int)
    parser.add_argument("--atom-naux", nargs="+", type=int)
    parser.add_argument("--atom-types", nargs="+", type=int, help="one-based type index for each producer atom")
    parser.add_argument("--basis-wfc", type=Path, help="ABACUS split basis_wfc_out")
    parser.add_argument("--basis-aux", type=Path, help="ABACUS split basis_aux_out")
    parser.add_argument("--provenance-file", type=Path, action="append", default=[], help="additional identity input; repeatable")
    parser.add_argument("--confirm-serial-isolated-gamma", action="store_true", required=True,
                        help="attest to serial, isolated free-space, real-AO Gamma producer identity")
    args = parser.parse_args(argv)
    if not args.atom_nao and not args.basis_wfc:
        parser.error("provide --atom-nao or --basis-wfc")
    if (args.basis_wfc or args.basis_aux) and not args.atom_types:
        parser.error("basis files require explicit --atom-types in producer atom order")
    return args


def run(args):
    producer = args.producer.resolve()
    if not producer.is_dir():
        raise ValueError("Producer directory does not exist")
    # Keep the final component unresolved so a dangling output symlink is also refused.
    output = Path(os.path.abspath(str(args.output)))
    for path in (output, Path(str(output) + ".json")):
        if os.path.lexists(str(path)):
            raise FileExistsError("Refusing to overwrite output: " + str(path))
    check_serial_directory(producer)
    inventory = InputInventory()
    v_path, cs_path = producer / "v1_coulomb_full_iq_1_rank0.dat", producer / "v1_Cs_data_0.txt"
    inventory.add(v_path, "full Coulomb V")
    inventory.add(cs_path, "localized RI Cs")
    inventory.add(Path(__file__), "exporter source")
    for path in (args.basis_wfc, args.basis_aux):
        if path is not None:
            inventory.add(path, "basis dimension metadata")
    for path in args.provenance_file:
        inventory.add(path, "explicit provenance")
    kpt_checked = collect_provenance(producer, inventory)
    coulomb = read_coulomb_v1(v_path)
    atom_nao = positive_counts(args.atom_nao, "atom_nao") if args.atom_nao else None
    if args.basis_wfc:
        basis_nao = read_basis_counts(args.basis_wfc, args.atom_types)
        if atom_nao is not None and atom_nao != basis_nao:
            raise ValueError("AO count arguments differ from basis file")
        atom_nao = basis_nao
    atom_naux = coulomb["atom_naux"]
    if len(atom_nao) != len(atom_naux):
        raise ValueError("AO/aux atom counts differ")
    if args.atom_types and len(args.atom_types) != len(atom_nao):
        raise ValueError("Atom type and AO metadata lengths differ")
    if args.atom_types:
        positive_counts(args.atom_types, "atom_types (one-based)")
    if args.atom_naux and positive_counts(args.atom_naux, "atom_naux") != atom_naux:
        raise ValueError("Auxiliary count arguments differ from Coulomb header")
    if args.basis_aux and read_basis_counts(args.basis_aux, args.atom_types) != atom_naux:
        raise ValueError("Auxiliary basis file differs from Coulomb header")
    coefficients = read_cs_v1(cs_path, atom_nao, atom_naux)
    potentials = contract(coefficients, coulomb["matrix"])
    metadata = {
        "schema_version": 1, "status": "exported_not_physically_validated",
        "external_identity_verified": False,
        "producer_kpt_gamma_checked": kpt_checked,
        "user_confirmed_serial_isolated_gamma": bool(args.confirm_serial_isolated_gamma),
        "producer_directory": str(producer), "python_version": platform.python_version(),
        "numpy_version": np.__version__, "argv": sys.argv[1:],
        "atom_nao": atom_nao, "atom_naux": atom_naux, "atom_types": args.atom_types,
        "formula": "T[mu,p,q] = sum_nu C[p,q,nu] * V[nu,mu]",
        "cs_convention": "each directed record plus AO-transposed record, including onsite",
        "units": "Hartree", "layout": "auxiliary-major, row-major AO; real imaginary",
        "hermitian_tolerance": HERMITIAN_TOLERANCE,
        "coulomb_hermitian_relative_error": coulomb["hermitian_relative_error"],
        "potentials_hermitian_relative_error": check_hermitian(potentials, "AO potentials"),
        "checks": ["iq=1", "complete unique serial atom pairs", "Cs ncell=0 and all R=0",
                   "exact contiguous offsets and byte counts", "finite Hermitian V and AO potentials"],
        "assumptions": [
            "User confirms serial isolated free-space Gamma producer with real AO/auxiliary functions.",
            "Producer full v1 Coulomb payload is in Hartree; its binary header has no unit tag.",
            "q index 1 alone does not encode or prove the physical Gamma coordinate or Coulomb kernel.",
            "V and Cs and declared atom/basis order refer to the same unchanged producer.",
            "No auxiliary rotation, shrink, whitening, cell or basis change is permitted at consumption.",
            "PP, orbitals, ABFS, geometry, kernel and executable identity must be checked externally.",
            "Hashes establish file identity only, not physical convergence or response accuracy."],
        "inputs": [inventory.files[key] for key in sorted(inventory.files)],
        "missing_provenance": sorted(inventory.missing)}
    manifest = write_outputs(output, potentials, metadata, inventory)
    return {"status": metadata["status"], "output": str(output), "manifest": str(manifest),
            "shape": list(potentials.shape)}


def main():
    args = parse_args()
    try:
        print(json.dumps(run(args), sort_keys=True))
    except (ValueError, OSError, MemoryError) as error:
        print("error: " + str(error), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
