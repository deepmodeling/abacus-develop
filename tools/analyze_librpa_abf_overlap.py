#!/usr/bin/env python3
"""Offline PSD analysis for ABACUS LibRPA v1 active-ABF diagnostics."""

import argparse
import glob
import json
import os
import struct
import sys

import numpy as np


RAW_MARKER = -40817329
RAW_VERSION = 1
RAW_KIND_ACTIVE = 1
COULOMB_MARKER = -20129433
RAW_HEADER = struct.Struct("<6i d 3d")
COULOMB_HEADER = struct.Struct("<6i")


class AnalysisError(RuntimeError):
    pass


def _read_exact(handle, size, label):
    data = handle.read(size)
    if len(data) != size:
        raise AnalysisError("{}: truncated {}".format(handle.name, label))
    return data


def _read_basis(path):
    try:
        with open(path, "r", encoding="utf-8") as handle:
            rows = [line.split() for line in handle if line.split()]
    except OSError as exc:
        raise AnalysisError("{}: cannot read basis: {}".format(path, exc))
    if not rows or len(rows[0]) != 3:
        raise AnalysisError("{}: invalid split basis header".format(path))
    try:
        ntypes, declared_total, label = int(rows[0][0]), int(rows[0][1]), rows[0][2]
    except ValueError:
        raise AnalysisError("{}: invalid split basis header".format(path))
    if ntypes <= 0 or declared_total <= 0:
        raise AnalysisError("{}: invalid split basis dimensions".format(path))
    cursor = 1
    type_sizes = [None] * ntypes
    for _ in range(ntypes):
        if cursor >= len(rows) or len(rows[cursor]) != 2:
            raise AnalysisError("{}: truncated per-type basis metadata".format(path))
        try:
            itype, size = int(rows[cursor][0]), int(rows[cursor][1])
        except ValueError:
            raise AnalysisError("{}: invalid per-type basis metadata".format(path))
        cursor += 1
        if not 1 <= itype <= ntypes or size <= 0 or type_sizes[itype - 1] is not None:
            raise AnalysisError("{}: invalid or duplicate per-type basis metadata".format(path))
        type_sizes[itype - 1] = size
    if any(size is None for size in type_sizes):
        raise AnalysisError("{}: incomplete per-type basis metadata".format(path))

    seen = [False] * ntypes
    for _ in range(ntypes):
        if cursor >= len(rows) or len(rows[cursor]) != 2:
            raise AnalysisError("{}: truncated shell-layout metadata".format(path))
        try:
            itype, nshell = int(rows[cursor][0]), int(rows[cursor][1])
        except ValueError:
            raise AnalysisError("{}: invalid shell-layout metadata".format(path))
        cursor += 1
        if not 1 <= itype <= ntypes or nshell < 0 or seen[itype - 1]:
            raise AnalysisError("{}: invalid or duplicate shell layout".format(path))
        seen[itype - 1] = True
        shell_size = 0
        for _ in range(nshell):
            if cursor >= len(rows) or len(rows[cursor]) != 1:
                raise AnalysisError("{}: truncated shell layout".format(path))
            try:
                l_value = int(rows[cursor][0])
            except ValueError:
                raise AnalysisError("{}: invalid angular momentum".format(path))
            cursor += 1
            if l_value < 0:
                raise AnalysisError("{}: negative angular momentum".format(path))
            shell_size += 2 * l_value + 1
        if shell_size != type_sizes[itype - 1]:
            raise AnalysisError("{}: shell layout does not match type size".format(path))
    if cursor != len(rows):
        raise AnalysisError("{}: trailing tokens after split basis layout".format(path))
    # The header total is sum(type_size * atom multiplicity). This file does not
    # contain atom-to-type mapping, so it is only a declared cross-check target.
    return {"ntypes": ntypes, "declared_total": declared_total, "label": label,
            "type_sizes": type_sizes}


def _read_raw(path):
    try:
        size = os.path.getsize(path)
        with open(path, "rb") as handle:
            fields = RAW_HEADER.unpack(_read_exact(handle, RAW_HEADER.size, "raw header"))
            marker, version, iq, kind, naux, natom = fields[:6]
            q_weight = fields[6]
            q = fields[7:10]
            if marker != RAW_MARKER or version != RAW_VERSION:
                raise AnalysisError("{}: bad raw marker/version".format(path))
            if iq <= 0 or kind != RAW_KIND_ACTIVE or naux <= 0 or natom <= 0:
                raise AnalysisError("{}: invalid raw metadata".format(path))
            if not np.isfinite(q_weight) or not np.all(np.isfinite(q)):
                raise AnalysisError("{}: non-finite raw q metadata".format(path))
            atom_naux = list(struct.unpack("<{}i".format(natom),
                                            _read_exact(handle, 4 * natom, "raw atom_naux")))
            if any(value <= 0 for value in atom_naux) or sum(atom_naux) != naux:
                raise AnalysisError("{}: inconsistent raw atom_naux".format(path))
            payload_bytes = 16 * naux * naux
            payload = np.frombuffer(_read_exact(handle, payload_bytes, "raw payload"),
                                    dtype="<c16").copy().reshape((naux, naux))
            if handle.tell() != size:
                raise AnalysisError("{}: trailing bytes after raw payload".format(path))
    except OSError as exc:
        raise AnalysisError("{}: cannot read raw overlap: {}".format(path, exc))
    if not np.all(np.isfinite(payload.real)) or not np.all(np.isfinite(payload.imag)):
        raise AnalysisError("{}: non-finite raw payload".format(path))
    hermitian_residual = float(np.max(np.abs(payload - payload.conj().T)))
    if not np.allclose(payload, payload.conj().T, rtol=1e-10, atol=1e-10):
        raise AnalysisError("{}: raw S is not Hermitian".format(path))
    return {"path": path, "iq": iq, "q_weight": q_weight, "q": q,
            "naux": naux, "natom": natom, "atom_naux": atom_naux, "S": payload,
            "hermitian_residual": hermitian_residual}


def _pair_from_index(index, natom):
    for i in range(natom):
        for j in range(i, natom):
            if index == 0:
                return i, j
            index -= 1
    raise AnalysisError("invalid Coulomb pair index")


def _read_coulomb(path):
    try:
        size = os.path.getsize(path)
        with open(path, "rb") as handle:
            marker, iq, naux, value_flag, natom, nblocks = COULOMB_HEADER.unpack(
                _read_exact(handle, COULOMB_HEADER.size, "Coulomb header"))
            if marker != COULOMB_MARKER or iq <= 0 or naux <= 0 or natom <= 0 or nblocks < 0:
                raise AnalysisError("{}: invalid Coulomb header".format(path))
            if value_flag != 1:
                raise AnalysisError("{}: expected complex-double Coulomb payload".format(path))
            atom_naux = list(struct.unpack("<{}i".format(natom),
                                            _read_exact(handle, 4 * natom, "Coulomb atom_naux")))
            if any(value <= 0 for value in atom_naux) or sum(atom_naux) != naux:
                raise AnalysisError("{}: inconsistent Coulomb atom_naux".format(path))
            header_size = COULOMB_HEADER.size + 4 * natom + 12 * nblocks
            records = []
            seen = set()
            ranges = []
            npairs = natom * (natom + 1) // 2
            for _ in range(nblocks):
                pair_index, offset = struct.unpack("<iq", _read_exact(handle, 12, "Coulomb block record"))
                if pair_index in seen or not 0 <= pair_index < npairs or offset < header_size:
                    raise AnalysisError("{}: invalid or duplicate Coulomb block".format(path))
                seen.add(pair_index)
                i, j = _pair_from_index(pair_index, natom)
                nbytes = 16 * atom_naux[i] * atom_naux[j]
                if offset + nbytes > size:
                    raise AnalysisError("{}: Coulomb payload exceeds file".format(path))
                ranges.append((offset, offset + nbytes))
                records.append((pair_index, i, j, offset, nbytes))
            sorted_ranges = sorted(ranges)
            if not sorted_ranges:
                if size != header_size:
                    raise AnalysisError("{}: trailing or unreferenced Coulomb bytes".format(path))
            else:
                if sorted_ranges[0][0] != header_size:
                    raise AnalysisError("{}: unreferenced gap before Coulomb payload".format(path))
                for (_, end), (begin, _) in zip(sorted_ranges, sorted_ranges[1:]):
                    if begin < end:
                        raise AnalysisError("{}: overlapping Coulomb payloads".format(path))
                    if begin != end:
                        raise AnalysisError("{}: unreferenced gap between Coulomb payloads".format(path))
                if sorted_ranges[-1][1] != size:
                    raise AnalysisError("{}: trailing or unreferenced Coulomb bytes".format(path))
            blocks = {}
            for pair_index, i, j, offset, nbytes in records:
                handle.seek(offset)
                blocks[pair_index] = np.frombuffer(
                    _read_exact(handle, nbytes, "Coulomb payload"), dtype="<c16").copy().reshape(
                        (atom_naux[i], atom_naux[j]))
    except OSError as exc:
        raise AnalysisError("{}: cannot read Coulomb shard: {}".format(path, exc))
    return {"path": path, "iq": iq, "naux": naux, "natom": natom,
            "atom_naux": atom_naux, "blocks": blocks}


def _assemble_coulomb(shards, raw):
    expected = raw["natom"] * (raw["natom"] + 1) // 2
    V = np.zeros((raw["naux"], raw["naux"]), dtype=np.complex128)
    shifts = np.cumsum([0] + raw["atom_naux"])
    blocks = {}
    for shard in shards:
        if (shard["iq"], shard["naux"], shard["natom"], shard["atom_naux"]) != (
                raw["iq"], raw["naux"], raw["natom"], raw["atom_naux"]):
            raise AnalysisError("{}: Coulomb/raw metadata mismatch".format(shard["path"]))
        for pair_index, payload in shard["blocks"].items():
            if pair_index in blocks:
                raise AnalysisError("duplicate Coulomb pair {} for iq {}".format(pair_index, raw["iq"]))
            i, j = _pair_from_index(pair_index, raw["natom"])
            blocks[pair_index] = payload
            i0, i1 = shifts[i], shifts[i + 1]
            j0, j1 = shifts[j], shifts[j + 1]
            V[i0:i1, j0:j1] = payload
            if i != j:
                V[j0:j1, i0:i1] = payload.conj().T
    if len(blocks) != expected:
        missing = sorted(set(range(expected)) - set(blocks))
        raise AnalysisError("missing Coulomb pairs for iq {}: {}".format(raw["iq"], missing))
    if not np.all(np.isfinite(V.real)) or not np.all(np.isfinite(V.imag)):
        raise AnalysisError("non-finite Coulomb payload for iq {}".format(raw["iq"]))
    if not np.allclose(V, V.conj().T, rtol=1e-10, atol=1e-10):
        raise AnalysisError("assembled Coulomb matrix is not Hermitian for iq {}".format(raw["iq"]))
    return V


def _spectrum(matrix):
    return np.linalg.eigvalsh(matrix).real


def analyze_directory(directory, basis_path=None, eig_abs=1e-10, eig_rel=1e-8,
                       psd_tol=1e-10):
    basis_path = basis_path or os.path.join(directory, "basis_aux_shrink_out")
    basis = _read_basis(basis_path)
    raw_paths = sorted(glob.glob(os.path.join(directory, "v1_abf_overlap_active_iq_*.dat")))
    if not raw_paths:
        raise AnalysisError("no raw active-ABF overlap files found")
    raw_by_iq = {}
    for path in raw_paths:
        raw = _read_raw(path)
        if raw["iq"] in raw_by_iq:
            raise AnalysisError("duplicate raw iq {}".format(raw["iq"]))
        raw_by_iq[raw["iq"]] = raw
        if raw["naux"] != basis["declared_total"]:
            raise AnalysisError("{}: raw naux does not match basis declared_total".format(path))
        if len(raw_by_iq) > 1:
            reference = raw_by_iq[min(raw_by_iq)]
            if (raw["natom"], raw["naux"], raw["atom_naux"]) != (
                    reference["natom"], reference["naux"], reference["atom_naux"]):
                raise AnalysisError("{}: raw q-point atom metadata mismatch".format(path))

    coulomb_paths = glob.glob(os.path.join(directory, "v1_coulomb_full_iq_*_rank*.dat"))
    coulomb_iqs = set()
    for path in coulomb_paths:
        name = os.path.basename(path)
        try:
            coulomb_iqs.add(int(name.split("_iq_", 1)[1].split("_rank", 1)[0]))
        except (IndexError, ValueError):
            raise AnalysisError("{}: invalid Coulomb shard filename".format(path))
    if coulomb_iqs != set(raw_by_iq):
        raise AnalysisError("raw/Coulomb q-point set mismatch: raw={}, Coulomb={}".format(
            sorted(raw_by_iq), sorted(coulomb_iqs)))

    results = []
    for iq in sorted(raw_by_iq):
        raw = raw_by_iq[iq]
        shard_paths = sorted(glob.glob(os.path.join(directory,
                                                     "v1_coulomb_full_iq_{}_rank*.dat".format(iq))))
        if not shard_paths:
            raise AnalysisError("missing Coulomb shards for iq {}".format(iq))
        shards = [_read_coulomb(path) for path in shard_paths]
        V = _assemble_coulomb(shards, raw)
        s_values = _spectrum(raw["S"])
        s_scale = max(1.0, float(np.max(np.abs(s_values))))
        cutoff = max(float(eig_abs), float(eig_rel) * s_scale)
        keep = s_values > cutoff
        if not np.any(keep):
            raise AnalysisError("no positive S eigenvalues remain for iq {}".format(iq))
        s_eigvals, s_vectors = np.linalg.eigh(raw["S"])
        keep = s_eigvals > cutoff
        X = s_vectors[:, keep] / np.sqrt(s_eigvals[keep])[None, :]
        whitened = X.conj().T.dot(V).dot(X)
        w_values = _spectrum(whitened)
        raw_psd = float(np.min(s_values)) >= -float(psd_tol) * s_scale
        w_scale = max(1.0, float(np.max(np.abs(w_values))))
        whitened_psd = float(np.min(w_values)) >= -float(psd_tol) * w_scale
        retained = s_eigvals[keep]
        results.append({
            "iq": iq,
            "q": list(raw["q"]),
            "q_weight": raw["q_weight"],
            "raw_s": {"min_eigenvalue": float(np.min(s_values)),
                      "max_eigenvalue": float(np.max(s_values)), "psd": raw_psd,
                      "hermitian_max_residual": raw["hermitian_residual"]},
            "whitening": {"cutoff": float(cutoff), "eig_abs": float(eig_abs),
                          "eig_rel": float(eig_rel), "psd_tol": float(psd_tol),
                          "condition_number": float(np.max(retained) / np.min(retained))},
            "whitened_v": {"min_eigenvalue": float(np.min(w_values)),
                           "max_eigenvalue": float(np.max(w_values)),
                           "psd": whitened_psd, "rank": int(np.count_nonzero(keep))},
        })
    return {"basis_kind": "active", "naux": basis["declared_total"], "natom": raw_by_iq[next(iter(raw_by_iq))]["natom"],
            "q_points": results}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory")
    parser.add_argument("--basis", dest="basis_path")
    parser.add_argument("--eig-abs", type=float, default=1e-10)
    parser.add_argument("--eig-rel", type=float, default=1e-8)
    parser.add_argument("--psd-tol", type=float, default=1e-10)
    args = parser.parse_args(argv)
    if args.eig_abs < 0.0 or args.eig_rel < 0.0 or args.psd_tol < 0.0:
        raise AnalysisError("eigenvalue and PSD tolerances must be non-negative")
    result = analyze_directory(args.directory, args.basis_path, args.eig_abs, args.eig_rel, args.psd_tol)
    print(json.dumps(result, sort_keys=True, indent=2))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (AnalysisError, OSError, ValueError, np.linalg.LinAlgError) as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        sys.exit(2)
