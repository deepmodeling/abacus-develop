#!/usr/bin/env python3
import argparse
import glob
import os
import struct
import sys


COULOMB_MARKER = -20129433
CS_MARKER = -10267453


def read_exact(handle, nbytes, label):
    data = handle.read(nbytes)
    if len(data) != nbytes:
        raise RuntimeError(f"truncated {label}")
    return data


def validate_coulomb(path):
    size = os.path.getsize(path)
    with open(path, "rb") as handle:
        marker, iq, naux, value_flag, natoms, nblocks = struct.unpack("<6i", read_exact(handle, 24, path))
        if marker != COULOMB_MARKER:
            raise RuntimeError(f"{path}: bad Coulomb marker {marker}")
        if iq <= 0 or naux <= 0 or natoms <= 0 or nblocks < 0:
            raise RuntimeError(f"{path}: invalid Coulomb header")
        if value_flag not in (0, 1):
            raise RuntimeError(f"{path}: invalid Coulomb value flag {value_flag}")
        atom_naux = list(struct.unpack(f"<{natoms}i", read_exact(handle, 4 * natoms, path)))
        if sum(atom_naux) != naux or any(v <= 0 for v in atom_naux):
            raise RuntimeError(f"{path}: inconsistent atom_naux")
        npairs = natoms * (natoms + 1) // 2
        if nblocks > npairs:
            raise RuntimeError(f"{path}: too many Coulomb blocks")
        header_size = 24 + 4 * natoms + 12 * nblocks
        seen = set()
        ranges = []
        pairs = [(i, j) for i in range(natoms) for j in range(i, natoms)]
        value_bytes = 16 if value_flag == 1 else 8
        for _ in range(nblocks):
            ipair, offset = struct.unpack("<iq", read_exact(handle, 12, path))
            if ipair in seen or ipair < 0 or ipair >= npairs:
                raise RuntimeError(f"{path}: bad Coulomb pair index {ipair}")
            seen.add(ipair)
            i, j = pairs[ipair]
            nbytes = atom_naux[i] * atom_naux[j] * value_bytes
            if offset < header_size or offset + nbytes > size:
                raise RuntimeError(f"{path}: bad Coulomb payload range")
            ranges.append((offset, offset + nbytes))
    for (_, end), (begin, _) in zip(sorted(ranges), sorted(ranges)[1:]):
        if begin < end:
            raise RuntimeError(f"{path}: overlapping Coulomb payload ranges")
    return nblocks


def validate_cs(path):
    size = os.path.getsize(path)
    with open(path, "rb") as handle:
        marker, natoms, ncell = struct.unpack("<3i", read_exact(handle, 12, path))
        nblocks, nblocks_max = struct.unpack("<2q", read_exact(handle, 16, path))
        if marker != CS_MARKER:
            raise RuntimeError(f"{path}: bad Cs marker {marker}")
        if natoms <= 0 or ncell < 0 or nblocks < 0 or nblocks_max < nblocks:
            raise RuntimeError(f"{path}: invalid Cs header")
        header_size = 28 + 36 * nblocks_max
        if header_size > size:
            raise RuntimeError(f"{path}: Cs header exceeds file size")
        ranges = []
        for iblock in range(nblocks_max):
            ia1, ia2, r0, r1, r2 = struct.unpack("<5i", read_exact(handle, 20, path))
            max_abs, offset = struct.unpack("<dq", read_exact(handle, 16, path))
            if iblock >= nblocks:
                if (ia1, ia2, r0, r1, r2, max_abs, offset) != (0, 0, 0, 0, 0, 0.0, 0):
                    raise RuntimeError(f"{path}: nonzero Cs padding record")
                continue
            if ia1 <= 0 or ia1 > natoms or ia2 <= 0 or ia2 > natoms:
                raise RuntimeError(f"{path}: invalid Cs atom index")
            if max_abs < 0.0 or offset < header_size or offset >= size:
                raise RuntimeError(f"{path}: invalid Cs block metadata")
            ranges.append((offset, size))
    return nblocks


def main():
    parser = argparse.ArgumentParser(description="Validate ABACUS direct LibRPA reader-v1 output headers.")
    parser.add_argument("directory")
    parser.add_argument("--coul-prefix", action="append")
    parser.add_argument("--cs-prefix", action="append")
    args = parser.parse_args()

    total = 0
    coul_prefixes = args.coul_prefix or ["v1_coulomb_full_iq_", "v1_coulomb_cut_iq_"]
    cs_prefixes = args.cs_prefix or ["v1_Cs_data_"]
    for prefix in coul_prefixes:
        for path in sorted(glob.glob(os.path.join(args.directory, prefix + "*"))):
            blocks = validate_coulomb(path)
            print(f"OK Coulomb {path}: blocks={blocks}")
            total += 1
    for prefix in cs_prefixes:
        for path in sorted(glob.glob(os.path.join(args.directory, prefix + "*"))):
            blocks = validate_cs(path)
            print(f"OK Cs {path}: blocks={blocks}")
            total += 1
    if total == 0:
        raise RuntimeError("no LibRPA v1 files found")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
