#!/usr/bin/env python3

import argparse
import csv
import hashlib
import math
import re
import struct
from pathlib import Path
from typing import NamedTuple


CHI0_V1_MARKER = -41073291
COMPLEX_FLAG = 1
TOTAL_FILE_PATTERN = re.compile(r"^v1_sternheimer_chi0_iq_(\d+)_ifreq_(\d+)_rank(\d+)\.dat$")


class ReaderV1Block(NamedTuple):
    pair_index: int
    offset: int
    count: int


class ReaderV1:
    def __init__(self, path, data, metadata, blocks, payload_start, values):
        self.path = path
        self.data = data
        self.metadata = metadata
        self.blocks = blocks
        self.payload_start = payload_start
        self.values = values

    @classmethod
    def read(cls, path):
        path = Path(path)
        data = path.read_bytes()
        position = 0

        def unpack(fmt, label):
            nonlocal position
            size = struct.calcsize("<" + fmt)
            if position + size > len(data):
                raise ValueError(f"truncated reader-v1 {label} in {path}")
            values = struct.unpack_from("<" + fmt, data, position)
            position += size
            return values[0] if len(values) == 1 else values

        marker, iq, ifrequency, naux, complex_flag, natom = unpack("6i", "header")
        omega, weight = unpack("2d", "frequency metadata")
        nblocks = unpack("i", "block count")
        if marker != CHI0_V1_MARKER:
            raise ValueError(f"invalid reader-v1 marker in {path}")
        if complex_flag != COMPLEX_FLAG:
            raise ValueError(f"reader-v1 file is not complex in {path}")
        if iq <= 0 or ifrequency <= 0 or naux <= 0 or natom <= 0:
            raise ValueError(f"invalid reader-v1 dimensions in {path}")
        if nblocks != natom * (natom + 1) // 2:
            raise ValueError(f"invalid reader-v1 block count in {path}")

        atom_naux = unpack(f"{natom}i", "atom auxiliary dimensions")
        if natom == 1:
            atom_naux = (atom_naux,)
        if any(size <= 0 for size in atom_naux) or sum(atom_naux) != naux:
            raise ValueError(f"invalid reader-v1 atom auxiliary dimensions in {path}")

        raw_blocks = tuple((unpack("i", "pair index"), unpack("q", "block offset")) for _ in range(nblocks))
        header_end = position
        expected_counts = []
        for iatom in range(natom):
            for jatom in range(iatom, natom):
                expected_counts.append(atom_naux[iatom] * atom_naux[jatom])

        blocks = []
        expected_offset = header_end
        for index, ((pair_index, offset), count) in enumerate(zip(raw_blocks, expected_counts)):
            if pair_index != index or offset != expected_offset:
                raise ValueError(f"invalid reader-v1 block layout in {path}")
            blocks.append(ReaderV1Block(pair_index, offset, count))
            expected_offset += count * 2 * struct.calcsize("<d")
        if expected_offset != len(data):
            raise ValueError(f"invalid reader-v1 payload size in {path}")

        values = []
        for block in blocks:
            values.extend(
                complex(*struct.unpack_from("<2d", data, offset))
                for offset in range(block.offset, block.offset + 16 * block.count, 16)
            )
        metadata = (marker, iq, ifrequency, naux, complex_flag, natom, omega, weight, nblocks, atom_naux)
        return cls(path, data, metadata, tuple(blocks), header_end, values)


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _relative_norm(difference, reference):
    difference_norm = math.sqrt(sum(abs(value) ** 2 for value in difference))
    reference_norm = math.sqrt(sum(abs(value) ** 2 for value in reference))
    return difference_norm / max(reference_norm, 1.0e-300)


def _ensure_compatible(reference, candidate):
    if reference.metadata != candidate.metadata:
        raise ValueError(f"reader-v1 metadata mismatch: {reference.path} vs {candidate.path}")
    if reference.blocks != candidate.blocks:
        raise ValueError(f"reader-v1 block layout mismatch: {reference.path} vs {candidate.path}")


def ensure_matching_basis_files(basis40, basis50):
    basis40 = Path(basis40)
    basis50 = Path(basis50)
    digest40 = _sha256(basis40)
    digest50 = _sha256(basis50)
    if digest40 != digest50:
        raise ValueError(f"auxiliary basis order mismatch: {basis40} vs {basis50}")
    return digest40


def build_hybrid_file(total40_path, component40_path, component50_path, output_path):
    total40 = ReaderV1.read(total40_path)
    component40 = ReaderV1.read(component40_path)
    component50 = ReaderV1.read(component50_path)
    _ensure_compatible(total40, component40)
    _ensure_compatible(total40, component50)

    output_path = Path(output_path)
    if output_path.exists():
        raise FileExistsError(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    component_change = [value50 - value40 for value40, value50 in zip(component40.values, component50.values)]
    if component40.data[component40.payload_start :] == component50.data[component50.payload_start :]:
        output_data = total40.data
    else:
        hybrid_values = [total + change for total, change in zip(total40.values, component_change)]
        payload = b"".join(struct.pack("<2d", value.real, value.imag) for value in hybrid_values)
        output_data = total40.data[: total40.payload_start] + payload
    output_path.write_bytes(output_data)
    return {
        "total40_sha256": _sha256(total40_path),
        "component40_sha256": _sha256(component40_path),
        "component50_sha256": _sha256(component50_path),
        "output_sha256": _sha256(output_path),
        "relative_component_change": _relative_norm(component_change, total40.values),
    }


def _indexed_files(directory, pattern):
    indexed = {}
    for path in Path(directory).iterdir():
        match = pattern.match(path.name)
        if match is None:
            continue
        key = tuple(int(value) for value in match.groups())
        if key in indexed:
            raise ValueError(f"duplicate reader-v1 key {key} in {directory}")
        indexed[key] = path
    return indexed


def build_hybrid_campaign(run40, run50, output_root, components=("sos", "pulay", "qspace")):
    run40 = Path(run40)
    run50 = Path(run50)
    output_root = Path(output_root)
    if output_root.exists():
        raise FileExistsError(output_root)
    basis_sha256 = ensure_matching_basis_files(run40 / "basis_aux_out", run50 / "basis_aux_out")
    total40_files = _indexed_files(run40, TOTAL_FILE_PATTERN)
    if not total40_files:
        raise ValueError(f"no total reader-v1 files in {run40}")

    campaign_files = {}
    for component in components:
        if component not in ("sos", "pulay", "qspace"):
            raise ValueError(f"unknown Sternheimer component {component}")
        component_pattern = re.compile(
            rf"^v1_sternheimer_component_{component}_iq_(\d+)_ifreq_(\d+)_rank(\d+)\.dat$"
        )
        component40_files = _indexed_files(run40, component_pattern)
        component50_files = _indexed_files(run50, component_pattern)
        expected_keys = set(total40_files)
        if set(component40_files) != expected_keys or set(component50_files) != expected_keys:
            raise ValueError(f"reader-v1 frequency/rank keys mismatch for component {component}")
        for key in expected_keys:
            total_reader = ReaderV1.read(total40_files[key])
            _ensure_compatible(total_reader, ReaderV1.read(component40_files[key]))
            _ensure_compatible(total_reader, ReaderV1.read(component50_files[key]))
        campaign_files[component] = (component40_files, component50_files)

    output_root.mkdir(parents=True)
    rows = []
    for component in components:
        component40_files, component50_files = campaign_files[component]
        expected_keys = set(total40_files)

        component_output = output_root / component
        component_output.mkdir()
        for iq, ifrequency, rank in sorted(expected_keys):
            key = (iq, ifrequency, rank)
            output_path = component_output / total40_files[key].name
            result = build_hybrid_file(
                total40_files[key], component40_files[key], component50_files[key], output_path
            )
            rows.append(
                {
                    "component": component,
                    "iq": iq,
                    "ifrequency": ifrequency,
                    "rank": rank,
                    "basis_sha256": basis_sha256,
                    "total40_path": str(total40_files[key]),
                    "component40_path": str(component40_files[key]),
                    "component50_path": str(component50_files[key]),
                    "output_path": str(output_path),
                    **result,
                }
            )

    manifest_path = output_root / "hybrid_manifest.csv"
    fieldnames = list(rows[0])
    with manifest_path.open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return manifest_path


def main(argv=None):
    parser = argparse.ArgumentParser(description="Build Delta-ST component hybrid reader-v1 matrices.")
    subparsers = parser.add_subparsers(dest="command")

    file_parser = subparsers.add_parser("file", help="Build one hybrid reader-v1 file.")
    file_parser.add_argument("total40", type=Path)
    file_parser.add_argument("component40", type=Path)
    file_parser.add_argument("component50", type=Path)
    file_parser.add_argument("output", type=Path)
    file_parser.add_argument("--basis40", type=Path)
    file_parser.add_argument("--basis50", type=Path)

    campaign_parser = subparsers.add_parser("campaign", help="Build all matching frequency files.")
    campaign_parser.add_argument("run40", type=Path)
    campaign_parser.add_argument("run50", type=Path)
    campaign_parser.add_argument("output_root", type=Path)
    campaign_parser.add_argument(
        "--components",
        nargs="+",
        choices=("sos", "pulay", "qspace"),
        default=("sos", "pulay", "qspace"),
    )

    args = parser.parse_args(argv)
    if args.command == "file":
        if (args.basis40 is None) != (args.basis50 is None):
            file_parser.error("--basis40 and --basis50 must be provided together")
        if args.basis40 is not None:
            ensure_matching_basis_files(args.basis40, args.basis50)
        result = build_hybrid_file(args.total40, args.component40, args.component50, args.output)
        print(",".join(result.keys()))
        print(",".join(str(value) for value in result.values()))
    elif args.command == "campaign":
        print(build_hybrid_campaign(args.run40, args.run50, args.output_root, tuple(args.components)))
    else:
        parser.error("choose either the file or campaign command")
    return 0


if __name__ == "__main__":
    main()
