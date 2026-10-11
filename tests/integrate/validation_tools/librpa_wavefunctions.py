"""Structural and same-run checks for gauge-dependent LibRPA producer data."""

import math
import struct


def read_ks_wfc(path):
    data = path.read_bytes()
    if len(data) < 24:
        raise ValueError("{} has a truncated KS header".format(path))
    marker, kind, nk, nspin, nbands, nbasis = struct.unpack_from("<6i", data)
    if marker != -12345679 or kind != 28 or min(nk, nspin, nbands, nbasis) <= 0:
        raise ValueError("{} has an invalid KS header".format(path))
    if nspin not in (1, 2):
        raise ValueError("{} has an invalid KS spin count".format(path))
    start = 24 + 12 * nk
    count = nspin * nbands * nbasis
    if len(data) != start + 16 * nk * count:
        raise ValueError("{} has an inconsistent KS payload size".format(path))
    records = []
    for ik in range(nk):
        index, offset = struct.unpack_from("<iq", data, 24 + 12 * ik)
        if index != ik + 1 or offset != start + 16 * ik * count:
            raise ValueError("{} has an invalid KS directory record".format(path))
        values = struct.unpack_from("<{}d".format(2 * count), data, offset)
        if not all(math.isfinite(value) for value in values):
            raise ValueError("{} has non-finite KS coefficients".format(path))
        records.append([complex(values[i], values[i + 1]) for i in range(0, len(values), 2)])
    return (nk, nspin, nbands, nbasis), records


def check_ks_nao(path, nao_dir, tolerance):
    """Compare the MPI-IO assembly to the independent same-run text writer."""
    (nk, nspin, nbands, nbasis), records = read_ks_wfc(path)
    for ik, record in enumerate(records, 1):
        for spin in range(1, nspin + 1):
            name = "wfk{}_nao.txt".format(ik) if nspin == 1 else "wfk{}s{}_nao.txt".format(ik, spin)
            source = nao_dir / name
            lines = [line.strip() for line in source.read_text().splitlines() if line.strip()]
            if len(lines) < 4 or "(index of k points)" not in lines[0]:
                raise ValueError("{} has an invalid text KS header".format(source))
            text_index = ik + (spin - 1) * nk
            if (int(lines[0].split()[0]), int(lines[2].split()[0]), int(lines[3].split()[0])) != (text_index, nbands, nbasis):
                raise ValueError("{} dimensions disagree with {}".format(path, source))
            bands = [i for i, line in enumerate(lines) if line.endswith("(band)")]
            if len(bands) != nbands:
                raise ValueError("{} has the wrong text band count".format(source))
            for band, begin in enumerate(bands):
                end = bands[band + 1] if band + 1 < nbands else len(lines)
                tokens = [float(token) for line in lines[begin + 3:end] for token in line.split()]
                if len(tokens) != 2 * nbasis:
                    raise ValueError("{} has the wrong text coefficient count".format(source))
                values = [complex(tokens[i], tokens[i + 1]) for i in range(0, len(tokens), 2)]
                offset = ((spin - 1) * nbands + band) * nbasis
                for actual, expected in zip(record[offset:offset + nbasis], values):
                    if not math.isfinite(expected.real) or not math.isfinite(expected.imag) or abs(actual - expected) > tolerance:
                        raise ValueError("{} KS payload disagrees with {}".format(path, source))


def check_velocity(path):
    lines = [line.split() for line in path.read_text().splitlines() if line.strip()]
    if len(lines) < 4 or any(len(line) != 1 for line in lines[:4]):
        raise ValueError("{} has a truncated velocity header".format(path))
    nk, nspin, nbands, nbasis = [int(line[0]) for line in lines[:4]]
    if min(nk, nspin, nbands, nbasis) <= 0 or nspin not in (1, 2):
        raise ValueError("{} has invalid velocity dimensions".format(path))
    if len(lines) != 4 + 3 * nk * nspin * (1 + nbands * nbands):
        raise ValueError("{} has an inconsistent velocity payload size".format(path))
    cursor = 4
    for spin in range(1, nspin + 1):
        for ik in range(1, nk + 1):
            for direction in range(1, 4):
                if [int(value) for value in lines[cursor]] != [direction, ik, spin]:
                    raise ValueError("{} has an invalid velocity block index".format(path))
                cursor += 1
                for _ in range(nbands * nbands):
                    row = lines[cursor]
                    if len(row) != 2 or not all(math.isfinite(float(value)) for value in row):
                        raise ValueError("{} has invalid velocity values".format(path))
                    cursor += 1
    return nk, nspin, nbands, nbasis
