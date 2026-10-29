#!/usr/bin/env python3
"""Compare complete norm-conserving PW eigenvectors up to one phase per band.

This deliberately rejects degenerate eigenvalues: an independently validated
subspace comparison is required for that situation. No coefficient is dropped
or renormalized. Both files must contain normalized, mutually orthogonal states.
"""
import argparse
import json
import math
from pathlib import Path
import re
import sys

HEADER = ('Kpoint nKpoint kv.x kv.y kv.z weight ngtot nband ecut lat0 2pi/lat0').split()


def finite(value, context):
    try:
        number = float(value)
    except ValueError as error:
        raise ValueError(f'{context}: invalid number {value!r}') from error
    if not math.isfinite(number):
        raise ValueError(f'{context}: nonfinite number')
    return number


def integer(value, context):
    if not re.fullmatch(r'[+-]?\d+', value):
        raise ValueError(f'{context}: invalid integer')
    return int(value)


def read_wavefunction(path):
    path = Path(path)
    print(f'Read {path.resolve()}')
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    cursor = 0

    def line():
        nonlocal cursor
        if cursor == len(lines):
            raise ValueError(f'{path}: truncated wavefunction')
        result = lines[cursor]
        cursor += 1
        return result

    def expect(text):
        if line() != text:
            raise ValueError(f'{path}: expected {text!r}')

    if line().split() != HEADER:
        raise ValueError(f'{path}: invalid header')
    fields = line().split()
    if len(fields) != len(HEADER):
        raise ValueError(f'{path}: invalid header field count')
    header = {key: finite(value, str(path)) for key, value in zip(HEADER, fields)}
    for key in ('Kpoint', 'nKpoint', 'ngtot', 'nband'):
        header[key] = integer(fields[HEADER.index(key)], str(path))
        if header[key] <= 0:
            raise ValueError(f'{path}: {key} must be positive')
    if header['Kpoint'] > header['nKpoint']:
        raise ValueError(f'{path}: invalid k-point index')
    if any(header[key] <= 0 for key in ('ecut', 'lat0', '2pi/lat0')):
        raise ValueError(f'{path}: nonpositive basis metadata')
    expect('<Reciprocal Lattice Vector>')
    lattice = []
    for _ in range(3):
        values = line().split()
        if len(values) != 3:
            raise ValueError(f'{path}: invalid reciprocal lattice')
        lattice.append(tuple(finite(x, str(path)) for x in values))
    expect('<Reciprocal Lattice Vector>')
    determinant = sum(lattice[0][i] * (lattice[1][(i+1)%3] * lattice[2][(i+2)%3]
                      - lattice[1][(i+2)%3] * lattice[2][(i+1)%3]) for i in range(3))
    if determinant == 0 or not math.isfinite(determinant):
        raise ValueError(f'{path}: singular reciprocal lattice')
    expect('<G vectors>')
    vectors = []
    for _ in range(header['ngtot']):
        values = line().split()
        if len(values) != 3:
            raise ValueError(f'{path}: invalid G vector')
        vector = tuple(integer(x, str(path)) for x in values)
        if any(x < 0 for x in vector):
            raise ValueError(f'{path}: negative FFT-grid index')
        vectors.append(vector)
    if len(vectors) != len(set(vectors)):
        raise ValueError(f'{path}: duplicate G vector')
    expect('<G vectors>')
    bands = []
    polarizations = None
    order = sorted(range(len(vectors)), key=vectors.__getitem__)
    for index in range(1, header['nband'] + 1):
        marker = f'< Band {index} >'
        expect(marker)
        values = []
        while True:
            record = line()
            if record == marker:
                break
            values.extend(finite(x, str(path)) for x in record.split())
        if len(values) not in (2 * header['ngtot'], 4 * header['ngtot']):
            raise ValueError(f'{path}: incorrect coefficient count in band {index}')
        npol = len(values) // (2 * header['ngtot'])
        if polarizations is not None and polarizations != npol:
            raise ValueError(f'{path}: inconsistent polarization count')
        polarizations = npol
        raw = [complex(values[i], values[i+1]) for i in range(0, len(values), 2)]
        bands.append([raw[spin * header['ngtot'] + i] for spin in range(npol) for i in order])
    if cursor != len(lines):
        raise ValueError(f'{path}: unexpected trailing data')
    return {'header': header, 'lattice': lattice, 'vectors': [vectors[i] for i in order],
            'npol': polarizations, 'bands': bands}


def inner(left, right):
    products = [a.conjugate() * b for a, b in zip(left, right)]
    return complex(math.fsum(z.real for z in products), math.fsum(z.imag for z in products))


def check_orthonormal(wave, tolerance):
    maximum = 0.0
    for i, left in enumerate(wave['bands']):
        for j in range(i + 1):
            value = inner(left, wave['bands'][j])
            expected = 1.0 if i == j else 0.0
            error = abs(value - expected)
            if not math.isfinite(error) or error > tolerance:
                raise ValueError(f'wavefunction normalization/orthogonality error: {error:.9e}')
            maximum = max(maximum, error)
    return maximum


def compare_wavefunctions(reference, actual, tolerance):
    if not math.isfinite(tolerance) or tolerance <= 0:
        raise ValueError('coefficient tolerance must be finite and positive')
    for key in ('header', 'lattice', 'vectors', 'npol'):
        if reference[key] != actual[key]:
            raise ValueError(f'wavefunction {key} mismatch')
    orthogonality = [check_orthonormal(wave, tolerance) for wave in (reference, actual)]
    errors = []
    for band, (left, right) in enumerate(zip(reference['bands'], actual['bands']), 1):
        overlap = inner(left, right)
        if abs(overlap) == 0 or not math.isfinite(abs(overlap)):
            raise ValueError(f'band {band}: zero/nonfinite phase overlap')
        phase = overlap / abs(overlap)
        differences = [abs(b - phase * a) for a, b in zip(left, right)]
        maximum = max(differences)
        distance = math.sqrt(math.fsum(x*x for x in differences))
        if not math.isfinite(maximum) or not math.isfinite(distance) or maximum > tolerance:
            raise ValueError(f'band {band}: full coefficient error {maximum:.9e} exceeds {tolerance:.9e}')
        errors.append({'band': band, 'max_coefficient_error': maximum, 'vector_l2_error': distance})
    return {'orthogonality_errors': orthogonality, 'bands': errors}


def read_bands(path, count, bands):
    path = Path(path)
    print(f'Read {path.resolve()}')
    lines = [line.split() for line in path.read_text().splitlines() if line.strip()]
    if len(lines) != count:
        raise ValueError(f'{path}: wrong number of k-points')
    result = []
    for index, fields in enumerate(lines, 1):
        if len(fields) != bands + 2 or integer(fields[0], str(path)) != index:
            raise ValueError(f'{path}: invalid band table shape/index')
        row = [finite(x, str(path)) for x in fields[1:]]
        result.append(row)
    return result


def read_directory(path):
    files = sorted(Path(path).glob('wfk*_pw.txt'))
    if not files:
        raise ValueError(f'{path}: no PW wavefunctions')
    result = {}
    for file in files:
        wave = read_wavefunction(file)
        index = wave['header']['Kpoint']
        if index in result:
            raise ValueError(f'{path}: duplicate k-point')
        result[index] = wave
    count = len(result)
    if sorted(result) != list(range(1, count+1)):
        raise ValueError(f'{path}: missing k-point')
    if any(wave['header']['nKpoint'] != count for wave in result.values()):
        raise ValueError(f'{path}: incomplete k-point set')
    if len({wave['header']['nband'] for wave in result.values()}) != 1:
        raise ValueError(f'{path}: inconsistent band count')
    return result


def compare_directories(reference_dir, actual_dir, reference_bands, actual_bands,
                        coefficient_tolerance, energy_tolerance):
    if not math.isfinite(energy_tolerance) or energy_tolerance <= 0:
        raise ValueError('energy tolerance must be finite and positive')
    reference = read_directory(reference_dir)
    actual = read_directory(actual_dir)
    if set(reference) != set(actual):
        raise ValueError('k-point sets differ')
    count = len(reference)
    bands = reference[1]['header']['nband']
    left = read_bands(reference_bands, count, bands)
    right = read_bands(actual_bands, count, bands)
    maximum_energy = 0.0
    for a, b in zip(left, right):
        if a[0] != b[0]:
            raise ValueError('band-path coordinates differ')
        for row in (a, b):
            if any(row[i+1] - row[i] <= 2 * energy_tolerance for i in range(1, bands)):
                raise ValueError('degenerate/unordered bands require a separate subspace comparison')
        maximum_energy = max(maximum_energy, max(abs(x-y) for x,y in zip(a[1:], b[1:])))
    if maximum_energy > energy_tolerance:
        raise ValueError(f'band energy error {maximum_energy:.9e} exceeds {energy_tolerance:.9e}')
    results = {str(index): compare_wavefunctions(reference[index], actual[index], coefficient_tolerance)
               for index in sorted(reference)}
    return {'max_energy_error_eV': maximum_energy, 'kpoints': results}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('reference_dir')
    parser.add_argument('actual_dir')
    parser.add_argument('reference_bands')
    parser.add_argument('actual_bands')
    parser.add_argument('--coefficient-tolerance', type=float, default=1e-7)
    parser.add_argument('--energy-tolerance', type=float, default=1e-8)
    args = parser.parse_args(argv)
    try:
        result = compare_directories(args.reference_dir, args.actual_dir, args.reference_bands,
                                     args.actual_bands, args.coefficient_tolerance, args.energy_tolerance)
    except (OSError, ValueError, OverflowError) as error:
        print(f'PW wavefunction comparison failed: {error}', file=sys.stderr)
        return 1
    print(json.dumps(result, allow_nan=False))
    return 0


if __name__ == '__main__':
    sys.exit(main())
