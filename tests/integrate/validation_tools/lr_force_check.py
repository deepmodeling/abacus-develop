#!/usr/bin/env python3
"""Compare signed single-point LR force components against a numerical baseline."""
import argparse
import math
from pathlib import Path
import re
import sys

CHANNEL = re.compile(r'Start to calculate excited-state force of (singlet|triplet|updown)\b')
ROW = re.compile(r'^\s*(?:(\d+)\s+)?(\d+)\s+force\s+(\S+)\s+(\S+)\s+(\S+)\s*$')


def add_row(rows, key, components):
    if key in rows:
        raise ValueError('duplicate LR force entry: {}'.format(key))
    values = tuple(float(value) for value in components)
    if len(values) != 3 or not all(math.isfinite(value) for value in values):
        raise ValueError('invalid LR force components: {}'.format(key))
    rows[key] = values


def read_forces(log):
    """Read only F_Omega tables, excluding ground-state and total excited-state forces."""
    rows = {}
    channel = None
    state = None
    active = False
    for line in log.splitlines():
        match = CHANNEL.search(line)
        if match:
            channel = match.group(1)
            active = False
            continue
        if 'Forces (-gradients) of each excited state:' in line:
            if channel is None:
                raise ValueError('LR force table has no spin channel')
            active = True
            state = None
            continue
        if not active:
            continue
        match = ROW.match(line)
        if match:
            if match.group(1) is not None:
                state = int(match.group(1))
            if state is None:
                raise ValueError('LR force row has no state index')
            key = (channel, state, int(match.group(2)))
            add_row(rows, key, match.group(3, 4, 5))
        elif not line.strip() or line.split() == ['state', 'atom', 'x', 'y', 'z']:
            continue
        elif re.search(r'\bforce\b', line):
            raise ValueError('malformed LR force row: ' + line)
        else:
            active = False
    if not rows:
        raise ValueError('no LR force components found')
    return rows


def read_reference(text):
    rows = {}
    for line in text.splitlines():
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        fields = line.split()
        if len(fields) != 6 or fields[0] not in ('singlet', 'triplet', 'updown'):
            raise ValueError('malformed LR force reference row: ' + line)
        key = (fields[0], int(fields[1]), int(fields[2]))
        add_row(rows, key, fields[3:])
    if not rows:
        raise ValueError('empty LR force reference')
    return rows


def compare(actual, expected, tolerance):
    if not math.isfinite(tolerance) or tolerance < 0:
        raise ValueError('LR force tolerance must be finite and nonnegative')
    if actual.keys() != expected.keys():
        raise ValueError('LR force entries differ: missing={}, extra={}'.format(
            sorted(expected.keys() - actual.keys()), sorted(actual.keys() - expected.keys())))
    for key in expected:
        for axis, value, reference in zip('xyz', actual[key], expected[key]):
            if abs(value - reference) > tolerance:
                raise ValueError('LR force {} {}: calculated={}, reference={}, tolerance={}'.format(
                    key, axis, value, reference, tolerance))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('log', type=Path)
    parser.add_argument('reference', type=Path)
    parser.add_argument('tolerance', type=float)
    args = parser.parse_args()
    try:
        actual = read_forces(args.log.read_text())
        expected = read_reference(args.reference.read_text())
        compare(actual, expected, args.tolerance)
    except (OSError, ValueError) as error:
        print('LR force comparison failed: {}'.format(error), file=sys.stderr)
        return 1
    print('Signed LR force comparison passed: {} components'.format(3 * len(expected)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
