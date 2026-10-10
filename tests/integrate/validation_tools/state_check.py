#!/usr/bin/env python3
"""Check two-step LR root following without depending on degenerate eigenvector gauges."""
import json
import math
from pathlib import Path
import re
import sys


def check(log, expected):
    steps = re.findall(r'EXCITED-STATE RELAX step (\d+): (.*)', log)
    if [int(index) for index, _ in steps] != list(range(expected['steps'])):
        raise ValueError('missing or repeated ionic steps')
    for _, line in steps:
        values = []
        for label in ['E_gs', 'Omega', 'E_exc', r'\|F_gs\|max', r'\|F_Omega\|max']:
            match = re.search(label + r'\s*=\s*(\S+)', line)
            if match is None:
                raise ValueError('missing energy or force summary')
            values.append(float(match.group(1).rstrip(',')))
        if not all(math.isfinite(value) for value in values):
            raise ValueError('nonfinite energy or force')
        if values[1] <= 0 or min(values[3:]) < 0:
            raise ValueError('unstable excitation or invalid force norm')
        if abs(values[0] + values[1] - values[2]) > 1e-6:
            raise ValueError('inconsistent excited-state energy')
    matches = re.findall(r'subspace dimensions (\d+) -> (\d+); (.*)', log)
    if len(matches) != expected['steps'] - 1:
        raise ValueError('missing cross-step subspace match')
    for old, new, line in matches:
        if [int(old), int(new)] != expected['dimensions']:
            raise ValueError('unexpected group split/merge')
        coverage = re.search(r'old/new coverage (\S+) / ([^;]+)', line)
        priority = re.search(r'reference-first ([01])', line)
        spectrum = re.search(r'singular values ([^;]+)', line)
        scores = re.search(r'best/runner-up score (\S+) / ([^;]+)', line)
        if not all([coverage, priority, spectrum, scores]):
            raise ValueError('missing matching diagnostics')
        if bool(int(priority.group(1))) != expected['reference_first']:
            raise ValueError('wrong group selection criterion')
        singular = [float(value) for value in spectrum.group(1).split()]
        retained = float(scores.group(1))
        if not expected['reference_first']:
            retained = min(float(coverage.group(1)), float(coverage.group(2)))
        diagnostics = singular + [retained, float(coverage.group(1)), float(coverage.group(2)),
                                  float(scores.group(1)), float(scores.group(2))]
        if not all(math.isfinite(value) for value in diagnostics):
            raise ValueError('nonfinite matching diagnostics')
        if any(value < 0 or value > 1.0 + 1e-6 for value in diagnostics):
            raise ValueError('invalid normalized overlap or coverage')
        if len(singular) != min(int(old), int(new)):
            raise ValueError('wrong singular spectrum dimension')
        if min(singular) < expected['min_overlap'] or retained < expected['min_overlap']:
            raise ValueError('lost tracked subspace or actual JT branch')
        if float(scores.group(1)) - float(scores.group(2)) < 0.5:
            raise ValueError('ambiguous group assignment')


if __name__ == '__main__':
    try:
        check(Path(sys.argv[1]).read_text(), json.loads(Path(sys.argv[2]).read_text()))
    except (ValueError, OSError, IndexError) as error:
        print('Root tracking check failed:', error)
        sys.exit(1)
    print('Root tracking check passed')
