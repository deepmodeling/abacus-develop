#!/usr/bin/env python3
"""Validate an existing onsite +U diagnostic run; does not launch calculations.

Usage: python check_dftu_scf_diagnostics.py CASE/OUT.autotest/scf_diagnostics.jsonl
       [--require-ramp] [--require-dmr] [--require-occupation]
"""
import argparse
import json
import math
from pathlib import Path


def close(a, b, tolerance=2e-10):
    return len(a) == len(b) and all(abs(x-y) < tolerance for x, y in zip(a, b))


def validate(path, require_ramp=False, require_dmr=False, require_occupation=False):
    records = [json.loads(line) for line in Path(path).read_text().splitlines()]
    assert records, 'Missing completed electronic steps'
    previous = None
    ramp_seen = False
    dmr_seen = False
    occupation_seen = False
    for record in records:
        assert record['schema'] == 4
        assert record['ionic_step'] >= 1
        if previous and previous['ionic_step'] == record['ionic_step']:
            assert record['iteration'] == previous['iteration'] + 1
        blocks = record['blocks']
        assert blocks
        correction = 0.0
        for index, block in enumerate(blocks):
            dim = block['dimension']
            for key in ['n_in', 'n_raw_out', 'n_before_occupation_mixing',
                        'n_after_occupation_mixing', 'n_after_dmr_mixing']:
                matrix = block[key]
                if matrix is not None:
                    assert len(matrix) == dim * dim
                    assert all(math.isfinite(x) for x in matrix)
                    assert all(abs(matrix[i*dim+j]-matrix[j*dim+i]) < 2e-10
                               for i in range(dim) for j in range(dim))
            mixed = block['n_after_occupation_mixing']
            if mixed is not None:
                occupation_seen = True
                assert close(mixed, block['n_in'])
                if record['basis_type'] == 'lcao' and previous:
                    old = previous['blocks'][index]['n_in']
                    raw = block['n_before_occupation_mixing']
                    if old is not None:
                        beta = record['mixing_beta']
                        expected = [beta*x+(1-beta)*y for x, y in zip(raw, old)]
                        assert close(mixed, expected), 'Incorrect onsite linear mixing'
            if previous and record['iteration'] > 1:
                prev = previous['blocks'][index]
                ramp_seen |= block['U_eff_eV'] != prev['U_eff_eV']
                if mixed is None and block['n_in'] is not None:
                    expected = prev['n_after_dmr_mixing']
                    if expected is None:
                        expected = prev['n_raw_out']
                    assert close(block['n_in'], expected), 'Hubbard input is not previous output'
            if block['n_in'] is not None:
                n = block['n_in']
                trace_n2 = sum(n[i*dim+j]*n[j*dim+i] for i in range(dim) for j in range(dim))
                correction += 0.5*block['U_eff_eV']*trace_n2
            assert (block['n_after_dmr_mixing'] is not None) == bool(record['dmr_mixing_executed'])
        assert abs(correction-record['source_band_correction_eV']) < 2e-9
        if record['dmr_mixing_executed']:
            dmr_seen = True
            assert record['mixing_restart_count'] > 0
            # No fresh charge coefficients on the step immediately before a restart.
            assert record['mixing_restart_scheduled_step'] != record['iteration'] + 1
        if record['scf_converged']:
            assert record['u_at_target']
        previous = record
    assert not require_ramp or ramp_seen, 'U never advanced'
    assert not require_dmr or dmr_seen, 'DMR mixing never executed'
    assert not require_occupation or occupation_seen, 'Occupation mixing never executed'
    return dict(steps=len(records), ramp_seen=ramp_seen, dmr_seen=dmr_seen,
                occupation_seen=occupation_seen, converged=bool(records[-1]['scf_converged']))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('path')
    parser.add_argument('--require-ramp', action='store_true')
    parser.add_argument('--require-dmr', action='store_true')
    parser.add_argument('--require-occupation', action='store_true')
    args = parser.parse_args()
    print(json.dumps(validate(args.path, args.require_ramp, args.require_dmr, args.require_occupation)))
