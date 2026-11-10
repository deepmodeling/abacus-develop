"""Verify fixed-input multi-GPU scaling results and report physical differences."""

import argparse
import json
from pathlib import Path
import re

from benchmark import digest, summarize
from run_multigpu import isolated_gpu_processes, parse_gpu_processes, summarize_gpu


ROOT = Path(__file__).resolve().parents[1]


def load_verified(folder: Path) -> dict:
    """Reparse one complete result and verify its input and sampling records.

    Parameters
    ----------
    folder : Path
        Directory produced by run_multigpu.py.

    Returns
    -------
    dict
        Verified result with k-point counts reported by the application.

    Raises
    ------
    ValueError
        If the calculation is incomplete, profiled, unconverged, inconsistent
        with its raw records, or lacks verified process/device placement.
    """
    result = json.loads((folder / 'result.json').read_text())
    for key, value in summarize(folder).items():
        if result[key] != value:
            raise ValueError(f'Reparsed output differs for {folder.name}: {key}')
    if result['return_code'] or not result['runtime_complete'] or not result['finite_outputs']:
        raise ValueError(f'Incomplete or nonfinite result: {folder}')
    if not result['native_scf_converged'] or result['step_limit_override'] is not None or result['profiled']:
        raise ValueError(f'Scaling requires a complete unprofiled converged run: {folder}')
    if not result['placement_verified']:
        raise ValueError(f'Unverified process/device placement: {folder}')
    observed = parse_gpu_processes((folder / 'processes.csv').read_text(), result['gpu_uuids'])
    if observed != result['observed_gpu_pids'] or not isolated_gpu_processes(observed, result['gpu_uuids']):
        raise ValueError(f'GPU processes do not match isolated device placement: {folder}')
    for key in ('input_sha256', 'asset_sha256'):
        for name, expected in result[key].items():
            if digest(folder / name) != expected:
                raise ValueError(f'Input or asset checksum mismatch: {folder / name}')
    measured = summarize_gpu(folder / 'gpu.csv', result['devices'])
    if measured != result['gpu_measurements']:
        raise ValueError(f'GPU sampling summary mismatch: {folder}')
    text = (folder / 'OUT.ABACUS/running_scf.log').read_text()
    total = re.findall(r'\bnkstot\s*=\s*(\d+)', text)
    local = re.findall(r'k-point number in this process\s*=\s*(\d+)', text)
    result['native_kpoints'] = {'total': int(total[-1]) if total else None,
                               'root_process': int(local[-1]) if local else None}
    return result


def physical_changes(reference: dict, candidate: dict) -> dict[str, float]:
    """Compare verified physical outputs with identical scientific inputs.

    Parameters
    ----------
    reference, candidate : dict
        Results with complete finite physical tables. MPI ranks may differ.

    Returns
    -------
    dict[str, float]
        Absolute energy difference and maximum absolute table differences in
        the units named by each result field.

    Raises
    ------
    ValueError
        If systems, input checksums, asset checksums or table lengths differ.
    """
    for key in ('case', 'input_sha256', 'asset_sha256'):
        if reference[key] != candidate[key]:
            raise ValueError(f'Comparison requires identical {key}')
    changes = {'energy_eV': abs(reference['energy_eV'] - candidate['energy_eV'])}
    for key in ('forces_eV_per_A', 'stress_kbar', 'moments_uB', 'magnetic_forces_eV_per_uB'):
        if len(reference[key]) != len(candidate[key]):
            raise ValueError(f'Physical table length mismatch: {key}')
        changes[key] = max(abs(a - b) for x, y in zip(reference[key], candidate[key]) for a, b in zip(x, y))
    return changes


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prefix', default='baseline')
    parser.add_argument('--threads', type=int, default=12)
    args = parser.parse_args()
    root = ROOT / 'multigpu/results'
    rows = []
    binary_sha = None
    for short in ('fe16', 'fe32'):
        reference = load_verified(root / f'{args.prefix}-{short}-g1-t{args.threads}-full')
        for count in (1, 2, 4):
            label = f'{args.prefix}-{short}-g{count}-t{args.threads}-full'
            result = load_verified(root / label)
            if result['gpu_count'] != count or result['ranks'] != count or result['threads'] != args.threads:
                raise ValueError(f'Parallel settings do not match the requested series: {label}')
            if binary_sha is None:
                binary_sha = result['binary_sha256']
            if result['binary_sha256'] != binary_sha:
                raise ValueError('Scaling series mixes application binaries')
            speedup = reference['elapsed_seconds'] / result['elapsed_seconds']
            rows.append({'label': label, 'case': result['case'], 'gpus': count,
                         'threads_per_rank': args.threads, 'total_omp_threads': count * args.threads,
                         'seconds': result['elapsed_seconds'], 'speedup_over_one_gpu': speedup,
                         'parallel_efficiency': speedup / count, 'scf_steps': result['scf_steps'],
                         'gpu_memory': result['gpu_measurements'], 'native_kpoints': result['native_kpoints'],
                         'maximum_absolute_changes': physical_changes(reference, result),
                         'spin_rms_uB': result['spin_rms_uB'], 'timers': result['timers']})
    output = {'scope': 'Same-node fixed-input full scaling with one MPI rank per GPU and fixed threads per rank. '
                        'Allocated CPU resources grow with GPU count. Per-device and concurrent-round total memory are distinct. '
                        'Native host timers overlap; magnetic diagnostics are separate from electronic convergence.',
              'binary_sha256': binary_sha, 'rows': rows}
    with (root / f'{args.prefix}-scaling-t{args.threads}.json').open('x') as stream:
        json.dump(output, stream, indent=2, allow_nan=False)
        stream.write('\n')
    for row in rows:
        print(json.dumps({key: row[key] for key in ('case', 'gpus', 'threads_per_rank', 'seconds',
                         'speedup_over_one_gpu', 'parallel_efficiency', 'scf_steps', 'maximum_absolute_changes')}))


if __name__ == '__main__':
    main()
