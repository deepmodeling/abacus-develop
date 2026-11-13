"""Run fixed-input one-rank-per-GPU cases with explicit CPU and device placement."""

import argparse
import csv
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone

from benchmark import digest, summarize


ROOT = Path(__file__).resolve().parents[1]


def gpu_inventory() -> dict[int, dict]:
    """Return physical GPU indices, UUIDs and total memory in MiB."""
    output = subprocess.check_output([
        'nvidia-smi', '--query-gpu=index,uuid,name,memory.total',
        '--format=csv,noheader,nounits'], text=True)
    return {int(row[0]): {'uuid': row[1].strip(), 'name': row[2].strip(),
                         'memory_MiB': float(row[3])}
            for row in csv.reader(output.splitlines())}


def parse_gpu_processes(text: str, uuids: list[str]) -> dict[str, list[int]]:
    """Parse observed driver PIDs for the requested GPU UUIDs.

    Driver PIDs can belong to a different namespace from MPI rank PIDs.
    Repeated sampling rows are deduplicated without equating those namespaces.
    Malformed nonempty CSV rows raise ValueError.
    """
    processes = {}
    for row in csv.reader(text.splitlines()):
        if not row:
            continue
        if len(row) != 3:
            raise ValueError(f'Unexpected GPU process sampling row: {row}')
        uuid = row[0].strip()
        if uuid in uuids:
            processes.setdefault(uuid, set()).add(int(row[1]))
    return {uuid: sorted(pids) for uuid, pids in processes.items()}


def isolated_gpu_processes(processes: dict[str, list[int]], uuids: list[str]) -> bool:
    """Return whether each requested GPU has one distinct observed driver PID."""
    return (set(processes) == set(uuids) and all(len(pids) == 1 for pids in processes.values())
            and len({pid for pids in processes.values() for pid in pids}) == len(uuids))


def summarize_gpu(path: Path, devices: list[int]) -> dict:
    """Summarize per-device memory and complete sampling rounds.

    Parameters
    ----------
    path : Path
        CSV output from the persistent nvidia-smi sampler.
    devices : list[int]
        Physical device indices present in every complete sampling round.

    Returns
    -------
    dict
        Per-device peaks and the largest sum within a complete round, in MiB.
        A trailing partial round contributes only to per-device statistics.

    Raises
    ------
    ValueError
        If rows are malformed or the device coverage is inconsistent.

    Notes
    -----
    Device readings within a round are not atomic. The aggregate statistic
    sums readings from each complete round, not independently observed peaks.
    """
    rows = []
    for values in csv.reader(path.read_text().splitlines()):
        if not values:
            continue
        if len(values) != 5:
            raise ValueError(f'Unexpected GPU sampling row in {path}: {values}')
        rows.append({'timestamp': values[0].strip(), 'index': int(values[1]),
                     'uuid': values[2].strip(), 'memory_MiB': float(values[3]),
                     'utilization_percent': float(values[4])})
    expected = set(devices)
    if {row['index'] for row in rows} != expected:
        raise ValueError(f'GPU samples do not cover requested devices {devices}')
    peaks = {str(device): max(row['memory_MiB'] for row in rows if row['index'] == device)
             for device in devices}
    utilization = {str(device): max(row['utilization_percent'] for row in rows if row['index'] == device)
                   for device in devices}
    totals = []
    for first in range(0, len(rows) - len(devices) + 1, len(devices)):
        group = rows[first:first + len(devices)]
        if {row['index'] for row in group} != expected:
            raise ValueError('GPU sample rounds have missing or repeated devices')
        totals.append(sum(row['memory_MiB'] for row in group))
    return {'peak_per_gpu_MiB': peaks, 'maximum_device_peak_MiB': max(peaks.values()),
            'peak_sampled_total_MiB': max(totals), 'maximum_utilization_percent': utilization,
            'complete_sampling_rounds': len(totals), 'trailing_sampling_rows': len(rows) % len(devices),
            'sample_interval_ms': 100}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--workspace', type=Path, default=ROOT)
    parser.add_argument('--case', required=True)
    parser.add_argument('--label', required=True)
    parser.add_argument('--devices', required=True, help='Comma-separated physical nvidia-smi indices')
    parser.add_argument('--threads', type=int, default=12)
    parser.add_argument('--cores-per-rank', type=int, default=6)
    parser.add_argument('--steps', type=int)
    parser.add_argument('--profile', action='store_true')
    args = parser.parse_args()
    if not re.fullmatch(r'[A-Za-z0-9_.-]+', args.label):
        parser.error('label must contain only letters, digits, underscores, periods, or hyphens')
    if not re.fullmatch(r'\d+(,\d+)*', args.devices):
        parser.error('devices must be comma-separated nonnegative GPU indices')
    devices = [int(value) for value in args.devices.split(',')]
    if len(set(devices)) != len(devices) or min(args.threads, args.cores_per_rank) < 1:
        parser.error('GPU indices must be unique; thread and core counts must be positive')
    if args.steps is not None and args.steps < 1:
        parser.error('steps must be positive')
    inventory = gpu_inventory()
    if any(device not in inventory for device in devices):
        parser.error(f'Requested devices {devices} are not present in {list(inventory)}')
    binary = args.binary.resolve(strict=True)
    workspace = args.workspace.resolve(strict=True)
    inputs = (workspace / 'inputs' / args.case).resolve(strict=True)
    if inputs.parent != (workspace / 'inputs').resolve():
        parser.error('case must name a direct child of the input directory')
    uuids = [inventory[device]['uuid'] for device in devices]
    active = subprocess.check_output([
        'nvidia-smi', '--query-compute-apps=gpu_uuid,pid,used_gpu_memory',
        '--format=csv,noheader,nounits'], text=True)
    busy = parse_gpu_processes(active, uuids)
    if busy:
        raise RuntimeError(f'Requested GPUs are already in use; wait or select idle devices: {busy}')
    directory = ROOT / 'multigpu/results' / args.label
    directory.mkdir(parents=True, exist_ok=False)
    for name in ('INPUT', 'STRU', 'KPT'):
        shutil.copy2(inputs / name, directory / name)
    for asset in (workspace / 'assets').iterdir():
        (directory / asset.name).symlink_to(asset)
    if args.steps is not None:
        path = directory / 'INPUT'
        content, count = re.subn(r'(?m)^scf_nmax\s+.*$', f'scf_nmax {args.steps}', path.read_text())
        if count != 1:
            raise ValueError('Expected exactly one scf_nmax entry')
        path.write_text(content)
    environment = dict(os.environ)
    environment.update(OMP_NUM_THREADS=str(args.threads), OPENBLAS_NUM_THREADS='1',
                       OMP_PROC_BIND='spread', OMP_PLACES='threads',
                       OMPI_ALLOW_RUN_AS_ROOT='1', OMPI_ALLOW_RUN_AS_ROOT_CONFIRM='1')
    command = ['mpirun', '-np', str(len(devices)), '--map-by', f'slot:PE={args.cores_per_rank}',
               '--bind-to', 'core', '--report-bindings', sys.executable,
               str(Path(__file__).with_name('rank_exec.py')), '--binary', str(binary),
               '--gpu-uuids', ','.join(uuids)]
    if args.profile:
        command = ['nsys', 'profile', '--trace=cuda,nvtx,mpi', '--mpi-impl=openmpi',
                   '--sample=none', '--cpuctxsw=none', '--output=profile'] + command
    metadata = {'case': args.case, 'label': args.label, 'binary': str(binary),
                'binary_sha256': digest(binary), 'started_at': datetime.now(timezone.utc).isoformat(),
                'command': command, 'threads': args.threads, 'cores_per_rank': args.cores_per_rank,
                'ranks': len(devices), 'gpu_count': len(devices), 'devices': devices,
                'gpu_uuids': uuids, 'gpu_inventory': inventory, 'workspace': str(workspace),
                'gpu_processes_before_run': busy,
                'expected_kpar_from_gpu_source': len(devices), 'step_limit_override': args.steps,
                'profiled': args.profile, 'omp_proc_bind': 'spread', 'omp_places': 'threads',
                'omp_wait_policy': environment.get('OMP_WAIT_POLICY', 'runtime_default'),
                'gomp_spincount': environment.get('GOMP_SPINCOUNT', 'runtime_default'),
                'runner_sha256': digest(Path(__file__)),
                'rank_wrapper_sha256': digest(Path(__file__).with_name('rank_exec.py')),
                'input_sha256': {name: digest(directory / name) for name in ('INPUT', 'STRU', 'KPT')},
                'asset_sha256': {path.name: digest(path) for path in (workspace / 'assets').iterdir()}}
    (directory / 'metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    with (directory / 'gpu.csv').open('x') as gpu_log, (directory / 'processes.csv').open('x') as process_log:
        samplers = [
            (['nvidia-smi', f'--id={args.devices}',
              '--query-gpu=timestamp,index,uuid,memory.used,utilization.gpu',
              '--format=csv,noheader,nounits', '--loop-ms=100'], gpu_log),
            (['nvidia-smi', f'--id={args.devices}', '--query-compute-apps=gpu_uuid,pid,used_gpu_memory',
              '--format=csv,noheader,nounits', '--loop-ms=200'], process_log)]
        monitors = []
        try:
            for sampler, log in samplers:
                monitors.append(subprocess.Popen(sampler, stdout=log, stderr=subprocess.STDOUT))
            start = time.monotonic()
            with (directory / 'run.log').open('x') as log:
                completed = subprocess.run(command, cwd=directory, env=environment, stdin=subprocess.DEVNULL,
                                           stdout=log, stderr=subprocess.STDOUT)
            elapsed = time.monotonic() - start
        finally:
            for monitor in monitors:
                monitor.terminate()
                monitor.wait(timeout=5)
    rank_files = sorted(directory.glob('rank-*.json'))
    placement = [json.loads(path.read_text()) for path in rank_files]
    placement_ok = len(placement) == len(devices) and {row['gpu_uuid'] for row in placement} == set(uuids)
    assigned_cpus = set()
    for row in placement:
        cpus = set(row['cpu_affinity'])
        placement_ok = placement_ok and len(cpus) >= args.threads and not assigned_cpus.intersection(cpus)
        assigned_cpus.update(cpus)
    observed_processes = parse_gpu_processes((directory / 'processes.csv').read_text(), uuids)
    placement_ok = placement_ok and isolated_gpu_processes(observed_processes, uuids)
    result = {**metadata, 'elapsed_seconds': elapsed, 'return_code': completed.returncode,
              **summarize(directory), 'placement_verified': bool(placement_ok), 'rank_placement': placement,
              'observed_gpu_pids': observed_processes,
              'gpu_measurements': summarize_gpu(directory / 'gpu.csv', devices)}
    (directory / 'result.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({key: result[key] for key in ('label', 'elapsed_seconds', 'return_code',
                     'native_scf_converged', 'scf_steps', 'placement_verified', 'gpu_measurements')}), flush=True)
    if completed.returncode or not result['runtime_complete'] or not placement_ok:
        raise SystemExit(completed.returncode or 1)


if __name__ == '__main__':
    main()
