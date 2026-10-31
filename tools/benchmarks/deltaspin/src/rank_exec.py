"""Assign one visible CUDA device to each local OpenMPI rank and execute ABACUS."""

import argparse
import json
import os
from pathlib import Path
import socket


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--gpu-uuids', required=True)
    args = parser.parse_args()
    devices = args.gpu_uuids.split(',')
    rank = int(os.environ['OMPI_COMM_WORLD_RANK'])
    local_rank = int(os.environ['OMPI_COMM_WORLD_LOCAL_RANK'])
    if local_rank >= len(devices):
        parser.error(f'Local rank {local_rank} has no assigned GPU among {len(devices)} devices')
    os.environ['CUDA_VISIBLE_DEVICES'] = devices[local_rank]
    record = {'rank': rank, 'local_rank': local_rank, 'pid': os.getpid(),
              'hostname': socket.gethostname(), 'gpu_uuid': devices[local_rank],
              'cpu_affinity': sorted(os.sched_getaffinity(0)),
              'threads': int(os.environ['OMP_NUM_THREADS'])}
    with Path(f'rank-{rank}.json').open('x') as stream:
        json.dump(record, stream, indent=2)
        stream.write('\n')
    os.execv(str(args.binary), [str(args.binary)])


if __name__ == '__main__':
    main()
