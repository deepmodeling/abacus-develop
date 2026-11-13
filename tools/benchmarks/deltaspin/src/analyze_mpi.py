"""Separate MPI call durations from late arrivals in a multi-rank Nsight trace."""

import argparse
import collections
import json
from pathlib import Path
import sqlite3


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('database', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    connection = sqlite3.connect(f'file:{args.database.resolve()}?mode=ro', uri=True)
    process_mask = ~((1 << 24) - 1)
    ranks = {tid & process_mask: rank for tid, rank in connection.execute('SELECT globalTid,rank FROM MPI_RANKS')}
    if len(ranks) < 2:
        raise ValueError('At least two application MPI ranks are required')
    phases = collections.defaultdict(list)
    for start, end, tid in connection.execute(
        'SELECT n.start,n.end,n.globalTid FROM NVTX_EVENTS n LEFT JOIN StringIds s ON n.textId=s.id '
        "WHERE coalesce(n.text,s.value)='ESolver_KS_PW/hamilt2density' AND n.end IS NOT NULL ORDER BY n.start"):
        phases[ranks[tid & process_mask]].append((start, end))
    if any(len(phases[rank]) < 3 for rank in ranks.values()):
        raise ValueError('Every rank requires at least three SCF ranges')
    windows = {rank: (values[2][0], values[-1][1]) for rank, values in phases.items()}
    worlds = {}
    comm_sizes = {}
    for rank, handle, parent, size in connection.execute(
            'SELECT rank,commHandle,parentHandle,size FROM MPI_COMMUNICATORS ORDER BY timestamp'):
        if (rank, handle) in comm_sizes and comm_sizes[(rank, handle)] != size:
            raise ValueError('Reused communicator handles with changed sizes require time-dependent analysis')
        comm_sizes[(rank, handle)] = size
        if parent == 0 and size == len(ranks):
            if rank in worlds and worlds[rank] != handle:
                raise ValueError(f'Ambiguous MPI world communicator for rank {rank}')
            worlds[rank] = handle
    if set(worlds) != set(ranks.values()):
        raise ValueError('MPI world communicator metadata is incomplete')
    counts, durations = collections.Counter(), collections.Counter()
    world_calls = collections.defaultdict(lambda: collections.defaultdict(list))
    for start, end, tid, name, handle, size, received in connection.execute(
        'SELECT m.start,m.end,m.globalTid,s.value,m.commHandle,m.size,m.recvSize '
        'FROM MPI_COLLECTIVES_EVENTS m JOIN StringIds s ON m.textId=s.id '
        'WHERE m.end IS NOT NULL ORDER BY m.start'):
        rank = ranks[tid & process_mask]
        steady = windows[rank][0] <= start and end <= windows[rank][1]
        comm_size = comm_sizes[(rank, handle)]
        for scope in (('all', 'steady') if steady else ('all',)):
            key = scope, rank, name, comm_size
            counts[key] += 1
            durations[key] += end - start
        if handle == worlds[rank] and name in ('MPI_Allgather', 'MPI_Allreduce'):
            world_calls[name][rank].append((start, end, size, received, steady))
    alignment = []
    for name, by_rank in world_calls.items():
        lengths = {len(values) for values in by_rank.values()}
        if set(by_rank) != set(worlds) or len(lengths) != 1:
            raise ValueError(f'Unmatched MPI world collective counts: {name}')
        selected = 0
        tail = spread = 0
        before_last_entry = collections.Counter()
        total_by_rank = collections.Counter()
        for group in zip(*(by_rank[rank] for rank in sorted(by_rank))):
            if len({row[2:4] for row in group}) != 1:
                raise ValueError(f'Collective order or traced payload mismatch: {name}')
            if not all(row[4] for row in group):
                continue
            selected += 1
            last_entry = max(row[0] for row in group)
            tail += max(row[1] for row in group) - last_entry
            spread += last_entry - min(row[0] for row in group)
            for rank, row in zip(sorted(by_rank), group):
                total_by_rank[rank] += row[1] - row[0]
                before_last_entry[rank] += max(0, min(row[1], last_entry) - row[0])
        alignment.append({'name': name, 'matched_steady_calls_per_rank': selected,
                          'sum_time_after_last_entry_seconds': tail / 1e9,
                          'sum_entry_spread_seconds': spread / 1e9,
                          'call_seconds_per_rank': {rank: value / 1e9 for rank, value in total_by_rank.items()},
                          'seconds_before_last_rank_entry': {rank: value / 1e9 for rank, value in before_last_entry.items()}})
    rows = [{'scope': key[0], 'rank': key[1], 'name': key[2], 'communicator_size': key[3],
             'calls': counts[key], 'seconds': value / 1e9} for key, value in durations.items()]
    result = {'scope': 'Per-rank collective call times include waiting for other ranks. '
                       'World collectives align by operation order and matching trace payloads. '
                       'Time after the last observed entry includes MPI progress and software overhead; '
                       'it is not pure network latency or a prediction of recoverable wall time. '
                       'Steady excludes the first two SCF ranges and final forces/stress.',
              'rank_processes': {rank: process for process, rank in ranks.items()},
              'calls': sorted(rows, key=lambda row: (row['scope'], row['rank'], -row['seconds'])),
              'world_collective_alignment': alignment}
    destination = args.output or args.database.with_name('mpi-attribution.json')
    with destination.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps(alignment, indent=2))


if __name__ == '__main__':
    main()
