"""Attribute device-idle intervals to the main thread's active NVTX range."""

import argparse
import collections
import json
import sqlite3
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("database", type=Path)
    parser.add_argument("--global-pid", type=int, help="Encoded Nsight process ID for one MPI rank")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    connection = sqlite3.connect(f"file:{args.database.resolve()}?mode=ro", uri=True)
    process_parameters = () if args.global_pid is None else (~((1 << 24) - 1), args.global_pid)
    nvtx_filter = "" if args.global_pid is None else " AND (n.globalTid & ?) = ?"
    ranges = list(connection.execute(
        "SELECT start, end, globalTid, coalesce(text, s.value) FROM NVTX_EVENTS n "
        "LEFT JOIN StringIds s ON n.textId=s.id WHERE end IS NOT NULL" + nvtx_filter + " ORDER BY start",
        process_parameters,
    ))
    phases = [entry for entry in ranges if entry[3] == "ESolver_KS_PW/hamilt2density"]
    if len({row[2] & ~((1 << 24) - 1) for row in phases}) > 1:
        raise ValueError("Multiple application processes require --global-pid")
    if len(phases) < 3:
        raise ValueError("At least three hamilt2density ranges are required")
    start, stop, tid = phases[2][0], phases[-1][1], phases[2][2]
    events = []
    for table in ("KERNEL", "MEMCPY", "MEMSET"):
        device_filter = "deviceId=0" if args.global_pid is None else "globalPid=?"
        parameters = (start, stop) if args.global_pid is None else (args.global_pid, start, stop)
        for left, right in connection.execute(
            f"SELECT start,end FROM CUPTI_ACTIVITY_KIND_{table} WHERE {device_filter} AND end>? AND start<?",
            parameters,
        ):
            events.extend(((max(start, left), 1, None), (min(stop, right), -1, None)))
    labels = {}
    for index, (left, right, thread, label) in enumerate(ranges):
        if thread == tid and right > start and left < stop:
            labels[index] = (left, label)
            events.extend(((max(start, left), 0, index), (min(stop, right), 0, -index - 1)))
    events.sort(key=lambda event: event[0])
    active = set()
    gpu_count = 0
    previous = start
    idle = collections.Counter()
    for timestamp, change, identifier in events:
        if gpu_count == 0 and timestamp > previous:
            label = max((labels[index] for index in active), default=(0, "outside_marked_ranges"))[1]
            idle[label] += timestamp - previous
        if identifier is None:
            gpu_count += change
        elif identifier >= 0:
            active.add(identifier)
        else:
            active.remove(-identifier - 1)
        previous = timestamp
    output = {
        "global_pid": args.global_pid,
        "scope": "Main-thread NVTX context during GPU-idle gaps in the steady SCF interval. "
                 "Default covers device 0; --global-pid covers that process's device activity. Not CPU stack samples.",
        "window_seconds": (stop - start) / 1e9,
        "idle_seconds": sum(idle.values()) / 1e9,
        "ranges": [{"range": label, "seconds": duration / 1e9} for label, duration in idle.most_common()],
    }
    destination = args.output or args.database.with_name("idle-attribution.json")
    with destination.open("x") as stream:
        stream.write(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
