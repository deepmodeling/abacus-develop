"""Attribute CUDA activity to the NVTX range containing its host launch."""

import argparse
import collections
import heapq
import json
import sqlite3
from pathlib import Path


def union_seconds(intervals: list[tuple[int, int]]) -> float:
    total, previous_end = 0, 0
    for start, end in sorted(intervals):
        total += max(0, end - max(start, previous_end))
        previous_end = max(previous_end, end)
    return total / 1e9


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("database", type=Path)
    args = parser.parse_args()
    connection = sqlite3.connect(f"file:{args.database.resolve()}?mode=ro", uri=True)
    strings = dict(connection.execute("SELECT id, value FROM StringIds"))
    ranges = collections.defaultdict(list)
    phase_windows = []
    for start, end, tid, label in connection.execute(
        "SELECT start, end, globalTid, coalesce(text, s.value) FROM NVTX_EVENTS n "
        "LEFT JOIN StringIds s ON n.textId=s.id WHERE end IS NOT NULL ORDER BY start"
    ):
        ranges[tid].append((start, end, label))
        if label == "ESolver_KS_PW/hamilt2density":
            phase_windows.append((start, end))
    if len(phase_windows) < 3:
        raise ValueError("At least three hamilt2density ranges are required for the steady-state view")
    steady_start = phase_windows[2][0]
    steady_end = phase_windows[-1][1]
    positions = collections.defaultdict(int)
    active = collections.defaultdict(dict)
    ends = collections.defaultdict(list)
    launch_info = {}
    api_times = collections.Counter()
    for start, end, tid, correlation, name_id in connection.execute(
        "SELECT start, end, globalTid, correlationId, nameId "
        "FROM CUPTI_ACTIVITY_KIND_RUNTIME ORDER BY globalTid, start"
    ):
        while positions[tid] < len(ranges[tid]) and ranges[tid][positions[tid]][0] <= start:
            index = positions[tid]
            entry = ranges[tid][index]
            active[tid][index] = entry
            heapq.heappush(ends[tid], (entry[1], index))
            positions[tid] += 1
        while ends[tid] and ends[tid][0][0] <= start:
            _, index = heapq.heappop(ends[tid])
            active[tid].pop(index)
        ancestors = tuple(value[2] for value in sorted(active[tid].values()))
        label = ancestors[-1] if ancestors else "outside_marked_ranges"
        steady = steady_start <= start < steady_end
        launch_info[(tid & ~((1 << 24) - 1), correlation)] = (label, ancestors, steady)
        api_times[("all", label, strings[name_id])] += end - start
        if steady:
            api_times[("steady", label, strings[name_id])] += end - start

    by_range = collections.Counter()
    inclusive = collections.Counter()
    kernels = collections.Counter()
    transfers = collections.Counter()
    transfer_counts = collections.Counter()
    intervals = collections.defaultdict(list)
    unknown = collections.Counter()
    for kind, table, extra in (
        ("kernel", "CUPTI_ACTIVITY_KIND_KERNEL", "demangledName, 0"),
        ("copy", "CUPTI_ACTIVITY_KIND_MEMCPY", "copyKind, bytes"),
        ("memset", "CUPTI_ACTIVITY_KIND_MEMSET", "value, bytes"),
    ):
        for start, end, pid, correlation, device, detail, size in connection.execute(
            f"SELECT start, end, globalPid, correlationId, deviceId, {extra} FROM {table}"
        ):
            info = launch_info.get((pid, correlation))
            if info is None:
                unknown[kind] += 1
                info = ("uncorrelated", (), False)
            label, ancestors, steady = info
            for scope in (("all", "steady") if steady else ("all",)):
                by_range[(scope, label, kind)] += end - start
                for ancestor in ancestors:
                    inclusive[(scope, ancestor, kind)] += end - start
                intervals[(scope, device)].append((start, end))
                if kind == "kernel":
                    kernels[(scope, strings[detail])] += end - start
                if kind == "copy":
                    transfers[(scope, label, detail, size)] += end - start
                    transfer_counts[(scope, label, detail, size)] += 1

    def rows(counter: collections.Counter, names: tuple[str, ...], limit: int = 40) -> list[dict]:
        result = []
        for scope in ("all", "steady"):
            selected = [(key, value) for key, value in counter.items() if key[0] == scope]
            for key, value in sorted(selected, key=lambda item: -item[1])[:limit]:
                result.append({**dict(zip(names, key)), "seconds": value / 1e9})
        return result

    output = {
        "scope": "Steady excludes the first two hamilt2density calls and final forces/stress. "
                 "GPU activity is assigned by CUDA launch correlation; API times include waits. "
                 "Inclusive range times overlap and must not be summed.",
        "hamilt2density_calls": len(phase_windows),
        "steady_host_window_seconds": (steady_end - steady_start) / 1e9,
        "gpu_activity_union_seconds": [
            {"scope": scope, "device": device, "seconds": union_seconds(events)}
            for (scope, device), events in intervals.items()
        ],
        "uncorrelated_event_counts": dict(unknown),
        "exclusive_range_activity": rows(by_range, ("scope", "range", "kind")),
        "inclusive_range_activity": rows(inclusive, ("scope", "range", "kind")),
        "kernels": rows(kernels, ("scope", "kernel"), 20),
        "api_times": rows(api_times, ("scope", "range", "api"), 25),
        "transfers": rows(transfers, ("scope", "range", "copy_kind", "bytes_per_copy"), 25),
    }
    for entry in output["transfers"]:
        key = tuple(entry[name] for name in ("scope", "range", "copy_kind", "bytes_per_copy"))
        entry["count"] = transfer_counts[key]
    destination = args.database.with_name("launch-attribution.json")
    destination.write_text(json.dumps(output, indent=2) + "\n")
    print(destination)


if __name__ == "__main__":
    main()
