"""Attribute application ZGEMMs to device kernels using CUDA launch correlation."""

import argparse
import collections
import heapq
import json
from pathlib import Path
import sqlite3


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("database", type=Path)
    args = parser.parse_args()
    connection = sqlite3.connect(f"file:{args.database.resolve()}?mode=ro", uri=True)
    ranges = collections.defaultdict(list)
    phases = []
    for start, end, tid, label in connection.execute(
        "SELECT start, end, globalTid, coalesce(text, s.value) FROM NVTX_EVENTS n "
        "LEFT JOIN StringIds s ON n.textId=s.id WHERE end IS NOT NULL ORDER BY start"
    ):
        ranges[tid].append((start, end, label))
        if label == "ESolver_KS_PW/hamilt2density":
            phases.append((start, end))
    if len(phases) < 3:
        raise ValueError("At least three hamilt2density ranges are required")
    steady_start, steady_end = phases[2][0], phases[-1][1]
    positions = collections.defaultdict(int)
    active = collections.defaultdict(dict)
    ends = collections.defaultdict(list)
    launch_info = {}
    calls = collections.defaultdict(set)
    for start, tid, correlation in connection.execute(
        "SELECT start, globalTid, correlationId FROM CUPTI_ACTIVITY_KIND_RUNTIME ORDER BY globalTid, start"
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
        ancestors = sorted(active[tid].values())
        shape = next((entry for entry in reversed(ancestors) if entry[2].startswith("ZGEMM/")), None)
        if shape is None:
            continue
        parents = [entry[2] for entry in ancestors if entry[0] <= shape[0]
                   and entry[1] >= shape[1] and not entry[2].startswith("ZGEMM/")]
        parent = parents[-1] if parents else "outside_marked_ranges"
        parts = shape[2].split("/")
        key = (parent, parts[1], *(int(value) for value in parts[2:]))
        scopes = ("all", "steady") if steady_start <= shape[0] < steady_end else ("all",)
        launch_info[(tid & ~((1 << 24) - 1), correlation)] = (key, scopes)
        for scope in scopes:
            calls[(scope, *key)].add((tid, shape[0], shape[1]))
    durations = collections.Counter()
    kernel_counts = collections.Counter()
    for start, end, pid, correlation in connection.execute(
        "SELECT start, end, globalPid, correlationId FROM CUPTI_ACTIVITY_KIND_KERNEL"
    ):
        info = launch_info.get((pid, correlation))
        if info is None:
            continue
        key, scopes = info
        for scope in scopes:
            durations[(scope, *key)] += end - start
            kernel_counts[(scope, *key)] += 1
    names = ("scope", "parent", "transpose", "m", "n", "k", "lda", "ldb", "ldc")
    rows = []
    for key, identifiers in calls.items():
        row = dict(zip(names, key))
        count = len(identifiers)
        seconds = durations[key] / 1e9
        flops = 8 * row["m"] * row["n"] * row["k"] * count
        rows.append({**row, "calls": count, "gpu_kernels": kernel_counts[key],
                     "gpu_seconds": seconds, "nominal_flops": flops,
                     "tflops_per_second": flops / seconds / 1e12 if seconds else None})
    summaries = []
    for scope in ("all", "steady"):
        selected = [row for row in rows if row["scope"] == scope]
        seconds = sum(row["gpu_seconds"] for row in selected)
        flops = sum(row["nominal_flops"] for row in selected)
        summaries.append({"scope": scope, "calls": sum(row["calls"] for row in selected),
                          "gpu_seconds": seconds, "nominal_flops": flops,
                          "tflops_per_second": flops / seconds / 1e12 if seconds else None})
    result = {"scope": "Application double-complex GEMM wrapper only; excludes internal cuSolver GEMMs. "
                       "Nominal work is 8*m*n*k per call, excluding alpha/beta operations. "
                       "GPU duration is the sum of correlated kernels, not host API duration. "
                       "Steady excludes the first two hamilt2density calls and final forces/stress.",
              "summary": summaries,
              "shapes": sorted(rows, key=lambda row: (row["scope"], -row["gpu_seconds"]))}
    destination = args.database.with_name("gemm-attribution.json")
    destination.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(summaries))


if __name__ == "__main__":
    main()
