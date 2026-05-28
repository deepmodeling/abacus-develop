#!/usr/bin/env python3
"""Summarize ABACUS TIME STATISTICS by function.

The benchmark runner writes one stdout file per case/config/repetition.  This
parser scans every TIME STATISTICS table in those files and reports how much
time each timer entry used in total.
"""

import argparse
import os
import re
from collections import defaultdict


RUN_PATTERN = re.compile(r"np(\d+)_omp(\d+)_r(\d+)\.stdout$")
SEPARATOR_PATTERN = re.compile(r"^-+$")
TOTAL_TIME_PATTERN = re.compile(r"TOTAL  Time  : (\d+)")


def empty_timer():
    return {
        "time_s": 0.0,
        "calls": 0,
        "per_pct_sum": 0.0,
        "entries": 0,
        "runs": set(),
    }


def parse_timer_line(line):
    """Parse one timer table row.

    ABACUS rows are whitespace separated.  The first row is special because the
    CLASS_NAME column is empty:

        total              5.07   17       0.30   100.00
        PW_Basis_K         real2recip      0.75   12091    0.00   14.88
    """
    parts = line.strip().split()
    if len(parts) < 5:
        return None

    try:
        if len(parts) == 5:
            class_name = ""
            name = parts[0]
            time_s = float(parts[1])
            calls = int(parts[2])
            avg_s = float(parts[3])
            per_pct = float(parts[4])
        else:
            class_name = parts[0]
            name = parts[1]
            time_s = float(parts[2])
            calls = int(parts[3])
            avg_s = float(parts[4])
            per_pct = float(parts[5])
    except ValueError:
        return None

    key = f"{class_name}::{name}" if class_name else name
    return key, time_s, calls, avg_s, per_pct


def parse_timer_tables(stdout_path):
    """Return aggregated timer entries from all TIME STATISTICS tables in a file."""
    with open(stdout_path, "r", encoding="utf-8", errors="replace") as f:
        lines = f.read().splitlines()

    timers = defaultdict(empty_timer)
    table_count = 0
    idx = 0

    while idx < len(lines):
        if lines[idx].strip() != "TIME STATISTICS":
            idx += 1
            continue

        separators = []
        scan = idx + 1
        while scan < len(lines):
            if SEPARATOR_PATTERN.fullmatch(lines[scan].strip()):
                separators.append(scan)
                if len(separators) == 3:
                    break
            scan += 1

        if len(separators) < 3:
            idx += 1
            continue

        table_count += 1
        data_start = separators[1] + 1
        data_end = separators[2]
        for line in lines[data_start:data_end]:
            parsed = parse_timer_line(line)
            if parsed is None:
                continue
            key, time_s, calls, _avg_s, per_pct = parsed
            entry = timers[key]
            entry["time_s"] += time_s
            entry["calls"] += calls
            entry["per_pct_sum"] += per_pct
            entry["entries"] += 1

        idx = separators[2] + 1

    if not timers:
        return None

    with open(stdout_path, "r", encoding="utf-8", errors="replace") as f:
        content = f.read()

    total_match = TOTAL_TIME_PATTERN.search(content)
    wall_time = int(total_match.group(1)) if total_match else None
    scf_iters = len(re.findall(r"#ELEC ITER#", content))

    return {
        "timers": timers,
        "wall_time": wall_time,
        "scf_iters": scf_iters,
        "table_count": table_count,
    }


def iter_benchmark_runs(results_dir):
    for case_name in sorted(os.listdir(results_dir)):
        case_dir = os.path.join(results_dir, case_name)
        if not os.path.isdir(case_dir):
            continue

        for fname in sorted(os.listdir(case_dir)):
            match = RUN_PATTERN.match(fname)
            if not match:
                continue

            np_count = int(match.group(1))
            omp_count = int(match.group(2))
            rep = int(match.group(3))
            stdout_path = os.path.join(case_dir, fname)
            parsed = parse_timer_tables(stdout_path)
            if parsed is None:
                continue

            yield case_name, np_count, omp_count, rep, stdout_path, parsed


def aggregate_runs(runs):
    grouped = {}

    for case_name, np_count, omp_count, rep, stdout_path, parsed in runs:
        group_key = (case_name, np_count, omp_count)
        group = grouped.setdefault(
            group_key,
            {
                "runs": [],
                "timers": defaultdict(empty_timer),
                "wall_time_s": 0.0,
                "wall_time_count": 0,
                "scf_iters": 0,
                "table_count": 0,
            },
        )

        run_id = f"r{rep}:{os.path.basename(stdout_path)}"
        group["runs"].append(run_id)
        group["scf_iters"] += parsed["scf_iters"]
        group["table_count"] += parsed["table_count"]
        if parsed["wall_time"] is not None:
            group["wall_time_s"] += parsed["wall_time"]
            group["wall_time_count"] += 1

        for timer_name, timer in parsed["timers"].items():
            entry = group["timers"][timer_name]
            entry["time_s"] += timer["time_s"]
            entry["calls"] += timer["calls"]
            entry["per_pct_sum"] += timer["per_pct_sum"]
            entry["entries"] += timer["entries"]
            entry["runs"].add(run_id)

    return grouped


def format_table(group, top=None, filter_text=None):
    n_runs = len(group["runs"])
    rows = []

    for timer_name, timer in group["timers"].items():
        if filter_text and filter_text.lower() not in timer_name.lower():
            continue
        total_s = timer["time_s"]
        total_calls = timer["calls"]
        avg_run_s = total_s / n_runs if n_runs else 0.0
        avg_call_s = total_s / total_calls if total_calls else 0.0
        avg_pct = timer["per_pct_sum"] / timer["entries"] if timer["entries"] else 0.0
        rows.append((total_s, timer_name, total_calls, avg_run_s, avg_call_s, avg_pct, len(timer["runs"])))

    rows.sort(key=lambda item: item[0], reverse=True)
    if top is not None:
        rows = rows[:top]

    output = []
    output.append(
        f"  {'Timer':<42} {'Total/s':>12} {'Calls':>10} {'AvgRun/s':>12} "
        f"{'AvgCall/ms':>13} {'Avg%':>8} {'Runs':>6}"
    )
    output.append(f"  {'-' * 42} {'-' * 12} {'-' * 10} {'-' * 12} {'-' * 13} {'-' * 8} {'-' * 6}")

    for total_s, timer_name, total_calls, avg_run_s, avg_call_s, avg_pct, run_count in rows:
        display_name = timer_name
        if len(display_name) > 42:
            display_name = display_name[:41]
        output.append(
            f"  {display_name:<42} {total_s:>12.6f} {total_calls:>10} "
            f"{avg_run_s:>12.6f} {avg_call_s * 1000.0:>13.6f} {avg_pct:>7.3f}% {run_count:>6}"
        )

    return output


def main():
    parser = argparse.ArgumentParser(description="Summarize ABACUS TIME STATISTICS by timer function.")
    parser.add_argument("results_dir", nargs="?", default="results", help="benchmark result directory")
    parser.add_argument("--top", type=int, default=None, help="show only the N most expensive timers per config")
    parser.add_argument("--filter", default=None, help="show only timers containing this text")
    args = parser.parse_args()

    grouped = aggregate_runs(iter_benchmark_runs(args.results_dir))
    if not grouped:
        print(f"No benchmark TIME STATISTICS found under: {args.results_dir}")
        return

    for (case_name, np_count, omp_count), group in sorted(grouped.items()):
        n_runs = len(group["runs"])
        print("")
        print("=" * 96)
        print(f"Case: {case_name} | np={np_count} | omp={omp_count} | runs={n_runs}")
        print("=" * 96)
        print("Function timer totals across all repetitions in this case/config.")

        if group["wall_time_count"]:
            avg_wall = group["wall_time_s"] / group["wall_time_count"]
            print(
                f"Wall total/s={group['wall_time_s']:.3f}, "
                f"wall avg/run/s={avg_wall:.3f}, "
                f"scf_iters={group['scf_iters']}, "
                f"timer_tables={group['table_count']}"
            )
        else:
            print(f"scf_iters={group['scf_iters']}, timer_tables={group['table_count']}")

        for line in format_table(group, top=args.top, filter_text=args.filter):
            print(line)


if __name__ == "__main__":
    main()
