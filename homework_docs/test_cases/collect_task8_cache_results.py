#!/usr/bin/env python3
"""Collect ABACUS Task 8 cache benchmark logs into CSV and Markdown summaries."""

from __future__ import annotations

import argparse
import csv
import math
import re
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional


TIMER_FIELD_SPECS = [
    {
        "output": "pw_init_s",
        "candidates": [
            "PW_Basis::setuptransform",
            "PW_Basis_K::setuptransform",
            "PW_Basis_Sup::setuptransform",
            "setuptransform",
        ],
        "note": "missing pw_init timer",
    },
    {
        "output": "plane_wave_generation_s",
        "candidates": [
            "PW_Basis::collect_local_pw",
            "PW_Basis_K::collect_local_pw",
            "PW_Basis::collect_uniqgg",
            "PW_Basis::distributeg",
            "PW_Basis_K::distributeg",
            "collect_local_pw",
            "collect_uniqgg",
            "distributeg",
        ],
        "note": "missing plane-wave generation timer",
    },
    {
        "output": "gcar_or_gk2_s",
        "candidates": [
            "PW_Basis_K::setupIndGk",
            "setupIndGk",
        ],
        "note": "missing gcar/gk2 timer",
    },
    {
        "output": "real2recip_s",
        "candidates": [
            "PW_Basis::real2recip",
            "PW_Basis_K::real2recip",
            "real2recip",
        ],
        "note": "missing real2recip timer",
    },
    {
        "output": "recip2real_s",
        "candidates": [
            "PW_Basis::recip2real",
            "PW_Basis_K::recip2real",
            "recip2real",
        ],
        "note": "missing recip2real timer",
    },
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        action="append",
        dest="input_dirs",
        required=True,
        help="Benchmark output directory produced by run_task8_cache_benchmark.sh. Repeat this option to merge multiple runs.",
    )
    parser.add_argument(
        "--out-dir",
        required=True,
        help="Directory for raw_results.csv, summary_results.csv, summary_table.md.",
    )
    return parser.parse_args()


def load_key_value_file(path: Path) -> Dict[str, str]:
    data: Dict[str, str] = {}
    if not path.exists():
        return data
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or "=" not in line or line.endswith("<<EOF"):
                continue
            key, value = line.split("=", 1)
            data[key] = value
    return data


def parse_time_statistics(content: str) -> Dict[str, Dict[str, float]]:
    timers: Dict[str, Dict[str, float]] = {}
    if "TIME STATISTICS" not in content:
        return timers

    lines = content.splitlines()
    start = None
    for idx, line in enumerate(lines):
        if "TIME STATISTICS" in line:
            start = idx + 4
            break
    if start is None:
        return timers

    for line in lines[start:]:
        if not line.strip():
            continue
        if set(line.strip()) == {"-"}:
            break
        parts = line.split()
        if len(parts) < 5:
            continue
        if len(parts) == 5:
            class_name = ""
            name, time_s, calls, avg_s, per_pct = parts
        else:
            class_name = parts[0]
            name = parts[1]
            time_s = parts[2]
            calls = parts[3]
            avg_s = parts[4]
            per_pct = parts[5]
        key = f"{class_name}::{name}" if class_name else name
        try:
            timers[key] = {
                "time_s": float(time_s),
                "calls": float(calls),
                "avg_s": float(avg_s),
                "per_pct": float(per_pct),
            }
        except ValueError:
            continue
    return timers


def parse_wall_time(content: str) -> Optional[float]:
    match = re.search(r"TOTAL\s+Time\s*:\s*([0-9]+(?:\.[0-9]+)?)", content)
    if match:
        return float(match.group(1))
    return None


def parse_stdout(stdout_path: Path) -> Dict[str, object]:
    try:
        content = stdout_path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return {
            "wall_time_s": None,
            "has_time_statistics": False,
            "timers": {},
        }
    return {
        "wall_time_s": parse_wall_time(content),
        "has_time_statistics": "TIME STATISTICS" in content,
        "timers": parse_time_statistics(content),
    }


def timer_value_from_candidates(
    timers: Dict[str, Dict[str, float]],
    candidates: Iterable[str],
) -> Optional[float]:
    matched_values: List[float] = []
    for candidate in candidates:
        value = timers.get(candidate, {}).get("time_s")
        if isinstance(value, (int, float)):
            matched_values.append(float(value))
    if matched_values:
        return sum(matched_values)
    return None


def fmt_float(value: Optional[float], digits: int = 3) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "NA"
    return f"{value:.{digits}f}"


def collect_rows(input_dir: Path) -> List[Dict[str, object]]:
    env_info = load_key_value_file(input_dir / "environment_info.txt")
    manifest_path = input_dir / "run_manifest.csv"
    rows: List[Dict[str, object]] = []
    if not manifest_path.exists():
        return rows

    with manifest_path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            stdout_path = Path(row["stdout_path"])
            parsed = parse_stdout(stdout_path)
            timers = parsed["timers"]

            record: Dict[str, object] = dict(row)
            record["input_dir"] = str(input_dir)
            record["hostname"] = env_info.get("hostname", "NA")
            record["timestamp"] = env_info.get("timestamp", "NA")
            record["compiler_cxx"] = env_info.get("compiler_cxx", "NA")
            record["compiler_cxx_id"] = env_info.get("compiler_cxx_id", "NA")
            record["compiler_cxx_version"] = env_info.get("compiler_cxx_version", "NA")
            record["mpirun_version"] = env_info.get("mpirun_version", "NA")
            record["omp_proc_bind"] = env_info.get("omp_proc_bind", "NA")
            record["omp_places"] = env_info.get("omp_places", "NA")
            record["omp_schedule"] = env_info.get("omp_schedule", "NA")
            record["wall_time_s"] = parsed["wall_time_s"]
            record["has_time_statistics"] = "yes" if parsed["has_time_statistics"] else "no"

            for spec in TIMER_FIELD_SPECS:
                record[spec["output"]] = timer_value_from_candidates(timers, spec["candidates"])

            rows.append(record)
    return rows


def write_raw_results(rows: List[Dict[str, object]], out_path: Path) -> None:
    fieldnames = [
        "input_dir",
        "label",
        "case_name",
        "git_branch",
        "git_commit",
        "hostname",
        "timestamp",
        "nproc",
        "threads",
        "run_kind",
        "run_index",
        "exit_code",
        "start_time",
        "end_time",
        "stdout_path",
        "stderr_path",
        "out_abacus_path",
        "wall_time_s",
        "has_time_statistics",
        "pw_init_s",
        "plane_wave_generation_s",
        "gcar_or_gk2_s",
        "real2recip_s",
        "recip2real_s",
        "compiler_cxx",
        "compiler_cxx_id",
        "compiler_cxx_version",
        "mpirun_version",
        "omp_proc_bind",
        "omp_places",
        "omp_schedule",
    ]
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            formatted = {}
            for field in fieldnames:
                value = row.get(field)
                if isinstance(value, float):
                    formatted[field] = fmt_float(value)
                elif value is None or value == "":
                    formatted[field] = "NA"
                else:
                    formatted[field] = value
            writer.writerow(formatted)


def summarize_rows(rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    grouped: Dict[tuple, List[Dict[str, object]]] = defaultdict(list)
    for row in rows:
        if row.get("run_kind") != "repeat":
            continue
        grouped[
            (
                row.get("label", "NA"),
                row.get("case_name", "NA"),
                row.get("git_commit", "NA"),
                row.get("hostname", "NA"),
                row.get("nproc", "NA"),
                row.get("threads", "NA"),
            )
        ].append(row)

    summary_rows: List[Dict[str, object]] = []
    for key, items in sorted(grouped.items()):
        label, case_name, git_commit, hostname, nproc, threads = key
        wall_times = [
            row["wall_time_s"]
            for row in items
            if isinstance(row.get("wall_time_s"), (int, float))
        ]
        timer_notes: List[str] = []
        timer_stats: Dict[str, Optional[float]] = {}
        for spec in TIMER_FIELD_SPECS:
            field_name = spec["output"]
            values = [
                row[field_name]
                for row in items
                if isinstance(row.get(field_name), (int, float))
            ]
            if values:
                timer_stats[f"{field_name}_median"] = statistics.median(values)
            else:
                timer_stats[f"{field_name}_median"] = None
                timer_notes.append(spec["note"])

        summary_rows.append(
            {
                "label": label,
                "case_name": case_name,
                "git_commit": git_commit,
                "hostname": hostname,
                "nproc": nproc,
                "threads": threads,
                "repeat_count": len(items),
                "wall_times_s": ";".join(fmt_float(v) for v in wall_times) if wall_times else "NA",
                "median_time_s": statistics.median(wall_times) if wall_times else None,
                "mean_time_s": statistics.mean(wall_times) if wall_times else None,
                "min_time_s": min(wall_times) if wall_times else None,
                "max_time_s": max(wall_times) if wall_times else None,
                "successful_runs": sum(1 for row in items if str(row.get("exit_code")) == "0"),
                "time_statistics_runs": sum(1 for row in items if row.get("has_time_statistics") == "yes"),
                "notes": "; ".join(timer_notes) if timer_notes else "",
                **timer_stats,
            }
        )
    return summary_rows


def add_speedup(summary_rows: List[Dict[str, object]]) -> None:
    baseline_by_config: Dict[tuple, Dict[str, object]] = {}
    for row in summary_rows:
        if str(row.get("label")) != "baseline":
            continue
        config_key = (
            row.get("case_name"),
            row.get("nproc"),
            row.get("threads"),
        )
        baseline_by_config[config_key] = row

    for row in summary_rows:
        row["speedup_vs_baseline"] = None
        row["baseline_median_time_s"] = None
        config_key = (
            row.get("case_name"),
            row.get("nproc"),
            row.get("threads"),
        )
        baseline_row = baseline_by_config.get(config_key)
        current_median = row.get("median_time_s")
        if baseline_row and isinstance(current_median, (int, float)):
            baseline_median = baseline_row.get("median_time_s")
            row["baseline_median_time_s"] = baseline_median
            if isinstance(baseline_median, (int, float)) and current_median > 0:
                row["speedup_vs_baseline"] = baseline_median / current_median


def write_summary_csv(summary_rows: List[Dict[str, object]], out_path: Path) -> None:
    fieldnames = [
        "label",
        "case_name",
        "git_commit",
        "hostname",
        "nproc",
        "threads",
        "repeat_count",
        "successful_runs",
        "time_statistics_runs",
        "wall_times_s",
        "median_time_s",
        "mean_time_s",
        "min_time_s",
        "max_time_s",
        "baseline_median_time_s",
        "speedup_vs_baseline",
        "pw_init_s_median",
        "plane_wave_generation_s_median",
        "gcar_or_gk2_s_median",
        "real2recip_s_median",
        "recip2real_s_median",
        "notes",
    ]
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in summary_rows:
            formatted = {}
            for field in fieldnames:
                value = row.get(field)
                if isinstance(value, float):
                    formatted[field] = fmt_float(value)
                elif value is None or value == "":
                    formatted[field] = "NA"
                else:
                    formatted[field] = value
            writer.writerow(formatted)


def write_markdown(summary_rows: List[Dict[str, object]], out_path: Path) -> None:
    lines = [
        "| label | case | commit | np | omp | repeats | median(s) | mean(s) | min(s) | max(s) | speedup vs baseline | pw_init(s) | pw_gen(s) | gcar/gk2(s) | real2recip(s) | recip2real(s) | notes |",
        "| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in summary_rows:
        lines.append(
            "| {label} | {case_name} | {commit} | {nproc} | {threads} | {repeat_count} | {median} | {mean} | {min_t} | {max_t} | {speedup} | {pw_init} | {pw_gen} | {gcar} | {r2c} | {r2r} | {notes} |".format(
                label=row["label"],
                case_name=row["case_name"],
                commit=str(row["git_commit"])[:12],
                nproc=row["nproc"],
                threads=row["threads"],
                repeat_count=row["repeat_count"],
                median=fmt_float(row.get("median_time_s")),
                mean=fmt_float(row.get("mean_time_s")),
                min_t=fmt_float(row.get("min_time_s")),
                max_t=fmt_float(row.get("max_time_s")),
                speedup=fmt_float(row.get("speedup_vs_baseline")),
                pw_init=fmt_float(row.get("pw_init_s_median")),
                pw_gen=fmt_float(row.get("plane_wave_generation_s_median")),
                gcar=fmt_float(row.get("gcar_or_gk2_s_median")),
                r2c=fmt_float(row.get("real2recip_s_median")),
                r2r=fmt_float(row.get("recip2real_s_median")),
                notes=row.get("notes") or "NA",
            )
        )
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    all_rows: List[Dict[str, object]] = []
    for input_dir_str in args.input_dirs:
        input_dir = Path(input_dir_str).resolve()
        all_rows.extend(collect_rows(input_dir))

    raw_csv_path = out_dir / "raw_results.csv"
    summary_csv_path = out_dir / "summary_results.csv"
    markdown_path = out_dir / "summary_table.md"

    write_raw_results(all_rows, raw_csv_path)
    summary_rows = summarize_rows(all_rows)
    add_speedup(summary_rows)
    write_summary_csv(summary_rows, summary_csv_path)
    write_markdown(summary_rows, markdown_path)

    print(f"Wrote {raw_csv_path}")
    print(f"Wrote {summary_csv_path}")
    print(f"Wrote {markdown_path}")


if __name__ == "__main__":
    main()
