#!/usr/bin/env python3
"""Generate the Phase 4 Markdown report from benchmark CSV files."""

import argparse
import csv
import json
import re
import statistics
from collections import defaultdict
from pathlib import Path


def read_csv(path):
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def median(rows, field):
    values = [float(row[field]) for row in rows if row.get(field) not in ("", None, "None")]
    return statistics.median(values) if values else None


def format_number(value, digits=3):
    return "N/A" if value is None else f"{value:.{digits}f}"


def grouped(rows, fields):
    result = defaultdict(list)
    for row in rows:
        result[tuple(row[field] for field in fields)].append(row)
    return result


def grid_path_section(rows):
    lines = [
        "## 三路径核心搜索对比",
        "",
        "| Cell | Radius | Atoms | beta.1 legacy total (ms) | Current Full27 total (ms) | Current Half14 total (ms) | beta.1/Half14 | Full27/Half14 | Pair hash |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    groups = grouped(
        [row for row in rows if row.get("status") == "ok"],
        ["cell", "radius", "atoms", "implementation"],
    )
    keys = sorted(
        {(row["cell"], row["radius"], row["atoms"]) for row in rows},
        key=lambda key: (key[0], float(key[1]), int(key[2])),
    )
    for cell, radius, atoms in keys:
        beta = groups.get((cell, radius, atoms, "beta1_legacy"), [])
        full = groups.get((cell, radius, atoms, "current_full27"), [])
        half = groups.get((cell, radius, atoms, "current_half14"), [])
        beta_time = median(beta, "total_us")
        full_time = median(full, "total_us")
        half_time = median(half, "total_us")
        beta_ratio = beta_time / half_time if beta_time and half_time else None
        full_ratio = full_time / half_time if full_time and half_time else None
        hashes = {row["pair_hash"] for row in beta + full + half}
        hash_status = "一致" if len(hashes) == 1 else "不一致"
        if not beta:
            hash_status += "（beta.1 大规模跳过）"
        lines.append(
            f"| {cell} | {float(radius):.1f} | {int(atoms):,} | "
            f"{format_number(beta_time / 1000 if beta_time else None)} | "
            f"{format_number(full_time / 1000 if full_time else None)} | "
            f"{format_number(half_time / 1000 if half_time else None)} | "
            f"{format_number(beta_ratio)} | {format_number(full_ratio)} | {hash_status} |"
        )
    return lines


def version_section(rows):
    lines = [
        "## 版本整体对比",
        "",
        "| Case | Ranks | beta.1 中位时间 (s) | 当前版中位时间 (s) | 加速比 | beta.1 RSS (MiB) | 当前版 RSS (MiB) | 能量差 (eV) |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    groups = grouped(rows, ["case", "ranks", "version"])
    keys = sorted({(row["case"], row["ranks"]) for row in rows})
    for case, ranks in keys:
        beta = groups.get((case, ranks, "beta1"), [])
        current = groups.get((case, ranks, "current"), [])
        beta_time = median(beta, "wall_s")
        current_time = median(current, "wall_s")
        beta_rss = median(beta, "peak_rss_kib")
        current_rss = median(current, "peak_rss_kib")
        beta_energy = median(beta, "final_energy_ev")
        current_energy = median(current, "final_energy_ev")
        speedup = beta_time / current_time if beta_time and current_time else None
        energy_difference = abs(beta_energy - current_energy) if beta_energy is not None and current_energy is not None else None
        lines.append(
            f"| {case} | {ranks} | {format_number(beta_time)} | {format_number(current_time)} | "
            f"{format_number(speedup)} | {format_number(beta_rss / 1024 if beta_rss else None)} | "
            f"{format_number(current_rss / 1024 if current_rss else None)} | "
            f"{format_number(energy_difference, 8)} |"
        )
    return lines


def kernel_section(rows):
    lines = [
        "## 当前版本核心路径消融",
        "",
        "| Atoms | Half14 search (ms) | Full27 search (ms) | Full27/Half14 | Half14+reference (ms) | reference/Half14 | Half14 RSS (MiB) | reference RSS (MiB) |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    groups = grouped(rows, ["atoms", "mode"])
    for atoms in sorted({row["atoms"] for row in rows}, key=int):
        half = groups.get((atoms, "half14"), [])
        full = groups.get((atoms, "full27"), [])
        reference = groups.get((atoms, "half14_ref"), [])
        half_time = median(half, "search_us")
        full_time = median(full, "search_us")
        ref_time = median(reference, "search_us")
        ratio = full_time / half_time if half_time and full_time else None
        reference_ratio = ref_time / half_time if half_time and ref_time else None
        half_rss = median(half, "peak_rss_kib")
        reference_rss = median(reference, "peak_rss_kib")
        lines.append(
            f"| {int(atoms):,} | {format_number(half_time / 1000 if half_time else None)} | "
            f"{format_number(full_time / 1000 if full_time else None)} | {format_number(ratio)} | "
            f"{format_number(ref_time / 1000 if ref_time else None)} | {format_number(reference_ratio)} | "
            f"{format_number(half_rss / 1024 if half_rss else None)} | "
            f"{format_number(reference_rss / 1024 if reference_rss else None)} |"
        )
    return lines


def mpi_section(rows):
    lines = ["## MPI 单节点扩展", ""]
    groups = grouped(rows, ["scaling", "cell", "ranks"])
    for scaling in ("strong", "weak"):
        for cell in ("orthogonal", "triclinic"):
            keys = [key for key in groups if key[0] == scaling and key[1] == cell]
            if not keys:
                continue
            title = "强扩展" if scaling == "strong" else "弱扩展"
            lines.extend(
                [
                    f"### {cell} {title}",
                    "",
                    "| Ranks | 总原子数 | Exchange max (ms) | Search max (ms) | 加速比/弱扩展比 | 效率 | Ghost 总数 | Payload sent 总数 |",
                    "|---:|---:|---:|---:|---:|---:|---:|---:|",
                ]
            )
            base_time = None
            for _, _, ranks in sorted(keys, key=lambda key: int(key[2])):
                values = groups[(scaling, cell, ranks)]
                exchange = median(values, "max_exchange_us")
                search = median(values, "max_search_us")
                total = (exchange or 0.0) + (search or 0.0)
                if int(ranks) == 1:
                    base_time = total
                ratio = base_time / total if base_time and total else None
                efficiency = ratio / int(ranks) if scaling == "strong" and ratio is not None else ratio
                ghosts = median(values, "sum_ghost_atoms")
                payload = median(values, "sum_sent_payload")
                global_atoms = int(float(values[0]["global_atoms"]))
                lines.append(
                    f"| {ranks} | {global_atoms:,} | "
                    f"{format_number(exchange / 1000 if exchange is not None else None)} | "
                    f"{format_number(search / 1000 if search is not None else None)} | "
                    f"{format_number(ratio)} | {format_number(efficiency)} | "
                    f"{format_number(ghosts, 0)} | {format_number(payload, 0)} |"
                )
            lines.append("")
    return lines


def verlet_section(rows):
    lines = [
        "## Verlet 增量重建消融",
        "",
        "| Case | skin=0 时间 (s) | skin=3 时间 (s) | 时间比 (skin0/skin3) | skin=0 重建/复用 | skin=3 重建/复用 | 能量差 (eV) |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    groups = grouped(rows, ["case", "skin_bohr"])
    for case in sorted({row["case"] for row in rows}):
        baseline = groups.get((case, "0.0"), [])
        reused = groups.get((case, "3.0"), [])
        baseline_time = median(baseline, "total_time_s")
        reused_time = median(reused, "total_time_s")
        ratio = baseline_time / reused_time if baseline_time and reused_time else None
        baseline_energy = median(baseline, "final_energy_ev")
        reused_energy = median(reused, "final_energy_ev")
        difference = (
            abs(baseline_energy - reused_energy)
            if baseline_energy is not None and reused_energy is not None
            else None
        )
        baseline_counts = (
            f"{baseline[0]['rebuild_count']}/{baseline[0]['reuse_count']}" if baseline else "N/A"
        )
        reused_counts = f"{reused[0]['rebuild_count']}/{reused[0]['reuse_count']}" if reused else "N/A"
        lines.append(
            f"| {case} | {format_number(baseline_time)} | {format_number(reused_time)} | "
            f"{format_number(ratio)} | {baseline_counts} | {reused_counts} | "
            f"{format_number(difference, 8)} |"
        )
    return lines


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    results = args.results.resolve()
    environment = json.loads((results / "environment.json").read_text())
    version_rows = read_csv(results / "version.csv")
    grid_path_rows = read_csv(results / "grid_paths.csv")
    kernel_rows = read_csv(results / "kernel.csv")
    mpi_rows = read_csv(results / "mpi.csv")
    verlet_rows = read_csv(results / "verlet" / "neighbor_verlet_validation.csv")

    cpu_match = re.search(r"Model name:\s+(.+)", environment.get("lscpu", ""))
    core_match = re.search(r"Core\(s\) per socket:\s+(\d+)", environment.get("lscpu", ""))
    thread_match = re.search(r"CPU\(s\):\s+(\d+)", environment.get("lscpu", ""))
    memory_lines = environment.get("memory", "").splitlines()
    memory = memory_lines[1].split()[1] if len(memory_lines) > 1 and len(memory_lines[1].split()) > 1 else "unknown"

    lines = [
        "# Phase 4 近邻搜索性能 Benchmark 报告",
        "",
        "## 测试口径",
        "",
        f"- 基线：ABACUS `3.11.0-beta.1`，提交 `{environment['beta_commit']}`。",
        f"- 优化版：`final-project`，提交 `{environment['current_commit']}`。",
        f"- 优化版工作树：{'包含待提交修改' if environment.get('current_worktree_dirty') else '干净'}；"
        f"neighbor diff SHA-256 为 `{environment.get('current_neighbor_diff_sha256', 'N/A')}`。",
        "- 两个版本使用相同编译器、依赖、Release 配置、算例和 CPU 绑定。",
        "- 当前版在解决上游 PR 冲突时包含少量 develop 接口适配、测试依赖和 mock 链接修复；这些修改未改变近邻搜索核心算法。",
        "- 本报告为单节点测试，不据此声称多节点扩展效率。",
        "- 时间采用正式重复运行的中位数；提交版原始数据位于 `docs/benchmark_results/phase4/`。",
        "- beta.1 与当前 conda ScaLAPACK/MKL 组合存在 ABI 崩溃，因此整体版本对比统一使用单 rank LAPACK；MPI 由独立 neighbor probe 测量。",
        "",
        "## 环境与方法",
        "",
        f"- CPU：{cpu_match.group(1) if cpu_match else 'unknown'}，"
        f"{core_match.group(1) if core_match else 'unknown'} 核/"
        f"{thread_match.group(1) if thread_match else 'unknown'} 线程；内存约 {memory}。",
        "- 环境：WSL 单节点，GCC 14.3.0、OpenMPI 5.0.10、Conda `bxcx`。",
        "- 构建：Release、MPI/LCAO/OpenMP 开启，ELPA 和 native optimization 关闭。",
        "- kernel：4,096/32,768/131,072 原子，每个模式预热一次并正式运行 10 次。",
        "- 三路径：beta.1 legacy、current Full27、current Half14 使用相同 PBC 晶胞、坐标和 cutoff；beta.1 因 O(N²) 仅运行到 32,768 原子，其中 4,096 原子为 10 次正式重复、32,768 原子为 3 次正式重复。",
        "- MPI：固定 65,536 原子的强扩展，以及每 rank 4,096 原子的弱扩展；覆盖 1/2/4/8/16 ranks 和 orthogonal/triclinic，每组运行 5 次。",
        "- 版本对比：Si16/Si64 LCAO SCF，预热一次、正式运行 3 次，使用中位数。",
        "- Verlet：Si2 relax 5 步、Si8 MD 10 步和 Si64 MD 20 步，比较 `neighbor_skin=0.0/3.0`。",
        "",
        "## 结论摘要",
        "",
        "- 完整 SCF 用于验证端到端物理结果；其总耗时会被 Hamiltonian、求解器等主流程主导，不能单独代表近邻搜索热点收益。",
        "- 三路径核心搜索显示 beta.1 legacy 在 32K 原子时约需 97-104 秒，而当前 Half14 为 83-181 毫秒；131K 下 beta.1 因 O(N²) 成本不再实际运行。",
        "- 重构后的 Half14 在 4K/32K/131K 原子上分别比 Full27 快约 1.68/2.05/2.24 倍，pair hash 完全一致，峰值 RSS 也未增加。",
        "- 取消默认 Full27 reference 继续提供明确收益：131,072 原子时避免约 81 MiB 峰值 RSS 和额外 reference 构建。",
        "- Verlet 显著降低重建次数；Si8 小体系时间持平，Si64 20 步 MD 从 301 秒降至 297 秒，收益约 1.3%。",
        "- 单节点 MPI strong/weak 在 1/2/4 ranks 表现较好，8/16 ranks 受到 ghost 增长、通信和 WSL 单节点资源竞争影响。",
        "",
    ]
    if grid_path_rows:
        lines.extend(grid_path_section(grid_path_rows))
        lines.append("")
    if version_rows:
        lines.extend(version_section(version_rows))
        lines.append("")
    if kernel_rows:
        lines.extend(kernel_section(kernel_rows))
        lines.append("")
    if mpi_rows:
        lines.extend(mpi_section(mpi_rows))
        lines.append("")
    if verlet_rows:
        lines.extend(verlet_section(verlet_rows))
        lines.append("")
    lines.extend(
        [
            "## 回归测试",
            "",
            "- 完整 `MODULE_CELL` CTest 36/36 通过，包括 MPI wrapper。",
            "- `MODULE_CELL_NEIGHBOR_mpi_grid` 在 1/2/4/8/16 ranks 下通过。",
            "- `MODULE_LCAO_record_adj_mpi` 在 1/2/4 ranks 下通过。",
            "- no-MPI 配置下 `cell`、`neighbor`、`hamilt_lcao` 和 `esolver` 构建通过。",
            "",
            "## 正确性与限制",
            "",
            "- Half14、Full27 和显式 Full27 reference 的 pair hash 必须逐规模一致。",
            "- 版本对比仅在两个版本均正常结束且物理结果满足容差时计入。",
            "- 峰值 RSS 是完整独立进程的最大常驻集，适合版本/模式相对比较，不等同于容器理论容量。",
            "- WSL 单节点上的 MPI 数据主要用于观察算法和通信趋势，不能替代集群多节点 Benchmark。",
            "",
            "## 原始结果",
            "",
            "- `docs/benchmark_results/phase4/environment.json`：硬件、软件、命令参数和提交号。",
            "- `docs/benchmark_results/phase4/summary.json`：各配置的中位数、最小值和最大值。",
            "- `docs/benchmark_results/phase4/version.csv`：beta.1 与当前版完整 ABACUS 运行结果。",
            "- `docs/benchmark_results/phase4/grid_paths.csv`：beta.1 legacy、current Full27 与 current Half14 的直接对比。",
            "- `docs/benchmark_results/phase4/kernel.csv`：Half14、Full27 与 reference 消融结果。",
            "- `docs/benchmark_results/phase4/mpi.csv`：MPI local/ghost、payload 和耗时结果。",
            "- `docs/benchmark_results/phase4/verlet.csv`：skin=0/3 的真实 relax/MD 对照。",
            "- 完整标准输出保存在本地 Benchmark 工作目录的 `logs/`，不纳入仓库。",
        ]
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
