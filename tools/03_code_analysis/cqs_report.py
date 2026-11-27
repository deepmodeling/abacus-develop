"""Report data model and renderers for the ABACUS code-quality scorer."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import Dict, List, Optional

from cqs_config import PASS_THRESHOLD


@dataclass
class Finding:
    rule: str
    line: Optional[int]
    reason: str
    deduction: int


@dataclass
class FileReport:
    path: str
    score: int
    findings: List[Finding] = field(default_factory=list)
    real_score: int = 100

    def to_dict(self) -> dict:
        return {
            "path": self.path,
            "score": self.score,
            "real_score": self.real_score,
            "findings": [
                {
                    "rule": f.rule,
                    "line": f.line,
                    "reason": f.reason,
                    "deduction": f.deduction,
                }
                for f in self.findings
            ],
        }


def get_module(path: str) -> str:
    """Group a file path into its module for rollup purposes.

    Returns the `source/source_*` prefix when present (e.g.
    `source/source_lcao`); otherwise falls back to the file's parent
    directory. Backslashes are normalised to forward slashes.
    """
    parts = path.replace("\\", "/").split("/")
    for i, p in enumerate(parts):
        if (
            p == "source"
            and i + 1 < len(parts)
            and parts[i + 1].startswith("source_")
        ):
            return f"source/{parts[i + 1]}"
    return "/".join(parts[:-1]) if len(parts) > 1 else path


def compute_module_rollup(
    reports: List[FileReport],
) -> List[Tuple[str, int, float, int]]:
    """Group reports by module and return rows sorted by avg score ascending.

    Each row is (module, file_count, average_score, passing_count).
    Worst module (lowest average) appears first.
    """
    groups: dict = {}
    for r in reports:
        m = get_module(r.path)
        groups.setdefault(m, []).append(r)
    rows: List[Tuple[str, int, float, int]] = []
    for module, rs in groups.items():
        n = len(rs)
        avg = sum(r.score for r in rs) / n if n else 0.0
        passing = sum(1 for r in rs if r.score >= PASS_THRESHOLD)
        rows.append((module, n, avg, passing))
    rows.sort(key=lambda x: x[2])
    return rows


def render_text(reports: List[FileReport], min_score: Optional[int]) -> str:
    """Render reports as plain text."""
    if not reports:
        return "No files to analyze.\n"

    visible = reports if min_score is None else [r for r in reports if r.score <= min_score]
    sorted_reports = sorted(visible, key=lambda r: (r.score, r.real_score))

    out: List[str] = []
    width = 70
    out.append("=" * width)
    out.append("Code Quality Score Report")
    out.append("=" * width)

    rollup = compute_module_rollup(reports)
    if rollup:
        out.append("")
        out.append(
            "Per-module rollup (sorted by average score, worst first):"
        )
        out.append(f"  {'Module':<40} {'Files':>5} {'Avg':>6} {'Pass':>12}")
        out.append("-" * width)
        for module, n, avg, passing in rollup:
            module_disp = module if len(module) <= 40 else module[-37:] + "..."
            pass_str = f"{passing}/{n}"
            out.append(
                f"  {module_disp:<40} {n:>5} {avg:>6.1f} {pass_str:>12}"
            )
        out.append("")

    out.append(f"{'Score':>6}  {'File'}")
    out.append("-" * width)
    for r in sorted_reports:
        suffix = f"  (real_score={r.real_score})" if r.real_score < 0 else ""
        out.append(f"{r.score:>6}  {r.path}{suffix}")
    out.append("")

    for r in sorted_reports:
        if not r.findings:
            continue
        out.append("-" * width)
        out.append(f"File: {r.path}  (score: {r.score})")
        out.append("-" * width)
        for f in r.findings:
            loc = f"line {f.line}" if f.line else "file"
            out.append(f"  [{f.rule}] {loc}: {f.reason} (-{f.deduction})")
        out.append("")

    total = len(reports)
    visible_n = len(visible)
    avg = sum(r.score for r in reports) / total if total else 0.0
    passing = sum(1 for r in reports if r.score >= PASS_THRESHOLD)
    out.append("=" * width)
    out.append(f"Files scanned:        {total}")
    out.append(f"Files shown:          {visible_n}")
    out.append(f"Average score:        {avg:.1f}")
    out.append(f"Passing (>= {PASS_THRESHOLD}):  {passing}/{total}")
    out.append("=" * width)
    return "\n".join(out) + "\n"


def render_json(reports: List[FileReport], min_score: Optional[int]) -> str:
    visible = reports if min_score is None else [r for r in reports if r.score <= min_score]
    rollup = compute_module_rollup(reports)
    return json.dumps({
        "summary": {
            "total_scanned": len(reports),
            "total_shown": len(visible),
            "average_score": (
                sum(r.score for r in reports) / len(reports) if reports else 0.0
            ),
            "passing": sum(1 for r in reports if r.score >= PASS_THRESHOLD),
            "pass_threshold": PASS_THRESHOLD,
            "by_module": [
                {
                    "module": module,
                    "file_count": n,
                    "average_score": avg,
                    "passing": passing,
                }
                for module, n, avg, passing in rollup
            ],
        },
        "files": [r.to_dict() for r in sorted(visible, key=lambda r: (r.score, r.real_score))],
    }, indent=2)

