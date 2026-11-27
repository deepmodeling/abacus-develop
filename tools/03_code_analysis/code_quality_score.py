#!/usr/bin/env python3
"""ABACUS code quality scoring tool.

Scans source files/directories and assigns each file a quality score
starting from 100, deducting points for each rule violation found.

General semantics:
  - Score is clamped to a minimum of 0 (a file cannot go negative even
    when post-C++11 + .hpp alone would force it below 0). The
    unclamped value is preserved as `real_score`: shown in the text
    report as `(real_score=N)` next to clamped-0 files, included in
    JSON output, and used to break ties among clamped-0 files.
  - All caps listed below are per-file limits on the cumulative
    deduction for that rule within a single file.
  - Files under any of these directories are skipped entirely: test/,
    tests/, test_serial/, test_parallel/, test_gpu/, unit_test/, unittest/,
    unittests/.
  - Pass threshold: score >= 60. Output is sorted by score ascending
    (worst files first); ties are broken by `real_score` ascending
    (more negative = more total deduction = ranked first). Text and
    JSON output also include a per-module rollup grouped by the
    `source/source_*` prefix (worst module first). Written to
    ./code_quality_score.txt by default when format is text and no -o
    is given.

Rules (per-file score starts at 100), grouped by concern and ordered
by per-rule severity within each group:

  1. Architecture & compliance (AGENTS.md red lines):
    - post-C++11 feature usage: -40 base + -8 per distinct feature
      (cap 5). Total max -80, matching the previous one-shot weight,
      but a single misuse no longer swamps the whole file. High weight
      because it violates the repo-wide C++11 baseline contract
      (AGENTS.md §7). Detects high-confidence C++14/17/20/23 tokens
      such as `std::make_unique`, `if constexpr`, `[[nodiscard]]`,
      `auto [...]` structured bindings, `concept`, `requires`,
      `consteval`, `co_await`, `std::optional`, `std::variant`,
      `std::any`, `std::span`, `std::expected`, `std::format`, etc.
    - filename extension is .hpp: -50 (AGENTS.md §4)
    - GlobalV::/GlobalC::/PARAM.* cross-layer dependency: -3 per
      occurrence (cap 10) (AGENTS.md §1)
    - function declaration with default parameter: -2 per occurrence
      (cap 5) (AGENTS.md §5). Detects multi-line declarations where
      the parameter list spans lines.
    - #include of .hpp implementation header: -2 per occurrence
      (cap 5) (AGENTS.md §3, §4)
    - function call passing a ternary `?:`, `&&` or `||` expression as
      a top-level argument: -1 per offending call (cap 5)
      (AGENTS.md §14). Assign the value to a named local variable
      first, then pass it. Operators nested inside an inner call's
      parentheses are judged on the inner call; declarations,
      definitions, control-statement headers and lambdas are excluded.
    - each public member variable in a class/struct: -1
    - each static member variable in a class/struct: -1 (cap 10)
      Excludes static_assert, static member functions, and static_cast.
    - `friend` keyword exposing internals: -1 per occurrence (cap 5).
      Excludes friend declarations to standard-library types (e.g.
      `friend class std::hash<T>`, `friend std::ostream& operator<<`).

  2. Size, complexity & memory (implementation-level maintainability):
    - function with more than 7 parameters: -1 per extra param
      (cap 30 per function; multiple over-parameterised functions in
      the same file accumulate)
    - function cyclomatic complexity > 10 (if/for/while/switch/case/
      &&/||): -1 per extra point (cap 30 per function; multiple
      complex functions in the same file accumulate)
    - file longer than 500 lines: -2 per additional 50-line block
      (counts all physical lines including blanks and comments; no cap)
    - .h header longer than 120 lines: -1 per additional 10-line block
      (cap 5; in addition to the generic 500-line file rule above)
    - non-inline, non-template function body in a .h header longer than
      20 lines: -1 per additional 10-line block (cap 5; `inline` and
      `template <...>` functions are exempt, short accessor bodies at or
      below 20 lines are exempt by construction)
    - member function longer than 50 lines: -1 per additional 50-line
      block (counts physical lines after comments are blanked, blanks
      preserved)
    - each `new` keyword usage: -1 per occurrence (no cap)
    - unpaired `new` without matching `delete`: -1 per occurrence
      (cap 5). NOTE: stacks with the `new` rule above, so an unpaired
      `new` costs -2 total (-1 from `new` + -1 from unpaired).
    - local variable shadowing a member variable: -1 per occurrence
      (cap 5)

  3. Naming & formatting (lightweight surface rules, regex-based):
    - filename stem longer than 20 chars: -1 (stem = filename without
      extension, e.g. `matrix_orbs11.cpp` -> stem `matrix_orbs11`)
    - filename contains uppercase letters: -1
    - UPPERCASE constant naming (>3 chars all caps): -1 per occurrence
      (cap 5). Excludes string literals and member accesses such as
      `X::UPPER` or `x.UPPER` (those are governed by the
      cross-layer-dependency rule or another file's class definition).
    - tab indentation: -1 per line (cap 5)
    - `using namespace <any>;` directive: -1 per occurrence (cap 5)
    - `goto` keyword: -2 per occurrence (cap 5)
    - `auto` keyword: -1 per occurrence (cap 10). Counted after
      comments and string literals are stripped; `auto_ptr` does not
      match. Prefer explicit types for readability.
    - line longer than 120 chars: -1 per line (cap 5)
    - Chinese characters in comments/code: -1 per line (cap 5)
    - preprocessor directive not at column 0 (indented #if/#include/...):
      -1 per line (cap 5)
    - `#if`/`#ifdef`/`#ifndef` block inside a function call's argument
      list (conditional compilation wrapping individual arguments):
      -1 per opening directive (cap 5). Select the argument (or the
      overload) before the call instead of wrapping arguments in the
      preprocessor.
    - comment block above the first include guard of a .h file: -1 per
      contiguous comment block (a file-header is one block, so typically -1)
    - `@file` name does not match the actual filename: -1
    - comment text duplicated across comment blocks (same >=40-char line
      appearing in two or more blocks): -1 per duplicated line (cap 5)
    - new `#define` macro in a .h/.cpp file: -3 per definition (cap 10).
      Include guards (macro name ending in `_H`/`_HPP` with no parameter
      list or replacement text) and `#ifdef`/`#ifndef` usages of
      existing macros are exempt (AGENTS.md rule 20).
    - member variable initialized from `PARAM` in a .h/.hpp file: -5 per
      declaration (cap 10). Pass the value through the constructor or a
      setter instead (AGENTS.md rule 21).

Usage:
    python3 code_quality_score.py source/source_base
    # writes report to ./code_quality_score.txt (text default output file)
    python3 code_quality_score.py source/ -o report.txt
    python3 code_quality_score.py --format json path/to/file.cpp
    python3 code_quality_score.py --format json source/ -o report.json
    python3 code_quality_score.py --min-score 80 source/
"""

"""Entry-point module: CLI parsing and orchestration only.

Implementation was split into:
  cqs_config.py  - weights/caps/thresholds/regexes/prefix helpers
  cqs_text.py    - comment/string stripping, comment-block extraction
  cqs_rules.py   - rule-detection (find_*/analyze_*) functions
  cqs_analyze.py - per-file analysis orchestration, file discovery
  cqs_report.py  - Finding/FileReport dataclasses, text/JSON renderers
"""

import argparse
import sys
from pathlib import Path
from typing import Optional, Sequence

from cqs_analyze import analyze_file, discover_files
from cqs_config import PASS_THRESHOLD
from cqs_report import render_json, render_text


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="ABACUS code quality scoring tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("paths", nargs="+", help="file or directory paths to scan")
    parser.add_argument(
        "--format", choices=["text", "json"], default="text",
        help="output format (default: text)",
    )
    parser.add_argument(
        "--min-score", type=int, default=None,
        help="only show files with score <= this value in output",
    )
    parser.add_argument(
        "--output", "-o", default=None,
        help="write report to this file instead of stdout "
             "(default: code_quality_score.txt when format is text and "
             "no explicit output is given)",
    )
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    files = discover_files(args.paths)
    if not files:
        print(f"No source files found in: {args.paths}", file=sys.stderr)
        return 1

    reports = [analyze_file(f) for f in files]

    if args.format == "json":
        rendered = render_json(reports, args.min_score)
    else:
        rendered = render_text(reports, args.min_score)

    out_path = args.output
    if out_path is None and args.format == "text":
        out_path = "code_quality_score.txt"

    if out_path:
        Path(out_path).write_text(rendered, encoding="utf-8")
        print(
            f"Report written to {out_path} "
            f"({len(reports)} files scanned, "
            f"{sum(1 for r in reports if r.score >= PASS_THRESHOLD)} passing)",
            file=sys.stderr,
        )
    else:
        print(rendered)
    return 0


if __name__ == "__main__":
    sys.exit(main())
