#!/usr/bin/env python3
"""Diff-oriented governance checks for ABACUS agent and PR review."""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


BLOCK = "error"
WARN = "warning"
INFO = "info"

TEXT_EXTENSIONS = {
    ".c",
    ".cc",
    ".cpp",
    ".cxx",
    ".cu",
    ".cuh",
    ".h",
    ".hh",
    ".hpp",
    ".hxx",
    ".md",
    ".rst",
    ".txt",
    ".yaml",
    ".yml",
    ".py",
    ".sh",
    ".cmake",
    ".json",
}
HEADER_EXTENSIONS = {".h", ".hh", ".hpp", ".hxx"}
SOURCE_EXTENSIONS = {".c", ".cc", ".cpp", ".cxx", ".cu"}
SOURCE_REVIEW_EXTENSIONS = SOURCE_EXTENSIONS | HEADER_EXTENSIONS | {".cuh"}
CODE_EXTENSIONS = SOURCE_EXTENSIONS | HEADER_EXTENSIONS | {".cuh", ".py", ".cmake"}
WINDOWS_SCRIPT_EXTENSIONS = {".bat", ".cmd"}
TEST_PATH_PREFIXES = ("tests/", "test/", "unit_test/", "examples/")
TEST_NAME_PATTERNS = ("test", "tests", "unittest", "pytest", "ctest", "case")
HETEROGENEOUS_MARKERS = ("/cuda/", "/rocm/", "/kernels/", "cuda/", "rocm/", "kernels/")
PR_SECTION_RE = re.compile(r"^###\s+(.+?)\s*$", re.MULTILINE)


@dataclass
class Finding:
    rule: str
    severity: str
    path: str
    line: Optional[int]
    reason: str
    suggestion: str
    allow_exception: bool


@dataclass
class DiffLine:
    path: str
    line: Optional[int]
    content: str


class GitError(RuntimeError):
    pass


def git(args: Sequence[str], cwd: Path, *, text: bool = True) -> subprocess.CompletedProcess:
    result = subprocess.run(
        ["git", *args],
        cwd=str(cwd),
        text=text,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if result.returncode != 0:
        raise GitError(result.stderr.strip() or "git command failed")
    return result


def repo_root() -> Path:
    result = subprocess.run(
        ["git", "rev-parse", "--show-toplevel"],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if result.returncode != 0:
        return Path.cwd()
    return Path(result.stdout.strip())


def parse_name_status(output: str) -> Tuple[Dict[str, str], List[str]]:
    statuses: Dict[str, str] = {}
    changed: List[str] = []
    for raw in output.splitlines():
        if not raw:
            continue
        parts = raw.split("\t")
        status = parts[0]
        path = parts[-1]
        statuses[path] = status
        changed.append(path)
    return statuses, changed


def changed_paths(root: Path, args: argparse.Namespace) -> Tuple[Dict[str, str], List[str]]:
    if args.staged:
        output = git(["diff", "--cached", "--name-status"], root).stdout
    elif args.base and args.head:
        output = git(["diff", "--name-status", args.base, args.head], root).stdout
    else:
        output = ""
    return parse_name_status(output)


def parse_changed_lines(diff_text: str) -> Tuple[List[DiffLine], List[DiffLine]]:
    added: List[DiffLine] = []
    removed: List[DiffLine] = []
    old_path = ""
    new_path = ""
    old_line: Optional[int] = None
    new_line: Optional[int] = None
    hunk_re = re.compile(r"@@ -(\d+)(?:,\d+)? \+(\d+)(?:,\d+)? @@")
    for raw in diff_text.splitlines():
        if raw.startswith("--- a/"):
            old_path = raw[6:]
            continue
        if raw.startswith("+++ b/"):
            new_path = raw[6:]
            continue
        if raw.startswith("--- ") or raw.startswith("+++ "):
            continue
        match = hunk_re.match(raw)
        if match:
            old_line = int(match.group(1))
            new_line = int(match.group(2))
            continue
        if old_line is None or new_line is None:
            continue
        if raw.startswith("\\"):
            continue
        if raw.startswith("+") and not raw.startswith("+++"):
            added.append(DiffLine(new_path, new_line, raw[1:]))
            new_line += 1
        elif raw.startswith("-") and not raw.startswith("---"):
            removed.append(DiffLine(old_path, old_line, raw[1:]))
            old_line += 1
        else:
            old_line += 1
            new_line += 1
    return added, removed


def changed_lines(root: Path, args: argparse.Namespace) -> Tuple[List[DiffLine], List[DiffLine]]:
    if args.staged:
        output = git(["diff", "--cached", "--ignore-cr-at-eol", "-U0"], root).stdout
    elif args.base and args.head:
        output = git(["diff", "--ignore-cr-at-eol", "-U0", args.base, args.head], root).stdout
    else:
        output = ""
    return parse_changed_lines(output)


def read_changed_file_bytes(root: Path, path: str, args: argparse.Namespace) -> bytes:
    if args.staged:
        result = subprocess.run(
            ["git", "show", f":{path}"],
            cwd=str(root),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if result.returncode == 0:
            return result.stdout
    if args.head:
        result = subprocess.run(
            ["git", "show", f"{args.head}:{path}"],
            cwd=str(root),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if result.returncode == 0:
            return result.stdout
    with open(root / path, "rb") as handle:
        return handle.read()


def is_text_path(path: str) -> bool:
    suffix = Path(path).suffix.lower()
    if suffix in WINDOWS_SCRIPT_EXTENSIONS:
        return False
    return suffix in TEXT_EXTENSIONS or path in {".gitattributes", ".gitignore"}


def has_crlf_line_endings(content: bytes) -> bool:
    return any(line.endswith(b"\r\n") for line in content.splitlines(keepends=True))


def add_finding(
    findings: List[Finding],
    rule: str,
    severity: str,
    path: str,
    line: Optional[int],
    reason: str,
    suggestion: str,
    allow_exception: bool = True,
) -> None:
    findings.append(Finding(rule, severity, path, line, reason, suggestion, allow_exception))


def check_line_endings(
    findings: List[Finding],
    root: Path,
    paths: Iterable[str],
    statuses: Dict[str, str],
    args: argparse.Namespace,
) -> None:
    for path in paths:
        if statuses.get(path, "").startswith("D") or not is_text_path(path):
            continue
        try:
            content = read_changed_file_bytes(root, path, args)
        except OSError:
            continue
        if has_crlf_line_endings(content):
            add_finding(
                findings,
                "LF line endings",
                BLOCK,
                path,
                None,
                "Changed text file contains CRLF line endings.",
                "Convert the file to LF. Windows .bat and .cmd scripts are the only CRLF exception.",
                allow_exception=False,
            )


GLOBAL_DEPENDENCY_RE = re.compile(r"\b(GlobalV::|GlobalC::|PARAM(?:\.|->|::|\b))")


def is_global_dependency_check_path(path: str) -> bool:
    if path.startswith("tools/03_code_analysis/"):
        return False
    return Path(path).suffix.lower() in CODE_EXTENSIONS


def global_dependency_hits(lines: Iterable[DiffLine]) -> List[Tuple[DiffLine, int]]:
    hits: List[Tuple[DiffLine, int]] = []
    for line in lines:
        if not is_global_dependency_check_path(line.path):
            continue
        count = len(GLOBAL_DEPENDENCY_RE.findall(line.content))
        if count:
            hits.append((line, count))
    return hits


def check_global_dependencies(
    findings: List[Finding],
    added_lines: Iterable[DiffLine],
    removed_lines: Iterable[DiffLine],
) -> None:
    added_hits = global_dependency_hits(added_lines)
    removed_hits = global_dependency_hits(removed_lines)
    added_count = sum(count for _, count in added_hits)
    removed_count = sum(count for _, count in removed_hits)
    delta = added_count - removed_count
    if added_count == 0:
        return

    severity = BLOCK if delta > 0 else WARN
    action = (
        "Reduce or explicitly pass dependencies so this PR does not increase global dependency usage."
        if delta > 0
        else "Confirm this is a migration-neutral move or partial cleanup, and explain the remaining global dependency rationale."
    )
    for line, count in added_hits:
        add_finding(
            findings,
            "Global dependency budget",
            severity,
            line.path,
            line.line,
            (
                f"Added line introduces {count} GlobalV/GlobalC/PARAM reference(s); "
                f"PR total added={added_count}, removed={removed_count}, net_delta={delta}."
            ),
            action,
        )


def check_default_parameters(findings: List[Finding], lines: Iterable[DiffLine]) -> None:
    default_arg = re.compile(r"[(,]\s*[^()=;,{}]+\b\w+\s*=\s*[^,);{}]+")
    control_flow = re.compile(r"^(for|if|while|switch|catch)\s*\(")
    for line in lines:
        if Path(line.path).suffix.lower() not in HEADER_EXTENSIONS:
            continue
        stripped = line.content.strip()
        if not stripped or stripped.startswith("//") or stripped.startswith("*"):
            continue
        if control_flow.match(stripped):
            continue
        if "(" in stripped and ")" in stripped and default_arg.search(stripped):
            add_finding(
                findings,
                "No new default parameters",
                BLOCK,
                line.path,
                line.line,
                "Header diff adds a function declaration with a default argument.",
                "Update call sites explicitly or introduce a clearer overload/configuration object.",
            )


def check_hpp_warnings(
    findings: List[Finding],
    statuses: Dict[str, str],
    lines: Iterable[DiffLine],
) -> None:
    for path, status in statuses.items():
        if status.startswith("A") and Path(path).suffix.lower() == ".hpp":
            add_finding(
                findings,
                "Avoid new .hpp propagation",
                WARN,
                path,
                None,
                "New .hpp files are discouraged unless they are narrowly justified.",
                "Prefer .h declarations with .cpp implementation, or explain why header-only implementation is needed.",
            )
    include_hpp = re.compile(r'^\s*#\s*include\s+[<"][^>"]+\.hpp[>"]')
    for line in lines:
        if Path(line.path).suffix.lower() in HEADER_EXTENSIONS and include_hpp.search(line.content):
            add_finding(
                findings,
                "Avoid new .hpp propagation",
                WARN,
                line.path,
                line.line,
                "Header diff includes a .hpp file.",
                "Avoid propagating implementation-heavy headers, or document why this include is necessary.",
            )


def check_header_include_warnings(findings: List[Finding], lines: Iterable[DiffLine]) -> None:
    include_re = re.compile(r"^\s*#\s*include\s+[<\"][^>\"]+[>\"]")
    for line in lines:
        if Path(line.path).suffix.lower() in HEADER_EXTENSIONS and include_re.search(line.content):
            add_finding(
                findings,
                "Header dependency review",
                WARN,
                line.path,
                line.line,
                "Header diff adds an include dependency.",
                "Confirm the declaration requires this include; prefer forward declarations where practical.",
            )


def is_source_under_source_tree(path: str) -> bool:
    p = Path(path)
    return path.startswith("source/") and p.suffix.lower() in (SOURCE_EXTENSIONS | {".cuh"})


def is_heterogeneous_path(path: str) -> bool:
    lowered = path.lower()
    suffix = Path(lowered).suffix
    return (
        suffix in {".cu", ".cuh"}
        or lowered.endswith(".hip.cu")
        or any(marker in lowered for marker in HETEROGENEOUS_MARKERS)
    )


def has_related_cmake_change(path: str, changed: Sequence[str]) -> bool:
    changed_cmake_dirs = {
        str(Path(changed_path).parent)
        for changed_path in changed
        if Path(changed_path).name == "CMakeLists.txt"
    }
    if not changed_cmake_dirs:
        return False
    parent = Path(path).parent
    parent_chain = [str(parent)]
    parent_chain.extend(str(p) for p in parent.parents if str(p) not in {".", ""})
    return any(directory in changed_cmake_dirs for directory in parent_chain)


def check_cmake_linkage(findings: List[Finding], statuses: Dict[str, str], changed: Sequence[str]) -> None:
    for path, status in statuses.items():
        if not status.startswith("A") or not is_source_under_source_tree(path):
            continue
        if is_heterogeneous_path(path):
            if not has_related_cmake_change(path, changed):
                add_finding(
                    findings,
                    "CMake linkage for heterogeneous sources",
                    BLOCK,
                    path,
                    None,
                    "New heterogeneous source has no related CMakeLists.txt change.",
                    "Update the same-directory or parent CMakeLists.txt, or explain generated/indirect inclusion in the PR.",
                )
            continue
        if not has_related_cmake_change(path, changed):
            add_finding(
                findings,
                "CMake linkage for new sources",
                BLOCK,
                path,
                None,
                "New source file under source/ has no related CMakeLists.txt change.",
                "Update the relevant CMakeLists.txt or explain why the file is generated or included indirectly.",
            )


def input_parameter_changed(paths: Sequence[str], lines: Sequence[DiffLine]) -> bool:
    parameter_paths = [
        path
        for path in paths
        if path.startswith("source/source_io/module_parameter/")
        and Path(path).suffix.lower() in {".cpp", ".h", ".hpp"}
    ]
    if parameter_paths:
        sensitive = re.compile(
            r"\b(Input_Item|add_item|default_value|description|category|availability|read_value|reset_value|check_value|type)\b"
        )
        return any(
            line.path in parameter_paths
            and not line.content.lstrip().startswith("//")
            and sensitive.search(line.content)
            for line in lines
        )
    return any(
        line.path.startswith("source/")
        and re.search(r"\bInput_Item\s+\w+|add_item\s*\(", line.content)
        for line in lines
    )


def pr_body_allows_no_input_doc_update(body: str) -> bool:
    lowered = body.lower()
    needles = [
        "input parameter documentation: not needed",
        "input docs: not needed",
        "no input documentation update required",
        "无需更新 input",
    ]
    return any(needle in lowered for needle in needles)


def check_input_parameter_docs(
    findings: List[Finding],
    changed: Sequence[str],
    statuses: Dict[str, str],
    lines: Sequence[DiffLine],
    pr_body: str,
) -> None:
    if not input_parameter_changed(changed, lines):
        return
    has_yaml = "docs/parameters.yaml" in changed and not statuses.get("docs/parameters.yaml", "").startswith("D")
    has_markdown = (
        "docs/advanced/input_files/input-main.md" in changed
        and not statuses.get("docs/advanced/input_files/input-main.md", "").startswith("D")
    )
    if has_yaml and has_markdown:
        return
    if pr_body and pr_body_allows_no_input_doc_update(pr_body):
        return
    add_finding(
        findings,
        "INPUT parameter documentation linkage",
        BLOCK,
        "source/source_io/module_parameter",
        None,
        "INPUT parameter behavior appears to change without both docs/parameters.yaml and input-main.md updates.",
        "Regenerate docs/parameters.yaml and docs/advanced/input_files/input-main.md, or state why no INPUT documentation update is required in the PR.",
    )


def read_pr_body(event_path: Optional[str]) -> Optional[str]:
    if not event_path:
        return None
    try:
        with open(event_path, "r", encoding="utf-8") as handle:
            payload = json.load(handle)
    except (OSError, json.JSONDecodeError):
        return None
    if "pull_request" not in payload or not isinstance(payload["pull_request"], dict):
        return None
    body = payload["pull_request"].get("body")
    if body is None:
        return ""
    return str(body)


def pr_sections(body: str) -> Dict[str, str]:
    matches = list(PR_SECTION_RE.finditer(body))
    sections: Dict[str, str] = {}
    for index, match in enumerate(matches):
        start = match.end()
        end = matches[index + 1].start() if index + 1 < len(matches) else len(body)
        sections[match.group(1).strip()] = body[start:end].strip()
    return sections


def section_is_placeholder(content: str) -> bool:
    stripped = content.strip()
    if not stripped:
        return True
    lowered = stripped.lower()
    placeholder_patterns = [
        r"fix #\.\.\.",
        r"example:",
        r"ignore if not applicable",
        r"\byes/no/not applicable\b",
        r"a unit test is added for each new feature or bug fix",
        r"my changes might affect",
        r"because \.\.\.",
    ]
    if any(re.search(pattern, lowered) for pattern in placeholder_patterns):
        return True
    meaningful = [
        line.strip()
        for line in stripped.splitlines()
        if line.strip() and not re.match(r"^[-*]\s*[^:]+:\s*$", line.strip())
    ]
    return not meaningful


def check_pr_metadata(findings: List[Finding], body: Optional[str]) -> None:
    if body is None:
        return
    required_sections = [
        "Linked Issue",
        "Unit Tests and/or Case Tests for my changes",
        "What's changed?",
    ]
    sections = pr_sections(body)
    missing = [section for section in required_sections if section not in sections]
    placeholders = [
        section
        for section in required_sections
        if section in sections and section_is_placeholder(sections[section])
    ]
    if missing or placeholders:
        reason_parts = []
        if missing:
            reason_parts.append("missing sections: " + ", ".join(missing))
        if placeholders:
            reason_parts.append("empty or placeholder sections: " + ", ".join(placeholders))
        add_finding(
            findings,
            "PR metadata completeness",
            WARN,
            "pull_request.body",
            None,
            "; ".join(reason_parts),
            "Fill the PR template with issue linkage, test evidence, and a concise change summary.",
        )


def pr_test_section_has_evidence(pr_body: str) -> bool:
    if not pr_body:
        return False
    sections = pr_sections(pr_body)
    content = sections.get("Unit Tests and/or Case Tests for my changes", "")
    if section_is_placeholder(content):
        return False
    lowered = content.lower()
    explicit_no_test_rationale = (
        re.search(r"\btests?\s+(?:are\s+)?not required\b.+\bbecause\b", lowered)
        or re.search(r"\bnot applicable\b\s*:?.*\b(?:docs|documentation) only\b", lowered)
        or "docs only" in lowered
        or "documentation only" in lowered
    )
    if explicit_no_test_rationale:
        return True
    missing_test_phrases = [
        r"\bno tests? (?:were )?run\b",
        r"\bno tests? (?:were )?added\b",
        r"\bno tests? added yet\b",
        r"\bnot run tests?\b",
        r"\btests? (?:were )?not run\b",
    ]
    if any(re.search(pattern, lowered) for pattern in missing_test_phrases):
        return False
    no_test_reason = (
        "not required" in lowered
        or "not applicable" in lowered
    )
    command_or_test = any(pattern in lowered for pattern in TEST_NAME_PATTERNS)
    return no_test_reason or command_or_test


def path_is_test(path: str) -> bool:
    lowered = path.lower()
    name = Path(lowered).name
    return lowered.startswith(TEST_PATH_PREFIXES) or any(pattern in name for pattern in TEST_NAME_PATTERNS)


def source_code_changed(changed: Sequence[str]) -> bool:
    return any(path.startswith("source/") and Path(path).suffix.lower() in SOURCE_REVIEW_EXTENSIONS for path in changed)


def check_test_evidence_warning(findings: List[Finding], changed: Sequence[str], pr_body: str) -> None:
    if not source_code_changed(changed):
        return
    if any(path_is_test(path) for path in changed) or pr_test_section_has_evidence(pr_body):
        return
    add_finding(
        findings,
        "Test evidence review",
        WARN,
        "pull_request.body",
        None,
        "Source code changed without test path changes or PR test evidence.",
        "Add focused tests, update a relevant case, or document why tests are not required.",
    )


def check_heterogeneous_test_warning(findings: List[Finding], changed: Sequence[str], pr_body: str) -> None:
    hetero_changed = any(path.startswith("source/") and is_heterogeneous_path(path) for path in changed)
    if not hetero_changed:
        return
    if any(path_is_test(path) for path in changed) or pr_test_section_has_evidence(pr_body):
        return
    add_finding(
        findings,
        "Heterogeneous test evidence review",
        WARN,
        "pull_request.body",
        None,
        "Heterogeneous source changed without test path changes or PR test evidence.",
        "Add backend-specific test evidence or explain why existing coverage is sufficient.",
    )


def check_documentation_warning(findings: List[Finding], changed: Sequence[str], pr_body: str) -> None:
    code_changed = source_code_changed(changed)
    docs_changed = any(path.startswith("docs/") for path in changed)
    if code_changed and not docs_changed and "no documentation update required" not in pr_body.lower():
        add_finding(
            findings,
            "Documentation sync review",
            WARN,
            "pull_request.body",
            None,
            "Source changes have no docs change or explicit no-docs-needed statement.",
            "Add documentation updates for behavior/interface changes, or state why documentation is not required.",
        )


def collect_findings(root: Path, args: argparse.Namespace) -> List[Finding]:
    findings: List[Finding] = []
    statuses, changed = changed_paths(root, args)
    lines, removed_lines = changed_lines(root, args)
    body = read_pr_body(args.event_path)
    body_text = body or ""

    check_line_endings(findings, root, changed, statuses, args)
    check_global_dependencies(findings, lines, removed_lines)
    check_default_parameters(findings, lines)
    check_hpp_warnings(findings, statuses, lines)
    check_header_include_warnings(findings, lines)
    check_cmake_linkage(findings, statuses, changed)
    check_input_parameter_docs(findings, changed, statuses, lines, body_text)
    check_pr_metadata(findings, body)
    check_test_evidence_warning(findings, changed, body_text)
    check_heterogeneous_test_warning(findings, changed, body_text)
    check_documentation_warning(findings, changed, body_text)
    return findings


def finding_location(finding: Finding) -> str:
    if finding.line is None:
        return finding.path
    return f"{finding.path}:{finding.line}"


def render_text(findings: Sequence[Finding]) -> str:
    if not findings:
        return "Agent governance check: no findings.\n"
    chunks = ["Agent governance check findings:"]
    for finding in findings:
        chunks.append(
            f"- [{finding.severity.upper()}] {finding.rule} at {finding_location(finding)}\n"
            f"  Reason: {finding.reason}\n"
            f"  Suggested action: {finding.suggestion}\n"
            f"  Exception allowed: {'yes' if finding.allow_exception else 'no'}"
        )
    return "\n".join(chunks) + "\n"


def render_markdown(findings: Sequence[Finding]) -> str:
    if not findings:
        return "## Agent Governance Check\n\nNo findings.\n"
    lines = [
        "## Agent Governance Check",
        "",
        "| Severity | Rule | Location | Reason | Suggested action | Exception |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for finding in findings:
        lines.append(
            "| {severity} | {rule} | `{location}` | {reason} | {suggestion} | {exception} |".format(
                severity=finding.severity,
                rule=finding.rule,
                location=finding_location(finding),
                reason=finding.reason.replace("|", "\\|"),
                suggestion=finding.suggestion.replace("|", "\\|"),
                exception="allowed" if finding.allow_exception else "not allowed",
            )
        )
    return "\n".join(lines) + "\n"


def render_json(findings: Sequence[Finding]) -> str:
    return json.dumps([asdict(finding) for finding in findings], indent=2, ensure_ascii=False) + "\n"


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    diff_group = parser.add_mutually_exclusive_group()
    diff_group.add_argument("--staged", action="store_true", help="Check staged changes.")
    parser.add_argument("--base", help="Base commit for diff checks.")
    parser.add_argument("--head", help="Head commit for diff checks.")
    parser.add_argument("--event-path", help="GitHub event JSON path for PR body checks.")
    parser.add_argument(
        "--format",
        choices=("text", "markdown", "json"),
        default="text",
        help="Output format.",
    )
    args = parser.parse_args(argv)
    if args.staged and (args.base or args.head):
        parser.error("--staged cannot be combined with --base/--head")
    if bool(args.base) ^ bool(args.head):
        parser.error("--base and --head must be provided together")
    return args


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    root = repo_root()
    try:
        findings = collect_findings(root, args)
    except GitError as exc:
        print(f"agent_governance_check.py: {exc}", file=sys.stderr)
        return 2

    if args.format == "markdown":
        sys.stdout.write(render_markdown(findings))
    elif args.format == "json":
        sys.stdout.write(render_json(findings))
    else:
        sys.stdout.write(render_text(findings))

    return 1 if any(finding.severity == BLOCK for finding in findings) else 0


if __name__ == "__main__":
    raise SystemExit(main())
