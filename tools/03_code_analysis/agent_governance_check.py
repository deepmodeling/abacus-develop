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
EXCEPTION_SECTION = "Governance Exception"
EXCEPTION_APPROVED_LABEL = "governance-exception-approved"
EXCEPTION_FIELDS = (
    "Rule",
    "Reason",
    "Scope",
    "User or maintenance risk",
    "Why the normal rule cannot be followed now",
    "Follow-up cleanup plan",
    "Requested approver",
)
EXCEPTION_FIELD_RE = re.compile(r"^[-*]\s*([^:]+):\s*(.*)$")


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

# `#define private public` / `#define protected public`. See AGENTS.md rule 10.
ACCESS_HACK_RE = re.compile(r"^\s*#\s*define\s+(?:private|protected)\s+public\b")

# `#pragma once`. See AGENTS.md rule 15.
PRAGMA_ONCE_RE = re.compile(r"^\s*#\s*pragma\s+once\b")


def is_global_dependency_check_path(path: str) -> bool:
    if path.startswith("tools/03_code_analysis/"):
        return False
    return Path(path).suffix.lower() in CODE_EXTENSIONS


def is_access_hack_check_path(path: str) -> bool:
    """Only C/C++ translation units can carry a real access-control hack.

    Without this filter a fenced `#define private public` inside AGENTS.md or
    docs/ counts against the budget, so documenting the anti-pattern would
    block CI and deleting that documentation would credit it.
    """
    if path.startswith("tools/03_code_analysis/"):
        return False
    return Path(path).suffix.lower() in SOURCE_REVIEW_EXTENSIONS


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


def check_access_hacks(
    findings: List[Finding],
    added_lines: Iterable[DiffLine],
    removed_lines: Iterable[DiffLine],
) -> None:
    """Ratchet on `#define private public` (AGENTS.md rule 10).

    Mirrors the global-dependency budget: removals are free, a net increase
    blocks. The remaining offenders therefore do not block unrelated work while
    they are being refactored away module by module.
    """
    added = [
        line
        for line in added_lines
        if is_access_hack_check_path(line.path) and ACCESS_HACK_RE.search(line.content)
    ]
    removed = [
        line
        for line in removed_lines
        if is_access_hack_check_path(line.path) and ACCESS_HACK_RE.search(line.content)
    ]
    if not added:
        return

    delta = len(added) - len(removed)
    severity = BLOCK if delta > 0 else WARN
    for line in added:
        add_finding(
            findings,
            "No access-control hacks",
            severity,
            line.path,
            line.line,
            (
                "Adds `#define private/protected public`, which reinterprets access "
                "control for the whole translation unit (standard library headers "
                "included) and makes this TU disagree with the rest of the build; "
                "PR total added={a}, removed={r}, net_delta={d}.".format(
                    a=len(added), r=len(removed), d=delta
                )
            ),
            (
                "Pass the INPUT values the code needs as explicit arguments instead of "
                "reading global PARAM inside it, so the test can drive it without "
                "touching PARAM at all; otherwise add a public const observer, or an "
                "explicit `friend class XxxTest;` on the class under test."
            ),
        )


def check_pragma_once(
    findings: List[Finding],
    added_lines: Iterable[DiffLine],
) -> None:
    """Block new `#pragma once` in headers (AGENTS.md rule 15).

    Use traditional `#ifndef`/`#define`/`#endif` include guards instead.
    """
    for line in added_lines:
        if Path(line.path).suffix.lower() not in HEADER_EXTENSIONS | {".cuh"}:
            continue
        if PRAGMA_ONCE_RE.search(line.content):
            add_finding(
                findings,
                "No #pragma once",
                BLOCK,
                line.path,
                line.line,
                "Adds `#pragma once` to a header file.",
                "Use traditional `#ifndef`/`#define`/`#endif` include guards instead.",
                allow_exception=False,
            )


def _has_default_arg_in_parens(stripped: str) -> bool:
    paren_depth = 0
    in_string = False
    string_char = ""
    in_line_comment = False
    in_block_comment = False
    i = 0
    while i < len(stripped):
        c = stripped[i]
        nxt = stripped[i + 1] if i + 1 < len(stripped) else ""

        if in_line_comment:
            break

        if in_block_comment:
            if c == "*" and nxt == "/":
                in_block_comment = False
                i += 2
                continue
            i += 1
            continue

        if not in_string:
            if c == "/" and nxt == "/":
                in_line_comment = True
                i += 2
                continue
            if c == "/" and nxt == "*":
                in_block_comment = True
                i += 2
                continue

        if in_string:
            if c == "\\":
                i += 2
                continue
            if c == string_char:
                in_string = False
            i += 1
            continue
        if c in ('"', "'"):
            in_string = True
            string_char = c
            i += 1
            continue
        if c == "(":
            paren_depth += 1
        elif c == ")":
            if paren_depth > 0:
                paren_depth -= 1
        elif c == "=" and paren_depth > 0:
            return True
        i += 1
    return False


def check_default_parameters(findings: List[Finding], lines: Iterable[DiffLine]) -> None:
    default_arg = re.compile(r"[(,]\s*[^()=;,{}]+\b\w+\s*=\s*[^,);{}]+")
    control_flow = re.compile(r"^(for|if|while|switch|catch)\s*\(")
    comment_strip_re = re.compile(r"//.*$|/\*.*?\*/")
    for line in lines:
        if Path(line.path).suffix.lower() not in HEADER_EXTENSIONS:
            continue
        stripped = line.content.strip()
        if not stripped or stripped.startswith("//") or stripped.startswith("*"):
            continue
        if control_flow.match(stripped):
            continue
        if "=" not in stripped:
            continue
        code_only = comment_strip_re.sub("", stripped)
        if "=" not in code_only:
            continue
        if not _has_default_arg_in_parens(code_only):
            continue
        if "(" in code_only and ")" in code_only and default_arg.search(code_only):
            add_finding(
                findings,
                "No new default parameters",
                WARN,
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
                "Confirm this include is needed by the header and does not introduce unnecessary or cyclic dependencies.",
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


# CMake governance checks. See AGENTS.md rule 24 and docs/developers_guide/cmake_guidelines.md.
CMAKE_MATH_LIBS_RE = re.compile(r"\$\{\s*math_libs\s*\}")
CMAKE_DEPRECATED_SYMBOL_PATTERNS = (
    (re.compile(r"\bUSE_LIBXC\b"), "USE_LIBXC", "Use `__LIBXC` instead of the renamed `USE_LIBXC` macro."),
    (
        re.compile(r"\b(?:ENABLE|USE)_RAPIDJSON\b"),
        "ENABLE_RAPIDJSON/USE_RAPIDJSON",
        "Use `ENABLE_JSON` (`__JSON`, backed by nlohmann_json) instead of removed RapidJSON options/macros.",
    ),
    (
        re.compile(r"\bENABLE_LIBCOMM\b"),
        "ENABLE_LIBCOMM",
        "Use `ENABLE_LIBRI` instead; LibComm is no longer a standalone ABACUS option.",
    ),
)
CMAKE_SOURCE_FORBIDDEN_COMMANDS = {
    "add_definitions", "add_compile_definitions", "add_compile_options",
    "include_directories", "link_directories", "link_libraries", "find_package", "find_library",
}
CMAKE_BRACKET_RE = re.compile(r"\[(=*)\[")


@dataclass
class CMakeToken:
    content: str
    line: int
    kind: str
    source: str


@dataclass
class CMakeCommand:
    path: str
    name: CMakeToken
    arguments: List[CMakeToken]


def is_cmake_file(path: str) -> bool:
    p = Path(path)
    return p.name == "CMakeLists.txt" or p.suffix.lower() == ".cmake"


def is_main_cmake_path(path: str) -> bool:
    """Main ABACUS build tree CMake files (excludes standalone subprojects like python/pyabacus)."""
    if not is_cmake_file(path):
        return False
    return path == "CMakeLists.txt" or path.startswith(("cmake/", "source/", "tests/"))


def is_source_cmake_path(path: str) -> bool:
    return is_cmake_file(path) and path.startswith("source/")


def is_source_unittest_cmake_path(path: str) -> bool:
    if not is_source_cmake_path(path):
        return False
    dirs = Path(path).parts[1:-1]
    return any(d == UNITTEST_DIR or is_legacy_test_dir(d) for d in dirs)


def cmake_tokens(content: str) -> Iterable[CMakeToken]:
    """Scan literal command syntax, retaining line spans and skipping all comments.

    Quoted/bracket arguments are single tokens, so their parentheses and hashes
    cannot become commands or comments. This does not evaluate CMake variables.
    """
    index = 0
    line = 1
    while index < len(content):
        start = index
        start_line = line
        ch = content[index]
        if ch.isspace():
            line += ch == "\n"
            index += 1
            continue
        is_comment = ch == "#"
        bracket = CMAKE_BRACKET_RE.match(content, index + int(is_comment))
        if bracket:
            close = "]" + bracket.group(1) + "]"
            end = content.find(close, bracket.end())
            stop = len(content) if end < 0 else end
            index = len(content) if end < 0 else end + len(close)
            line += content[start:index].count("\n")
            if not is_comment:
                yield CMakeToken(content[bracket.end():stop], start_line, "argument", content[start:index])
            continue
        if is_comment:
            end = content.find("\n", index)
            index = len(content) if end < 0 else end
            continue
        if ch in "()":
            index += 1
            yield CMakeToken(ch, line, ch, ch)
            continue
        # Quotes may also appear within legacy unquoted arguments. Escaped
        # delimiters and line continuations belong to the same argument.
        value = []
        quoted = False
        while index < len(content):
            ch = content[index]
            if ch == "\\" and index + 1 < len(content):
                escaped = content[index + 1]
                line += escaped == "\n"
                if escaped != "\n":
                    value.append(escaped)
                index += 2
            elif ch == '"':
                quoted = not quoted
                index += 1
            elif not quoted and (ch.isspace() or ch in "()#"):
                break
            else:
                value.append(ch)
                line += ch == "\n"
                index += 1
        yield CMakeToken("".join(value), start_line, "argument", content[start:index])


def cmake_commands(path: str, content: str) -> List[CMakeCommand]:
    tokens = list(cmake_tokens(content))
    commands = []
    index = 0
    while index + 1 < len(tokens):
        name = tokens[index]
        if name.kind != "argument" or tokens[index + 1].kind != "(":
            index += 1
            continue
        index += 2
        arguments = []
        depth = 1
        while index < len(tokens) and depth:
            token = tokens[index]
            if token.kind == "(":
                depth += 1
            elif token.kind == ")":
                depth -= 1
            if depth:
                arguments.append(token)
            index += 1
        commands.append(CMakeCommand(path, name, arguments))
    return commands


def changed_cmake_commands(
    root: Path, args: argparse.Namespace, added_lines: Sequence[DiffLine]
) -> List[Tuple[CMakeCommand, List[DiffLine]]]:
    """Read the index/head, not the worktree, before selecting added code tokens."""
    added_by_path: Dict[str, set] = {}
    for line in added_lines:
        if is_main_cmake_path(line.path) and line.line is not None:
            added_by_path.setdefault(line.path, set()).add(line.line)
    changed = []
    for path, added in added_by_path.items():
        content = read_changed_file_bytes(root, path, args).decode("utf-8")
        for command in cmake_commands(path, content):
            tokens = [command.name, *command.arguments]
            # Keep only added physical lines even when a quoted/bracket argument
            # spans both old and new text. Searching the whole argument would
            # report unchanged legacy symbols and point outside the diff.
            touched = [
                DiffLine(path, token.line + offset, text)
                for token in tokens
                for offset, text in enumerate(token.source.splitlines())
                if token.line + offset in added
            ]
            if touched:
                changed.append((command, touched))
    return changed


def check_cmake_deprecated_patterns(
    findings: List[Finding], commands: Sequence[Tuple[CMakeCommand, List[DiffLine]]]
) -> None:
    for command, touched in commands:
        if command.name.content.lower() == "message":
            continue
        if command.name.content.lower() == "remove_definitions":
            add_finding(
                findings, "CMake deprecated patterns", BLOCK, command.path, touched[0].line,
                "Added CMake code uses `remove_definitions()`, which does not strip target-scoped feature definitions.",
                "Use `abacus_disable_feature_definitions(...)` instead (see docs/developers_guide/cmake_guidelines.md).",
            )
        for token in touched:
            if CMAKE_MATH_LIBS_RE.search(token.content):
                add_finding(
                    findings, "CMake deprecated patterns", BLOCK, command.path, token.line,
                    "Added CMake argument references removed `${math_libs}` variable.",
                    "Remove `${math_libs}`; `AddTest()` links `abacus::linalg_libs` automatically and "
                    "`source/CMakeLists.txt` owns the executable link closure.",
                )
            for pattern, label, suggestion in CMAKE_DEPRECATED_SYMBOL_PATTERNS:
                if pattern.search(token.content):
                    add_finding(
                        findings, "CMake deprecated patterns", BLOCK, command.path, token.line,
                        f"Added CMake argument references deprecated or removed symbol `{label}`.", suggestion,
                    )


def check_cmake_source_commands(
    findings: List[Finding], commands: Sequence[Tuple[CMakeCommand, List[DiffLine]]]
) -> None:
    for command, touched in commands:
        cmd = command.name.content.lower()
        if not is_source_cmake_path(command.path) or cmd not in CMAKE_SOURCE_FORBIDDEN_COMMANDS:
            continue
        if cmd in {"find_package", "find_library"}:
            suggestion = (
                "Move package/library discovery to root `CMakeLists.txt` and `cmake/modules/Find<Pkg>.cmake` "
                "(exporting an IMPORTED target), or `cmake/Testing.cmake` for test-only packages."
            )
        else:
            suggestion = (
                "Use target-scoped `target_*()` commands, `cmake/CompilerConfiguration.cmake`, or "
                "`_abacus_feature_definitions` / `abacus_compile_requirements` in `source/CMakeLists.txt`."
            )
        add_finding(
            findings, "CMake directory-scoped commands under source/", BLOCK, command.path, touched[0].line,
            f"Added code in `{cmd}()` under `source/`.", suggestion,
        )


def check_cmake_test_file_preparation(
    findings: List[Finding], commands: Sequence[Tuple[CMakeCommand, List[DiffLine]]]
) -> None:
    for command, touched in commands:
        if not is_source_unittest_cmake_path(command.path):
            continue
        name = command.name.content.lower()
        if name == "install":
            label = "install()"
        elif name == "file" and command.arguments and command.arguments[0].content.upper() == "COPY":
            label = "file(COPY ...)"
        else:
            continue
        add_finding(
            findings, "CMake test file preparation", BLOCK, command.path, touched[0].line,
            f"Unit-test CMake file uses `{label}` to stage test files.",
            "Use `abacus_copy_test_files(...)` from `cmake/Testing.cmake` so test scripts and data are "
            "prepared at configure time.",
        )


def check_cmake_option_and_dep_sync(
    findings: List[Finding],
    root: Path,
    args: argparse.Namespace,
    statuses: Dict[str, str],
    changed: Sequence[str],
    commands: Sequence[Tuple[CMakeCommand, List[DiffLine]]],
) -> None:
    if "cmake/CollectBuildInfoVars.cmake" in changed:
        return
    option_commands = [
        (command, touched) for command, touched in commands
        if command.path == "CMakeLists.txt"
        and command.name.content.lower() in {"option", "cmake_dependent_option"}
        and command.arguments
    ]
    if option_commands:
        base = "HEAD" if args.staged else args.base
        old_names = set()
        if git(["ls-tree", "--name-only", base, "--", "CMakeLists.txt"], root).stdout.strip():
            before = git(["show", f"{base}:CMakeLists.txt"], root).stdout
            old_names = {
                command.arguments[0].content for command in cmake_commands("CMakeLists.txt", before)
                if command.name.content.lower() in {"option", "cmake_dependent_option"} and command.arguments
            }
        for command, touched in option_commands:
            if command.arguments[0].content in old_names:
                continue
            add_finding(
                findings, "CMake option and dependency sync", WARN, command.path, touched[0].line,
                "Root `CMakeLists.txt` adds an option without updating `cmake/CollectBuildInfoVars.cmake`.",
                "Update `cmake/CollectBuildInfoVars.cmake` (and `_abacus_feature_definitions` / "
                "`abacus_add_target_compile_requirements()` in `source/CMakeLists.txt` for external "
                "dependencies), or explain why the build summary does not need this option.",
            )
            break
    for path, status in statuses.items():
        if status.startswith("A") and path.startswith("cmake/modules/") and Path(path).suffix.lower() == ".cmake":
            add_finding(
                findings, "CMake option and dependency sync", WARN, path, None,
                "New CMake discovery module added without updating `cmake/CollectBuildInfoVars.cmake`.",
                "Record the new dependency in `cmake/CollectBuildInfoVars.cmake` and wire its IMPORTED "
                "target into `source/CMakeLists.txt` (`_abacus_feature_libs` and "
                "`abacus_add_target_compile_requirements()`), or explain why no summary update is needed.",
            )




# Unit test layout under source/. See AGENTS.md rule 11.
UNITTEST_DIR = "unittests"
LEGACY_TEST_DIR_NAMES = {"test", "tests", "unittest", "unit_test", "unit_tests"}
UNITTEST_CODE_EXTENSIONS = SOURCE_EXTENSIONS | HEADER_EXTENSIONS | {".cuh"}


def is_legacy_test_dir(name: str) -> bool:
    """`test/`, `test_serial/`, `UnitTests/`, ... -- anything test-like but `unittests/`."""
    if name == UNITTEST_DIR:
        return False
    lowered = name.lower()
    return lowered in LEGACY_TEST_DIR_NAMES or lowered == UNITTEST_DIR or lowered.startswith("test_")


def check_unittest_layout(findings: List[Finding], statuses: Dict[str, str]) -> None:
    """Enforce the `unittests/` layout for unit tests under source/.

    Added, renamed, and copied files are checked at their new path; modified
    files are checked too, so touching a test means bringing it into the
    layout. Deletions are always allowed, which is how legacy directories
    shrink.
    """
    for path, status in statuses.items():
        if not path.startswith("source/") or status.startswith("D"):
            continue
        added = status[:1] in {"A", "R", "C"}
        p = Path(path)
        suffix = p.suffix.lower()
        dirs = p.parts[1:-1]
        # Data subdirectories below unittests/ may use any name, `test/` included.
        owners = dirs[:dirs.index(UNITTEST_DIR)] if UNITTEST_DIR in dirs else dirs
        legacy = next((d for d in owners if is_legacy_test_dir(d)), None)
        if legacy is not None:
            if added:
                reason = f"Adds a file to the legacy test directory `{legacy}/`."
            elif suffix in SOURCE_EXTENSIONS:
                reason = f"Modifies a unit test source in the legacy test directory `{legacy}/`."
            else:
                continue
            add_finding(
                findings,
                "Unit test directory layout",
                BLOCK,
                path,
                None,
                reason,
                "Move the test into a sibling `unittests/` directory (named test_<file>.cpp) "
                "and wire it into that directory's CMakeLists.txt; legacy test directories "
                "may only have files deleted or migrated out.",
            )
            continue
        if UNITTEST_DIR not in dirs:
            continue
        below = dirs[dirs.index(UNITTEST_DIR) + 1:]
        if below:
            if suffix in UNITTEST_CODE_EXTENSIONS:
                add_finding(
                    findings,
                    "Unit test directory layout",
                    BLOCK,
                    path,
                    None,
                    f"Code file inside the `unittests/{'/'.join(below)}/` subdirectory.",
                    "Keep test sources and headers directly in `unittests/`; its subdirectories "
                    "may hold test data only.",
                )
            continue
        if suffix in SOURCE_EXTENSIONS and not p.name.startswith("test_"):
            add_finding(
                findings,
                "Unit test file naming",
                BLOCK,
                path,
                None,
                f"Source file `{p.name}` in `unittests/` does not start with `test_`.",
                "Name the test test_<file under test>.cpp; move mocks and shared helpers into a .h header.",
            )



def read_pull_request(event_path: Optional[str]) -> Optional[dict]:
    if not event_path:
        return None
    try:
        with open(event_path, "r", encoding="utf-8") as handle:
            payload = json.load(handle)
    except (OSError, json.JSONDecodeError):
        return None
    if "pull_request" not in payload or not isinstance(payload["pull_request"], dict):
        return None
    return payload["pull_request"]


def read_pr_body(event_path: Optional[str]) -> Optional[str]:
    pull_request = read_pull_request(event_path)
    if pull_request is None:
        return None
    body = pull_request.get("body")
    if body is None:
        return ""
    return str(body)


def read_pr_labels(event_path: Optional[str]) -> List[str]:
    pull_request = read_pull_request(event_path)
    if pull_request is None:
        return []
    labels = pull_request.get("labels") or []
    return [
        str(label.get("name", ""))
        for label in labels
        if isinstance(label, dict)
    ]


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


def parse_governance_exception(body: str) -> Optional[Dict[str, str]]:
    sections = pr_sections(body)
    content = sections.get(EXCEPTION_SECTION)
    if content is None:
        return None
    fields: Dict[str, str] = {}
    for raw in content.splitlines():
        match = EXCEPTION_FIELD_RE.match(raw.strip())
        if match:
            fields[match.group(1).strip().lower()] = match.group(2).strip()
    return fields


def apply_governance_exception(
    findings: List[Finding],
    body: Optional[str],
    labels: Sequence[str],
) -> None:
    """Downgrade blocking findings covered by an approved Governance Exception.

    The PR body records the exception; the maintainer-applied label is the
    approval. Findings with allow_exception=False are never downgraded.
    """
    if body is None:
        return
    approved = EXCEPTION_APPROVED_LABEL in labels
    fields = parse_governance_exception(body)
    if fields is None or not any(fields.values()):
        if approved:
            add_finding(
                findings,
                "Governance exception record",
                WARN,
                "pull_request.body",
                None,
                f"PR carries the `{EXCEPTION_APPROVED_LABEL}` label but has no filled `### {EXCEPTION_SECTION}` section.",
                "Record the exception with the template in docs/developers_guide/agent_governance.md, or remove the label.",
            )
        return
    missing = [field for field in EXCEPTION_FIELDS if not fields.get(field.lower())]
    if missing:
        add_finding(
            findings,
            "Governance exception record",
            WARN,
            "pull_request.body",
            None,
            f"`### {EXCEPTION_SECTION}` section is incomplete; missing fields: " + ", ".join(missing) + ".",
            "Fill every field of the exception template; incomplete exceptions do not downgrade blocking findings.",
        )
        return
    if not approved:
        add_finding(
            findings,
            "Governance exception record",
            WARN,
            "pull_request.body",
            None,
            f"Governance exception requested for: {fields['rule']}; awaiting maintainer approval.",
            f"A maintainer approves the exception by applying the `{EXCEPTION_APPROVED_LABEL}` label.",
        )
        return
    excepted_rules = fields["rule"].lower()
    for finding in findings:
        if finding.severity != BLOCK or not finding.allow_exception:
            continue
        if finding.rule.lower() not in excepted_rules:
            continue
        finding.severity = WARN
        finding.reason += " Downgraded from error by the approved governance exception."


def collect_findings(root: Path, args: argparse.Namespace) -> List[Finding]:
    findings: List[Finding] = []
    statuses, changed = changed_paths(root, args)
    lines, removed_lines = changed_lines(root, args)
    body = read_pr_body(args.event_path)
    labels = read_pr_labels(args.event_path)
    body_text = body or ""

    check_line_endings(findings, root, changed, statuses, args)
    check_global_dependencies(findings, lines, removed_lines)
    check_access_hacks(findings, lines, removed_lines)
    check_pragma_once(findings, lines)
    check_default_parameters(findings, lines)
    check_hpp_warnings(findings, statuses, lines)
    check_header_include_warnings(findings, lines)
    check_cmake_linkage(findings, statuses, changed)
    cmake_changes = changed_cmake_commands(root, args, lines)
    check_cmake_deprecated_patterns(findings, cmake_changes)
    check_cmake_source_commands(findings, cmake_changes)
    check_cmake_test_file_preparation(findings, cmake_changes)
    check_cmake_option_and_dep_sync(findings, root, args, statuses, changed, cmake_changes)
    check_unittest_layout(findings, statuses)
    check_pr_metadata(findings, body)
    check_test_evidence_warning(findings, changed, body_text)
    check_heterogeneous_test_warning(findings, changed, body_text)
    check_documentation_warning(findings, changed, body_text)
    apply_governance_exception(findings, body, labels)
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
