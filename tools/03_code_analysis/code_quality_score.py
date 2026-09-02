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
    tests/, test_serial/, test_parallel/, test_gpu/, unit_test/, unittest/.
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
    - `using namespace std;`: -1 per occurrence (cap 5)
    - `goto` keyword: -2 per occurrence (cap 5)
    - line longer than 120 chars: -1 per line (cap 5)
    - Chinese characters in comments/code: -1 per line (cap 5)

Usage:
    python3 code_quality_score.py source/source_base
    # writes report to ./code_quality_score.txt (text default output file)
    python3 code_quality_score.py source/ -o report.txt
    python3 code_quality_score.py --format json path/to/file.cpp
    python3 code_quality_score.py --format json source/ -o report.json
    python3 code_quality_score.py --min-score 80 source/
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Sequence, Tuple


SOURCE_EXTENSIONS = {".c", ".cc", ".cpp", ".cxx", ".cu", ".h", ".hh", ".hpp", ".hxx", ".cuh"}

SKIP_DIRS = {
    ".git", "build", "__pycache__", "node_modules", ".cache",
    "third_party", "thirdparty", ".vscode", ".idea", ".trae-cn",
    "Dependencies",
    "test", "tests", "test_serial", "test_parallel", "test_gpu", "unit_test", "unittest",
}

CAPS = {
    "tab_indentation": 5,
    "using_namespace_std": 5,
    "line_too_long": 5,
    "chinese_comment": 5,
    "uppercase_constant": 5,
    "default_parameter": 5,
    "global_dependency": 10,
    "hpp_include": 5,
    "friend_keyword": 5,
    "unpaired_new_delete": 5,
    "member_local_name_conflict": 5,
}

WEIGHTS = {
    "filename_too_long": 1,
    "filename_uppercase": 1,
    "hpp_implementation": 50,
    "public_member_variable": 1,
    "member_function_too_long": 1,
    "tab_indentation": 1,
    "using_namespace_std": 1,
    "line_too_long": 1,
    "chinese_comment": 1,
    "uppercase_constant": 1,
    "default_parameter": 2,
    "global_dependency": 3,
    "hpp_include": 2,
    "friend_keyword": 1,
    "unpaired_new_delete": 1,
    "member_local_name_conflict": 1,
    "file_too_long": 2,
    "raw_new_keyword": 1,
    "too_many_parameters": 1,
    "high_cyclomatic_complexity": 1,
    "post_cpp11_feature": 40,
    "post_cpp11_per_feature": 8,
    "goto_keyword": 2,
    "static_member_variable": 1,
}

CAPS = {
    "tab_indentation": 5,
    "using_namespace_std": 5,
    "line_too_long": 5,
    "chinese_comment": 5,
    "uppercase_constant": 5,
    "default_parameter": 5,
    "global_dependency": 10,
    "hpp_include": 5,
    "friend_keyword": 5,
    "unpaired_new_delete": 5,
    "member_local_name_conflict": 5,
    "too_many_parameters": 30,
    "high_cyclomatic_complexity": 30,
    "post_cpp11_per_feature": 5,
    "goto_keyword": 5,
    "static_member_variable": 10,
}

FUNCTION_LENGTH_THRESHOLD = 50
FUNCTION_LENGTH_STEP = 50
FILE_LENGTH_THRESHOLD = 500
FILE_LENGTH_STEP = 50
FUNCTION_PARAM_THRESHOLD = 7
CYCLO_THRESHOLD = 10
LINE_LENGTH_LIMIT = 120
FILENAME_LENGTH_LIMIT = 20
PASS_THRESHOLD = 60

CHINESE_RE = re.compile("[\u4e00-\u9fff]")
USING_NS_STD_RE = re.compile(r"\busing\s+namespace\s+std\b")
UPPERCASE_CONST_RE = re.compile(r"(?<![.:>])\b[A-Z][A-Z0-9_]{2,}\b(?![.:])")
ACCESS_RE = re.compile(r"^\s*(public|private|protected)\s*:")
CLASS_OPEN_RE = re.compile(r"\b(class|struct)\s+(\w+)\b")

GLOBAL_DEPENDENCY_RE = re.compile(r"\b(?:GlobalV::|GlobalC::|PARAM(?:\.|->|::))")
HPP_INCLUDE_RE = re.compile(r'^\s*#\s*include\s+[<"][^>"]+\.hpp[>"]')
FRIEND_RE = re.compile(r"\bfriend\b(?!\s+(?:class\s+|struct\s+)?std::)")
NEW_EXPR_RE = re.compile(r"\bnew\s+\w")
DELETE_EXPR_RE = re.compile(r"\bdelete\s*\[\s*\]?\s+\w")

# Occurrences of `new T(...)` or `new T[...]` that are clearly owned by a
# smart-pointer wrapper and therefore should NOT count as unpaired
# new/delete leaks. Two idiomatic C++11 shapes are recognised:
#   1. `ptr.reset(new T(args))`  /  `ptr->reset(new T(args))`
#      — reset() transfers ownership into an existing smart pointer;
#      most commonly used with unique_ptr/shared_ptr when make_unique/
#      make_shared are not available (pre-C++14) or unsuitable.
#   2. `unique_ptr<T> p(new T(args))` / `shared_ptr<T> p(new T(args))`
#      plus the allocator-taking `allocate_shared` form — direct
#      construction of a smart pointer that owns the freshly-allocated
#      object from day one.
# The pattern does NOT try to skip new in `raw_new_count` (that rule
# penalises *any* direct use of `new` vs. make_unique/make_shared); it
# only affects the ownership heuristic `unpaired_new`.
_OWNED_NEW_SHAPES = [
    # shape 1:  .reset(new T  or  ->reset(new T   (optionally with whitespace)
    r"\.(?:reset|operator\s*=)\s*\(\s*new\s+\w",
    r"->(?:reset|operator\s*=)\s*\(\s*new\s+\w",
    # shape 2:  unique_ptr<T[const & ...]>   ident ( new T
    #           shared_ptr<T[const & ...]>   ident ( new T
    #           allocate_shared<T,...> ( alloc , new T   -> not exact, but
    #           the simpler `allocate_shared<...> ( ... new T` heuristic
    r"\b(?:unique_ptr|shared_ptr|scoped_ptr|auto_ptr)\s*<[^>]*>\s*[A-Za-z_]\w*\s*\(\s*new\s+\w",
    r"\b(?:make_unique|make_shared|allocate_shared)\s*<[^>]*>\s*\([^;()]*?\bnew\s+\w",
]
OWNED_NEW_RE = re.compile("|".join(f"(?:{p})" for p in _OWNED_NEW_SHAPES))
GOTO_RE = re.compile(r"\bgoto\b")

DEFAULT_PARAM_RE = re.compile(
    r"\b\w+\s*\([^();]*\b\w+\s*=(?![=>])[^();]*\)\s*"
    r"(?:const\s*|noexcept\s*|override\s*|final\s*)*[;{]"
)

# Strict declaration regex: requires `type ... name ;` or `type ... name = ...;`.
# Excludes pure assignments like `counter = counter + 1;` which lack a leading type.
DECL_RE = re.compile(r"^\s*(?:\w+[\s\*]+)+(\w+)\s*(?:=[^;]*)?;\s*$")

BLOCK_KEYWORD_PREFIXES = (
    "public:", "private:", "protected:",
    "typedef", "using", "static_assert",
    "class", "struct", "friend",
    "//", "/*", "*",
    "return", "if ", "for ", "while ",
    "switch ", "case ", "template",
    "break", "continue", "goto", "default:",
    "else", "do", "try", "catch", "throw",
    "namespace", "enum", "union",
)

NON_FUNCTION_KEYWORDS = {
    "if", "for", "while", "switch", "sizeof", "return", "throw",
    "catch", "class", "struct", "namespace", "enum", "union",
    "template", "static_assert", "typedef", "using", "new",
    "delete", "operator", "do", "else", "goto", "continue",
    "break", "try",
}

# Keywords that, when they appear as the last identifier in the prefix
# before a `name(args)` candidate followed by `;` or `=`, strongly
# indicate the candidate is a function call in statement/expression
# context rather than a declaration.
STATEMENT_CONTEXT_KEYWORDS = NON_FUNCTION_KEYWORDS | {
    "co_await", "co_return", "co_yield",
    "typeid", "noexcept", "alignof", "alignas", "decltype",
    "static_cast", "const_cast", "dynamic_cast", "reinterpret_cast",
    "and", "or", "not", "xor", "bitor", "bitand", "compl",
}

# Punctuation characters that, when they are the last non-whitespace
# char before `name(args)` followed by `;` or `=`, indicate an
# expression / argument-list context (i.e. a call, not a declaration).
# `*` and `&` are intentionally NOT in this set because they also act
# as valid pointer/reference type qualifiers in a return type. The
# double-character operator forms (`&&`, `||`, `**`, etc.) are handled
# by the two-character check inside _is_declaration_prefix.
# Examples: `foo(1, bar(x), 3)`    ->  prefix ends with `,`
#           `if (cond && ok(x))`    ->  prefix ends with `&&`
#           `v.push_back(fn())`      ->  prefix ends with `.`
DECL_PREFIX_BAD_TAIL = set("=+-/%|^~!<>?([{,.;:")

# Last-identifier regex: extracts the trailing identifier word from a
# prefix string (strips trailing punctuation and whitespace first).
_LAST_IDENT_RE = re.compile(r"([A-Za-z_]\w*)\s*$")


def _is_declaration_prefix(prefix_stripped: str) -> bool:
    """Return True when `prefix_stripped` (the right-stripped text between
    the start of the line and the candidate `name(args)` opener) is
    consistent with a function declaration/definition context rather
    than a call-site expression.

    A prefix is *not* a declaration prefix when:
      * it is empty;
      * its last non-whitespace character is punctuation that typically
        appears inside expressions / argument lists (see
        DECL_PREFIX_BAD_TAIL), or it is the double-character operator
        form of `&&`, `||`, `**`, `->` that commonly appears inside
        expressions;
      * its trailing identifier token is a statement/expression
        keyword (see STATEMENT_CONTEXT_KEYWORDS) such as `return`,
        `throw`, `if`, `sizeof`, `static_cast`, etc.

    Constructors/destructors that open with `{` instead of `;` or `=`
    are handled separately by the caller and do not go through this
    check.
    """
    if not prefix_stripped:
        return False
    tail_char = prefix_stripped[-1]
    # Two-character expression operators that include `*` or `&` (which
    # are not declared in the single-char bad-tail set because they are
    # also valid type qualifiers).
    if len(prefix_stripped) >= 2:
        prev_char = prefix_stripped[-2]
        double = prev_char + tail_char
        if double in {"&&", "||", "**", "*&", "&*", "->"}:
            return False
    if tail_char in DECL_PREFIX_BAD_TAIL:
        return False
    m = _LAST_IDENT_RE.search(prefix_stripped)
    if m and m.group(1) in STATEMENT_CONTEXT_KEYWORDS:
        return False
    return True


FUNC_NAME_RE = re.compile(r"\b([A-Za-z_]\w*)\s*\(")
QUALIFIER_RE = re.compile(r"\b(?:const|override|final|noexcept)\b")
CYCLO_KEYWORDS_RE = re.compile(r"\b(?:if|for|while|switch|case)\b|&&|\|\|")

# Post-C++11 features with low false-positive rates.
# Each entry is (label, compiled_regex) — the label appears in the finding.
POST_CPP11_PATTERNS: List[Tuple[str, "re.Pattern[str]"]] = [
    # C++14
    (
        "std::make_unique",
        re.compile(r"\bstd::make_unique\s*<"),
    ),
    (
        "digit-separator in numeric literal (C++14)",
        re.compile(r"(?<!')\b\d[\d']*'[\d']*\b"),
    ),
    # C++17
    (
        "if constexpr (C++17)",
        re.compile(r"\bif\s+constexpr\b"),
    ),
    (
        "structured binding auto [...] (C++17)",
        re.compile(r"\bauto\s*&?\s*\[[^\[\]]+\]\s*="),
    ),
    (
        "fold expression (C++17)",
        re.compile(
            r"\(\s*[A-Za-z_]\w*\s*(?:\+\+|\&\&|\|\||[+\-*/%^&|<>]=?|[.<>=])"
            r"\s*\.\.\.\s*\)"
            r"|\(\s*\.\.\.\s*(?:\+\+|\&\&|\|\||[+\-*/%^&|<>]=?|[.<>=])"
            r"\s*[A-Za-z_]\w*\s*\)"
            r"|\(\s*[A-Za-z_]\w*\s*(?:\+\+|\&\&|\|\||[+\-*/%^&|<>]=?|[.<>=])"
            r"\s*\.\.\.\s*(?:\+\+|\&\&|\|\||[+\-*/%^&|<>]=?|[.<>=])"
            r"\s*[^,)]+\s*\)"
        ),
    ),
    (
        "std::optional (C++17)",
        re.compile(r"\bstd::optional\b"),
    ),
    (
        "std::variant (C++17)",
        re.compile(r"\bstd::variant\b"),
    ),
    (
        "std::any (C++17)",
        re.compile(r"\bstd::any\b"),
    ),
    (
        "[[nodiscard]] (C++17)",
        re.compile(r"\[\[nodiscard\b"),
    ),
    (
        "[[maybe_unused]] (C++17)",
        re.compile(r"\[\[maybe_unused\b"),
    ),
    # C++20
    (
        "concept (C++20)",
        re.compile(r"\bconcept\b"),
    ),
    (
        "requires (C++20)",
        re.compile(r"\brequires\b"),
    ),
    (
        "consteval (C++20)",
        re.compile(r"\bconsteval\b"),
    ),
    (
        "constinit (C++20)",
        re.compile(r"\bconstinit\b"),
    ),
    (
        "coroutine co_await / co_yield / co_return (C++20)",
        re.compile(r"\bco_(?:await|yield|return)\b"),
    ),
    (
        "std::span (C++20)",
        re.compile(r"\bstd::span\b"),
    ),
    (
        "std::ranges (C++20)",
        re.compile(r"\bstd::ranges::"),
    ),
    (
        "std::format (C++20)",
        re.compile(r"\bstd::format\b"),
    ),
    # C++23
    (
        "std::expected (C++23)",
        re.compile(r"\bstd::expected\b"),
    ),
    (
        "std::print / std::println (C++23)",
        re.compile(r"\bstd::print(?:ln)?\s*\("),
    ),
]


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


def discover_files(paths: Sequence[str]) -> List[Path]:
    """Walk input paths and return source files, skipping generated dirs."""
    result: List[Path] = []
    for p in paths:
        path = Path(p)
        if path.is_file():
            if path.suffix in SOURCE_EXTENSIONS:
                result.append(path)
        elif path.is_dir():
            for sub in path.rglob("*"):
                if not sub.is_file():
                    continue
                if any(part in SKIP_DIRS for part in sub.parts):
                    continue
                if sub.suffix in SOURCE_EXTENSIONS:
                    result.append(sub)
    return result


def _is_digit_separator(content: str, quote_pos: int) -> bool:
    """Return True if the single-quote at `quote_pos` is a C++14 digit
    separator (e.g. the middle `'` in `1'000'000`).

    Heuristic: a digit separator must be preceded and followed by a
    character that can belong to a numeric literal — digits,
    hexadecimal letters (a–f, A–F) for hex literals, or the base/type
    suffix letters (uUlLbB) that typically attach to numbers after the
    body. Floating-point '.' can also appear on either side (e.g.
    `1'000.5` or `1.000'5`). Any other previous/next char means the
    quote starts a character literal.
    """
    n = len(content)
    if quote_pos <= 0 or quote_pos >= n - 1:
        return False
    prev = content[quote_pos - 1]
    nxt = content[quote_pos + 1]
    numeric_chars = set("0123456789abcdefABCDEFuUlLbBxXoO.")
    return prev in numeric_chars and nxt in numeric_chars


def strip_comments(content: str) -> str:
    """Return content with comments replaced by spaces, preserving line numbers.

    Handles // line comments and /* */ block comments. String literals are
    respected so that '/' inside strings is not mistaken for a comment.

    C++14 digit separators such as `1'000'000` are NOT treated as character
    literals: the parser inspects the characters on both sides of a `'` and
    only enters character-literal mode when the quote is genuinely opening
    a char literal like `'x'` or `'\\n'`.
    """
    out = []
    i = 0
    n = len(content)
    in_string = False
    in_char = False
    while i < n:
        c = content[i]
        if in_string:
            out.append(c)
            if c == "\\" and i + 1 < n:
                out.append(content[i + 1])
                i += 2
                continue
            if c == '"':
                in_string = False
            i += 1
            continue
        if in_char:
            out.append(c)
            if c == "\\" and i + 1 < n:
                out.append(content[i + 1])
                i += 2
                continue
            if c == "'":
                in_char = False
            i += 1
            continue
        if c == '"':
            in_string = True
            out.append(c)
            i += 1
            continue
        if c == "'":
            if _is_digit_separator(content, i):
                out.append(c)
                i += 1
                continue
            in_char = True
            out.append(c)
            i += 1
            continue
        if c == "/" and i + 1 < n:
            if content[i + 1] == "/":
                # line comment until newline
                j = content.find("\n", i)
                if j < 0:
                    j = n
                out.append(" " * (j - i))
                i = j
                continue
            if content[i + 1] == "*":
                # block comment until */
                j = content.find("*/", i + 2)
                if j < 0:
                    j = n
                else:
                    j += 2
                block = content[i:j]
                # preserve newlines
                out.append(re.sub(r"[^\n]", " ", block))
                i = j
                continue
        out.append(c)
        i += 1
    return "".join(out)


def strip_strings(content: str) -> str:
    """Return content with string and character literals erased to spaces,
    preserving line numbers.

    The result still has the same number of characters and lines, but no
    letter/word survives inside "..." or '...'. This prevents false matches
    on keywords that happen to appear inside string literals.

    C++14 digit separators such as `1'000'000` are left alone: the parser
    inspects the characters on both sides of a `'` and only enters
    character-literal mode when the quote genuinely opens a char literal
    like `'x'` or `'\\n'`.
    """
    out = []
    i = 0
    n = len(content)
    in_string = False
    in_char = False
    while i < n:
        c = content[i]
        if in_string:
            if c == "\\" and i + 1 < n:
                # keep the escape sequence length but blank it out
                out.append(" ")
                if content[i + 1] == "\n":
                    out.append("\n")
                else:
                    out.append(" ")
                i += 2
                continue
            if c == '"':
                in_string = False
                out.append('"')  # keep the boundary marker (safest for line numbers)
                i += 1
                continue
            if c == "\n":
                out.append("\n")
            else:
                out.append(" ")
            i += 1
            continue
        if in_char:
            if c == "\\" and i + 1 < n:
                out.append(" ")
                if content[i + 1] == "\n":
                    out.append("\n")
                else:
                    out.append(" ")
                i += 2
                continue
            if c == "'":
                in_char = False
                out.append("'")
                i += 1
                continue
            if c == "\n":
                out.append("\n")
            else:
                out.append(" ")
            i += 1
            continue
        if c == '"':
            in_string = True
            out.append('"')
            i += 1
            continue
        if c == "'":
            if _is_digit_separator(content, i):
                out.append(c)
                i += 1
                continue
            in_char = True
            out.append("'")
            i += 1
            continue
        out.append(c)
        i += 1
    return "".join(out)


def find_class_blocks(code: str) -> List[Tuple[int, int, str, str]]:
    """Find top-level class/struct blocks. Returns list of
    (start_line_1indexed, end_line_1indexed, kind, name).

    Walks brace matching starting from each `class X {` / `struct X {` opener.
    Excludes `template <class T>` type parameters and `enum class` scoped enums.
    """
    blocks: List[Tuple[int, int, str, str]] = []
    i = 0
    n = len(code)
    while i < n:
        m = CLASS_OPEN_RE.search(code, i)
        if not m:
            break
        # reject template type parameters: `template <class T>` or `<..., class T, ...>`
        # look at non-whitespace chars before the match
        pre_start = m.start() - 1
        while pre_start >= 0 and code[pre_start] in " \t\n\r":
            pre_start -= 1
        if pre_start >= 0 and code[pre_start] in "<,":
            i = m.end()
            continue
        # reject `enum class` / `enum struct` scoped enumerations
        # look backwards for `enum` keyword (with optional whitespace)
        enum_check_pos = pre_start
        while enum_check_pos >= 0 and code[enum_check_pos] in " \t\n\r":
            enum_check_pos -= 1
        if (
            enum_check_pos >= 3
            and code[enum_check_pos - 3:enum_check_pos + 1].lower() == "enum"
            and (enum_check_pos == 3 or not code[enum_check_pos - 4].isalnum())
        ):
            i = m.end()
            continue
        # find the opening brace after the class/struct header
        # allow inheritance clauses: class X : public Y { ... }
        brace_pos = code.find("{", m.end())
        if brace_pos < 0:
            break
        # reject if there's a `;` before the brace (forward declaration)
        if ";" in code[m.end():brace_pos]:
            i = m.end()
            continue
        depth = 1
        j = brace_pos + 1
        in_string = False
        in_char = False
        while j < n and depth > 0:
            c = code[j]
            if in_string:
                if c == "\\":
                    j += 2
                    continue
                if c == '"':
                    in_string = False
                j += 1
                continue
            if in_char:
                if c == "\\":
                    j += 2
                    continue
                if c == "'":
                    in_char = False
                j += 1
                continue
            if c == '"':
                in_string = True
                j += 1
                continue
            if c == "'":
                in_char = True
                j += 1
                continue
            if c == "{":
                depth += 1
            elif c == "}":
                depth -= 1
                if depth == 0:
                    break
            j += 1
        if depth == 0:
            start_line = code[:m.start()].count("\n") + 1
            end_line = code[:j].count("\n") + 1
            blocks.append((start_line, end_line, m.group(1), m.group(2)))
            i = j + 1
        else:
            i = m.end()
    return blocks


def detect_access_at_line(class_lines: List[str], target_idx: int, default_access: str) -> str:
    """Determine the active access specifier for line `target_idx` inside class."""
    access = default_access
    for k in range(target_idx + 1):
        m = ACCESS_RE.match(class_lines[k])
        if m:
            access = m.group(1)
    return access


def is_public_member_var(line: str) -> bool:
    """Heuristic: does this stripped code line look like a public member variable
    declaration (not a function, not a typedef, not a nested class)?"""
    s = line.strip()
    if not s.endswith(";"):
        return False
    if "(" in s or ")" in s:
        return False
    if "{" in s or "}" in s:
        return False
    if "=" in s and "(" in s:
        return False
    bad_prefixes = (
        "public:", "private:", "protected:",
        "typedef", "using", "static_assert",
        "class", "struct", "friend",
        "//", "/*", "*",
        "return", "if ", "for ", "while ", "switch ", "case ",
        "template",
    )
    if s.startswith(bad_prefixes):
        return False
    # skip macro-like lines (all caps with parens already excluded above)
    # require at least one identifier character
    if not re.search(r"[A-Za-z_]", s):
        return False
    return True


def is_static_member_var(line: str) -> bool:
    """Heuristic: does this stripped code line look like a static member
    variable declaration inside a class/struct body?

    Targets `static <type> <name>;` and `static <type> <name> = ...;`
    while excluding static member functions, static_assert, and
    static_cast (the latter two carry parentheses).
    """
    s = line.strip()
    if not s.endswith(";"):
        return False
    if not s.startswith("static "):
        return False
    if s.startswith("static_assert"):
        return False
    if "(" in s or ")" in s:
        return False
    if "{" in s or "}" in s:
        return False
    return True


def _match_var_decl(stripped_line: str) -> Optional[str]:
    """If line looks like a variable declaration `type name;` or
    `type name = ...;`, return the variable name; else None.

    Excludes pure assignments (e.g. `counter = counter + 1;`) which lack a
    leading type token. This is a best-effort heuristic, not a full parser.
    """
    s = stripped_line.strip()
    if not s or s.startswith(BLOCK_KEYWORD_PREFIXES):
        return None
    m = DECL_RE.match(s)
    return m.group(1) if m else None


# Keywords that, when they appear on a class-body line that does not yet
# end with ';', mark the start of a declaration that *may* span several
# source lines. Following lines (up to the terminating ';') should not be
# mistaken for individual member-variable declarations.  This catches
# multi-line `using`/`typedef`/`template`/namespace-qualified type aliases
# such as:
#     using value_type
#         = std::vector<int>;
# Without the lookahead, the second line is inspected in isolation and
# reads as `= std::vector<int>;` — a bogus public-static member hit.
_MULTILINE_DECL_STARTERS = (
    # Both the space-separated form (`using value_type = ...`) and the
    # bare-keyword-alone form (`typedef` at the end of a line with a
    # line-break before the type) are valid.  Using tuples of both
    # spellings catches multi-line aliases regardless of where the
    # author wraps the line.
    "using", "typedef", "template", "typename",
    "namespace", "extern", "friend",
)


def _mark_continued_decl_lines(class_lines: List[str]) -> "set":
    """Return the set of 0-based line indices inside a class body that
    belong to a multi-line declaration started by a previous line (and
    therefore must not be treated as standalone member declarations).

    A line is considered a "continuation line" when:
      * there exists a strictly earlier non-empty line at the same or
        smaller brace-nesting depth that begins with a declaration
        starter keyword (`using`, `typedef`, `template`, `typename`,
        `namespace`, `extern`, `friend`);
      * the earlier starter line does NOT terminate with ';' (i.e. the
        declaration has not been closed yet);
      * no ';' has been seen on any intervening line at the starter's
        brace-nesting depth since the starter was opened.

    The heuristic is deliberately conservative: lines are only marked as
    continuations inside the same brace-depth run so that real member
    declarations following a `using` block are never incorrectly
    suppressed.
    """
    marked: set = set()
    depth = 0
    # stack of (depth_at_starter, starter_line_idx) for currently-open
    # multi-line declaration starters. Multiple starters may coexist at
    # different brace-nesting levels; a semicolon at depth D closes all
    # starters that were opened at exactly D (conservative, but this
    # matches how class bodies declare one item per terminating ';').
    open_starters: List[Tuple[int, int]] = []
    for idx, line in enumerate(class_lines):
        prev_depth = depth
        stripped_l = line.strip()
        # (1) Mark as continuation BEFORE processing any ';' on this line:
        # the terminating-semicolon line is still part of the multi-line
        # declaration and must be suppressed in the public/static-member
        # heuristics. The starter line itself is never marked; every
        # later line that runs while any starter is open is marked.
        if open_starters and any(idx > sidx for _, sidx in open_starters):
            marked.add(idx)
        # (2) Now close starters at prev_depth if the current line
        # contains a ';' at that depth.  Closure happens AFTER the
        # marking step so the line that carries the terminating ';' is
        # still counted as the final continuation of the declaration.
        if stripped_l and ";" in stripped_l:
            closed_depths = {d for d, _ in open_starters if d == prev_depth}
            if closed_depths:
                open_starters = [
                    (d, sidx) for (d, sidx) in open_starters
                    if d not in closed_depths
                ]
        # (3) Advance brace-depth for the current line (in/out of
        # functions, nested classes, etc.).  Starter-opening happens
        # after depth tracking so the starter is recorded at its
        # prev_depth, which is the same depth at which the terminating
        # ';' will be processed.
        depth += line.count("{") - line.count("}")
        # (4) Lines that begin with declaration-starter keywords but do
        # NOT contain ';' will (in a well-formed class body) continue
        # across at least one more line before reaching the terminating
        # ';'. Record them so subsequent lines can be marked.
        if prev_depth >= 1 and stripped_l and ";" not in stripped_l \
                and stripped_l.startswith(_MULTILINE_DECL_STARTERS):
            open_starters.append((prev_depth, idx))
    return marked


def _skip_trailing_return_type(s: str, start: int) -> int:
    """Advance past `-> ReturnType` starting at position `start`.

    `start` must point at the `-` in `->`. Scans over the return type by
    tracking bracket nesting; stops and returns the position of the first
    post-type character. Recognises two stop conditions at depth 0:

      (a) current char is a direct terminator (`;`, `{`, `=`) or a
          qualifier keyword starts here — stop on the spot so the caller
          consumes it on the next loop iteration.
      (b) current char is whitespace followed by a terminator or qualifier
          keyword — stop at the whitespace position so the caller's own
          whitespace-skip advances past it cleanly.

    `:` is NOT a terminator because qualified return types like
    `std::vector<int>` contain `::` scope-resolution colons, and in any
    case a trailing-return-type function cannot also be a constructor
    (which is the only case that introduces a member initializer list).

    Caller re-runs qualifier / suffix scanning after this helper returns.
    """
    i = start + 2  # skip past `->`
    n = len(s)
    # skip whitespace between `->` and the type
    while i < n and s[i] in " \t\n":
        i += 1
    depth = 0
    while i < n:
        c = s[i]
        if c in "([{<":
            depth += 1
        elif c in ")]}>":
            depth = max(0, depth - 1)
        if depth == 0:
            # (a) terminator/qualifier attached directly to the type, e.g. `int;`
            if c in ";{=":
                break
            if QUALIFIER_RE.match(s, i):
                break
            # (b) terminator/qualifier after whitespace, e.g. `int const` / `int ;`
            if c in " \t\n":
                next_word_start = i + 1
                while next_word_start < n and s[next_word_start] in " \t\n":
                    next_word_start += 1
                if next_word_start >= n:
                    break
                nc = s[next_word_start]
                if nc in ";{=":
                    break
                if QUALIFIER_RE.match(s, next_word_start):
                    break
        i += 1
    return i


def _skip_member_initializer_list(s: str, start: int) -> int:
    """Advance past `: member(args), member{args}, ...` constructor initializer
    list starting at position `start`.

    `start` must point at the `:` that introduces the list (and be preceded by
    `)` — caller ensures this is not a `::` scope-resolution colon). Scans
    commas between members while respecting nesting; returns the position of
    the opening `{` (or `;`, `=`, or end-of-string if body-less).

    Heuristic: at depth 0, `{` can either open a brace-initializer for a
    member (e.g. `a{1,2}`) or open the function body itself. The function
    body `{` is distinguished by looking at the previous non-whitespace
    character: if it is `)` or `}` (a previous initializer ended normally)
    or `,` / `:` (list-level punctuation), then `{` belongs to the function
    body and scanning stops. Otherwise the `{` is consumed as a member
    brace-initializer (depth increases normally).
    """
    i = start + 1  # skip past `:`
    n = len(s)
    depth = 0
    while i < n:
        c = s[i]
        # Decide `{` at depth 0 before bumping depth.
        if depth == 0 and c == "{":
            prev = i - 1
            while prev >= 0 and s[prev] in " \t\n":
                prev -= 1
            prev_c = s[prev] if prev >= 0 else ""
            # `{` after full initializer / separator = function body opener
            if prev_c in ")}:,":
                break
        if c in "([{<":
            depth += 1
        elif c in ")]}>":
            depth = max(0, depth - 1)
        if depth == 0 and c in ";=":
            break
        i += 1
    return i


def count_function_params(params_str: str) -> int:
    """Count top-level parameters by counting commas at bracket depth 0.

    Tracks (), [] {} and <> as nesting (so commas inside std::map<K, V> or
    function-pointer params are not counted as param separators).
    """
    s = params_str.strip()
    if not s or s == "void":
        return 0
    depth = 0
    count = 1
    for c in s:
        if c in "([{<":
            depth += 1
        elif c in ")]}>":
            depth = max(0, depth - 1)
        elif c == "," and depth == 0:
            count += 1
    return count


def find_long_function_signatures(
    content: str, threshold: int
) -> List[Tuple[int, str, int]]:
    """Find function declarations/definitions whose parameter count exceeds threshold.

    Returns list of (line_no, function_name, param_count). Best-effort:
    requires either a return-type prefix before the function name (for `;`
    and `=` endings) or a `{` ending (for constructors/destructors which
    have no return type).
    """
    stripped = strip_comments(content)
    n = len(stripped)
    findings: List[Tuple[int, str, int]] = []

    i = 0
    while i < n:
        m = FUNC_NAME_RE.search(stripped, i)
        if not m:
            break
        name = m.group(1)
        if name in NON_FUNCTION_KEYWORDS:
            i = m.end()
            continue

        # inspect the prefix before `name(` to reject calls/lambdas/macros
        line_start = stripped.rfind("\n", 0, m.start()) + 1
        prefix = stripped[line_start:m.start()]
        prefix_stripped = prefix.rstrip()
        prefix_lstripped = prefix.lstrip()

        # skip preprocessor lines
        if prefix_lstripped.startswith("#"):
            i = m.end()
            continue
        # skip lambdas: prefix ends with `]`
        if prefix_stripped.endswith("]"):
            i = m.end()
            continue
        # skip method calls: prefix ends with `.` or `->`
        if prefix_stripped.endswith(".") or prefix_stripped.endswith("->"):
            i = m.end()
            continue
        # skip assignment-result calls: prefix ends with `=` (but not `==`)
        if prefix_stripped.endswith("=") and not prefix_stripped.endswith("=="):
            i = m.end()
            continue
        # skip function-pointer typedefs
        if "typedef" in prefix_stripped:
            i = m.end()
            continue

        # match parens to find params string
        paren_open = m.end() - 1
        depth = 1
        j = paren_open + 1
        while j < n and depth > 0:
            c = stripped[j]
            if c == "(":
                depth += 1
            elif c == ")":
                depth -= 1
                if depth == 0:
                    break
            j += 1
        if depth != 0:
            i = m.end()
            continue
        params_str = stripped[paren_open + 1:j]

        # after `)`: optional qualifiers (const/override/final/noexcept),
        # optional trailing-return type (`-> ReturnType`), optional
        # constructor member initializer list (`: a(x), b{y}`), then ; { or =
        pos = j + 1
        while pos < n and stripped[pos] in " \t\n":
            pos += 1
        while True:
            # consume qualifiers
            mq = QUALIFIER_RE.match(stripped, pos)
            if mq:
                pos = mq.end()
                while pos < n and stripped[pos] in " \t\n":
                    pos += 1
                continue
            # consume trailing return type: `-> ReturnType`
            if (
                pos + 1 < n
                and stripped[pos] == "-"
                and stripped[pos + 1] == ">"
            ):
                pos = _skip_trailing_return_type(stripped, pos)
                while pos < n and stripped[pos] in " \t\n":
                    pos += 1
                continue
            # consume constructor member initializer list: `: members...`
            # but avoid scope-resolution `::` — the prior non-whitespace
            # char must be `)` (the closing param paren we just matched)
            if pos < n and stripped[pos] == ":":
                pre_colon = pos - 1
                while pre_colon >= 0 and stripped[pre_colon] in " \t\n":
                    pre_colon -= 1
                if pre_colon == j:
                    pos = _skip_member_initializer_list(stripped, pos)
                    while pos < n and stripped[pos] in " \t\n":
                        pos += 1
                    continue
            break
        if pos >= n or stripped[pos] not in ";{=":
            i = j + 1
            continue

        # For `;` and `=`, the preceding prefix must look like a genuine
        # declaration context (return type + optional specifiers) rather
        # than the tail of an expression / argument list.
        #   - `int f(int);`              prefix `int`  -> declaration ✓
        #   - `return f(1, 2, 3, ...);`  prefix `return` -> call ✗
        #   - `foo(1, bar(x), 3);`       prefix `foo(1, ` -> call ✗ (ends with `,(`)
        #   - Constructor/destructor bodies terminate with `{` and bypass
        #     this check entirely.
        if stripped[pos] in ";=" and not _is_declaration_prefix(prefix_stripped):
            i = j + 1
            continue

        param_count = count_function_params(params_str)
        if param_count > threshold:
            line_no = stripped[:paren_open].count("\n") + 1
            findings.append((line_no, name, param_count))

        i = j + 1

    return findings


def find_function_bodies(content: str) -> List[Tuple[int, str, int, int]]:
    """Find function definitions (with body), not just declarations.

    Returns list of (signature_line_no, function_name, body_start_pos,
    body_end_pos) where positions are absolute offsets in stripped content.
    Reuses the prefix/reject logic from find_long_function_signatures.
    """
    stripped = strip_comments(content)
    n = len(stripped)
    bodies: List[Tuple[int, str, int, int]] = []

    i = 0
    while i < n:
        m = FUNC_NAME_RE.search(stripped, i)
        if not m:
            break
        name = m.group(1)
        if name in NON_FUNCTION_KEYWORDS:
            i = m.end()
            continue

        line_start = stripped.rfind("\n", 0, m.start()) + 1
        prefix = stripped[line_start:m.start()]
        prefix_stripped = prefix.rstrip()
        prefix_lstripped = prefix.lstrip()

        if prefix_lstripped.startswith("#"):
            i = m.end()
            continue
        if prefix_stripped.endswith("]"):
            i = m.end()
            continue
        if prefix_stripped.endswith(".") or prefix_stripped.endswith("->"):
            i = m.end()
            continue
        if prefix_stripped.endswith("=") and not prefix_stripped.endswith("=="):
            i = m.end()
            continue
        if "typedef" in prefix_stripped:
            i = m.end()
            continue

        paren_open = m.end() - 1
        depth = 1
        j = paren_open + 1
        while j < n and depth > 0:
            c = stripped[j]
            if c == "(":
                depth += 1
            elif c == ")":
                depth -= 1
                if depth == 0:
                    break
            j += 1
        if depth != 0:
            i = m.end()
            continue

        # find what comes after `)`: skip whitespace, qualifiers, optional
        # trailing-return type (`-> ReturnType`), and optional constructor
        # member initializer list (`: a(x), b{y}`).
        pos = j + 1
        while pos < n and stripped[pos] in " \t\n":
            pos += 1
        while True:
            # consume qualifiers
            mq = QUALIFIER_RE.match(stripped, pos)
            if mq:
                pos = mq.end()
                while pos < n and stripped[pos] in " \t\n":
                    pos += 1
                continue
            # consume trailing return type: `-> ReturnType`
            if (
                pos + 1 < n
                and stripped[pos] == "-"
                and stripped[pos + 1] == ">"
            ):
                pos = _skip_trailing_return_type(stripped, pos)
                while pos < n and stripped[pos] in " \t\n":
                    pos += 1
                continue
            # consume constructor member initializer list: `: members...`
            # but avoid scope-resolution `::` — the prior non-whitespace
            # char must be `)` (the closing param paren we just matched)
            if pos < n and stripped[pos] == ":":
                pre_colon = pos - 1
                while pre_colon >= 0 and stripped[pre_colon] in " \t\n":
                    pre_colon -= 1
                if pre_colon == j:
                    pos = _skip_member_initializer_list(stripped, pos)
                    while pos < n and stripped[pos] in " \t\n":
                        pos += 1
                    continue
            break

        # we need a `{` body (not `;` declaration, not `= 0` pure virtual)
        if pos >= n or stripped[pos] != "{":
            i = j + 1
            continue

        # match braces to find body end
        body_open = pos
        depth = 1
        body_close = body_open + 1
        while body_close < n and depth > 0:
            c = stripped[body_close]
            if c == "{":
                depth += 1
            elif c == "}":
                depth -= 1
                if depth == 0:
                    break
            body_close += 1
        if depth != 0:
            i = j + 1
            continue

        sig_line_no = stripped[:paren_open].count("\n") + 1
        bodies.append((sig_line_no, name, body_open, body_close))
        i = body_close + 1

    return bodies


def find_high_complexity_functions(
    content: str, threshold: int
) -> List[Tuple[int, str, int]]:
    """Find functions whose cyclomatic complexity exceeds threshold.

    Cyclomatic complexity counts: if, for, while, switch, case, &&, ||.
    Returns list of (line_no, function_name, complexity).

    String and character literals are blanked before the regex runs so
    that control-flow words appearing inside messages (e.g.
    `const char* msg = "if not ok return";`) are not counted as real
    control flow statements.
    """
    stripped = strip_strings(strip_comments(content))
    findings: List[Tuple[int, str, int]] = []
    for sig_line, name, body_open, body_close in find_function_bodies(content):
        body = stripped[body_open + 1:body_close]
        complexity = len(CYCLO_KEYWORDS_RE.findall(body))
        if complexity > threshold:
            findings.append((sig_line, name, complexity))
    return findings


def find_post_cpp11_features(
    content: str,
) -> List[Tuple[int, str]]:
    """Find the first occurrence of each post-C++11 feature.

    Uses the patterns in POST_CPP11_PATTERNS against content with both
    comments AND string literals blanked out. String blanking prevents
    false matches on keywords embedded in user-facing messages (e.g.
    `"X requires Y"` in a WARNING_QUIT call).

    Returns a list of (line_no, feature_label) — one entry per distinct
    detected feature so the finding message is informative (the deduction
    is one-shot -80 per file regardless of how many features hit).
    """
    stripped = strip_strings(strip_comments(content))
    hits: List[Tuple[int, str]] = []
    for label, regex in POST_CPP11_PATTERNS:
        m = regex.search(stripped)
        if m:
            line_no = stripped[:m.start()].count("\n") + 1
            hits.append((line_no, label))
    return hits


def analyze_class_blocks(
    content: str,
) -> Tuple[List[Finding], List[Finding], List[Finding], List[Finding]]:
    """Analyze class/struct blocks for public member variables, static member
    variables, long member functions, and member/local name conflicts.

    Returns (public_member_findings, static_member_findings,
              long_function_findings, name_conflict_findings).
    """
    pub_findings: List[Finding] = []
    static_findings: List[Finding] = []
    long_func_findings: List[Finding] = []
    conflict_findings: List[Finding] = []

    stripped = strip_comments(content)
    lines = stripped.split("\n")

    for start_line, end_line, kind, name in find_class_blocks(stripped):
        default_access = "private" if kind == "class" else "public"
        body_start = start_line - 1
        body_end = end_line - 1
        class_lines = lines[body_start:body_end + 1]

        # pre-compute access per line (for public member variable rule)
        access_per_line: List[str] = []
        cur_access = default_access
        for k in range(len(class_lines)):
            m = ACCESS_RE.match(class_lines[k])
            if m:
                cur_access = m.group(1)
            access_per_line.append(cur_access)

        # pass 0: identify lines belonging to multi-line declarations that
        # start with keywords such as `using`, `typedef`, `template`, etc.
        # Continuation lines (e.g. the "= std::vector<int>;" half of a
        # two-line `using value_type = ...;` alias) would otherwise be
        # treated as standalone member declarations when the per-line
        # heuristics run below.
        continued_idxs = _mark_continued_decl_lines(class_lines)

        # pass 1: collect all member variable names (any access section).
        # depth starts at 0; the class header line brings it to 1. Only lines
        # that begin AND end at depth 1 are class-body member declarations
        # (lines that open a function body bring depth from 1 to 2).
        member_var_names: set = set()
        depth = 0
        for idx, line in enumerate(class_lines):
            prev_depth = depth
            depth += line.count("{") - line.count("}")
            if prev_depth == 1 and depth == 1 and idx not in continued_idxs:
                var_name = _match_var_decl(line)
                if var_name:
                    member_var_names.add(var_name)

        # pass 2: walk class body for findings
        depth = 0
        func_start_abs = -1
        for idx in range(len(class_lines)):
            line = class_lines[idx]
            abs_line = body_start + idx + 1  # 1-indexed absolute line
            prev_depth = depth

            # public member variable: only at depth 1 (class body, not in a
            # function). Lines flagged as continuations of an earlier
            # multi-line declaration are skipped regardless of how they
            # look in isolation — e.g. `    = std::vector<int>;` after a
            # two-line `using value_type` alias must not count as a member.
            # Only `class` bodies are penalised: in C++ a `struct` has
            # public access by default and public data members are a
            # legitimate, intended usage (POD aggregate, value types,
            # mixin tags). Treating them as a code-quality issue would
            # penalise perfectly idiomatic C++.
            if prev_depth == 1 and access_per_line[idx] == "public" \
                    and kind == "class" \
                    and idx not in continued_idxs:
                if is_public_member_var(line):
                    pub_findings.append(Finding(
                        rule="public_member_variable",
                        line=abs_line,
                        reason=f"public member in {kind} {name}: {line.strip()}",
                        deduction=WEIGHTS["public_member_variable"],
                    ))

            # static member variable: depth 1, any access section.
            # Excludes static_assert, static member functions, and
            # static_cast via is_static_member_var's heuristic. Continued
            # multi-line declarations are skipped for the same reason
            # documented above for the public-member rule.
            if prev_depth == 1 and idx not in continued_idxs \
                    and is_static_member_var(line):
                static_findings.append(Finding(
                    rule="static_member_variable",
                    line=abs_line,
                    reason=f"static member in {kind} {name}: {line.strip()}",
                    deduction=WEIGHTS["static_member_variable"],
                ))

            # member/local name conflict: only at depth >= 2 (function body)
            if prev_depth >= 2:
                var_name = _match_var_decl(line)
                if var_name and var_name in member_var_names:
                    conflict_findings.append(Finding(
                        rule="member_local_name_conflict",
                        line=abs_line,
                        reason=(
                            f"local variable '{var_name}' in {kind} {name} "
                            f"shadows a member variable"
                        ),
                        deduction=WEIGHTS["member_local_name_conflict"],
                    ))

            # apply brace counting for this line
            depth = prev_depth + line.count("{") - line.count("}")

            # function body start: transition from depth 1 to >=2
            if prev_depth == 1 and depth >= 2 and func_start_abs < 0:
                func_start_abs = abs_line
            # function body end: transition from >=2 back to 1
            if prev_depth >= 2 and depth == 1 and func_start_abs > 0:
                func_end_abs = abs_line
                func_lines = func_end_abs - func_start_abs + 1
                if func_lines > FUNCTION_LENGTH_THRESHOLD:
                    excess = func_lines - FUNCTION_LENGTH_THRESHOLD
                    blocks = (excess + FUNCTION_LENGTH_STEP - 1) // FUNCTION_LENGTH_STEP
                    long_func_findings.append(Finding(
                        rule="member_function_too_long",
                        line=func_start_abs,
                        reason=(
                            f"member function in {kind} {name} spans {func_lines} lines "
                            f"(exceeds {FUNCTION_LENGTH_THRESHOLD})"
                        ),
                        deduction=blocks * WEIGHTS["member_function_too_long"],
                    ))
                func_start_abs = -1

    return pub_findings, static_findings, long_func_findings, conflict_findings


def analyze_file(path: Path) -> FileReport:
    """Analyze one file and return its FileReport."""
    findings: List[Finding] = []
    name = path.name
    stem = path.stem
    suffix = path.suffix

    # filename rules
    if len(stem) > FILENAME_LENGTH_LIMIT:
        findings.append(Finding(
            rule="filename_too_long",
            line=None,
            reason=(
                f"filename '{name}' stem has {len(stem)} chars "
                f"(limit {FILENAME_LENGTH_LIMIT})"
            ),
            deduction=WEIGHTS["filename_too_long"],
        ))
    if any(c.isupper() for c in stem):
        findings.append(Finding(
            rule="filename_uppercase",
            line=None,
            reason=f"filename '{name}' contains uppercase letters",
            deduction=WEIGHTS["filename_uppercase"],
        ))
    if suffix == ".hpp":
        findings.append(Finding(
            rule="hpp_implementation",
            line=None,
            reason=".hpp implementation header prohibited (use .cpp + .h split)",
            deduction=WEIGHTS["hpp_implementation"],
        ))

    # read content
    try:
        content = path.read_text(encoding="utf-8", errors="replace")
    except OSError as e:
        findings.append(Finding(
            rule="read_error",
            line=None,
            reason=f"failed to read file: {e}",
            deduction=0,
        ))
        return FileReport(path=str(path), score=100, findings=findings)

    lines = content.split("\n")

    # file length rule (no cap: large files legitimately scale deduction)
    total_lines = len(lines)
    if total_lines > FILE_LENGTH_THRESHOLD:
        excess = total_lines - FILE_LENGTH_THRESHOLD
        blocks = (excess + FILE_LENGTH_STEP - 1) // FILE_LENGTH_STEP
        findings.append(Finding(
            rule="file_too_long",
            line=None,
            reason=(
                f"file has {total_lines} lines "
                f"(exceeds {FILE_LENGTH_THRESHOLD} by {excess})"
            ),
            deduction=blocks * WEIGHTS["file_too_long"],
        ))

    # line-based rules with caps
    tab_count = sum(1 for l in lines if "\t" in l)
    long_line_count = sum(1 for l in lines if len(l) > LINE_LENGTH_LIMIT)
    using_ns_count = sum(
        1 for l in lines if USING_NS_STD_RE.search(strip_comments(l))
    )
    chinese_count = sum(1 for l in lines if CHINESE_RE.search(l))

    stripped_content = strip_comments(content)
    # uppercase constants: strip strings too (so default-value strings like
    # "ABACUS" do not trigger), and the regex excludes member accesses
    # (e.g. GlobalV::MY_RANK, x.UPPER) via a negative lookbehind.
    stripped_for_upper = strip_strings(stripped_content)
    upper_const_count = sum(
        len(UPPERCASE_CONST_RE.findall(l))
        for l in stripped_for_upper.split("\n")
    )

    # interface & dependency rules
    global_dep_count = len(GLOBAL_DEPENDENCY_RE.findall(stripped_content))
    hpp_include_count = sum(1 for l in lines if HPP_INCLUDE_RE.match(l))
    friend_count = len(FRIEND_RE.findall(stripped_content))
    goto_count = len(GOTO_RE.findall(stripped_content))

    # default parameter: scan the full stripped content so multi-line
    # declarations (parameter list spanning multiple lines) are detected;
    # `[^();]*` in the regex already matches newlines.
    default_param_count = len(DEFAULT_PARAM_RE.findall(stripped_content))

    # unpaired new/delete (file-level): rough heuristic, cap applied below.
    # All `new` expressions are counted first. `delete` expressions are the
    # obvious symmetric ownership cancel, but there is a second class of
    # cancellations that never produces a `delete` keyword at all: a
    # C++11 (or later) smart pointer that owns the object from construction.
    # Recognised ownership-transfer patterns are aggregated in OWNED_NEW_RE
    # and subtracted before comparing with `delete` so the leak heuristic
    # stays meaningful on modern code that prefers std::unique_ptr /
    # std::shared_ptr to manual delete. Examples that now cancel:
    #     ptr.reset(new Foo());
    #     std::unique_ptr<Foo> p(new Foo(arg));
    # These cancellations deliberately do NOT touch `raw_new_count` below
    # because the "raw_new_keyword" rule charges per direct `new` regardless
    # of ownership — it is the stylistic counterpart to std::make_unique /
    # std::make_shared, so `reset(new Foo)` and `unique_ptr<T>(new T)`
    # must still contribute to that count.
    new_count = len(NEW_EXPR_RE.findall(stripped_content))
    delete_count = len(DELETE_EXPR_RE.findall(stripped_content))
    owned_new_count = len(OWNED_NEW_RE.findall(stripped_content))
    cancelled_new = delete_count + owned_new_count
    unpaired_new = max(0, new_count - cancelled_new)

    # raw `new` keyword usage: each occurrence costs 1 (no cap)
    raw_new_count = new_count

    def append_capped(rule: str, count: int) -> None:
        if count == 0:
            return
        cap = CAPS.get(rule)
        actual = min(count, cap) if cap else count
        cap_text = f" (capped at {cap})" if cap and count > cap else ""
        findings.append(Finding(
            rule=rule,
            line=None,
            reason=f"{count} occurrence(s){cap_text}",
            deduction=actual * WEIGHTS[rule],
        ))

    append_capped("tab_indentation", tab_count)
    append_capped("line_too_long", long_line_count)
    append_capped("using_namespace_std", using_ns_count)
    append_capped("chinese_comment", chinese_count)
    append_capped("uppercase_constant", upper_const_count)
    append_capped("default_parameter", default_param_count)
    append_capped("global_dependency", global_dep_count)
    append_capped("hpp_include", hpp_include_count)
    append_capped("friend_keyword", friend_count)
    append_capped("unpaired_new_delete", unpaired_new)
    append_capped("raw_new_keyword", raw_new_count)
    append_capped("goto_keyword", goto_count)

    # too-many-parameters rule: each function's deduction is capped
    # individually (per-function cap, not per-file). A file with many
    # over-parameterised functions accumulates deductions across all of
    # them, instead of one bad function hiding the rest.
    long_param_funcs = find_long_function_signatures(content, FUNCTION_PARAM_THRESHOLD)
    cap_params = CAPS.get("too_many_parameters")
    for line_no, fname, pcount in long_param_funcs:
        excess = pcount - FUNCTION_PARAM_THRESHOLD
        raw_deduction = excess * WEIGHTS["too_many_parameters"]
        if cap_params is not None and raw_deduction > cap_params:
            per_deduction = cap_params
            cap_text = f" (capped at {cap_params})"
        else:
            per_deduction = raw_deduction
            cap_text = ""
        findings.append(Finding(
            rule="too_many_parameters",
            line=line_no,
            reason=(
                f"function '{fname}' has {pcount} parameters "
                f"(exceeds {FUNCTION_PARAM_THRESHOLD} by {excess}){cap_text}"
            ),
            deduction=per_deduction,
        ))

    # high cyclomatic complexity rule: each function's deduction is capped
    # individually (per-function cap, not per-file).
    high_cyclo_funcs = find_high_complexity_functions(content, CYCLO_THRESHOLD)
    cap_cyclo = CAPS.get("high_cyclomatic_complexity")
    for line_no, fname, cyclo in high_cyclo_funcs:
        excess = cyclo - CYCLO_THRESHOLD
        raw_deduction = excess * WEIGHTS["high_cyclomatic_complexity"]
        if cap_cyclo is not None and raw_deduction > cap_cyclo:
            per_deduction = cap_cyclo
            cap_text = f" (capped at {cap_cyclo})"
        else:
            per_deduction = raw_deduction
            cap_text = ""
        findings.append(Finding(
            rule="high_cyclomatic_complexity",
            line=line_no,
            reason=(
                f"function '{fname}' has cyclomatic complexity {cyclo} "
                f"(exceeds {CYCLO_THRESHOLD} by {excess}){cap_text}"
            ),
            deduction=per_deduction,
        ))

    # post-C++11 feature rule: base deduction for any usage + per-feature
    # incremental deduction (capped). Total cap matches the previous
    # one-shot -80, but a single misuse no longer swamps the whole file.
    post_cpp11_hits = find_post_cpp11_features(content)
    if post_cpp11_hits:
        feature_list = ", ".join(
            f"'{label}' (line {ln})" for ln, label in post_cpp11_hits
        )
        findings.append(Finding(
            rule="post_cpp11_feature",
            line=post_cpp11_hits[0][0],
            reason=(
                f"uses post-C++11 feature(s): {feature_list} "
                f"(repo baseline is C++11)"
            ),
            deduction=WEIGHTS["post_cpp11_feature"],
        ))
        num_features = len(post_cpp11_hits)
        cap_features = CAPS.get("post_cpp11_per_feature", num_features)
        capped = min(num_features, cap_features)
        if capped > 0:
            cap_text = (
                f" (capped at {cap_features})" if num_features > cap_features else ""
            )
            findings.append(Finding(
                rule="post_cpp11_per_feature",
                line=post_cpp11_hits[0][0],
                reason=(
                    f"{num_features} distinct post-C++11 feature(s) detected"
                    f"{cap_text}"
                ),
                deduction=capped * WEIGHTS["post_cpp11_per_feature"],
            ))

    # class-based rules
    (
        pub_findings,
        static_findings,
        long_func_findings,
        conflict_findings,
    ) = analyze_class_blocks(content)
    findings.extend(pub_findings)
    # static member findings: cap per file
    if len(static_findings) > CAPS.get("static_member_variable", len(static_findings)):
        cap_static = CAPS["static_member_variable"]
        static_findings = static_findings[:cap_static]
    findings.extend(static_findings)
    findings.extend(long_func_findings)
    # name conflict findings already respect cap via per-file cap on the rule
    if len(conflict_findings) > CAPS.get("member_local_name_conflict", len(conflict_findings)):
        cap = CAPS["member_local_name_conflict"]
        conflict_findings = conflict_findings[:cap]
    findings.extend(conflict_findings)

    # compute score: real_score is the unclamped value (may be negative);
    # score is clamped to a minimum of 0 for display.
    real_score = 100
    for f in findings:
        real_score -= f.deduction
    score = max(0, real_score)

    return FileReport(
        path=str(path),
        score=score,
        findings=findings,
        real_score=real_score,
    )


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
