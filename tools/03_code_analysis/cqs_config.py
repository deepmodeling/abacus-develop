"""Shared configuration for the ABACUS code-quality scorer.

Single source of truth for: source extensions, skip dirs, per-rule
weights/caps, numeric thresholds, all compiled regexes, and the two
prefix-classification helpers used by multiple rule detectors.

Rules documentation lives in the entry script's module docstring
(`code_quality_score.py`) — keep both in sync when adding a rule.
"""

from __future__ import annotations

import re
from typing import List, Tuple

SOURCE_EXTENSIONS = {
    ".c", ".cc", ".cpp", ".cxx", ".cu", ".h", ".hh", ".hpp", ".hxx", ".cuh",
}

SKIP_DIRS = {
    ".git", "build", "__pycache__", "node_modules", ".cache",
    "third_party", "thirdparty", ".vscode", ".idea", ".trae-cn",
    "Dependencies",
    "test", "tests", "test_serial", "test_parallel", "test_gpu",
    "unit_test", "unittest", "unittests",
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
    "auto_keyword": 1,
    "static_member_variable": 1,
    "indented_preprocessor": 1,
    "complex_argument_expression": 1,
    "ifdef_in_call_args": 1,
    "comment_above_guard": 1,
    "doc_file_mismatch": 1,
    "duplicate_doc_block": 1,
    "header_too_long": 1,
    "header_function_body": 1,
    "new_macro_definition": 3,
    "param_member_init": 5,
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
    "auto_keyword": 10,
    "static_member_variable": 10,
    "indented_preprocessor": 5,
    "complex_argument_expression": 5,
    "ifdef_in_call_args": 5,
    "comment_above_guard": 3,
    "duplicate_doc_block": 5,
    "header_too_long": 5,
    "header_function_body": 5,
    "new_macro_definition": 10,
    "param_member_init": 10,
}

FUNCTION_LENGTH_THRESHOLD = 50
FUNCTION_LENGTH_STEP = 50
FILE_LENGTH_THRESHOLD = 500
FILE_LENGTH_STEP = 50
HEADER_LENGTH_THRESHOLD = 120
HEADER_LENGTH_STEP = 10
HEADER_FUNC_BODY_THRESHOLD = 20
HEADER_FUNC_BODY_STEP = 10
FUNCTION_PARAM_THRESHOLD = 7
CYCLO_THRESHOLD = 10
LINE_LENGTH_LIMIT = 120
FILENAME_LENGTH_LIMIT = 20
PASS_THRESHOLD = 60

CHINESE_RE = re.compile("[一-鿿]")
USING_NS_STD_RE = re.compile(r"\busing\s+namespace\s+\w[\w:]*\s*;")
UPPERCASE_CONST_RE = re.compile(r"(?<![.:>])\b[A-Z][A-Z0-9_]{2,}\b(?![.:])")
ACCESS_RE = re.compile(r"^\s*(public|private|protected)\s*:")
CLASS_OPEN_RE = re.compile(r"\b(class|struct)\s+(\w+)\b")

GLOBAL_DEPENDENCY_RE = re.compile(r"\b(?:GlobalV::|GlobalC::|PARAM(?:\.|->|::))")
# Member-variable declaration initialized from PARAM inside a class body
# (`.h`/`.hpp` rule 21). Matches `Type name = PARAM...;`,
# `Type name = <expr> PARAM...;` and `Type name{PARAM...};`, also when
# the declaration is split across lines (`\s*` spans newlines). The gap
# between the `=`/`{` and `PARAM` excludes `;{}()` so a match can never
# cross a statement boundary or dive into another brace/paren level;
# class-scope filtering (member function bodies, parameter lists) is
# done by the caller via brace/paren depth.
PARAM_MEMBER_INIT_RE = re.compile(
    r"(?:static\s+)?(?:const\s+)?[\w:<>,\s*&]+?\s+\w+\s*"
    r"[={]\s*[^;={}()]*?PARAM\s*(?:\.|->|::)"
)
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
# `auto` type deduction keyword. `\bauto\b` does not match `auto_ptr`
# (`_` is a word character, so there is no boundary before it).
AUTO_RE = re.compile(r"\bauto\b")

# Preprocessor directive that does not start at column 0 (leading spaces/tabs
# before the '#'). Matched against raw lines; string literals virtually never
# contain a line-start '#', so no string stripping is needed here.
INDENTED_PREPROC_RE = re.compile(
    r"^[ \t]+#\s*(?:if|ifdef|ifndef|define|endif|include|undef|pragma|else|elif)\b"
)
# First include guard / pragma-once line of a header.
GUARD_RE = re.compile(r"^\s*#\s*(?:ifndef|pragma\s+once)\b")
# New `#define NAME ...` directive; NAME is captured so guard-style
# macros (NAME ending in _H / _HPP with no parameter list or value) can
# be exempted by the caller.
MACRO_DEFINE_RE = re.compile(r"^\s*#\s*define\s+([A-Za-z_]\w*)")
# `@file name` inside a doxygen comment.
DOC_FILE_RE = re.compile(r"@file\s+(\S+)")
# Comment-only line (/// ... or a line inside a /* */ block is not fully
# distinguishable line-wise, so we only treat leading // or /* as comment).
COMMENT_LINE_RE = re.compile(r"^\s*(?://|/\*|\*)")

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

# Preprocessor conditional opener / closer, matched against raw lines.
PREPROC_IF_RE = re.compile(r"^\s*#\s*(?:if|ifdef|ifndef)\b")
PREPROC_ENDIF_RE = re.compile(r"^\s*#\s*endif\b")

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
