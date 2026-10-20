"""Rule-detection functions for the ABACUS code-quality scorer.

Every public `find_*`/`analyze_*` function scans source text and returns
violation locations (line numbers plus optional detail strings). Rule
weights, caps and regexes live in `cqs_config`.
"""

from __future__ import annotations

import re
from typing import Dict, List, Optional, Set, Tuple

from cqs_config import (
    ACCESS_RE, AUTO_RE, BLOCK_KEYWORD_PREFIXES, CLASS_OPEN_RE,
    CYCLO_KEYWORDS_RE, CYCLO_THRESHOLD, DECL_RE, DEFAULT_PARAM_RE,
    FUNC_NAME_RE, FUNCTION_LENGTH_STEP, FUNCTION_LENGTH_THRESHOLD,
    FUNCTION_PARAM_THRESHOLD,
    GLOBAL_DEPENDENCY_RE, HEADER_FUNC_BODY_THRESHOLD, HPP_INCLUDE_RE,
    INDENTED_PREPROC_RE,
    NEW_EXPR_RE, NON_FUNCTION_KEYWORDS, PARAM_MEMBER_INIT_RE,
    POST_CPP11_PATTERNS,
    PREPROC_IF_RE, QUALIFIER_RE, STATEMENT_CONTEXT_KEYWORDS, WEIGHTS,
    _is_declaration_prefix, _LAST_IDENT_RE,
)
from cqs_report import Finding
from cqs_text import strip_comments, strip_strings


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


_TEMPLATE_DECL_RE = re.compile(r"\btemplate\b\s*<")
_INLINE_RE = re.compile(r"\binline\b")
# how many consecutive non-blank lines above a signature are examined
# when looking for a (possibly multi-line) `template <...>` header
_TEMPLATE_HEADER_SCAN_LINES = 16


def find_header_function_bodies(
    content: str,
) -> List[Tuple[int, str, int]]:
    """Find non-inline, non-template function bodies in a .h header that
    exceed HEADER_FUNC_BODY_THRESHOLD lines.

    Exempts functions declared `inline` and template functions (they
    must live in the header): the `template <...>` header is recognised
    on the signature line itself, on the consecutive non-blank lines
    above the signature (in any line-splitting), and for member
    functions of a templated class. Short getter/setter-style bodies at
    or below the threshold are exempt by construction. Returns a list of
    (signature_line_no, function_name, body_line_count).
    """
    stripped = strip_comments(content)
    lines = stripped.split("\n")
    class_blocks = find_class_blocks(stripped)
    results: List[Tuple[int, str, int]] = []
    for sig_line, name, body_open, body_close in find_function_bodies(
        content
    ):
        body_lines = stripped[body_open:body_close].count("\n")
        if body_lines <= HEADER_FUNC_BODY_THRESHOLD:
            continue
        # exempt `inline ... name(...) {`: check the declaration text
        decl_text = _declaration_prefix_before(stripped, body_open)
        if _INLINE_RE.search(decl_text):
            continue
        if _is_template_context(lines, sig_line, decl_text, class_blocks):
            continue
        results.append((sig_line, name, body_lines))
    return results


def _is_template_context(
    lines: List[str],
    sig_line: int,
    decl_text: str,
    class_blocks: List[Tuple[int, int, str, str]],
) -> bool:
    """Return True when the function whose signature starts on 1-based
    `sig_line` is introduced by a `template <...>` header.

    Recognised forms (formatting must not change the result):
      * the `template <...>` header shares the signature line:
        `template <class T> void f(T& v) { ... }`;
      * it sits on the consecutive non-blank lines above the signature,
        split across any number of lines — as long as no `;` appears
        between the `template` keyword and the signature (a `;` would
        close an unrelated declaration in between);
      * the function is a member of a class introduced by a
        `template <...>` header.
    """
    if _TEMPLATE_DECL_RE.search(decl_text):
        return True
    # consecutive non-blank lines above the signature; only the text
    # after the last `;` belongs to this declaration's header
    buf = ""
    idx = sig_line - 2  # line index just above the signature
    steps = 0
    while (
        idx >= 0 and lines[idx].strip()
        and steps < _TEMPLATE_HEADER_SCAN_LINES
    ):
        buf = lines[idx] + "\n" + buf
        if _TEMPLATE_DECL_RE.search(buf.rsplit(";", 1)[-1]):
            return True
        idx -= 1
        steps += 1
    # member of a templated class: check the innermost enclosing class
    for start_line, end_line, _kind, _name in class_blocks:
        if start_line <= sig_line <= end_line:
            cls_decl = lines[start_line - 1]
            jdx = start_line - 2
            steps = 0
            while (
                jdx >= 0 and lines[jdx].strip()
                and steps < _TEMPLATE_HEADER_SCAN_LINES
            ):
                cls_decl = lines[jdx] + "\n" + cls_decl
                jdx -= 1
                steps += 1
            if _TEMPLATE_DECL_RE.search(cls_decl.rsplit(";", 1)[-1]):
                return True
            break
    return False


def find_param_member_inits(content: str) -> List[int]:
    """Find class member variables initialized from PARAM in a header
    (AGENTS.md rule 17 / 21).

    Scans only inside top-level class/struct bodies and only at class
    scope: a match must sit at the class body's own brace level (not
    inside a member function body or a nested class) and outside any
    parentheses (not a function default argument). The whole file is
    comment- AND string-stripped first so declarations inside comment
    blocks or string literals are never matched.

    Returns a list of 1-indexed line numbers of offending declarations.
    """
    stripped = strip_strings(strip_comments(content))
    line_starts = [0]
    for line in stripped.split("\n"):
        line_starts.append(line_starts[-1] + len(line) + 1)
    hits: List[int] = []
    for start_line, end_line, _kind, _name in find_class_blocks(stripped):
        body_lo = line_starts[start_line - 1]
        body_hi = line_starts[end_line]
        body = stripped[body_lo:body_hi]
        # the class body opens with the first `{` of the block (an
        # inheritance clause carries no braces)
        brace_pos = body.find("{")
        if brace_pos < 0:
            continue
        for m in PARAM_MEMBER_INIT_RE.finditer(body):
            seg = body[brace_pos + 1:m.start()]
            brace_depth = 1 + seg.count("{") - seg.count("}")
            paren_depth = seg.count("(") - seg.count(")")
            if brace_depth != 1 or paren_depth != 0:
                continue
            # the type part of the pattern absorbs the whitespace and
            # newline preceding the declaration, so skip it to report
            # the line where the declaration text starts
            decl_start = m.start()
            while decl_start < len(body) and body[decl_start] in " \t\n":
                decl_start += 1
            line_no = stripped[:body_lo + decl_start].count("\n") + 1
            hits.append(line_no)
    return hits


def _declaration_prefix_before(stripped: str, body_open: int) -> str:
    """Return the declaration text preceding a function body `{`, walking
    back over the parameter list to the start of the signature line."""
    depth = 0
    i = body_open
    while i > 0:
        c = stripped[i - 1]
        if c == ")":
            depth += 1
        elif c == "(":
            depth -= 1
            if depth == 0:
                break
        i -= 1
    line_start = stripped.rfind("\n", 0, i) + 1
    return stripped[line_start:body_open]


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


# Preprocessor conditional opener / closer, matched against raw lines.
PREPROC_IF_RE = re.compile(r"^\s*#\s*(?:if|ifdef|ifndef)\b")
PREPROC_ENDIF_RE = re.compile(r"^\s*#\s*endif\b")


def find_complex_call_arguments(content: str) -> List[Tuple[int, str]]:
    """Find function calls that embed a ternary `?:`, `&&` or `||`
    expression in a top-level argument (AGENTS.md §14).

    Reuses the prefix/reject logic of find_long_function_signatures to
    skip declarations, definitions, control-statement headers, lambdas
    and method-call chains. After matching a call's parentheses, only
    `?`, `&&`, `||` at the call's own paren depth are attributed to
    this call; operators nested inside an inner call's parentheses are
    judged when the inner call is visited (the scan resumes inside the
    outer call's argument list so nested callees are never skipped).
    Angle brackets `<>` of template argument lists inside the argument
    list are treated as nesting only in their "glued" form (`Foo<int>`);
    a spaced `<`/`>` is a comparison (`a < b ? c : d`) and does not
    nest, so the ternary is still attributed to the call.

    Returns a list of (line_no, callee_name) for each offending call.
    """
    stripped = strip_strings(strip_comments(content))
    n = len(stripped)
    hits: List[Tuple[int, str]] = []

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
        # NOTE: do NOT skip on a trailing `=` — `int b = add(...)` is a
        # call site, not a declaration; declarations are filtered later
        # by the return-type-prefix check against _is_declaration_prefix.
        if "typedef" in prefix_stripped:
            i = m.end()
            continue

        # match parens to find the argument list
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
        args_str = stripped[paren_open + 1:j]

        # skip declaration/definition contexts: the prefix (text between
        # the line start and the callee name) must look like a return
        # type, e.g. `void f(...)` / `int add(...)` / `unsigned long
        # f(...)`. A prefix ending in `)` (`(cast) f(...)`, `if (x) f(...)`)
        # or `=` (`int b = add(...)`) or a statement keyword (`return
        # f(...)`) is a call site and is NOT skipped — _is_declaration_prefix
        # already encodes the keyword/punctuation part of this.
        if not prefix_stripped.endswith(")") \
                and _is_declaration_prefix(prefix_stripped):
            i = j + 1
            continue

        # scan arguments for `?`, `&&`, `||` at this call's own depth:
        # parens/brackets/braces dive one level; glued template angle
        # brackets (`Foo<int>`) dive too, while spaced comparison
        # operators (`a < b`) do not
        arg_depth = 0
        complex_hit = False
        k = 0
        while k < len(args_str):
            c = args_str[k]
            if c == "<":
                # template opener only when glued to the preceding word
                # (`Foo<int>`); a spaced `<` is a comparison (`a < b`)
                # and `<<` is a shift operator
                if (
                    arg_depth == 0
                    and k > 0
                    and (args_str[k - 1].isalnum()
                         or args_str[k - 1] in "_<")
                    and k + 1 < len(args_str)
                    and args_str[k + 1] not in " \t\n<"
                ):
                    arg_depth += 1
            elif c == ">":
                # template closer only when glued to the preceding word;
                # a spaced `>` is a comparison (`a > b`), `->` is a
                # member access and `>=`/`>>=m` are operators
                if (
                    arg_depth > 0
                    and k > 0
                    and (args_str[k - 1].isalnum()
                         or args_str[k - 1] in "_<>")
                    and (k + 1 >= len(args_str)
                         or args_str[k + 1] in "(),;&*> \t\n")
                ):
                    arg_depth -= 1
            elif c in "([{":
                arg_depth += 1
            elif c in ")]}":
                arg_depth = max(0, arg_depth - 1)
            elif arg_depth == 0:
                if c == "?":
                    complex_hit = True
                    break
                if (c == "&" or c == "|") and k + 1 < len(args_str) \
                        and args_str[k + 1] == c:
                    complex_hit = True
                    break
            k += 1

        if complex_hit:
            line_no = stripped[:paren_open].count("\n") + 1
            hits.append((line_no, name))

        # resume INSIDE the argument list so nested callees (e.g. the
        # `g` of `f(g(a ? b : c))`) are visited as their own calls
        i = paren_open + 1

    return hits


def find_ifdef_in_call_args(content: str) -> List[int]:
    """Find `#if`/`#ifdef`/`#ifndef` directives that sit inside a
    function call's argument list (conditional compilation wrapping
    individual arguments).

    A single pass tracks a paren stack: each `(` opener is classified
    as a *call* context or an *other* context (control-statement header
    `if (/for (/while (/switch (`, declaration parameter list, etc.).
    When a `#if`-family directive line is reached while the innermost
    open paren is a call context, the directive is inside that call's
    argument list. `#endif` lines are not penalised (counting the
    opener once is enough) and carry no parens, so they are skipped.

    Runs on comment-stripped (but NOT string-stripped) content so
    preprocessor lines survive; string literals virtually never contain
    a line-start `#`.

    Returns a list of 1-indexed line numbers of offending directives.
    """
    stripped = strip_comments(content)
    hits: List[int] = []
    # stack entries: True  -> paren opened by a function call
    #                False -> paren opened by something else
    paren_stack: List[bool] = []

    def _classify_opener(text_before: str) -> bool:
        """Decide whether `(` following `text_before` opens a call's
        argument list."""
        s = text_before.rstrip()
        if not s:
            return False
        m = _LAST_IDENT_RE.search(s)
        if not m:
            return False
        # control-statement header: `if (`, `for (`, `while (`,
        # `switch (`, `sizeof(` ... — the keyword itself is the last
        # identifier, with no callee between it and `(`
        if m.group(1) in STATEMENT_CONTEXT_KEYWORDS:
            return False
        # strip the callee name and inspect what precedes it
        before = s[:m.start()].rstrip()
        if not before:
            # bare callee in statement position: `f(` — a call (a
            # declaration would carry a return type before the name)
            return True
        # method call: `obj.m(`, `ptr->m(`
        if before.endswith(".") or before.endswith("->"):
            return True
        # assignment-init: `... = name(` — `name` is the CALLEE and the
        # whole expression is a call, regardless of any preceding
        # declaration-looking tokens. Without this, a cross-line tail
        # like `int f = add` would be misread as the declaration prefix
        # `type=int f =` + `name=add`.
        if before.endswith("="):
            return True
        mb = _LAST_IDENT_RE.search(before)
        if mb and mb.group(1) in STATEMENT_CONTEXT_KEYWORDS:
            # expression context: `return f(`, `else f(`, `throw f(`
            return True
        # `type name(` — declaration/definition parameter list
        return False

    line_no = 0
    # rolling window of recent non-preprocessor code so that a callee
    # name ending a previous line (`int f = add(`) is still visible when
    # a subsequent `(` or `#if` line is classified. Updated at the end
    # of EVERY non-preprocessor line, NOT just after the paren scan.
    prev_line_tail = ""
    for line in stripped.split("\n"):
        line_no += 1
        if PREPROC_IF_RE.match(line):
            if paren_stack and paren_stack[-1]:
                hits.append(line_no)
            continue
        if PREPROC_ENDIF_RE.match(line):
            continue
        # other preprocessor lines (`#else`, `#elif`, `#include` ...)
        # carry no argument-list parens; skip them so they never
        # disturb the stack. Do NOT touch prev_line_tail here: the code
        # line before the `#if`-block is what we still need to see.
        if line.lstrip().startswith("#"):
            continue
        for idx, c in enumerate(line):
            if c == "(":
                paren_stack.append(
                    _classify_opener(prev_line_tail + line[:idx]))
            elif c == ")":
                if paren_stack:
                    paren_stack.pop()
        prev_line_tail = (prev_line_tail + "\n" + line).rstrip()[-40:]
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
