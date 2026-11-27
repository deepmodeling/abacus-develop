"""Text-level helpers for the ABACUS code-quality scorer.

Comment/string stripping that preserves line structure, comment-block
extraction for documentation rules, and include-guard helpers.
"""

import re
from typing import List, Tuple


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


def extract_comment_blocks(content: str) -> List[List[str]]:
    """Extract normalised comment lines, grouped by comment block.

    Returns a list of blocks; each block is a list of normalised comment
    lines (comment markers stripped, whitespace collapsed, empty lines
    dropped). A 'block' is a maximal run of `///` lines, or one `/* ... */`
    region. String literals are respected so that `//` inside a string is
    not treated as a comment opener.
    """
    blocks: List[List[str]] = []
    current: List[str] = []
    i = 0
    n = len(content)
    in_string = False
    in_char = False

    def flush() -> None:
        nonlocal current
        if current:
            blocks.append(current)
            current = []

    while i < n:
        c = content[i]
        if in_string:
            if c == "\\" and i + 1 < n:
                i += 2
                continue
            if c == '"':
                in_string = False
            i += 1
            continue
        if in_char:
            if c == "\\" and i + 1 < n:
                i += 2
                continue
            if c == "'":
                in_char = False
            i += 1
            continue
        if c == '"':
            in_string = True
            i += 1
            continue
        if c == "'":
            if _is_digit_separator(content, i):
                i += 1
                continue
            in_char = True
            i += 1
            continue
        if c == "/" and i + 1 < n and content[i + 1] == "/":
            # line comment: collect the text after // (or /// for doxygen)
            j = content.find("\n", i)
            if j < 0:
                j = n
            text = content[i + 2:j].lstrip("/ \t").strip()
            if text:
                current.append(" ".join(text.split()))
            i = j
            continue
        if c == "/" and i + 1 < n and content[i + 1] == "*":
            # block comment: collect until */
            j = content.find("*/", i + 2)
            if j < 0:
                j = n - 2
            body = content[i + 2:j]
            block: List[str] = []
            for ln in body.split("\n"):
                # strip leading comment markers: '*', '/', whitespace, and
                # any leftover '///' from doxygen-style lines
                ln = re.sub(r"^[\s*/]+", "", ln).strip()
                if ln:
                    block.append(" ".join(ln.split()))
            if block:
                flush()
                blocks.append(block)
            i = j + 2
            continue
        if c == "\n":
            # a blank/non-comment line ends a /// run; check next non-empty
            # cheaply by flushing on newline only when current is non-empty
            # and the following line does not start with optional ws + //
            rest = content[i + 1:i + 200]
            m = re.match(r"\s*//", rest)
            if not m:
                flush()
            i += 1
            continue
        if not c.isspace():
            flush()
        i += 1
    flush()
    return blocks
