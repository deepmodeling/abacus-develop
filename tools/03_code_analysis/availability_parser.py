#!/usr/bin/env python3
"""Structured parser for ABACUS INPUT parameter ``availability`` values.

The ``availability`` field of an ``Input_Item`` states under which condition a
parameter is applicable. Historically this field is free text mixing prose and
ad-hoc conditions. This module provides a small, dependency-free parser that
turns such strings into a lightweight structured form so that documentation,
validation and tooling can consume the actual condition.

Kinds of parsed results (see :class:`Availability`):

* ``Expression`` -- a machine-readable condition such as ``basis_type==pw``
  or ``calculation==nscf``, optionally combined with ``and``/``or``/``,``.
* ``Label`` -- a bare applicability tag such as ``Numerical atomic orbital
  basis`` or ``TDOFDFT`` (no boolean operators).
* ``Unstructured`` -- free text that cannot be reliably structured without
  human review (kept verbatim).

The parser is deliberately tolerant: it normalises the several historical
spellings (``=`` vs ``==``, ``is set to``, ``is``, ``contains``) onto a single
canonical form so that rules can be compared and documented consistently.
"""

from __future__ import annotations

import re

# ---------------------------------------------------------------------------
# Result model
# ---------------------------------------------------------------------------


class Condition:
    """A single condition ``param op values``."""

    __slots__ = ("param", "op", "values")

    def __init__(self, param, op, values):
        self.param = param            # str, validated parameter name
        self.op = op                  # one of "==", "in"
        self.values = values          # list[str], canonical values

    def __repr__(self):
        vals = ", ".join(self.values)
        return f"Condition({self.param!r}, {self.op!r}, [{vals}])"


class Availability:
    """Parsed availability value."""

    __slots__ = ("kind", "expr", "label", "text")

    def __init__(self, kind, expr=None, label=None, text=None):
        self.kind = kind              # "Expression" | "Label" | "Unstructured"
        self.expr = expr              # Expression tree (see below) for Expression
        self.label = label            # str for Label
        self.text = text              # original text always kept

    def __repr__(self):
        return f"Availability({self.kind!r}, text={self.text!r})"


class Expr:
    """A boolean expression tree node: ('and'|'or', [children]) or leaf Condition."""

    __slots__ = ("op", "children")

    def __init__(self, op, children):
        self.op = op                  # 'and' | 'or' | None (leaf->single Condition)
        self.children = children      # list[Expr] when composite, or [Condition]

    def __repr__(self):
        return f"Expr({self.op!r}, {self.children!r})"


# ---------------------------------------------------------------------------
# Normalisation helpers
# ---------------------------------------------------------------------------

# Known boolean-ish / keyword markers that turn prose into operators.
_IS_SETTO = re.compile(r"is set to", re.IGNORECASE)
_IS_EQ = re.compile(r"\bis\b", re.IGNORECASE)
_CONTAINS = re.compile(r"\bcontains\b", re.IGNORECASE)
_IN = re.compile(r"\bin\b", re.IGNORECASE)

_PARAM_LIKE = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")


def _canonical_value(v):
    """Strip trailing punctuation / quotes and lowercase booleans."""
    v = v.strip().strip('"').strip("'").rstrip('.').strip()
    if v.lower() in ("true", "false"):
        return v.lower()
    return v


def _split_values_tokens(tokens):
    """Group candidate value tokens, dropping connecting words."""
    out = []
    seen = 0
    for t in tokens:
        tv = _canonical_value(t)
        if tv.lower() in ("or", ",", ",)"):
            continue
        out.append(tv)
        seen += 1
    return out


# ---------------------------------------------------------------------------
# Main parser
# ---------------------------------------------------------------------------


def _parse_single_condition(text, param_regex):
    """Try to parse ``text`` (a single atom) as ``param <op> <values>``.

    Returns a Condition or None. Handles the canonical grammar (==, >=, <=,
    !=, >, <, ``in [a, b]``) plus the historical spellings (``=``,
    ``is set to``, ``is``, ``contains``).
    """
    text = text.strip().strip('"').strip()
    if not text:
        return None

    # canonical comparison operators (longest/most specific first)
    for op in ("==", ">=", "<=", "!=", ">", "<"):
        pos = text.find(op)
        if pos > 0:
            param = text[:pos].strip()
            rhs = text[pos + len(op):].strip()
            if not (param_regex(param) and rhs):
                return None
            if op == "==":
                # equality may carry a slash/comma-separated value list
                values = _split_values_tokens([t for t in re.split(r"[/,]", rhs)])
                values = [v for v in values if v]
                if not values:
                    return None
                return Condition(param, "==", values)
            return Condition(param, op, [rhs])

    # legacy single '=' as equality
    pos = text.find("=")
    if pos > 0:
        param = text[:pos].strip()
        rhs = text[pos + 1:].strip()
        if param_regex(param) and rhs:
            values = _split_values_tokens([t for t in re.split(r"[/,]", rhs)])
            values = [v for v in values if v]
            if values:
                return Condition(param, "==", values)
        return None

    # canonical in-list: param in [v1, v2]
    m_in = re.search(r"\bin\s*\[", text, re.IGNORECASE)
    if m_in:
        param = text[:m_in.start()].strip()
        rest = text[m_in.end():].strip()
        if rest.endswith("]") and param_regex(param):
            values = [v.strip().strip('"').strip("'").rstrip(".") for v in rest[:-1].split(",")]
            values = [v for v in values if v]
            if values:
                return Condition(param, "in", values)
        return None

    # historical: is set to / contains / is
    for marker, op in ((_IS_SETTO, "=="), (_CONTAINS, "contains")):
        m = marker.search(text)
        if m:
            param = text[: m.start()].strip()
            rhs = text[m.end():].strip()
            if param_regex(param) and rhs:
                tokens = re.split(r"[,]", rhs)
                values = _split_values_tokens(tokens)
                values = [v for v in values if v]
                if values:
                    return Condition(param, op, values)
            return None

    m = _IS_EQ.search(text)
    if m:
        param = text[: m.start()].strip()
        rhs = text[m.end():].strip()
        if param_regex(param) and rhs:
            values = _split_values_tokens([t for t in re.split(r"[,/]", rhs)])
            values = [v for v in values if v]
            if values:
                op = "==" if len(values) == 1 else "in"
                return Condition(param, op, values)
        return None

    return None


def _split_top_level(text, keywords):
    """Split ``text`` on the given keywords at bracket/paren depth 0.

    Comma counts as a separator; word keywords (and/or) require word
    boundaries. Everything inside ``(...)`` or ``[...]`` is kept intact.
    """
    keywords = sorted(keywords, key=len, reverse=True)
    parts = []
    depth = 0
    start = 0
    n = len(text)
    i = 0
    while i < n:
        c = text[i]
        if c in "([":
            depth += 1
            i += 1
            continue
        if c in ")]":
            depth = max(0, depth - 1)
            i += 1
            continue
        if depth == 0:
            matched = False
            for kw in keywords:
                if kw == ",":
                    if c == ",":
                        seg = text[start:i].strip()
                        if seg:
                            parts.append(seg)
                        start = i + 1
                        i = start
                        matched = True
                        break
                else:
                    prev_bound = (i == 0) or text[i - 1].isspace()
                    next_bound = (i + len(kw) >= n) or text[i + len(kw)].isspace() \
                        or text[i + len(kw)] in "(["
                    if prev_bound and next_bound and text.startswith(kw, i):
                        seg = text[start:i].strip()
                        if seg:
                            parts.append(seg)
                        start = i + len(kw)
                        i = start
                        matched = True
                        break
            if matched:
                continue
        i += 1
    seg = text[start:].strip()
    if seg:
        parts.append(seg)
    return parts


def _build_expr(text, param_regex):
    """Recursively build an Expr/leaf from a boolean availability string.

    Returns a Condition/Expr or None if any part is not parseable.
    """
    text = text.strip()
    if not text:
        return None

    # Unwrap a parenthesized group that wraps the entire expression.
    if text.startswith("(") and text.endswith(")"):
        depth = 0
        wraps_all = True
        for k, ch in enumerate(text):
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
                if depth == 0 and k != len(text) - 1:
                    wraps_all = False
                    break
        if wraps_all:
            return _build_expr(text[1:-1], param_regex)

    # OR (top level)
    or_parts = _split_top_level(text, ["or"])
    if len(or_parts) > 1:
        children = [_build_expr(p, param_regex) for p in or_parts]
        if all(children):
            return Expr("or", children)
        return None

    # AND (top level, including comma)
    and_parts = _split_top_level(text, ["and", ","])
    if len(and_parts) > 1:
        children = [_build_expr(p, param_regex) for p in and_parts]
        if all(children):
            return Expr("and", children)
        return None

    # A single atom (which may itself be a parenthesized group).
    atom = text.strip()
    if atom.startswith("(") and atom.endswith(")"):
        return _build_expr(atom[1:-1], param_regex)
    return _parse_single_condition(atom, param_regex)


def parse_availability(text, param_regex=_PARAM_LIKE.match):
    """Parse an availability string into an :class:`Availability`.

    :param text: raw availability string (may be empty).
    :param param_regex: callable ``(str) -> bool`` used to decide whether a
        leading token is a plausible parameter name. Defaults to a loose
        identifier check.
    """
    text = (text or "").strip()
    if not text:
        return Availability("Label", text=text)

    # Canonical label form: "label: <name>"
    if len(text) >= 6 and text[:6].lower() == "label:":
        return Availability("Label", label=text[6:].strip(), text=text)

    # Cheap rejection of obvious prose / bare labels (no operator present).
    has_operator = re.search(
        r"==|=| in | in\[| in$| is set to|\bis\b|\bcontains\b|>=|<=|!=|\b>\b|\b<\b",
        text, re.IGNORECASE,
    )
    if not has_operator:
        return Availability("Label", label=text.strip('"').rstrip('.'), text=text)

    expr = _build_expr(text, param_regex)
    if expr is None:
        return Availability("Unstructured", text=text)
    # Always expose a top-level Expr node so that consumers can iterate
    # ``expr.children`` (matching the historical flat-AND shape for simple
    # conditions).
    if isinstance(expr, Condition):
        expr = Expr("and", [expr])
    return Availability("Expression", expr=expr, text=text)
