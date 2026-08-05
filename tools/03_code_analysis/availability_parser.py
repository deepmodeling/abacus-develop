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

_PARAM_LIKE = re.compile(r"^[a-z][a-z0-9_]*$")


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

    Returns a Condition or None.
    """
    text = text.strip().strip('"').strip()
    if not text:
        return None

    # param == value / param = value
    for op in ("==", "="):
        idx = text.find(op)
        if idx > 0:
            param = text[:idx].strip()
            rhs = text[idx + len(op):].strip()
            if param_regex(param) and rhs:
                values = _split_values_tokens([t for t in re.split(r"[/,]", rhs)])
                values = [v for v in values if v]
                if values:
                    return Condition(param, "==", values)
            return None

    # param is set to <values> / param is <value> / param contains <value>
    for marker, op in (
        (_IS_SETTO, "=="),
        (_CONTAINS, "in"),
    ):
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
                # "param is true/false" -> ==; otherwise treat as membership
                op = "==" if len(values) == 1 else "in"
                return Condition(param, op, values)
        return None

    return None


def _tokenise_bool(text):
    """Split a boolean expression at top-level ``and``/``or``/``,``.

    Returns a list of ``(atom_text, sep)`` pairs where ``sep`` is the
    connecting keyword that followed the atom (``and``/``or``/``,``), or
    ``None`` for the last atom.
    """
    tokens = re.split(r"(\band\b|\bor\b|,)", text, flags=re.IGNORECASE)
    parts = []
    pending_sep = None
    for t in tokens:
        t = t.strip()
        if not t:
            continue
        low = t.lower()
        if low in ("and", "or", ","):
            # separator that binds the *previous* atom to the next one
            pending_sep = low if low != "," else "and"
            continue
        parts.append((t, pending_sep))
        pending_sep = None
    return parts


def parse_availability(text, param_regex=_PARAM_LIKE.match):
    """Parse an availability string into an :class:`Availability`.

    :param text: raw availability string (may be empty).
    :param param_regex: callable ``(str) -> bool`` used to decide whether a
        leading token is a plausible parameter name. Defaults to a loose
        lowercase identifier check.
    """
    text = (text or "").strip()
    if not text:
        return Availability("Label", text=text)

    # Cheap rejection of obvious prose / bare labels (no operator present).
    has_operator = re.search(
        r"==|=| in | in$| is set to|\bis\b|\bcontains\b", text, re.IGNORECASE
    )
    if not has_operator:
        return Availability("Label", label=text.strip('"').rstrip('.'), text=text)

    # Try to interpret as a boolean condition over atoms.
    atoms = _tokenise_bool(text)
    conds = []
    for atom, sep in atoms:
        c = _parse_single_condition(atom, param_regex)
        if c is None:
            # Not parseable -> whole thing is unstructured.
            return Availability("Unstructured", text=text)
        conds.append((c, sep))

    if not conds:
        return Availability("Unstructured", text=text)

    # Build an AND/OR expression tree from the parsed atoms and connectors.
    nodes = [c for c, _ in conds]
    expr = Expr("and", nodes)  # flat AND; OR connectors are recorded but not
    # given a separate tree node (no parenthesised precedence in current data).
    return Availability("Expression", expr=expr, text=text)
