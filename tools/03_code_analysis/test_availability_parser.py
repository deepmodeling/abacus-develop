#!/usr/bin/env python3
"""Unit tests for :mod:`availability_parser`."""

import os
import sys
import unittest

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "tools" / "03_code_analysis"))

from availability_parser import (
    Availability,
    Condition,
    Expr,
    parse_availability,
)


def _conds(availability):
    """Return the flat list of Condition leaves of an Expression result."""
    assert availability.kind == "Expression"
    return [c for c in availability.expr.children if isinstance(c, Condition)]


class AvailabilityParserTest(unittest.TestCase):
    def test_empty_is_label(self):
        r = parse_availability("")
        self.assertEqual(r.kind, "Label")

    def test_bare_label(self):
        r = parse_availability("Numerical atomic orbital basis")
        self.assertEqual(r.kind, "Label")
        self.assertEqual(r.label, "Numerical atomic orbital basis")

    def test_label_with_period(self):
        r = parse_availability("TDOFDFT.")
        self.assertEqual(r.kind, "Label")

    def test_double_equal(self):
        r = parse_availability("basis_type==lcao")
        self.assertEqual(r.kind, "Expression")
        c = _conds(r)[0]
        self.assertEqual((c.param, c.op, c.values), ("basis_type", "==", ["lcao"]))

    def test_single_equal_and_trailing_period(self):
        r = parse_availability("esolver_type = dp.")
        self.assertEqual(r.kind, "Expression")
        c = _conds(r)[0]
        self.assertEqual((c.param, c.op, c.values), ("esolver_type", "==", ["dp"]))

    def test_slash_value_list(self):
        r = parse_availability("basis_type==pw, ks_solver==cg/dav/dav_subspace/bpcg")
        self.assertEqual(r.kind, "Expression")
        conds = _conds(r)
        vals = {c.param: c.values for c in conds}
        self.assertEqual(vals["basis_type"], ["pw"])
        self.assertEqual(vals["ks_solver"], ["cg", "dav", "dav_subspace", "bpcg"])

    def test_and_combination(self):
        r = parse_availability("basis_type==lcao and esolver_type==tddft")
        self.assertEqual(r.kind, "Expression")
        self.assertEqual(len(_conds(r)), 2)

    def test_is_true_boolean(self):
        r = parse_availability("imp_sol is true")
        self.assertEqual(r.kind, "Expression")
        c = _conds(r)[0]
        self.assertEqual((c.param, c.op, c.values), ("imp_sol", "==", ["true"]))

    def test_is_set_to(self):
        r = parse_availability("vdw_method is set to d2")
        self.assertEqual(r.kind, "Expression")
        c = _conds(r)[0]
        self.assertEqual((c.param, c.op, c.values), ("vdw_method", "==", ["d2"]))

    def test_unstructured_prose(self):
        r = parse_availability("Only used when relax_method is bfgs or cg_bfgs")
        self.assertEqual(r.kind, "Unstructured")
        self.assertEqual(r.text, "Only used when relax_method is bfgs or cg_bfgs")

    def test_text_preserved(self):
        r = parse_availability("basis_type==lcao")
        self.assertEqual(r.text, "basis_type==lcao")


if __name__ == "__main__":
    unittest.main()
