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

    def test_in_list(self):
        r = parse_availability("vdw_method in [d2, d3_0, d3_bj]")
        self.assertEqual(r.kind, "Expression")
        c = _conds(r)[0]
        self.assertEqual((c.param, c.op, c.values),
                         ("vdw_method", "in", ["d2", "d3_0", "d3_bj"]))

    def test_comparison_operator(self):
        r = parse_availability("vdw_C6_file!=default")
        self.assertEqual(r.kind, "Expression")
        c = _conds(r)[0]
        self.assertEqual((c.param, c.op, c.values), ("vdw_C6_file", "!=", ["default"]))

    def test_label_prefix(self):
        r = parse_availability("label: Numerical atomic orbital basis")
        self.assertEqual(r.kind, "Label")
        self.assertEqual(r.label, "Numerical atomic orbital basis")

    def test_parens_or_group(self):
        r = parse_availability(
            "symmetry==1 and (dft_functional in [hse, hf, pbe0, scan0] or rpa==true)")
        self.assertEqual(r.kind, "Expression")
        root = r.expr
        # "and" binds looser than "or"; the grouped "(A or B)" must stay a
        # separate child so precedence is preserved.
        self.assertIsInstance(root, Expr)
        self.assertEqual(root.op, "and")
        self.assertEqual(len(root.children), 2)
        c0, c1 = root.children
        self.assertIsInstance(c0, Condition)
        self.assertEqual((c0.param, c0.op, c0.values), ("symmetry", "==", ["1"]))
        self.assertIsInstance(c1, Expr)
        self.assertEqual(c1.op, "or")
        self.assertEqual(len(c1.children), 2)
        d, rpa = c1.children
        self.assertEqual((d.param, d.op, d.values),
                         ("dft_functional", "in", ["hse", "hf", "pbe0", "scan0"]))
        self.assertEqual((rpa.param, rpa.op, rpa.values), ("rpa", "==", ["true"]))

    def test_contains_vector_semantics(self):
        # td_ttype is a Vector: "contains 2" is containment, distinct from
        # scalar membership ("in [2]"), so the operator must be preserved.
        r = parse_availability("td_ttype contains 2")
        self.assertEqual(r.kind, "Expression")
        c = _conds(r)[0]
        self.assertEqual((c.param, c.op, c.values), ("td_ttype", "contains", ["2"]))

    def test_parameters_yaml_availability_are_all_expression(self):
        """Every non-empty availability in the canonical parameters.yaml must be
        a concrete boolean Expression. Any Label/Unstructured (prose) non-empty
        value fails the PR build: we do not silently accept non-expressions.
        """
        try:
            import yaml
        except ImportError:
            self.skipTest("PyYAML not available")
        yaml_path = REPO_ROOT / "docs" / "parameters.yaml"
        if not yaml_path.exists():
            self.skipTest("docs/parameters.yaml not present")
        data = yaml.safe_load(yaml_path.read_text())
        params = data.get("parameters", [])
        bad = []
        for p in params:
            avail = (p.get("availability") or "").strip()
            if not avail:
                continue
            r = parse_availability(avail)
            if r.kind != "Expression":
                bad.append((p.get("name"), avail, r.kind))
        self.assertEqual(bad, [],
                         "non-Expression availability present (not all "
                         "availability are machine-evaluable conditions): "
                         + repr(bad))


if __name__ == "__main__":
    unittest.main()
