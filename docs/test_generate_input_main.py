"""Checks for legal-value metadata in the generated INPUT reference."""

import unittest

from generate_input_main import generate_parameter_markdown


class InputDomainDocumentationTest(unittest.TestCase):
    def test_integer_choices_include_zero(self):
        rendered = generate_parameter_markdown(
            {"name": "gga_grad", "type": "Integer", "allowed_values": [0, 1, 2]}, {}
        )
        self.assertIn("- **Allowed values**: 0, 1, 2", rendered)

    def test_undeclared_domain_is_not_an_empty_restriction(self):
        rendered = generate_parameter_markdown({"name": "legacy", "type": "Real"}, {})
        self.assertNotIn("Allowed values", rendered)


if __name__ == "__main__":
    unittest.main()
