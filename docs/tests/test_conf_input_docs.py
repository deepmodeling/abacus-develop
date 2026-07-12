import stat
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import conf


PARAMETER_YAML = """\
parameters:
  - key: ecutwfc
    name: ecutwfc
    category: Plane wave related variables
    type: double
    default_value: 100
    unit: Ry
    description: Plane wave kinetic energy cutoff.
"""


class InputDocsGenerationTest(unittest.TestCase):
    def test_refreshes_input_docs_from_abacus_binary_without_writing_yaml_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            docs_dir = tmp_path / "docs"
            docs_dir.mkdir()
            fake_abacus = tmp_path / "abacus"
            fake_abacus.write_text(
                "#!/bin/sh\n"
                "if [ \"$1\" = \"--generate-parameters-yaml\" ]; then\n"
                f"  cat <<'EOF'\n{PARAMETER_YAML}EOF\n"
                "  exit 0\n"
                "fi\n"
                "exit 1\n"
            )
            fake_abacus.chmod(fake_abacus.stat().st_mode | stat.S_IXUSR)

            warnings = []
            refreshed = conf.refresh_input_docs(
                docs_dir=docs_dir,
                abacus_binary=fake_abacus,
                warn=warnings.append,
            )

            output = docs_dir / "advanced" / "input_files" / "input-main.md"
            self.assertTrue(refreshed)
            self.assertTrue(output.exists())
            self.assertIn("ecutwfc", output.read_text())
            self.assertFalse((docs_dir / "parameters.yaml").exists())
            self.assertEqual([], warnings)

    def test_warns_when_abacus_binary_is_not_available(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            docs_dir = Path(tmpdir) / "docs"
            docs_dir.mkdir()
            warnings = []

            refreshed = conf.refresh_input_docs(
                docs_dir=docs_dir,
                abacus_binary=docs_dir / "missing-abacus",
                warn=warnings.append,
            )

            self.assertFalse(refreshed)
            self.assertEqual([], list(docs_dir.rglob("input-main.md")))
            self.assertEqual(1, len(warnings))
            self.assertIn("may not be up to date", warnings[0])


if __name__ == "__main__":
    unittest.main()
