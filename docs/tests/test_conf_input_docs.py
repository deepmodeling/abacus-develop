import stat
import sys
import tempfile
import unittest
from unittest import mock
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
    def test_readthedocs_environment_requires_fresh_input_docs(self):
        self.assertTrue(conf.input_docs_refresh_required({"READTHEDOCS": "True"}))
        self.assertFalse(conf.input_docs_refresh_required({"READTHEDOCS": "False"}))
        self.assertFalse(conf.input_docs_refresh_required({}))

    def test_explicit_abacus_binary_path_disables_fallback_search(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = Path(tmpdir)
            fallback_binary = repo_root / "build-rtd-docs" / "source" / "abacus_pw_ser"
            fallback_binary.parent.mkdir(parents=True)
            fallback_binary.write_text("#!/bin/sh\nexit 0\n")
            fallback_binary.chmod(fallback_binary.stat().st_mode | stat.S_IXUSR)

            with mock.patch.dict("os.environ", {"ABACUS_BINARY": str(repo_root / "missing-abacus")}):
                self.assertIsNone(conf.find_abacus_binary(repo_root))

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

    def test_required_refresh_raises_when_abacus_binary_is_not_available(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            docs_dir = Path(tmpdir) / "docs"
            docs_dir.mkdir()
            warnings = []

            with self.assertRaises(RuntimeError) as caught:
                conf.refresh_input_docs(
                    docs_dir=docs_dir,
                    abacus_binary=docs_dir / "missing-abacus",
                    warn=warnings.append,
                    fail_on_error=True,
                )

            self.assertIn("may not be up to date", str(caught.exception))
            self.assertEqual([], warnings)

    def test_required_refresh_raises_when_abacus_metadata_generation_fails(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            docs_dir = tmp_path / "docs"
            docs_dir.mkdir()
            fake_abacus = tmp_path / "abacus"
            fake_abacus.write_text(
                "#!/bin/sh\n"
                "echo 'metadata generation failed' >&2\n"
                "exit 2\n"
            )
            fake_abacus.chmod(fake_abacus.stat().st_mode | stat.S_IXUSR)
            warnings = []

            with self.assertRaises(RuntimeError) as caught:
                conf.refresh_input_docs(
                    docs_dir=docs_dir,
                    abacus_binary=fake_abacus,
                    warn=warnings.append,
                    fail_on_error=True,
                )

            self.assertIn("could not generate INPUT parameter metadata", str(caught.exception))
            self.assertIn("metadata generation failed", str(caught.exception))
            self.assertEqual([], warnings)


if __name__ == "__main__":
    unittest.main()
