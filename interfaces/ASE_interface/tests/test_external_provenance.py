import tempfile
import unittest
from pathlib import Path

from verify_external import verify_locations


class ExternalProvenanceTest(unittest.TestCase):
    def test_installed_wheel_location(self):
        with tempfile.TemporaryDirectory() as root:
            root = Path(root)
            site = root / 'venv' / 'site-packages' / 'abacuslite' / '__init__.py'
            source = root / 'checkout' / 'interfaces' / 'ASE_interface'
            site.parent.mkdir(parents=True)
            source.mkdir(parents=True)
            site.touch()
            verify_locations(site, site, source)

    def test_source_shadow_is_rejected(self):
        with tempfile.TemporaryDirectory() as root:
            root = Path(root)
            source = root / 'checkout' / 'interfaces' / 'ASE_interface'
            shadow = source / 'abacuslite' / '__init__.py'
            installed = root / 'venv' / 'site-packages' / 'abacuslite' / '__init__.py'
            shadow.parent.mkdir(parents=True)
            installed.parent.mkdir(parents=True)
            shadow.touch()
            installed.touch()
            with self.assertRaisesRegex(AssertionError, 'source shadow'):
                verify_locations(shadow, installed, root / 'checkout')

    def test_wrong_distribution_is_rejected(self):
        with tempfile.TemporaryDirectory() as root:
            root = Path(root)
            imported = root / 'one' / 'abacuslite' / '__init__.py'
            installed = root / 'two' / 'abacuslite' / '__init__.py'
            imported.parent.mkdir(parents=True)
            installed.parent.mkdir(parents=True)
            imported.touch()
            installed.touch()
            with self.assertRaisesRegex(AssertionError, 'installed distribution'):
                verify_locations(imported, installed, root / 'checkout')
