# Tests extracted without semantic changes from the attributed upstream snapshot.
import os
import unittest
from unittest.mock import patch
from pathlib import Path

from abacuslite.core import AbacusProfile, SOCKET_METADATA_SCHEMA


class TestAbacusProfileAssetDefaults(unittest.TestCase):

    def test_explicit_directories_take_precedence(self):
        with patch.dict(os.environ, {
            "ABACUS_PP_PATH": "/env/pp",
            "ABACUS_ORB_PATH": "/env/orb",
        }):
            profile = AbacusProfile(
                "abacus",
                pseudo_dir="/explicit/pp",
                orbital_dir="/explicit/orb",
            )
        self.assertEqual(profile.pseudo_dir, "/explicit/pp")
        self.assertEqual(profile.orbital_dir, "/explicit/orb")

    def test_community_environment_defaults_are_used(self):
        with patch.dict(os.environ, {
            "ABACUS_PP_PATH": "/env/pp",
            "ABACUS_ORB_PATH": "/env/orb",
        }):
            profile = AbacusProfile("abacus")
        self.assertEqual(profile.pseudo_dir, "/env/pp")
        self.assertEqual(profile.orbital_dir, "/env/orb")

    def test_explicit_none_argument_uses_environment(self):
        with patch.dict(os.environ, {
            "ABACUS_PP_PATH": "/env/pp",
            "ABACUS_ORB_PATH": "/env/orb",
        }):
            profile = AbacusProfile("abacus", pseudo_dir=None, orbital_dir=None)
        self.assertEqual(profile.pseudo_dir, "/env/pp")
        self.assertEqual(profile.orbital_dir, "/env/orb")

    def test_explicit_empty_directory_is_preserved(self):
        with patch.dict(os.environ, {
            "ABACUS_PP_PATH": "/env/pp",
            "ABACUS_ORB_PATH": "/env/orb",
        }):
            profile = AbacusProfile("abacus", pseudo_dir="", orbital_dir="")
        self.assertEqual(profile.pseudo_dir, "")
        self.assertEqual(profile.orbital_dir, "")

    def test_unset_or_empty_environment_leaves_directories_unset(self):
        with patch.dict(os.environ, {
            "ABACUS_PP_PATH": "",
            "ABACUS_ORB_PATH": "",
        }, clear=True):
            profile = AbacusProfile("abacus")
        self.assertIsNone(profile.pseudo_dir)
        self.assertIsNone(profile.orbital_dir)

    def test_other_asset_variables_are_ignored(self):
        with patch.dict(os.environ, {
            "ABACUS_ORBITAL_PATH": "/legacy/orb",
            "ABACUS_PP_ORB": "/both",
        }, clear=True):
            profile = AbacusProfile("abacus")
        self.assertIsNone(profile.pseudo_dir)
        self.assertIsNone(profile.orbital_dir)


class TestAbacusCalculator(unittest.TestCase):

    here = Path(__file__).resolve().parents[2] / 'abacuslite/io'
    pporb = Path(os.environ.get('ABACUS_PP_ORB', 'missing-PP_ORB'))

    def test_socket_metadata_schema_is_v1(self):
        self.assertEqual(SOCKET_METADATA_SCHEMA, 'abacus.socket.properties.v1')
