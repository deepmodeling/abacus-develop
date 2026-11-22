import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from abacusw90.interface import ABACUSWannier90


class TestSpinOutput(unittest.TestCase):
    def test_spn_default(self):
        job = ABACUSWannier90()
        for overrides, expected in (
            ({}, True), ({"out_wannier_spn": "0"}, False),
            ({"towannier90": "no"}, False), ({"calculation": "scf"}, False),
            ({"basis_type": "pw"}, False), ({"wannier_method": 1}, False), ({"nspin": 2}, False),
        ):
            with self.subTest(overrides=overrides):
                params = {"basis_type": "lcao", **overrides}
                job.set_abacus_parameters(ecutwfc=20, nbands=14, **params)
                self.assertEqual(job._spn_enabled(), expected)

    def test_disabled_companion_rejected_before_nscf(self):
        with tempfile.TemporaryDirectory() as scratch:
            work = Path(scratch)
            out = work / "OUT.ABACUS"
            out.mkdir()
            for ext in ("mmn", "amn", "eig", "spn"):
                (out / f"wannier90.{ext}").write_text("previous NSCF run")

            job = ABACUSWannier90(work_dir=work)
            with patch("abacusw90.interface.runner.run_command") as launch:
                for ext, disabled in (("mmn", 0), ("amn", " .FALSE. "), ("eig", False)):
                    flag = f"out_wannier_{ext}"
                    with self.subTest(flag=flag):
                        job.set_abacus_parameters(
                            ecutwfc=20, nbands=14, out_wannier_spn="YES", **{flag: disabled}
                        )
                        with self.assertRaisesRegex(ValueError, flag):
                            job.step3_run_abacus()
                        launch.assert_not_called()
                        self.assertEqual(list(work.glob("wannier90.*")), [])
