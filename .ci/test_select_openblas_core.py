"""Policy and environment propagation tests; real-library checks run separately."""

from contextlib import redirect_stdout
import io
import math
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import Mock, patch

import select_openblas_core as selector


AVX2 = {"avx", "avx2", "fma"}
AVX512 = AVX2 | {"avx512f", "avx512cd", "avx512dq", "avx512bw", "avx512vl"}


def result(core, error):
    return {"core": core, "config": "OpenBLAS test fixture", "max_abs_error": error}


class SelectionTests(unittest.TestCase):
    def choose(self, machine, flags, probe):
        with redirect_stdout(io.StringIO()):
            return selector.select_core("/library.so", machine, flags, probe)

    def test_healthy_autodetection_keeps_every_architecture_unchanged(self):
        for machine, core in [("x86_64", "Cooperlake"), ("aarch64", "ARMV8")]:
            with self.subTest(machine=machine):
                probe = Mock(return_value=result(core, 4e-16))
                self.assertIsNone(self.choose(machine, set(), probe))
                probe.assert_called_once_with("/library.so", None)

    def test_broken_cooperlake_uses_validated_avx512_fallback(self):
        probe = Mock(side_effect=[result("Cooperlake", 2.786), result("SkylakeX", 2e-16)])
        self.assertEqual(self.choose("x86_64", AVX512, probe), "SkylakeX")
        self.assertEqual([call.args[1] for call in probe.call_args_list], [None, "SkylakeX"])

    def test_avx2_host_never_attempts_avx512(self):
        probe = Mock(side_effect=[result("Cooperlake", 2.786), result("Haswell", 4e-16)])
        self.assertEqual(self.choose("x86_64", AVX2, probe), "Haswell")
        self.assertEqual([call.args[1] for call in probe.call_args_list], [None, "Haswell"])

    def test_broken_avx512_fallback_can_continue_to_haswell(self):
        probe = Mock(side_effect=[result("Cooperlake", 2.786), result("SkylakeX", 0.1),
                                  result("Haswell", 4e-16)])
        self.assertEqual(self.choose("x86_64", AVX512, probe), "Haswell")

    def test_missing_requested_kernel_is_not_accepted_as_that_kernel(self):
        probe = Mock(side_effect=[result("Cooperlake", 2.786), result("Barcelona", 4e-16),
                                  result("Haswell", 4e-16)])
        self.assertEqual(self.choose("x86_64", AVX512, probe), "Haswell")

    def test_no_x86_override_on_arm_or_without_required_features(self):
        for machine, flags in [("aarch64", AVX512), ("x86_64", {"avx", "avx2"}),
                               ("x86_64", set())]:
            with self.subTest(machine=machine, flags=flags):
                probe = Mock(return_value=result("broken", 2.786))
                with self.assertRaisesRegex(RuntimeError, "No CPU-compatible"):
                    self.choose(machine, flags, probe)
                probe.assert_called_once_with("/library.so", None)

    def test_signal_and_failed_fallbacks_cannot_approve_a_kernel(self):
        probe = Mock(side_effect=[{"error": "probe exit -4"}, result("SkylakeX", math.nan),
                                  result("Haswell", math.inf)])
        with self.assertRaisesRegex(RuntimeError, "No CPU-compatible"):
            self.choose("x86_64", AVX512, probe)

    def test_cpu_flags_require_support_on_every_cpu(self):
        cpuinfo = "processor : 0\nflags : avx avx2 fma avx512f\n\nprocessor : 1\nflags : avx avx2 fma\n"
        flags = selector.common_cpu_flags(cpuinfo)
        self.assertEqual(flags, AVX2)
        self.assertEqual(selector.fallback_cores("x86_64", flags), ["Haswell"])

    def test_nonfinite_missing_and_large_errors_fail(self):
        for error in [math.nan, math.inf, None, -1.0, 2.786]:
            with self.subTest(error=error):
                self.assertFalse(selector.probe_passed(result("test", error)))


class LibraryTests(unittest.TestCase):
    def test_blas_alias_is_resolved_and_deduplicated(self):
        with tempfile.TemporaryDirectory() as directory:
            library = Path(directory) / "libopenblasp-r0.3.20.so"
            library.touch()
            alias = Path(directory) / "libblas.so.3"
            alias.symlink_to(library)
            listing = f"libblas.so.3 => {alias} (0x1)\nlibopenblas.so.0 => {library} (0x2)"
            self.assertEqual(selector.linked_openblas(listing), str(library.resolve()))

    def test_other_blas_is_not_subject_to_an_openblas_override(self):
        self.assertIsNone(selector.linked_openblas("libmkl_rt.so => /opt/mkl/libmkl_rt.so (0x1)"))

    def test_missing_library_is_an_error(self):
        with self.assertRaisesRegex(RuntimeError, "Unresolved"):
            selector.linked_openblas("libopenblas.so.0 => not found")

    def test_multiple_libraries_are_not_silently_reduced_to_one(self):
        with self.assertRaisesRegex(RuntimeError, "Multiple"):
            selector.linked_openblas("libopenblas.so.0 => /a/libopenblas.so.0 (0x1)\n"
                                     "libopenblas_other.so => /b/libopenblas_other.so (0x2)")


class EnvironmentTests(unittest.TestCase):
    def run_main(self, output, probe):
        ldd = subprocess.CompletedProcess([], 0, "libopenblas.so.0 => /lib/libopenblas.so.0 (0x1)")
        with patch.object(selector.shutil, "which", return_value="/abacus"), \
                patch.object(selector.subprocess, "run", return_value=ldd), \
                patch.object(selector.platform, "machine", return_value="x86_64"), \
                patch.object(selector, "common_cpu_flags", return_value=AVX512), \
                patch.object(selector, "run_probe", probe), redirect_stdout(io.StringIO()):
            return selector.main(["--github-env", str(output)])

    def test_only_a_validated_fallback_is_appended_to_github_env(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "github-env"
            output.write_text("EXISTING=value\n")
            probe = Mock(side_effect=[result("Cooperlake", 2.786), result("SkylakeX", 2e-16)])
            self.assertEqual(self.run_main(output, probe), 0)
            self.assertEqual(output.read_text(), "EXISTING=value\nOPENBLAS_CORETYPE=SkylakeX\n")

    def test_healthy_or_failed_checks_do_not_write_an_override(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "github-env"
            output.write_text("EXISTING=value\n")
            self.assertEqual(self.run_main(output, Mock(return_value=result("Cooperlake", 2e-16))), 0)
            self.assertEqual(output.read_text(), "EXISTING=value\n")
            with self.assertRaises(RuntimeError):
                self.run_main(output, Mock(return_value={"error": "library cannot initialize"}))
            self.assertEqual(output.read_text(), "EXISTING=value\n")


if __name__ == "__main__":
    unittest.main()
