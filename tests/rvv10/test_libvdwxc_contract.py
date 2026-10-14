"""Exercise the optional-backend build report without external dependencies."""

import os
from pathlib import Path
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]


class LibvdwxcContract(unittest.TestCase):
    def test_makefile_reports_only_compiled_libvdwxc_support(self):
        environment = os.environ.copy()
        environment["ENABLE_LIBVDWXC"] = "ON"
        with tempfile.TemporaryDirectory() as directory:
            header = Path(directory) / "build_info.h"
            subprocess.run(
                ["bash", str(ROOT / "generate_build_info.sh"), str(header)],
                env=environment, check=True, capture_output=True, text=True,
            )
            self.assertIn('#define ABACUS_LIBVDWXC_VERSION "no"', header.read_text())


if __name__ == "__main__":
    unittest.main()
