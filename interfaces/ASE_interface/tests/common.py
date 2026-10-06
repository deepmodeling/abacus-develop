"""Shared helpers for the ASE interface integration tests."""

import contextlib
import os
import shutil
import tempfile
from pathlib import Path


@contextlib.contextmanager
def preserved_tmpdir(prefix):
    """Temporary ABACUS run directory kept on failure for post-mortem CI.

    On success the directory is removed like tempfile.TemporaryDirectory.
    On exception it is copied into ASE_ABACUS_KEEP_DIR (default
    /tmp/abacus_ase_failure) so remote-CI-only failures (e.g. issue #7794)
    can be root-caused from the uploaded artifact.
    """
    tmpdir = tempfile.mkdtemp(prefix=prefix)
    try:
        yield tmpdir
    except Exception:
        keep_root = Path(os.environ.get(
            'ASE_ABACUS_KEEP_DIR', '/tmp/abacus_ase_failure'))
        keep_root.mkdir(parents=True, exist_ok=True)
        keep = keep_root / Path(tmpdir).name
        shutil.copytree(tmpdir, keep, dirs_exist_ok=True)
        print(f'[ase-test] failure, run dir preserved at {keep}')
        raise
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
