#!/usr/bin/env python3
"""Check ABACUS's OpenBLAS DSYMM and select a CPU-compatible CI fallback.

Keep the current kernel when the regression probe passes. On failure, probe
compatible x86 alternatives in fresh processes before writing GITHUB_ENV.
This checks the known Si failure, not every possible BLAS operation/shape.
Only the Python standard library is needed; no NumPy BLAS is loaded.
"""

import argparse
import ctypes as ct
import json
import math
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys


def linked_openblas(ldd_output):
    """Resolve aliases so the probe loads the library linked to ABACUS."""
    libraries = set()
    for name, target in re.findall(r"^\s*(\S+)\s+=>\s+(\S+)", ldd_output, re.MULTILINE):
        if not name.startswith(("libopenblas", "libblas.")):
            continue
        if not target.startswith("/"):
            raise RuntimeError(f"Unresolved BLAS dependency: {name}")
        path = Path(target).resolve()
        if name.startswith("libopenblas") or path.name.startswith("libopenblas"):
            libraries.add(str(path))
    if len(libraries) > 1:
        raise RuntimeError(f"Multiple OpenBLAS libraries linked: {sorted(libraries)}")
    return next(iter(libraries), None)


def common_cpu_flags(cpuinfo):
    """Use features present on every listed CPU, including heterogeneous hosts."""
    flags = [set(line.split(":", 1)[1].split())
             for line in cpuinfo.splitlines() if line.split(":", 1)[0].strip() == "flags"]
    return set.intersection(*flags) if flags else set()


def fallback_cores(machine, flags):
    if machine.lower() not in ("x86_64", "amd64"):
        return []
    cores = []
    avx2 = {"avx", "avx2", "fma"}
    avx512 = avx2 | {"avx512f", "avx512cd", "avx512dq", "avx512bw", "avx512vl"}
    if avx512 <= flags:
        cores.append("SkylakeX")
    if avx2 <= flags:
        cores.append("Haswell")
    return cores


def probe_library(library):
    blas = ct.CDLL(library)
    blas.openblas_get_corename.restype = ct.c_char_p
    blas.openblas_get_config.restype = ct.c_char_p
    blas.dsymm_.restype = None
    config = blas.openblas_get_config().decode()
    if "USE64BITINT" in config:
        raise RuntimeError("This regression probe requires the LP64 BLAS ABI used by ABACUS")

    # A = 2I makes the expected result exactly 2B. This 13-orbital Si block
    # with strided columns exposes the packaged Cooperlake error (16.0).
    m = 13
    n = 3
    stride = 2 * m
    a = (ct.c_double * (m * m))()
    for i in range(m):
        a[i + m * i] = 2.0
    b = (ct.c_double * (stride * n))(
        *((i + 3 * j) % 13 - 6 for j in range(n) for i in range(stride)))
    c = (ct.c_double * (stride * n))()
    cm = ct.c_int(m)
    cn = ct.c_int(n)
    ld = ct.c_int(stride)
    alpha = ct.c_double(1.0)
    beta = ct.c_double(1.0)
    blas.dsymm_(b"L", b"U", ct.byref(cm), ct.byref(cn), ct.byref(alpha),
                a, ct.byref(cm), b, ct.byref(ld), ct.byref(beta), c, ct.byref(ld))
    errors = (abs(c[i + stride * j] - 2 * b[i + stride * j])
              for j in range(n) for i in range(m))
    worst = max(error if math.isfinite(error) else math.inf for error in errors)
    return {"core": blas.openblas_get_corename().decode(), "config": config,
            "max_abs_error": worst}


def run_probe(library, core):
    # OpenBLAS chooses its kernel when loaded; changing the environment in an
    # already initialized process would not validate the requested alternative.
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    if core is not None:
        env["OPENBLAS_CORETYPE"] = core
    try:
        result = subprocess.run([sys.executable, str(Path(__file__).resolve()), "--probe", library],
                                env=env, capture_output=True, text=True, timeout=30)
        if result.returncode != 0:
            return {"error": f"probe exit {result.returncode}: {result.stderr.strip()}"}
        record = json.loads(result.stdout)
        if not isinstance(record, dict) or not isinstance(record.get("core"), str):
            raise ValueError("missing OpenBLAS kernel identity")
        return record
    except (OSError, ValueError, subprocess.TimeoutExpired) as error:
        return {"error": str(error)}


def probe_passed(record):
    error = record.get("max_abs_error")
    return (not record.get("error") and isinstance(error, (float, int))
            and math.isfinite(error) and 0.0 <= error <= 1e-12)


def select_core(library, machine, flags, probe):
    for core in [None] + fallback_cores(machine, flags):
        result = probe(library, core)
        print(f"OpenBLAS {core or 'current'}: " + json.dumps(result), flush=True)
        if probe_passed(result) and (core is None or result["core"].lower() == core.lower()):
            return core
    raise RuntimeError("No CPU-compatible OpenBLAS kernel passed the DSYMM regression probe. "
                       "Use a corrected BLAS library; the ASE tests have not been started.")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--abacus", default="abacus")
    parser.add_argument("--github-env", type=Path, default=os.environ.get("GITHUB_ENV"))
    parser.add_argument("--probe", help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    if args.probe:
        print(json.dumps(probe_library(args.probe)))
        return 0
    if args.github_env is None:
        parser.error("--github-env or GITHUB_ENV is required")
    executable = shutil.which(args.abacus)
    if executable is None:
        raise RuntimeError(f"ABACUS executable not found: {args.abacus}")
    dependencies = subprocess.run(["ldd", executable], check=True, capture_output=True, text=True)
    library = linked_openblas(dependencies.stdout)
    if library is None:
        print("ABACUS does not link OpenBLAS; no kernel override is needed.")
        return 0
    print(f"Checking ABACUS OpenBLAS library: {library}")
    cpuinfo = Path("/proc/cpuinfo")
    try:
        flags = common_cpu_flags(cpuinfo.read_text())
    except OSError:
        flags = set()
    core = select_core(library, platform.machine(), flags, run_probe)
    if core is not None:
        with args.github_env.open("a") as output:
            output.write(f"OPENBLAS_CORETYPE={core}\n")
        print(f"Selected OPENBLAS_CORETYPE={core} for subsequent CI steps.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"OpenBLAS check failed: {error}", file=sys.stderr)
        raise SystemExit(1)
