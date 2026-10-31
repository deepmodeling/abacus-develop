"""Link portable GEMM checks to a configured production build."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expect-disabled", action="store_true", help="Expect an eligible product to bypass the 3M API")
    args = parser.parse_args()
    build = args.build.resolve(strict=True)
    output = args.output.resolve()
    output.mkdir(exist_ok=False)
    cache = (build / "CMakeCache.txt").read_text().splitlines()
    ninja = next(line.split("=", 1)[1] for line in cache if line.startswith("CMAKE_MAKE_PROGRAM:"))
    root = Path(next(line.split("=", 1)[1] for line in cache if line.startswith("CMAKE_HOME_DIRECTORY:")))
    commands = subprocess.check_output([ninja, "-t", "commands", "abacus"], cwd=build, text=True).splitlines()
    compile_template = shlex.split(next(line for line in commands if " -c " in line and "/source/main.cpp" in line))
    old_main = compile_template[compile_template.index("-o") + 1]
    link_template = shlex.split(next(line for line in reversed(commands) if " -o abacus " in line))
    if link_template[:2] == [":", "&&"]:
        link_template = link_template[2:]
    if link_template[-2:] == ["&&", ":"]:
        link_template = link_template[:-2]
    metadata = {"application_sha256": digest(build / "abacus"),
                "math_sha256": digest(root / "source/module_hsolver/kernels/cuda/math_kernel_op.cu"),
                "checks": [], "expected_disabled_dispatch": args.expect_disabled,
                "runtime": subprocess.check_output(["ldd", str(build / "abacus")], text=True)}
    names = ("check_gemm_dispatch",) if args.expect_disabled else ("check_gemm_dispatch", "check_large_gemm", "check_davidson", "check_davidson_panels")
    for name in names:
        source = root / "tools/benchmarks/deltaspin/src" / (name + ".cpp")
        saved = output / source.name
        saved.write_bytes(source.read_bytes())
        obj = saved.with_suffix(".o")
        command = compile_template.copy()
        for option, value in (("-c", saved), ("-o", obj), ("-MT", obj), ("-MF", obj.with_suffix(".d"))):
            command[command.index(option) + 1] = str(value)
        command += ["-std=c++17", "-UNDEBUG"]
        binary = output / name
        link = link_template.copy()
        link[link.index(old_main)] = str(obj)
        link[link.index("-o") + 1] = str(binary)
        if name == "check_gemm_dispatch":
            link += ["-Wl,--wrap=cublasZgemm3m"]
        with (output / (name + ".build.log")).open("x") as log:
            for invocation in (command, link):
                subprocess.run(invocation, cwd=build, stdout=log, stderr=subprocess.STDOUT, check=True)
        with (output / (name + ".log")).open("x") as log:
            result = subprocess.run([str(binary)], cwd=output, stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT)
        expected = 2 if args.expect_disabled else 0
        if result.returncode != expected:
            raise RuntimeError(f"{name}: expected exit {expected}, got {result.returncode}")
        log_path = output / (name + ".log")
        if args.expect_disabled:
            assert "rows=129 expected_3m_calls=1 actual_3m_calls=0" in log_path.read_text()
        metadata["checks"].append({"name": name, "source_sha256": digest(source), "binary_sha256": digest(binary),
                                   "log_sha256": digest(log_path), "return_code": result.returncode, "commands": [command, link]})
        print(name + (" expected disabled dispatch" if args.expect_disabled else " PASS"), flush=True)
    (output / "validation.json").write_text(json.dumps(metadata, indent=2) + "\n")


if __name__ == "__main__":
    main()
