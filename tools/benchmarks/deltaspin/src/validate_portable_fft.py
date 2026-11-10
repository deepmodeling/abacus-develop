"""Link production FFT checks to a configured ABACUS build and retain evidence."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    names = ("check_fused_fft_dispatch", "check_spinor_potential", "check_spinor_density",
             "check_density_bands", "check_z_fft_fallback")
    parser.add_argument("--checks", nargs="+", choices=names, default=names)
    args = parser.parse_args()
    build = args.build.resolve(strict=True)
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    sources = Path(__file__).resolve().parent
    cache = (build / "CMakeCache.txt").read_text().splitlines()
    ninja = next(line.split("=", 1)[1] for line in cache if line.startswith("CMAKE_MAKE_PROGRAM:"))
    root = Path(next(line.split("=", 1)[1] for line in cache if line.startswith("CMAKE_HOME_DIRECTORY:")))
    native = subprocess.check_output([ninja, "-t", "commands", "abacus"], cwd=build, text=True).splitlines()
    compile_template = shlex.split(next(line for line in native if " -c " in line and "/source/main.cpp" in line))
    old_main = compile_template[compile_template.index("-o") + 1]
    link_template = shlex.split(next(line for line in reversed(native) if " -o abacus " in line))
    if link_template[:2] == [":", "&&"]:
        link_template = link_template[2:]
    if link_template[-2:] == ["&&", ":"]:
        link_template = link_template[:-2]
    source_paths = [root / "source/module_basis/module_pw/CMakeLists.txt",
                    root / "source/module_basis/module_pw/kernels/fused_veff.h"]
    source_paths.extend(sorted((root / "source/module_basis/module_pw/kernels/cuda").glob("fused*")))
    source_hashes = {str(path.relative_to(root)): digest(path) for path in source_paths}
    metadata = {"build": str(build), "application_sha256": digest(build / "abacus"),
                "sources": source_hashes,
                "runtime": subprocess.check_output(["ldd", str(build / "abacus")], text=True),
                "checks": []}
    environment = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    shutil.copy2(sources / "check_spinor_potential.cpp", output / "check_spinor_potential.cpp")
    for name in args.checks:
        source = output / (name + ".cpp")
        shutil.copy2(sources / source.name, source)
        obj = source.with_suffix(".o")
        compile_command = compile_template.copy()
        for flag, value in (("-c", source), ("-o", obj), ("-MT", obj), ("-MF", obj.with_suffix(".d"))):
            compile_command[compile_command.index(flag) + 1] = str(value)
        compile_command += ["-std=c++17", "-UNDEBUG"]
        binary = output / name
        link_command = link_template.copy()
        link_command[link_command.index(old_main)] = str(obj)
        link_command[link_command.index("-o") + 1] = str(binary)
        if name == "check_z_fft_fallback":
            link_command += ["-Wl,--wrap=cudaMalloc", "-Wl,--wrap=cudaFree"]
        (output / (name + ".commands.json")).write_text(json.dumps([compile_command, link_command], indent=2) + "\n")
        with (output / (name + ".build.log")).open("x") as log:
            for command in (compile_command, link_command):
                subprocess.run(command, cwd=build, env=environment, stdin=subprocess.DEVNULL,
                               stdout=log, stderr=subprocess.STDOUT, check=True)
        log_path = output / (name + ".log")
        with log_path.open("x") as log:
            subprocess.run([str(binary)], cwd=output, env=environment, stdin=subprocess.DEVNULL,
                           stdout=log, stderr=subprocess.STDOUT, check=True)
        record = {"name": name, "source_sha256": digest(source), "binary_sha256": digest(binary),
                  "log_sha256": digest(log_path), "commands": [compile_command, link_command]}
        metadata["checks"].append(record)
        (output / (name + ".json")).write_text(json.dumps(record, indent=2) + "\n")
        print(name + " PASS", flush=True)
    for name, expected in source_hashes.items():
        if digest(root / name) != expected:
            raise RuntimeError(f"Production source changed during validation: {name}")
    (output / "validation.json").write_text(json.dumps(metadata, indent=2) + "\n")


if __name__ == "__main__":
    main()
