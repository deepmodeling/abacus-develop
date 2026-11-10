"""Build a focused executable against the current ABACUS objects."""

import argparse
import os
import shlex
import subprocess
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--cxxflag", action="append", default=[])
    parser.add_argument("--ldflag", action="append", default=[])
    parser.add_argument("--library", type=Path, action="append", default=[])
    args = parser.parse_args()
    source = args.source.resolve(strict=True)
    root = source.parents[1]
    build = Path(os.environ["ABACUS_HOME"]) / "build"
    cache = (build / "CMakeCache.txt").read_text().splitlines()
    ninja = next(line.split("=", 1)[1] for line in cache if line.startswith("CMAKE_MAKE_PROGRAM:"))
    commands = subprocess.check_output([ninja, "-t", "commands", "abacus"], cwd=build, text=True).splitlines()
    compile_line = next(line for line in commands if " -c " in line and "/source/main.cpp" in line)
    link_line = next(line for line in reversed(commands) if " -o abacus " in line)
    object_path = root / "results" / (source.stem + ".o")
    executable = root / "results" / source.stem
    compile_command = shlex.split(compile_line)
    original_object = compile_command[compile_command.index("-o") + 1]
    for option, value in (("-o", str(object_path)), ("-MT", str(object_path)),
                          ("-MF", str(object_path.with_suffix(".d"))), ("-c", str(source))):
        compile_command[compile_command.index(option) + 1] = value
    compile_command.extend(args.cxxflag)
    link_command = shlex.split(link_line)
    if link_command[:2] == [":", "&&"]:
        link_command = link_command[2:]
    if link_command[-2:] == ["&&", ":"]:
        link_command = link_command[:-2]
    link_command[link_command.index(original_object)] = str(object_path)
    link_command[link_command.index("-o") + 1] = str(executable)
    link_command.extend(str(path.resolve(strict=True)) for path in args.library)
    link_command.extend(args.ldflag)
    log_path = executable.with_suffix(".build.log")
    with log_path.open("w") as log:
        for command in (compile_command, link_command):
            completed = subprocess.run(command, cwd=build, stdout=log, stderr=subprocess.STDOUT)
            if completed.returncode:
                raise SystemExit(f"Build failed; see {log_path}")
    print(executable)


if __name__ == "__main__":
    main()
