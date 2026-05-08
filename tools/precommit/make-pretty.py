#!/usr/bin/env python3

import json
import os
import shutil
import subprocess
from pathlib import Path

root = Path("/repo")
build_dir = Path(os.environ.get("BUILD_DIR", "build"))

cpp_format_exts = {
    ".c", ".cc", ".cpp", ".cxx", ".c++",
    ".h", ".hh", ".hpp", ".hxx",
    ".ipp", ".tpp",
    ".cu", ".cuh",
}

tidy_source_exts = {
    ".cc", ".cpp", ".cxx", ".c++",
}

exclude_names = {
    ".git",
    "build",
    build_dir.name,
    ".cache",
}

def is_excluded(path: Path) -> bool:
    return any(part in exclude_names or part.startswith("cmake-build-") for part in path.parts)

def run(cmd, check=True):
    print("+ " + " ".join(cmd), flush=True)
    return subprocess.run(cmd, check=check)

os.chdir(root)

compile_db = root / build_dir / "compile_commands.json"

if not compile_db.exists():
    cmake_args = os.environ.get("CMAKE_ARGS", "").split()
    run([
        "cmake",
        "-S", ".",
        "-B", str(build_dir),
        "-G", "Ninja",
        "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
        *cmake_args,
    ])

if not compile_db.exists():
    raise SystemExit(f"ERROR: {compile_db} was not generated")

with compile_db.open("r", encoding="utf-8") as f:
    db = json.load(f)

tidy_files = []
seen = set()

for entry in db:
    filename = entry.get("file")
    if not filename:
        continue

    path = Path(filename)
    if not path.is_absolute():
        path = Path(entry.get("directory", root)) / path

    try:
        rel = path.resolve().relative_to(root.resolve())
    except ValueError:
        continue

    if is_excluded(rel):
        continue

    if rel.suffix in tidy_source_exts and rel.exists() and rel not in seen:
        tidy_files.append(rel)
        seen.add(rel)

print(f"==> clang-tidy translation units: {len(tidy_files)}", flush=True)

tidy_failures = []
extra = os.environ.get("CLANG_TIDY_EXTRA_ARGS", "").split()
strict = os.environ.get("STRICT_CLANG_TIDY", "0") == "1"

for rel in tidy_files:
    cmd = [
        "clang-tidy",
        str(rel),
        f"-p={build_dir}",
        "--fix-errors",
        *extra,
    ]
    print("+ " + " ".join(cmd), flush=True)
    ret = subprocess.run(cmd).returncode
    if ret != 0:
        tidy_failures.append((str(rel), ret))
        print(f"WARNING: clang-tidy failed for {rel} with exit code {ret}", flush=True)

if tidy_failures:
    print("==> clang-tidy failures:", flush=True)
    for filename, ret in tidy_failures:
        print(f"  {ret}: {filename}", flush=True)
    if strict:
        raise SystemExit("ERROR: clang-tidy failed and STRICT_CLANG_TIDY=1")

format_files = []

for path in root.rglob("*"):
    if not path.is_file():
        continue

    rel = path.relative_to(root)
    if is_excluded(rel):
        continue

    if path.suffix in cpp_format_exts:
        format_files.append(rel)

print(f"==> clang-format files: {len(format_files)}", flush=True)

batch = []
for rel in format_files:
    batch.append(str(rel))
    if len(batch) >= 100:
        run(["clang-format", "-i", "-style=file", "--fallback-style=GNU", *batch])
        batch.clear()

if batch:
    run(["clang-format", "-i", "-style=file", "--fallback-style=GNU", *batch])

out = Path("/out")
if out.exists():
    shutil.rmtree(out)
out.mkdir(parents=True)

for rel in format_files:
    src = root / rel
    dst = out / rel
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)

print("==> Exported C/C++ files only to /out", flush=True)
