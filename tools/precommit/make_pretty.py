#!/usr/bin/env python3

import json
import os
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

root = Path("/repo").resolve()
build_dir = root / "build"
out = Path("/out").resolve()

jobs = int(os.environ.get("JOBS", "1"))

precommit_exts = {
    ".c", ".cc", ".cpp", ".cxx", ".c++",
    ".h", ".hh", ".hpp", ".hxx",
    ".ipp", ".tpp",
    ".cu", ".cuh",
}

exclude_names = {
    ".git",
    "build",
    ".cache",
    "__pycache__",
}

def is_excluded(path: Path) -> bool:
    return any(
        part in exclude_names
        or part.startswith("cmake-build-")
        or part.startswith("build-")
        for part in path.parts
    )

def run(cmd, check=True):
    print("+ " + " ".join(str(x) for x in cmd), flush=True)
    return subprocess.run(cmd, check=check)

def collect_precommit_files():
    files = []

    for path in root.rglob("*"):
        if not path.is_file():
            continue

        rel = path.relative_to(root)

        if is_excluded(rel):
            continue

        if path.suffix in precommit_exts:
            files.append(rel)

    return sorted(files)

def language_for_file(path: Path):
    suffix = path.suffix

    if suffix == ".c":
        return "clang", "c", "-std=c11"

    if suffix in {".cu", ".cuh"}:
        return "clang++", "cuda", "-std=c++17"

    if suffix in {".h", ".hh", ".hpp", ".hxx", ".ipp", ".tpp"}:
        return "clang++", "c++-header", "-std=c++17"

    return "clang++", "c++", "-std=c++17"

def collect_include_dirs(files):
    include_dirs = {
        root,
        root / "source",
        root / "source" / "source_base",
        root / "python" / "pyabacus" / "src",
    }

    for rel in files:
        include_dirs.add((root / rel).parent.resolve())

    return sorted(include_dirs)

def generate_compile_commands(files):
    if build_dir.exists():
        shutil.rmtree(build_dir)

    build_dir.mkdir(parents=True, exist_ok=True)

    include_dirs = collect_include_dirs(files)
    commands = []

    for rel in files:
        path = (root / rel).resolve()
        compiler, lang, std_flag = language_for_file(path)

        args = [
            compiler,
            "-x", lang,
            std_flag,
            "-fsyntax-only",
            "-Wno-unknown-warning-option",
            "-Wno-unused-command-line-argument",
            *[f"-I{inc}" for inc in include_dirs],
            str(path),
        ]

        commands.append({
            "directory": str(root),
            "file": str(path),
            "arguments": args,
        })

    compile_db = build_dir / "compile_commands.json"
    compile_db.write_text(json.dumps(commands, indent=2), encoding="utf-8")

    print(f"==> Generated synthetic compile database: {compile_db}", flush=True)
    print(f"==> Synthetic entries: {len(commands)}", flush=True)

def run_clang_format(files):
    print(f"==> clang-format files: {len(files)}", flush=True)

    batch = []

    for rel in files:
        batch.append(str(rel))

        if len(batch) >= 100:
            run(["clang-format", "-i", "-style=file", "--fallback-style=GNU", *batch])
            batch.clear()

    if batch:
        run(["clang-format", "-i", "-style=file", "--fallback-style=GNU", *batch])

def run_tidy(rel):
    abs_file = str((root / rel).resolve())

    line_filter = json.dumps([
        {"name": abs_file, "lines": [[1, 100000000]]},
        {"name": str(rel), "lines": [[1, 100000000]]},
    ])

    cmd = [
        "clang-tidy",
        abs_file,
        f"-p={build_dir}",
        f"-line-filter={line_filter}",
        "--fix",
        "--quiet",
    ]

    print("+ " + " ".join(cmd), flush=True)
    ret = subprocess.run(cmd).returncode
    return str(rel), ret

def run_clang_tidy(files):
    print(f"==> clang-tidy files: {len(files)}", flush=True)

    failures = []

    with ThreadPoolExecutor(max_workers=max(1, jobs)) as pool:
        futures = [pool.submit(run_tidy, rel) for rel in files]

        for fut in as_completed(futures):
            filename, ret = fut.result()

            if ret != 0:
                failures.append((filename, ret))
                print(
                    f"WARNING: clang-tidy failed for {filename} with exit code {ret}",
                    flush=True,
                )

    if failures:
        print("==> clang-tidy failures:", flush=True)
        for filename, ret in failures:
            print(f"  {ret}: {filename}", flush=True)

def export_files(files):
    if out.exists():
        shutil.rmtree(out)

    out.mkdir(parents=True)

    for rel in files:
        src = root / rel
        dst = out / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)

    print("==> Exported C/C++ files only to /out", flush=True)

def main():
    os.chdir(root)

    files = collect_precommit_files()
    print(f"==> C/C++ files discovered: {len(files)}", flush=True)

    generate_compile_commands(files)

    run_clang_tidy(files)
    run_clang_format(files)

    export_files(files)

if __name__ == "__main__":
    main()
