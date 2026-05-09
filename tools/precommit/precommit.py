#!/usr/bin/env python3

# Based on CP2K's tools/precommit/precommit.py, adapted for ABACUS.

import ast
import argparse
import concurrent.futures
from difflib import unified_diff
from http.client import HTTPResponse
import json
import os
from pathlib import Path
import re
from subprocess import PIPE, STDOUT
import subprocess
import shutil
import sys
from time import sleep, time
from typing import cast, Iterator
from urllib.error import HTTPError
from urllib.request import Request, urlopen
import uuid

SCRATCH_DIR = Path("./obj/precommit")
CACHE_FILE = SCRATCH_DIR / "cache.json"
SERVER = os.environ.get("ABACUS_PRECOMMIT_SERVER", "https://precommit.abacus.org")

C_CPP_EXT_RE = re.compile(r".*\.(c|cc|cpp|cxx|h|hh|hpp|hxx|cu|cuh|hip|hip.hpp|cl)$")


# ======================================================================================
def main() -> None:
    parser = argparse.ArgumentParser(
        description="Check source code for formatting and linter problems."
    )
    parser.add_argument(
        "-a",
        "--no-cache",
        action="store_true",
        help="ignore the cache and check all files",
    )
    parser.add_argument(
        "-m",
        "--allow-modifications",
        action="store_true",
        help="allow the tools to modify files",
    )
    parser.add_argument(
        "-j",
        "--num_workers",
        type=int,
        default=min(16, (os.cpu_count() or 0) + 2),
        help="number of parallel workers",
    )
    parser.add_argument(
        "--progressbar-wait",
        metavar="SECONDS",
        type=int,
        default=1,
        help="number seconds in between progressbar updates",
    )
    parser.add_argument(
        "-t",
        "--timeout",
        type=int,
        default=10,
        help="server wait timeout in seconds",
    )
    parser.add_argument("files", metavar="FILE", nargs="*", help="files to process")
    args = parser.parse_args()

    print(
        f"Running precommit checks using {args.num_workers} workers and server: {SERVER} "
        f"(server wait timeout: {args.timeout} seconds)"
    )
    server_hello = (
        urlopen(Request(SERVER + "/"), timeout=args.timeout).read().decode("utf8")
    )
    assert server_hello.startswith("abacus precommit server")

    # Store candidates before changing base directory and creating scratch dir.
    file_list = [os.path.abspath(fn) for fn in args.files]
    base_dir = Path(__file__).resolve().parent.parent.parent
    os.chdir(base_dir)
    SCRATCH_DIR.mkdir(parents=True, exist_ok=True)

    # Collect candidate files.  This intentionally follows CP2K's behavior:
    # without explicit file arguments, check git-tracked files first and fall
    # back to walking the working tree only if git is not available.
    if not file_list:
        sys.stdout.write("Searching for files...\r")
        sys.stdout.flush()

        walk_list = []
        try:
            output = subprocess.check_output(["git", "ls-files"], encoding="utf8")
            for line in filter(None, output.split("\n")):
                dir_name, file_name = os.path.split(line)
                walk_list.append((os.path.join(".", dir_name), [dir_name], [file_name]))
        except Exception:
            walk_list = []

        for root, dirs, files in walk_list if walk_list else os.walk("."):
            if should_skip_root(root):
                continue
            file_list += [os.path.join(root, fn) for fn in files]

    # Filter symlinks, backup copies, logs, hidden files, and very large files.
    file_list = [fn for fn in file_list if os.path.exists(fn)]
    file_list = [fn for fn in file_list if not os.path.isdir(fn)]
    file_list = [fn for fn in file_list if not os.path.islink(fn)]
    file_list = [fn for fn in file_list if not fn[-1] in ("~", "#")]
    file_list = [fn for fn in file_list if not fn.endswith(".log")]
    file_list = [fn for fn in file_list if not os.path.basename(fn).startswith(".")]

    # Sort files by size as larger ones will take longer to process.
    file_list.sort(reverse=True, key=lambda fn: os.path.getsize(fn))

    # Load cache.  Keep CP2K's mtime-based cache semantics.
    should_load_cache = CACHE_FILE.exists() and not args.no_cache
    cache = json.loads(CACHE_FILE.read_text()) if should_load_cache else {}

    # Launch async processing of files.
    futures = {}
    executor = concurrent.futures.ThreadPoolExecutor(max_workers=args.num_workers)
    for fn in file_list:
        if os.path.getmtime(fn) != cache.get(fn, -1):
            futures[fn] = executor.submit(process_file, fn, args.allow_modifications)
    num_skipped = len(file_list) - len(futures)

    # Continuously update progressbar, save cache file, and print errors.
    failed_files = set()
    while True:
        num_done = num_skipped
        for fn, f in futures.items():
            if f.done():
                num_done += 1
                if not f.exception():
                    cache[fn] = os.path.getmtime(fn)
                elif fn not in failed_files:
                    failed_files.add(fn)
                    print_box(fn, str(f.exception()))
        CACHE_FILE.write_text(json.dumps(cache))
        if file_list:
            progressbar = "=" * int(60 * num_done / len(file_list))
            sys.stdout.write(
                f"[{progressbar:60s}] {num_done} / {len(file_list)} files processed\r"
            )
            sys.stdout.flush()
        if num_done == len(file_list) or len(failed_files) >= 10:
            executor.shutdown(wait=False)
            break
        sleep(args.progressbar_wait)

    print(
        f"Summary: Found {len(file_list)}, "
        f"skipped {num_skipped}, "
        f"checked {num_done - num_skipped}, "
        f"and failed {len(failed_files)} files." + (" " * 50)
    )
    print("Status: " + ("FAILED" if failed_files else "OK"))
    sys.exit(len(failed_files))


# ======================================================================================
def should_skip_root(root: str) -> bool:
    # Keep the CP2K style: explicit project-local exclusions in the driver.
    skipped_prefixes = (
        "./.git",
        "./obj",
        "./build",
        "./cmake-build",
        "./CMakeFiles",
        "./external",
        "./third_party",
        "./deps",
        "./lib",
        "./bin",
        "./docs/_build",
        "./__pycache__",
    )
    if root.startswith(skipped_prefixes):
        return True
    skipped_fragments = (
        "/.mypy_cache/",
        "/.pytest_cache/",
        "/__pycache__/",
        "/CMakeFiles/",
    )
    return any(fragment in root for fragment in skipped_fragments)


# ======================================================================================
def print_box(fn: str, message: str) -> None:
    print("+" + "-" * 160 + "+")
    print(f"| {fn:^158s} |")
    print("+" + "-" * 160 + "+")
    for line in message.strip().split("\n"):
        print(f"| {line:<158s} |")
    print("+" + "-" * 160 + "+\n\n")


# ======================================================================================
def process_file(fn: str, allow_modifications: bool) -> None:
    basename = Path(fn).name
    orig_content = Path(fn).read_bytes()
    bak_fn = SCRATCH_DIR / f"{basename}_{time()}.bak"
    shutil.copy2(fn, bak_fn)

    # C, C++, CUDA, HIP, OpenCL: ABACUS is C++-first-class, unlike CP2K.
    if C_CPP_EXT_RE.match(fn):
        run_remote_tool("clangformat", fn)

    # Python.
    if re.match(r".*\.py$", fn):
        ast.parse(orig_content, filename=fn)
        run_remote_tool("black", fn)

    # Shell. (TODO: Re-enable shellcheck in the future)
    #    if re.match(r".*\.sh$", fn):
    #        run_remote_tool("shfmt", fn)
    #        run_remote_tool("shellcheck", fn)

    # Markdown.
    if re.match(r".*\.md$", fn):
        run_remote_tool("mdformat", fn)

    # CMake.
    if re.match(r"(.*/CMakeLists.txt)|(.*\.cmake)$", fn):
        run_remote_tool("cmakeformat", fn)

    # Makefile.
    if re.match(r".*/Makefile", fn) or basename == "Makefile":
        run_local_tool("./tools/precommit/format_makefile.py", fn)

    # Always run local project convention checks after formatters.
    run_check_file_properties(fn)

    new_content = Path(fn).read_bytes()
    if new_content == orig_content:
        bak_fn.unlink()
    elif not allow_modifications:
        bak_fn.replace(fn)
        diff: Iterator[str]
        try:
            orig_lines = orig_content.decode("utf8").split("\n")
            new_lines = new_content.decode("utf8").split("\n")
            diff = unified_diff(orig_lines, new_lines, "before", "after", lineterm="")
        except Exception:
            diff = iter([])
        raise Exception("File modified:\n" + "\n".join(diff))


# ======================================================================================
def run_check_file_properties(fn: str) -> None:
    run_local_tool("./tools/precommit/check_file_properties.py", fn)


# ======================================================================================
def run_local_tool(*cmd: str, timeout: int = 30) -> None:
    p = subprocess.run(cmd, timeout=timeout, stdout=PIPE, stderr=STDOUT)
    if p.returncode != 0:
        raise Exception(p.stdout.decode("utf8"))


# ======================================================================================
def run_remote_tool(tool: str, fn: str) -> None:
    if os.path.getsize(fn) > 2**20:
        return  # skip files larger than 1MiB, matching CP2K behavior

    url = f"{SERVER}/{tool}"
    r = http_post(url, fn)
    if r.status == 304:
        pass
    elif r.status == 200:
        Path(fn).write_bytes(r.read())
    else:
        raise Exception(r.read().decode("utf8"))


# ======================================================================================
def http_post(url: str, fn: str) -> HTTPResponse:
    boundary = uuid.uuid1().hex
    name = Path(fn).name
    data = b"".join(
        [
            f"--{boundary}\r\nContent-Disposition: ".encode("utf8"),
            f'form-data; name="{name}"; filename="{name}"\r\n\r\n'.encode("utf8"),
            Path(fn).read_bytes(),
            f"\r\n--{boundary}--\r\n".encode("utf8"),
        ]
    )
    headers = {
        "Content-Length": f"{len(data)}",
        "Content-Type": f"multipart/form-data; boundary={boundary}",
    }
    try:
        response = urlopen(Request(url, data=data, headers=headers), timeout=60)
        return cast(HTTPResponse, response)
    except HTTPError as err:
        return cast(HTTPResponse, err)


# ======================================================================================
if __name__ == "__main__":
    main()

# EOF
