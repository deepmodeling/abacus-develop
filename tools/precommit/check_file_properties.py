#!/usr/bin/env python3

# Based on CP2K's tools/precommit/check_file_properties.py, adapted for ABACUS.

import argparse
from functools import lru_cache
import os
import pathlib
import re
import sys
from typing import Iterable, List, Set, Tuple, TypeVar

T = TypeVar("T")

# We assume this script is in tools/precommit/.
ABACUS_DIR = pathlib.Path(__file__).resolve().parents[2]

PORTABLE_FILENAME_RE = re.compile(r"^[a-zA-Z0-9._/#~=+-]*$")
OP_RE = re.compile(r"[\\|()!&><=*/+-]")
NUM_RE = re.compile(r"[0-9]+[ulUL]*")
CMAKE_OPTION_RE = re.compile(r"option\(\s*(\w+)", re.DOTALL)
MERGE_CONFLICT_RE = re.compile(
    r"^(<<<<<<<(?: .*)?|\|\|\|\|\|\|\|(?: .*)?|=======$|>>>>>>>(?: .*)?)$",
    re.MULTILINE,
)

C_EXTENSIONS = (
    ".c",
    ".cc",
    ".cpp",
    ".cxx",
    ".h",
    ".hh",
    ".hpp",
    ".hxx",
    ".cu",
    ".cuh",
    ".hip",
    ".cl",
)
HEADER_EXTENSIONS = (".h", ".hh", ".hpp", ".hxx", ".cuh")
TEXT_EXTENSIONS = C_EXTENSIONS + (
    ".cmake",
    ".md",
    ".py",
    ".sh",
    ".txt",
    ".yml",
    ".yaml",
    ".json",
)

# Keep this list intentionally broad. Project-specific additions should be made
# here rather than suppressing warnings at call sites.
FLAG_EXCEPTIONS = (
    r"\$\{.+\}\$",
    r"__.+__",
    r"_M_.+",
    r"_WIN32",
    r"_OPENMP",
    r"__cplusplus",
    r"__GNUC__",
    r"__GNUC_MINOR__",
    r"__GNUC_PATCHLEVEL__",
    r"__clang__",
    r"__clang_major__",
    r"__INTEL_COMPILER",
    r"__INTEL_LLVM_COMPILER",
    r"__NVCC__",
    r"__CUDACC__",
    r"__HIPCC__",
    r"CUDA_VERSION",
    r"HIP_VERSION",
    r"NDEBUG",
    r"M_PI",
    r"EIGEN_.+",
    r"ELPA_.+",
    r"FFTW_.+",
    r"MPI_.+",
    r"MKL_.+",
    r"LAPACK_.+",
    r"BLAS_.+",
    r"OPENBLAS_.+",
    r"LIBINT.+",
    r"LIBXC.+",
    r"CEREAL_.+",
    r"GTEST_.+",
)
FLAG_EXCEPTIONS_RE = re.compile(r"|".join(FLAG_EXCEPTIONS))


@lru_cache(maxsize=None)
def get_cmake_sources() -> str:
    files: List[pathlib.Path] = []
    root_cmake = ABACUS_DIR / "CMakeLists.txt"
    if root_cmake.exists():
        files.append(root_cmake)
    for directory in ("cmake", "toolchain", "tools"):
        base = ABACUS_DIR / directory
        if base.exists():
            files.extend(base.glob("**/*.cmake"))
            files.extend(base.glob("**/CMakeLists.txt"))
    return "\n".join(read_text_safely(fn) for fn in sorted(set(files)))


@lru_cache(maxsize=None)
def get_build_docs() -> str:
    docs = ABACUS_DIR / "docs"
    if not docs.exists():
        return ""
    files = list(docs.glob("**/*.md")) + list(docs.glob("**/*.rst"))
    return "\n".join(read_text_safely(fn) for fn in sorted(set(files)))


def read_text_safely(path: pathlib.Path) -> str:
    try:
        return path.read_text(encoding="utf8")
    except Exception:
        return ""


def check_file(path: pathlib.Path) -> List[str]:
    """Check one file for ABACUS source-tree convention violations."""
    warnings: List[str] = []

    fn_ext = path.suffix
    abspath = path.resolve()
    basefn = path.name
    is_executable = os.access(abspath, os.X_OK)

    if not PORTABLE_FILENAME_RE.match(str(path)):
        warnings += [f"Filename '{path}' not portable"]

    if not abspath.exists():
        return warnings

    raw_content = abspath.read_bytes()
    if b"\0" in raw_content:
        return warnings

    try:
        content = raw_content.decode("utf8")
    except UnicodeDecodeError:
        if fn_ext in TEXT_EXTENSIONS or not fn_ext:
            warnings += [f"{path}: is not valid UTF-8"]
        return warnings

    if "\r\n" in content:
        warnings += [f"{path}: contains DOS linebreaks"]

    if (
        fn_ext not in (".pot", ".patch", ".diff")
        and basefn not in ("Makefile",)
        and "\t" in content
    ):
        warnings += [f"{path}: contains tab character"]

    if MERGE_CONFLICT_RE.search(content):
        warnings += [f"{path}: contains merge conflict marker"]

    # check shebang, matching CP2K's behavior for executable Python scripts
    PY_SHEBANG = "#!/usr/bin/env python3"
    if fn_ext == ".py" and is_executable and not content.startswith(f"{PY_SHEBANG}\n"):
        warnings += [f"{path}: Wrong shebang, please use '{PY_SHEBANG}'"]

    # C++ headers should not introduce namespace pollution into downstream TUs.
    # Be conservative: only reject a clearly global using-directive.  Local
    # using-directives inside inline/template functions are not ideal, but they
    # do not leak into downstream translation units and should not be a first
    # version hard failure.
    if fn_ext in HEADER_EXTENSIONS:
        for i, line in enumerate(content.splitlines(), start=1):
            stripped = line.strip()
            if line == stripped and stripped.startswith("using namespace "):
                warnings += [
                    f"{path}:{i} Do not use '{stripped}' at global scope in a header file"
                ]

    # Find likely ABACUS build-controlled preprocessor condition flags and check
    # whether ABACUS' CMake files know about them.  Unlike CP2K, ABACUS has many
    # C/C++ header guards and third-party compatibility macros, so this must not
    # treat every uppercase token in #if/#ifdef/#ifndef as a project flag.
    if fn_ext in C_EXTENSIONS:
        flags = find_preprocessor_flags(content, fn_ext, basefn, path)
        cmake_sources = get_cmake_sources()
        for flag in sorted(flags):
            if flag not in cmake_sources:
                warnings += [f"{path}: Flag '{flag}' not mentioned in CMake files"]

    # Check CMake options against docs if docs exist.
    if "cmake" in str(path).lower() or basefn == "CMakeLists.txt":
        build_docs = get_build_docs()
        if build_docs:
            options = CMAKE_OPTION_RE.findall(content)
            for opt in options:
                if opt not in build_docs:
                    warnings += [f"{path}: CMake option {opt} not mentioned in docs"]

    return warnings


BUILD_FLAG_RE = re.compile(
    r"^(ENABLE_.+|USE_.+|WITH_.+|ABACUS_.+|__UT_.+|__MPI|__CUDA|__ROCM|__EXX|__LCAO|__OPENMP)$"
)
HEADER_GUARD_RE = re.compile(r"^[A-Z0-9_]+_(H|HH|HPP|HXX|CUH)_?$")
SYSTEM_OR_COMPILER_MACRO_RE = re.compile(r"^__[A-Z0-9_]+$|^__[A-Z0-9_]+__$")


def find_preprocessor_flags(
    content: str, fn_ext: str, basefn: str, path: pathlib.Path
) -> Set[str]:
    flags: Set[str] = set()
    line_continuation = False
    include_guards = find_include_guards(content)
    include_guard = basefn.upper().replace(".", "_")

    for line in content.splitlines():
        line = line.lstrip()

        if not line_continuation:
            if not line or line[0] != "#":
                continue
            if line.split()[0] not in ("#if", "#ifdef", "#ifndef", "#elif"):
                continue

        line = line.split("/*", 1)[0]
        line = line.split("//", 1)[0]
        line_continuation = line.rstrip().endswith("\\")
        line = OP_RE.sub(" ", line)
        line = line.replace("defined", " ")

        for word in line.split()[1:]:
            if NUM_RE.match(word):
                continue
            if fn_ext in HEADER_EXTENSIONS and (
                word == include_guard
                or word in include_guards
                or HEADER_GUARD_RE.fullmatch(word)
            ):
                continue
            if FLAG_EXCEPTIONS_RE.fullmatch(word):
                continue
            if SYSTEM_OR_COMPILER_MACRO_RE.fullmatch(word):
                continue
            if not BUILD_FLAG_RE.fullmatch(word):
                continue
            flags.add(word)

    return {flag for flag in flags if not FLAG_EXCEPTIONS_RE.fullmatch(flag)}


def find_include_guards(content: str) -> Set[str]:
    """Return classic #ifndef/#define include guards found near file starts.

    This catches guards such as BASE_THIRD_PARTY_CUSOLVER_H_, not just the
    basename-derived CUSOLVER_H used by the first draft.
    """
    guards: Set[str] = set()
    lines = content.splitlines()[:50]
    for idx, line in enumerate(lines[:-1]):
        m_ifndef = re.match(r"\s*#\s*ifndef\s+([A-Za-z_][A-Za-z0-9_]*)\b", line)
        if not m_ifndef:
            continue
        guard = m_ifndef.group(1)
        for next_line in lines[idx + 1 : idx + 6]:
            m_define = re.match(
                r"\s*#\s*define\s+([A-Za-z_][A-Za-z0-9_]*)\b", next_line
            )
            if m_define:
                if m_define.group(1) == guard:
                    guards.add(guard)
                break
    return guards


def pairwise(iterable: Iterable[T]) -> Iterable[Tuple[T, T]]:
    """itertools.pairwise is not available before Python 3.10."""
    a, b = iter(iterable), iter(iterable)
    next(b, None)
    return zip(a, b)


# ======================================================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Check the given FILENAME for ABACUS source-tree conventions"
    )
    parser.add_argument("files", metavar="FILENAME", type=pathlib.Path, nargs="+")
    args = parser.parse_args()

    all_warnings = []
    for fpath in args.files:
        all_warnings += check_file(fpath)

    for warning in all_warnings:
        print(warning)

    if all_warnings:
        sys.exit(1)

# EOF
