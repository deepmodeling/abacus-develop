#!/usr/bin/env python3
"""Reject direct legacy-global access from module_hs C++ sources."""

from pathlib import Path
import re
import sys


FORBIDDEN = (
    re.compile(r"\bPARAM\s*\."),
    re.compile(r"\bGlobalV\s*::"),
    re.compile(r"\bGlobalC\s*::"),
    re.compile(r'#\s*include\s*[<\"][^\">]*(?:parameter\.h|global_variable\.h|exx_info\.h)[\">]'),
)
SOURCE_SUFFIXES = {".h", ".hpp", ".cpp", ".cc"}


def main() -> int:
    module_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent
    violations = []
    for source in sorted(module_dir.iterdir()):
        if source.suffix not in SOURCE_SUFFIXES:
            continue
        for line_number, line in enumerate(source.read_text(encoding="utf-8").splitlines(), 1):
            if any(pattern.search(line) for pattern in FORBIDDEN):
                violations.append(f"{source}:{line_number}:{line.strip()}")

    if violations:
        print("module_hs directly depends on legacy global state:", file=sys.stderr)
        print("\n".join(violations), file=sys.stderr)
        return 1

    print("module_hs legacy-global scan: 0 violations")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
