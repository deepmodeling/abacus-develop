"""Build an isolated NVTX-instrumented executable from the current ABACUS objects."""

import argparse
import hashlib
import json
import os
import re
import shlex
import subprocess
from pathlib import Path


def instrument_gemms(source: str) -> str:
    """Mark complete conventional and three-real-product cuBLAS calls."""
    pattern = r"(?m)^([ \t]*)cublasErrcheck\((cublasZgemm(?:3m)?)\([^;]+;"

    def mark_call(match: re.Match) -> str:
        indent, routine = match.group(1, 2)
        method = "ZGEMM3M" if routine.endswith("3m") else "ZGEMM"
        return (
            f"{indent}char shape[160];\n"
            f'{indent}std::snprintf(shape, sizeof(shape), "{method}/%c%c/%d/%d/%d/%d/%d/%d",\n'
            f"{indent}              transa, transb, m, n, k, lda, ldb, ldc);\n"
            f"{indent}nvtxRangePushA(shape);\n"
            + match.group(0) + f"\n{indent}nvtxRangePop();"
        )

    modified, count = re.subn(pattern, mark_call, source)
    if count == 0:
        raise ValueError("No double-complex cuBLAS GEMM calls found in the CUDA wrapper")
    return '#include <nvtx3/nvToolsExt.h>\n#include <cstdio>\n' + modified


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("label")
    parser.add_argument("--gemm-shapes", action="store_true",
                        help="Record dimensions and transpose modes of double-complex GEMMs")
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", args.label):
        parser.error("label must contain only letters, digits, underscores, periods, or hyphens")
    root = Path(os.environ["ABACUS_HOME"])
    output = Path(__file__).resolve().parents[1] / "results"
    build = root / "build"
    source = root / "source/source/module_base/timer.cpp"
    instrumented = output / (args.label + "-profile-timer.cpp")
    object_path = instrumented.with_suffix(".o")
    executable = output / ("abacus-" + args.label + "-profile")
    if instrumented.exists() or executable.exists():
        raise FileExistsError(f"Profile artifacts already exist for {args.label}")
    original = source.read_text()
    modified = original.replace('#include "timer.h"', '#include "module_base/timer.h"\n#include <nvtx3/nvToolsExt.h>', 1)
    marker = 'Timer_One &timer_one = timer_pool[class_name][name];'
    addition = '''
        const bool traced = class_name == "ESolver_KS_PW" || class_name == "HSolverPW"
                            || class_name == "SpinConstrain" || class_name == "Charge"
                            || class_name == "XC_Functional" || class_name == "Diago_DavSubspace"
                            || class_name == "Operator" || class_name == "FS_Nonlocal_tools"
                            || class_name == "ElecStatePW" || class_name == "Potential"
                            || class_name == "H_Hartree_pw" || class_name == "Charge_Mixing"
                            || class_name == "PW_Basis" || class_name == "PW_Basis_K"
                            || class_name == "PW_Basis_Sup" || class_name == "OnsiteProj"
                            || class_name == "Nonlocal" || class_name == "pp_cell_vnl"
                            || class_name.compare(0, 6, "Stress") == 0 || class_name == "Forces";
        if (traced && !(class_name == "OnsiteProj" && name == "cubspl_tabulate"))
        {
            if (timer_one.start_flag) { nvtxRangePushA((class_name + "/" + name).c_str()); }
            else { nvtxRangePop(); }
        }
'''
    assert modified.count(marker) == 1
    instrumented.write_text(modified.replace(marker, marker + addition, 1))
    cache = (build / "CMakeCache.txt").read_text().splitlines()
    ninja = next(line.split("=", 1)[1] for line in cache if line.startswith("CMAKE_MAKE_PROGRAM:"))
    commands = subprocess.check_output([ninja, "-t", "commands", "abacus"], cwd=build, text=True).splitlines()
    compile_line = next(line for line in commands if " -c " in line and "/module_base/timer.cpp" in line)
    link_line = next(line for line in reversed(commands) if " -o abacus " in line)
    compile_command = shlex.split(compile_line)
    original_object = compile_command[compile_command.index("-o") + 1]
    for option, value in (("-o", str(object_path)), ("-MT", str(object_path)),
                          ("-MF", str(object_path.with_suffix(".d"))), ("-c", str(instrumented))):
        compile_command[compile_command.index(option) + 1] = value
    link_command = shlex.split(link_line)
    if link_command[:2] == [":", "&&"]:
        link_command = link_command[2:]
    if link_command[-2:] == ["&&", ":"]:
        link_command = link_command[:-2]
    link_command[link_command.index(original_object)] = str(object_path)
    link_command[link_command.index("-o") + 1] = str(executable)
    link_command.append("-ldl")
    compile_commands = [compile_command]
    shape_metadata = {}
    if args.gemm_shapes:
        gemm_source = root / "source/source/module_hsolver/kernels/cuda/math_kernel_op.cu"
        gemm_original = gemm_source.read_text()
        gemm_variant = output / (args.label + "-profile-gemm.cu")
        gemm_object = gemm_variant.with_suffix(".o")
        gemm_modified = instrument_gemms(gemm_original)
        gemm_variant.write_text(gemm_modified)
        gemm_command = shlex.split(next(line for line in commands if " -c " in line and str(gemm_source) in line))
        gemm_old_object = gemm_command[gemm_command.index("-o") + 1]
        for option, value in (("-o", str(gemm_object)), ("-MT", str(gemm_object)),
                              ("-MF", str(gemm_object.with_suffix(".d"))), ("-c", str(gemm_variant))):
            gemm_command[gemm_command.index(option) + 1] = value
        compile_commands.append(gemm_command)
        link_command[link_command.index(gemm_old_object)] = str(gemm_object)
        shape_metadata = {"gemm_source_sha256": hashlib.sha256(gemm_original.encode()).hexdigest(),
                          "instrumented_gemm_sha256": hashlib.sha256(gemm_modified.encode()).hexdigest()}
    with (output / (args.label + "-profile-build.log")).open("w") as log:
        for command in compile_commands + [link_command]:
            subprocess.run(command, cwd=build, stdout=log, stderr=subprocess.STDOUT, check=True)
    metadata = {"timer_source_sha256": hashlib.sha256(original.encode()).hexdigest(),
                "instrumented_source_sha256": hashlib.sha256(instrumented.read_bytes()).hexdigest(),
                "binary_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
                "scope": "Selected host timer ranges; GPU kernels must be attributed by CUDA launch correlation.",
                **shape_metadata}
    (output / (args.label + "-profile-build.json")).write_text(json.dumps(metadata, indent=2) + "\n")
    print(executable)


if __name__ == "__main__":
    main()
