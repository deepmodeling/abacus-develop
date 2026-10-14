"""Check real full/irreducible-grid producer handoffs using temporary cases."""

import argparse
from collections import Counter
import json
import math
import os
from pathlib import Path
import re
import shutil
import struct
import subprocess
import tempfile

from librpa_wavefunctions import check_ks_nao, read_ks_wfc
from run_librpa_mpi import ignore_outputs, mpi_command, use_mpi_solver


def check_output(case, mesh, symmetry, nspin, wfc_nspin=None,
                 require_spin_symmetry=False, spinor=False):
    out = case / "OUT.librpa"
    rows = [line.split() for line in (out / "bz_sample.txt").read_text().splitlines()]
    assert tuple(map(int, rows[0])) == mesh
    nk, nq = map(int, rows[1])
    full = mesh[0] * mesh[1] * mesh[2]
    assert 0 < nk <= full and nk == nq
    weights = [float(row[1]) for row in rows[2:2 + nk]]
    assert abs(sum(weights) - 1) < 1e-12
    for ik, row in enumerate(rows[2:2 + nk], 1):
        assert len(row) == 10 and int(row[0]) == ik
        assert list(map(int, row[-2:])) == [ik, ik]
    if nk < full:
        tail = rows[2 + nk:]
        assert tail[0] == ["full_kmap", str(full)] and len(tail) == full + 1
        reps = Counter()
        coordinates = set()
        for ik, row in enumerate(tail[1:], 1):
            assert len(row) == 5 and int(row[0]) == ik
            rep = int(row[4])
            assert 1 <= rep <= nk
            point = tuple(float(x) for x in row[1:4])
            assert all(math.isfinite(x) for x in point) and point not in coordinates
            coordinates.add(point)
            reps[rep] += 1
        for ik, weight in enumerate(weights, 1):
            assert abs(weight - reps[ik] / full) < 1e-12
    else:
        assert len(rows) == nk + 2
    if symmetry == -1:
        assert nk == full
    if symmetry == 1 or (symmetry == 0 and mesh == (1, 1, 3)):
        assert nk < full, "The test must exercise an actual reduction"
    if wfc_nspin is None:
        wfc_nspin = nspin
    dims, _ = read_ks_wfc(out / "KS_wfc_0.dat")
    bands = (out / "band_out.txt").read_text().split()
    assert dims[:2] == (nk, wfc_nspin)
    assert tuple(map(int, bands[:4])) == (nk, nspin, dims[2], dims[3])
    manifest = json.loads((case / "librpa_producer_manifest.json").read_text())
    tolerance = next(entry["abs_tol"] for entry in manifest["required"] if entry["kind"] == "ks_wfc_v1")
    check_ks_nao(out / "KS_wfc_0.dat", case / "OUT.autotest", tolerance, spinor=spinor)
    for iq in range(1, nq + 1):
        for rank in (0, 1):
            data = (out / "V_full_{}_r{}.dat".format(iq, rank)).read_bytes()
            header = struct.unpack_from("<6i", data)
            assert header[0] == -20129433 and header[1] == iq and header[3] == 1
    assert (out / "Cs_1.dat").is_file()
    stru = (out / "stru_out.txt").read_text().splitlines()
    if symmetry == 1:
        begin = 7 + int(stru[6])
        count, convention = stru[begin].split()
        assert int(count) > 1 and convention == "row"
        row_end = begin + 1 + int(count)
        assert len(stru) >= row_end
        assert all(len(row.split()) == 12 for row in stru[begin + 1:row_end])
        for row in stru[begin + 1:row_end]:
            rotation = list(map(int, row.split()[:9]))
            for i in range(3):
                for j in range(3):
                    step = rotation[3 * i + j] * mesh[i] / mesh[j]
                    assert abs(step - round(step)) < 1e-12, "Operation does not preserve the BvK mesh"
        spin_start = row_end
        if require_spin_symmetry or spin_start < len(stru):
            assert spin_start < len(stru), "SOC stru_out has no spin_symmetry block"
            header = stru[spin_start].split()
            assert header[:1] == ["spin_symmetry"] and len(header) == 3
            assert int(header[1]) in (0, 1)
            spin_source = int(header[2])
            assert spin_source in (0, 1)
            spin_rows = stru[spin_start + 1:]
            assert len(spin_rows) == int(count)
            expected_fields = 9 if spin_source == 1 else 1
            for row in spin_rows:
                fields = row.split()
                assert len(fields) == expected_fields and int(fields[0]) in (0, 1)
                for value in fields[1:]:
                    assert math.isfinite(float(value))
    return {"case": str(case), "mesh": mesh, "symmetry": symmetry,
            "nspin": nspin, "nspin_wfc": wfc_nspin, "nk_full": full, "nk_scf": nk, "nq": nq}


def run(args, root):
    source = args.cases.resolve()
    work = root / "08_RI"
    work.mkdir()
    (root / "PP_ORB").symlink_to(source.parent / "PP_ORB", target_is_directory=True)
    results = []
    env = dict(os.environ,
               OMP_NUM_THREADS="1",
               MKL_NUM_THREADS="1",
               ABACUS_MPIEXEC=args.mpirun,
               ABACUS_MPIEXEC_NUMPROC_FLAG=args.mpi_np_flag,
               ABACUS_MPIEXEC_PREFLAGS=" ".join(args.mpi_preflag),
               ABACUS_MPIEXEC_POSTFLAGS=" ".join(args.mpi_postflag))
    cases = (source / "CASES_LIBRPA_PRODUCER.txt").read_text().split()
    for name, nspin in zip(cases[:2], (1, 2)):
        for mesh in ((1, 1, 3), (2, 2, 2)):
            for symmetry in (-1, 0, 1):
                case = work / "{}_{}_s{}".format(name, "".join(map(str, mesh)), symmetry)
                shutil.copytree(source / name, case, ignore=ignore_outputs)
                use_mpi_solver(case)
                inp = case / "INPUT"
                inp.write_text(re.sub(r"(?m)^symmetry\s+.*$", "symmetry {}".format(symmetry), inp.read_text()))
                with inp.open("a") as stream:
                    stream.write("scf_thr 1e-10\n")
                if nspin == 2:
                    # Fix the spin populations in these comparisons. Otherwise
                    # the metallic full-grid run can spontaneously magnetize.
                    with inp.open("a") as stream:
                        stream.write("nupdown 1e-12\n")
                (case / "KPT").write_text("K_POINTS\n0\nGamma\n{} {} {} 0 0 0\n".format(*mesh))
                with (case / "producer.log").open("w") as log:
                    subprocess.run(mpi_command(args, 2), cwd=case, env=env,
                                   stdout=log, stderr=subprocess.STDOUT, check=True, timeout=120)
                result = check_output(case, mesh, symmetry, nspin)
                results.append(result)
                print(json.dumps(result), flush=True)
    soc_case = work / "63_LIBRPA_SOC_SYMMETRY_222_s1"
    shutil.copytree(source / "63_LIBRPA_SOC_SYMMETRY", soc_case, ignore=ignore_outputs)
    use_mpi_solver(soc_case)
    soc_input = soc_case / "INPUT"
    soc_input.write_text(re.sub(r"(?m)^symmetry\s+.*$", "symmetry 1", soc_input.read_text()))
    with soc_input.open("a") as stream:
        stream.write("scf_thr 1e-10\n")
    (soc_case / "KPT").write_text("K_POINTS\n0\nGamma\n2 2 2 0 0 0\n")
    with (soc_case / "producer.log").open("w") as log:
        subprocess.run(mpi_command(args, 2), cwd=soc_case, env=env,
                       stdout=log, stderr=subprocess.STDOUT, check=True, timeout=180)
    soc_result = check_output(soc_case, (2, 2, 2), 1, 1,
                              require_spin_symmetry=True, spinor=True)
    results.append(soc_result)
    print(json.dumps(soc_result), flush=True)
    (root / "summary.json").write_text(json.dumps(results, indent=2) + "\n")
    print("Symmetry producer checks: PASS (12 two-rank runs and one SOC run)")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--abacus", required=True)
    parser.add_argument("--cases", required=True, type=Path)
    parser.add_argument("--mpirun", default="mpirun")
    parser.add_argument("--mpi-np-flag", default="-np")
    parser.add_argument("--mpi-preflag", action="append", default=[])
    parser.add_argument("--mpi-postflag", action="append", default=[])
    parser.add_argument("--work", type=Path, help="Retain run artifacts in a new directory")
    args = parser.parse_args()
    if args.work:
        args.work.mkdir(parents=True)
        run(args, args.work.resolve())
    else:
        with tempfile.TemporaryDirectory(prefix="abacus-librpa-sym-") as directory:
            run(args, Path(directory))


if __name__ == "__main__":
    main()
