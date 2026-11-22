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
from librpa_symmetry import preserves_bvk_mesh
from run_librpa_mpi import ignore_outputs, mpi_command, use_mpi_solver


def read_symrot_operations(text):
    """Return the (unitary, antiunitary) spatial operations of symrot_R.txt.

    Each operation is a flattened row-major integer rotation. The antiunitary
    block is the optional ``Antiunitary operations: N`` trailer that ABACUS
    writes for magnetic nspin=4 runs.
    """
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    lmax = [int(re.search(r"Lmax of " + label + r":\s*(\d+)", text).group(1))
            for label in ("AOs", "ABFs")]
    orbital_rows = (max(lmax) + 1) ** 2
    start = next(i for i, line in enumerate(lines) if line == "0")

    def one_block(cursor):
        # ABACUS prints these through std::to_string, so "-1" arrives as
        # "-1.000000"; the rotations are integers only up to the print format.
        values = [float(x) for row in lines[cursor + 1:cursor + 4] for x in row.split()]
        assert len(values) == 9, "Expected a 3x3 rotation in {}".format(lines[cursor:cursor + 5])
        rotation = []
        for value in values:
            assert abs(value - round(value)) < 1e-8, "Noninteger symmetry rotation: {}".format(value)
            rotation.append(int(round(value)))
        translation = lines[cursor + 4]
        assert translation.startswith("(") and translation.endswith(")")
        return rotation

    unitary = []
    pos = start
    while pos < len(lines) and not lines[pos].startswith("Antiunitary operations:"):
        unitary.append(one_block(pos))
        pos += 5 + orbital_rows
    antiunitary = []
    if pos < len(lines):
        count = int(lines[pos].split(":", 1)[1])
        cursor = pos + 1
        for _ in range(count):
            antiunitary.append(one_block(cursor))
            cursor += 5
    return unitary, antiunitary


def check_output(case, mesh, symmetry, nspin, wfc_nspin=None,
                 require_spin_symmetry=False, require_antiunitary=False, spinor=False,
                 offset=(0.0, 0.0, 0.0)):
    out = case / "OUT.librpa"
    rows = [line.split() for line in (out / "librpa_bz_sample.txt").read_text().splitlines()]
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
    dims, _ = read_ks_wfc(out / "librpa_KS_wfc_0.dat")
    bands = (out / "librpa_band_out.txt").read_text().split()
    assert dims[:2] == (nk, wfc_nspin)
    assert tuple(map(int, bands[:4])) == (nk, nspin, dims[2], dims[3])
    manifest = json.loads((case / "librpa_producer_manifest.json").read_text())
    tolerance = next(entry["abs_tol"] for entry in manifest["required"] if entry["kind"] == "ks_wfc_v1")
    check_ks_nao(out / "librpa_KS_wfc_0.dat", case / "OUT.autotest", tolerance, spinor=spinor)
    for iq in range(1, nq + 1):
        for rank in (0, 1):
            data = (out / "librpa_V_full_{}_r{}.dat".format(iq, rank)).read_bytes()
            header = struct.unpack_from("<6i", data)
            assert header[0] == -20129433 and header[1] == iq and header[3] == 1
    assert (out / "librpa_Cs_1.dat").is_file()
    # ABACUS no longer writes stru_out.txt: taesc-prep builds it from the run
    # directory. The symmetry data this check needs now comes from
    # OUT.autotest/symrot_R.txt, which exports the unitary operations plus an
    # antiunitary trailer for magnetic nspin=4 runs.
    if symmetry == 1:
        symrot = case / "OUT.autotest" / "symrot_R.txt"
        assert symrot.is_file(), "symmetry=1 run did not export symrot_R.txt"
        unitary, antiunitary = read_symrot_operations(symrot.read_text())
        assert len(unitary) > 1, "symmetry=1 run exported a single operation"
        for row_index, rotation in enumerate(unitary):
            assert preserves_bvk_mesh(rotation, mesh, offset, False), (
                "Operation {} does not preserve the BvK mesh: {}".format(row_index, rotation)
            )
        for row_index, rotation in enumerate(antiunitary):
            assert preserves_bvk_mesh(rotation, mesh, offset, True), (
                "Antiunitary operation {} does not preserve the BvK mesh: {}".format(
                    row_index, rotation)
            )
        if require_antiunitary:
            assert antiunitary, "shifted magnetic case exported no antiunitary operation"
    return {"case": str(case), "mesh": mesh, "symmetry": symmetry,
            "nspin": nspin, "nspin_wfc": wfc_nspin, "nk_full": full, "nk_scf": nk, "nq": nq}


def mesh_variants(input_text):
    """Mesh variants for a producer case, honoring its dimensionality.

    The strict-2D Coulomb path (exx_ewald_dimension 2) requires a single kz
    layer, so its symmetry variants must stay in-plane; the generic 3D
    variants would trip the nkz = 1 guard in Ewald_Vq::init.
    """
    if re.search(r"(?m)^exx_ewald_dimension\s+2\b", input_text):
        return ((1, 3, 1), (2, 2, 1))
    return ((1, 1, 3), (2, 2, 2))


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
    for name in cases:
        input_text = (source / name / "INPUT").read_text()
        meshes = mesh_variants(input_text)
        nspin_match = re.search(r"(?m)^nspin\s+(\d+)", input_text)
        nspin = int(nspin_match.group(1)) if nspin_match else 1
        for mesh in meshes:
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
    soc_case = work / "63_librpa_soc_symmetry_222_s1"
    shutil.copytree(source / "63_librpa_soc_symmetry", soc_case, ignore=ignore_outputs)
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

    magnetic_case = work / "63_librpa_soc_symmetry_shifted_magnetic"
    shutil.copytree(source / "63_librpa_soc_symmetry", magnetic_case, ignore=ignore_outputs)
    use_mpi_solver(magnetic_case)
    magnetic_stru = []
    diagonal_vectors = iter(("3.000000 0.000000 0.000000",
                             "0.000000 3.000000 0.000000",
                             "0.000000 0.000000 3.000000"))
    replace_vectors = False
    vector_index = 0
    for line in (magnetic_case / "STRU").read_text().splitlines():
        if line == "LATTICE_VECTORS":
            replace_vectors = True
            magnetic_stru.append(line)
            continue
        if replace_vectors:
            line = next(diagonal_vectors)
            vector_index += 1
            if vector_index == 3:
                replace_vectors = False
        fields = line.split()
        if len(fields) == 6 and fields[3:] == ["0", "0", "0"]:
            line += " mag 0 0 1"
        magnetic_stru.append(line)
    (magnetic_case / "STRU").write_text("\n".join(magnetic_stru) + "\n")
    magnetic_input = magnetic_case / "INPUT"
    magnetic_input.write_text(re.sub(r"(?m)^symmetry\s+.*$", "symmetry 1", magnetic_input.read_text()))
    with magnetic_input.open("a") as stream:
        stream.write("scf_thr 1e-10\n")
    (magnetic_case / "KPT").write_text("K_POINTS\n0\nGamma\n2 3 3 0.25 0 0\n")
    with (magnetic_case / "producer.log").open("w") as log:
        subprocess.run(mpi_command(args, 2), cwd=magnetic_case, env=env,
                       stdout=log, stderr=subprocess.STDOUT, check=True, timeout=180)
    magnetic_result = check_output(magnetic_case, (2, 3, 3), 1, 1,
                                   require_spin_symmetry=True,
                                   require_antiunitary=True,
                                   spinor=True,
                                   offset=(0.25, 0.0, 0.0))
    results.append(magnetic_result)
    print(json.dumps(magnetic_result), flush=True)
    (root / "summary.json").write_text(json.dumps(results, indent=2) + "\n")
    print("Symmetry producer checks: PASS ({} baseline runs and two SOC runs)".format(len(results) - 2))


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
