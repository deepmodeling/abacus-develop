"""
This example validates and benchmarks ABACUS socket I/O from ASE.

ASE runs as the i-PI socket server and ABACUS runs as the socket client.
ABACUS keeps calculation=scf and enables socket_driver internally.

The script checks two PR-review relevant points:
1. socket SCF gives the same energy and forces as a normal non-socket SCF;
2. repeated socket calculations avoid relaunching ABACUS and are faster than
   the normal FileIO calculator for a sequence of SCF force evaluations.

The i-PI protocol does not carry element symbols. AbacusSocketIO handles
the required STRU/socket atom-order alignment internally and returns forces
in the original ASE Atoms order.
"""
import os
import shutil
import time
from pathlib import Path

import numpy as np
from ase import Atoms
from abacuslite import Abacus, AbacusProfile, AbacusSocketIO

here = Path(__file__).parent
pporb = here.parent.parent.parent / 'tests' / 'PP_ORB'

aprof = AbacusProfile(
    command=os.environ.get('ABACUS_COMMAND', 'mpirun -np 4 abacus'),
    pseudo_dir=pporb,
    orbital_dir=pporb,
    omp_num_threads=1,
)

common_kwargs = {
    'pseudopotentials': {'Si': 'Si_ONCV_PBE-1.0.upf'},
    'basissets': {'Si': 'Si_gga_8au_100Ry_2s2p1d.orb'},
    'inp': {
        'calculation': 'scf',
        'nspin': 1,
        'basis_type': 'lcao',
        'ks_solver': 'scalapack_gvx',
        'ecutwfc': 30,
        'symmetry': 0,
        'kspacing': 0.5,
        'scf_thr': 1e-8,
        'scf_nmax': 40,
        'chg_extrap': 'atomic',
        'cal_force': 1,
    },
}

base_atoms = Atoms(
    'Si2',
    positions=[[0.0, 0.0, 0.0], [1.25, 1.25, 1.25]],
    cell=[5.43, 5.43, 5.43],
    pbc=True,
)


def clean(directory):
    shutil.rmtree(directory, ignore_errors=True)


def run_fileio(atoms, directory):
    clean(directory)
    calc = Abacus(profile=aprof, directory=str(directory), **common_kwargs)
    atoms = atoms.copy()
    atoms.calc = calc
    forces = atoms.get_forces()
    energy = atoms.get_potential_energy()
    return energy, forces


def run_socketio(atoms, directory, socket_name):
    clean(directory)
    calc = AbacusSocketIO(
        profile=aprof,
        directory=str(directory),
        unixsocket=socket_name,
        timeout=120,
        **common_kwargs,
    )
    atoms = atoms.copy()
    with calc:
        atoms.calc = calc
        energy = atoms.get_potential_energy()
        forces = atoms.get_forces()
    return energy, forces


def displaced_structures():
    structures = []
    for scale in (0.00, 0.03, -0.02, 0.05):
        atoms = base_atoms.copy()
        atoms.positions[1] += scale
        structures.append(atoms)
    return structures


fileio_dir = here / 'socketio_fileio_scf'
socket_dir = here / 'socketio_socket_scf'
bench_fileio_dir = here / 'socketio_bench_fileio'
bench_socket_dir = here / 'socketio_bench_socket'

try:
    reference_energy, reference_forces = run_fileio(base_atoms, fileio_dir)
    socket_energy, socket_forces = run_socketio(
        base_atoms, socket_dir, 'abacus_si_check')

    energy_diff = abs(socket_energy - reference_energy)
    force_diff = np.max(np.abs(socket_forces - reference_forces))
    print(f'FileIO SCF energy: {reference_energy:.12f} eV')
    print(f'Socket SCF energy: {socket_energy:.12f} eV')
    print(f'|dE|: {energy_diff:.3e} eV')
    print(f'max |dF|: {force_diff:.3e} eV/Angstrom')
    assert energy_diff < 1e-4
    assert force_diff < 1e-5

    structures = displaced_structures()

    clean(bench_fileio_dir)
    fileio_calc = Abacus(
        profile=aprof,
        directory=str(bench_fileio_dir),
        **common_kwargs,
    )
    t0 = time.perf_counter()
    for atoms in structures:
        atoms = atoms.copy()
        atoms.calc = fileio_calc
        atoms.get_forces()
    fileio_seconds = time.perf_counter() - t0

    clean(bench_socket_dir)
    socket_calc = AbacusSocketIO(
        profile=aprof,
        directory=str(bench_socket_dir),
        unixsocket='abacus_si_bench',
        timeout=120,
        **common_kwargs,
    )
    t0 = time.perf_counter()
    with socket_calc as calc:
        for atoms in structures:
            atoms = atoms.copy()
            atoms.calc = calc
            atoms.get_forces()
    socket_seconds = time.perf_counter() - t0

    speedup = fileio_seconds / socket_seconds
    print(f'FileIO repeated SCF force time: {fileio_seconds:.2f} s')
    print(f'Socket repeated SCF force time: {socket_seconds:.2f} s')
    print(f'Socket speedup vs FileIO: {speedup:.2f}x')
finally:
    for directory in (fileio_dir, socket_dir, bench_fileio_dir, bench_socket_dir):
        clean(directory)
