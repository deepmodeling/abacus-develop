#!/usr/bin/env python3
"""Compare two RI backends using identical GPU LCAO HSE inputs."""
import argparse
import json
import os
import re
import subprocess
import time
from pathlib import Path


def derivatives(log, size):
    """Read the final total force (eV/Angstrom) and stress (kbar)."""
    force_block = log.rsplit('#TOTAL-FORCE (eV/Angstrom)#', 1)[1]
    forces = re.findall(r'^\s*Si\d+\s+([-+\d.eE]+)\s+([-+\d.eE]+)\s+([-+\d.eE]+)', force_block, re.M)
    stress_block = log.rsplit('#TOTAL-STRESS (kbar)#', 1)[1]
    stress = re.findall(r'^\s*([-+\d.eE]+)\s+([-+\d.eE]+)\s+([-+\d.eE]+)\s*$', stress_block, re.M)[:3]
    if len(forces) != size or len(stress) != 3:
        raise RuntimeError('Missing force/stress output')
    return {key: [float(x) for row in values for x in row]
            for key, values in [('force', forces), ('stress', stress)]}


def write_inputs(work, source, size, orbital):
    points = [(0, 0, 0), (0, .5, .5), (.5, 0, .5), (.5, .5, 0),
              (.25, .25, .25), (.25, .75, .75), (.75, .25, .75), (.75, .75, .25)]
    nx = 2 if size == 32 else 1
    positions = [((x+i)/nx, (y+j)/nx, z) for i in range(nx)
                 for j in range(nx) for x, y, z in points]
    work.mkdir()
    (work / 'STRU').write_text(
        'ATOMIC_SPECIES\nSi 28.0855 Si_dojo_nsoc.upf\nNUMERICAL_ORBITAL\n'
        f'{orbital}\nLATTICE_CONSTANT\n10.26\nLATTICE_VECTORS\n'
        f'{nx} 0 0\n0 {nx} 0\n0 0 1\nATOMIC_POSITIONS\nDirect\nSi\n0\n{size}\n'
        + ''.join('%g %g %g 0 0 0\n' % p for p in positions))
    (work / 'KPT').write_text('K_POINTS\n0\nGamma\n1 1 1 0 0 0\n')
    values = dict(
        suffix='bench', calculation='scf', basis_type='lcao', device='gpu',
        ks_solver='cusolver', gamma_only=1, nspin=1, ntype=1, nbands=int(size*2.5),
        pseudo_dir=str(source/'tests/PP_ORB'), orbital_dir=str(source/'tests/PP_ORB'),
        ecutwfc=30, scf_thr='1e-9', scf_ene_thr='1e-8', scf_nmax=200, symmetry=0,
        smearing_method='gaussian', smearing_sigma=.002, mixing_beta=.1,
        dft_functional='hse', exx_real_number=1, exx_pca_threshold='1e-3',
        exx_ccp_rmesh_times=1, exx_separate_loop=1, exx_hybrid_step=50,
        cal_force=1, cal_stress=1)
    (work/'INPUT').write_text('INPUT_PARAMETERS\n' + ''.join(f'{k} {v}\n' for k, v in values.items()))


def run_case(binary, work, size):
    with (work/'gpu-memory.csv').open('w') as memory:
        monitor = subprocess.Popen(
            ['nvidia-smi', '--query-compute-apps=pid,used_gpu_memory',
             '--format=csv,noheader,nounits', '-lms', '100'],
            stdout=memory, stderr=subprocess.STDOUT)
        start = time.monotonic()
        try:
            with (work/'stdout.log').open('w') as stdout:
                result = subprocess.run(
                    ['timeout', '--kill-after=30s', '1200s', 'mpirun', '-np', '1', str(binary)],
                    cwd=work, stdout=stdout, stderr=subprocess.STDOUT, check=False)
        finally:
            monitor.terminate()
            monitor.wait()
    seconds = time.monotonic() - start
    if result.returncode:
        raise RuntimeError(f'{work}: exit code {result.returncode}')
    log = (work/'OUT.bench/running_scf.log').read_text()
    stdout = (work/'stdout.log').read_text()
    energy = re.findall(r'!FINAL_ETOT_IS\s+([-+\d.eE]+)', log+'\n'+stdout)
    deltas = [float(x) for x in re.findall(r'deltaE \(eV\) from outer loop:\s*([-+\d.eE]+)', stdout)]
    # An inner SCF marker alone does not establish outer EXX convergence.
    if not energy or '#SCF IS CONVERGED#' not in log or not deltas or deltas[-1] >= 1e-8 or len(deltas) >= 50:
        raise RuntimeError(f'{work}: outer energy convergence not established')
    timing = {cls+'.'+fn: float(value) for cls, fn, value in re.findall(
        r'^\s*(Exx_LRI|RI_2D_Comm|HSolverLCAO)\s+(\w+)\s+([\d.]+)', stdout, re.M)}
    memory = [int(x) for x in re.findall(r'^\s*\d+,\s*(\d+)\s*$', (work/'gpu-memory.csv').read_text(), re.M)]
    return dict(seconds=seconds, energy_ev=float(energy[-1]), outer_steps=len(deltas),
                outer_delta_ev=deltas[-1], timing=timing,
                sampled_max_process_mib=max(memory) if memory else None,
                **derivatives(log, size))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--abacus', type=Path, required=True, help='GPU RI executable')
    parser.add_argument('--baseline', type=Path, default=os.environ.get('ABACUS_CPU_RI'),
                        help='CPU RI executable, built with the same LibRI and GPU solver')
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--orbital', default='Si_dojo_6au_sz.orb')
    parser.add_argument('--repeats', type=int, default=2)
    args = parser.parse_args()
    if args.baseline is None or args.repeats < 1:
        parser.error('A baseline executable and positive repeat count are required')
    args.output.mkdir(parents=True, exist_ok=True)
    rows = []
    comparisons = []
    for size in (8, 32):
        for repeat in range(args.repeats):
            pair = {}
            order = ('cpu_ri', 'gpu_ri') if repeat % 2 == 0 else ('gpu_ri', 'cpu_ri')
            for backend in order:
                work = args.output/f'si{size}_{backend}_{repeat}'
                write_inputs(work, args.source, size, args.orbital)
                binary = args.baseline if backend == 'cpu_ri' else args.abacus
                row = dict(size=size, repeat=repeat, backend=backend,
                           orbital=args.orbital, **run_case(binary, work, size))
                rows.append(row)
                pair[backend] = row
                (args.output/'measurements.json').write_text(json.dumps(rows, indent=2)+'\n')
                print(json.dumps({k: v for k, v in row.items() if k not in ('force', 'stress')}), flush=True)
            cpu, gpu = pair['cpu_ri'], pair['gpu_ri']
            comparison = dict(size=size, repeat=repeat,
                              energy_ev_per_atom=abs(cpu['energy_ev']-gpu['energy_ev'])/size,
                              force_ev_angstrom=max(abs(x-y) for x, y in zip(cpu['force'], gpu['force'])),
                              stress_kbar=max(abs(x-y) for x, y in zip(cpu['stress'], gpu['stress'])),
                              same_outer_steps=cpu['outer_steps'] == gpu['outer_steps'])
            comparisons.append(comparison)
            (args.output/'comparison.json').write_text(json.dumps(comparisons, indent=2)+'\n')
            if (comparison['energy_ev_per_atom'] > 1e-6 or comparison['force_ev_angstrom'] > 1e-4
                    or comparison['stress_kbar'] > 1e-3 or not comparison['same_outer_steps']):
                raise RuntimeError(f'Backend comparison failed: {comparison}')


if __name__ == '__main__':
    main()
