# TDOFDFT Al timestep validation

This case evolves to **0.5 fs** using `md_dt 0.001` fs and `md_nstep 501`.
The first MD output is the initial state, so there are 500 advances. This
preserves the physical end time of the previous `0.25 fs`, three-output case.

The electronic solver uses the MD timestep for one RK4 advance without
electronic substeps. The old trajectory had about 23.73 Ry of MD total-energy
drift and a final stress absolute sum near 194932 kbar. Regenerating that
trajectory's reference, or loosening its tolerance, would not address its
timestep error. See [issue #8272](https://github.com/deepmodeling/abacus-develop/issues/8272).

## Step selection

Compare the same initial conditions, basis, pseudopotential, and initial
optimization settings at the same physical times, not just the last output
of runs with the same number of steps:

| dt (fs) | MD outputs | End time (fs) |
| --- | --- | --- |
| 0.001 | 501 | 0.5 |
| 0.0005 | 1001 | 0.5 |
| 0.00025 | 2001 | 0.5 |

`0.001 fs` is the least expensive member of this verified small-step set.
It is not asserted to be the largest stable timestep for arbitrary inputs.
Both finer-step final properties are within this case's unchanged Autotest
thresholds relative to the new reference. The entire aligned stress trajectory
also approaches a common result as the timestep is refined.

## Local measurements

The reference was generated with a fresh PW-only CPU/MPI build from upstream
commit `f4c3dd3a135ae8bd09ac7e7d09134fd5cdfa1912`, GCC 13.3, FFTW3,
OpenBLAS 0.3.26 and OpenMPI 4.1.6, with native CPU optimization disabled.
No production source was changed. Reference generation used 4 MPI processes,
`OMP_NUM_THREADS=1` and `OPENBLAS_NUM_THREADS=1`.

| dt (fs) | Final MD energy (Ry), 1 MPI | Max energy drift (meV/atom) | Final stress absolute sum (kbar) |
| --- | --- | --- | --- |
| 0.001 | -6.848336496 | 0.716680 | 11368.026722 |
| 0.0005 | -6.848338942 | 0.705707 | 11368.024609 |
| 0.00025 | -6.848340186 | 0.700265 | 11368.023706 |

Maximum stress-component differences over **all aligned physical times**:

- dt versus dt/2: `0.005632303 kbar`.
- dt/2 versus dt/4: `0.002937635 kbar`.
- 1 versus 4 MPI at dt: `1.3002e-9 kbar`.
- 4 MPI with 1 versus 2 OpenMP threads at dt: `1.9000e-9 kbar`.

The refinement differences decrease, but this does not establish fourth-order
accuracy of the coupled electronic/ionic algorithm. The small-step energy
drift has a remaining floor; this change does not fix the propagator or prove
exact energy conservation. The separate TDOFDFT `setfullpw()` dispatch issue
is also outside this case correction.

The reference is the final electronic energy/force/stress output at the chosen
timestep, not an independently exact physical solution. The MD total energy
reported above additionally includes ionic kinetic energy.

## Reproduce and check

Run the three inputs in independent directories, changing only `md_dt` and
`md_nstep` as in the first table. Also run the selected step with 1 and 4 MPI
processes, and with the CI layout of 4 MPI / 2 OpenMP threads. Preserve all
logs. Use the cluster's supported launcher and run within allocated resources.

The checker reads logs only and has no MPI implementation dependency:

```bash
python3 check_timestep.py \
  --coarse /absolute/path/dt001-mpi1/OUT.autotest/running_md.log \
  --fine /absolute/path/dt0005-mpi1/OUT.autotest/running_md.log \
  --finer /absolute/path/dt00025-mpi1/OUT.autotest/running_md.log \
  --parallel /absolute/path/dt001-mpi4/OUT.autotest/running_md.log \
  --output /absolute/path/validation.json
```

It rejects unconverged SCF output, non-finite/missing tensors, and unequal end
times. The case-specific numerical acceptance conditions are:

- All trajectories finish at 0.5 fs, with exact time alignment using rational
  timestep values.
- Maximum MD total-energy departure from the initial value is below
  **1 meV/atom** for these two Al atoms. This trajectory-quality budget is
  separate from Autotest's final-energy comparison tolerance.
- The dt versus dt/2 trajectory stress-component difference is below
  0.01 kbar, and decreases upon the next refinement.
- Adjacent final MD energies differ by less than 1e-5 Ry.
- When supplied, the parallel layout agrees in every aligned stress component
  within 1e-6 kbar and meets the same energy-drift budget.

Require a zero checker exit status and all JSON checks true. Also run the
repository's original `tests/integrate/Autotest.sh` against the committed
reference, including the two finer-step controls at the same end time. Its
energy, force and stress comparison thresholds remain unchanged. Never use
`Autotest.sh -g` to replace a failing cross-platform result during validation.

After checking the timestep controls, the neighboring cases
`31_TDOFDFT_Al_CD`, `32_TDOFDFT_Al_mCD`, and `35_OF_LPS_bspline` should retain
their original inputs and references for a focused regression check.

For every build, retain the exact source commit, binary SHA256, compiler and
linked-library identities, CMakeCache and build flags. Retain per-run inputs,
pseudopotential hashes, launcher commands, process/thread counts, running logs
and property output. The local numbers above were regenerated locally; the
independent Intel/MKL Beijixing results are reported in the linked issue and
must still be rechecked with the final PR checkout.
