# MgO64 Full-q Cache GPU Timing - 2026-06-14

## Case

- System: 64 atom MgO supercell, Mg32 O32
- K grid: 2 2 2 Gamma-centered, symmetry-reduced to 4 k-points
- Functional: HSE
- Basis: PW
- `ecutwfc`: 50 Ry
- `nbands`: 320
- `device`: gpu
- `precision`: single
- `exxace`: 1
- `exx_separate_loop`: 1
- `exx_batch_fft_size`: 8
- `exx_symmetry_realspace`: 0
- `OMP_NUM_THREADS`: 1
- Binary: `/tmp/abacus-batch-exx-pr-p0p1-gpu-build/abacus_pw_gpu`

## Run Directories

- Full-q cache on: `/tmp/abacus-mgo64-fullq-cache-on-20260614`
- Full-q cache off: `/tmp/abacus-mgo64-fullq-cache-off-20260614`

## Results

| Metric | `exx_full_q_cache 1` | `exx_full_q_cache 0` | Delta |
|---|---:|---:|---:|
| Total wall time | 343.94 s | 381.28 s | -37.34 s |
| `OperatorEXXPW construct_ace` | 307.12 s | 349.73 s | -42.61 s |
| `OperatorEXXPW act_op_batch` | 306.90 s | 349.55 s | -42.65 s |
| Avg ACE construction | 76.78 s | 87.43 s | -10.65 s |
| Final ETOT | -61302.76734772023 eV | -61302.76744314336 eV | +0.00009542313 eV |
| Final `E_exx` | -1963.6731245711 eV | -1963.6733132720 eV | +0.0001887009 eV |
| Final gap | 6.3900381201 eV | 6.3900551504 eV | -0.0000170303 eV |

## Cache Diagnostics

Cache on:

```text
EXX full-q cache = on, effective = on, ownership = local, reduced k = 4,
full q = 8, cached local q = 8, full-q npwk_max = 179156,
estimated memory = 3499.14 MB
```

Cache off:

```text
EXX full-q cache = off, effective = off, ownership = local, reduced k = 4,
full q = 8, cached local q = 0, full-q npwk_max = 0,
estimated memory = 0 MB
```

Observed active GPU memory during cache-on EXX construction was about 5705 MiB from `nvidia-smi`.

## Timer Breakdown

Cache on:

```text
OperatorEXXPW       build_full_q_cache          5.43   5        1.09
OperatorEXXPW       construct_ace               307.12 4        76.78
OperatorEXXPW       act_op_batch                306.90 16       19.18
act_op_batch        prepare_batch               93.75  1310720
PW_Basis_K          recip_to_real_batch gpu     93.55  1310720
act_op_batch        cal_density_recip_batch     83.04  1310720
act_op_batch        recip_to_real_batch_density 77.60  1310720
```

Cache off:

```text
OperatorEXXPW       construct_ace               349.73 4        87.43
OperatorEXXPW       act_op_batch                349.55 16       21.85
act_op_batch        prepare_batch               151.29 1310720
act_op_batch        recip_to_real_batch         47.92  655360
PW_Basis_K          recip2real_remapped_batch   48.29  655360
act_op_batch        cal_density_recip_batch     33.53  1310720
act_op_batch        recip_to_real_batch_density 76.73  1310720
```

## Interpretation

The full-q cache removes repeated symmetry-remapped q-state transforms during ACE construction. In this workload it reduces ACE construction by 42.61 s despite adding 5.43 s total cache-build time, and reduces end-to-end wall time by 37.34 s. Numerical differences between cache-on and cache-off are at the single-precision noise level for this run.
