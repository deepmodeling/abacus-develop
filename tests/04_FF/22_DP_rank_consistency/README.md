# DeePMD MPI rank comparison

Use the `frozen_model_op.pth` attached to [issue #8310](https://github.com/deepmodeling/abacus-develop/issues/8310). The model is not copied into this repository.

From the repository root, run:

```sh
python3 tests/integrate/validation_tools/check_dp_rank_consistency.py \
  --abacus /absolute/path/to/abacus \
  --mpi /absolute/path/to/mpirun \
  --model /absolute/path/to/frozen_model_op.pth \
  --fixture tests/04_FF/22_DP_rank_consistency/STRU
```

The script runs one and four ranks for the original two-atom cell, a displaced atom, and an expanded lattice. It then compares consecutive evaluations in two-step `relax` and `cell-relax` runs on the displaced cell; these keep the same DP instance while coordinates or lattice geometry change. Four ranks guarantee at least two ranks without owned atoms. It compares energy, every force component, and every stress component using the integration-suite tolerances of `1e-5` eV, `1e-4` eV/Angstrom, and `1e-3` kbar. It also requires the coordinate displacement to produce nonzero force and each geometry update to change the energy.

This is a manual model-backed check because the issue attachment is not vendored as a test fixture. The separate `04_FF_md_output_cadence` CTest check uses a fully committed LJ fixture.
