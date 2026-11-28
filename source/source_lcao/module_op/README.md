# module_op

Operator matrices of the LCAO (NAO) basis: position `r`, momentum `p` and
angular momentum `L`.

Everything in this module only needs the orbital tables, the two-center
integrators and, for the momentum operator, the pseudopotential projectors.
No SCF solution is required, so the same code path serves `calculation = scf`
and `calculation = get_s`.

## Operators and output

| Operator | INPUT keyword | Entry point | Output files |
| --- | --- | --- | --- |
| `r(R)` | `out_mat_r` | `Position_op` (`pos_op_mat.h`) | `rrx_nao.csr`, `rry_nao.csr`, `rrz_nao.csr` |
| `p = -i grad + i [V_NL, r]` | `out_mat_p` | `Momentum_op::calculate` (`mom_op_mat.h`) | `prx_nao.csr`, `pry_nao.csr`, `prz_nao.csr` |
| `Lx`, `Ly`, `Lz` | `out_mat_l` | `ModuleIO::Angmom_op` (`angmom_op_mat.h`) | `lrx_nao.csr`, `lry_nao.csr`, `lrz_nao.csr` |

The CSR files use the same self-describing header and block layout as the
H(R)/S(R) output of `module_hs`. When `out_app_flag = 2` the one-based ionic
step is inserted before `_nao` (`prxg1_nao.csr`, `lrxg1_nao.csr`, ...); with
`out_app_flag = 0` the fixed file is overwritten and with `1` it is appended.

## Layout

| File | Role |
| --- | --- |
| `pos_op_mat.{h,cpp}` | `Position_op` facade: owns the `PosOpBasis`/`PosOpCalc`/`PosOpWriter` pipeline |
| `pos_op_basis.{h,cpp}` | two-center `r` integrals on the orbital pairs |
| `pos_op_calc.{h,cpp}` | assembles `r(R)` matrix elements |
| `pos_op_writer.{h,cpp}` | writes `rrx/rry/rrz_nao.csr` |
| `pos_op_csr.{h,cpp}` | small CSR helpers shared with the writer (`lat_r_nonempty`, `assemble_csr`) |
| `velocity_op.{h,cpp}` | momentum/velocity kernel: `calculate_grad_term()` and `calculate_vcomm_r()` |
| `mom_op_mat.{h,cpp}` | `Momentum_op::calculate`: builds `p` from `velocity_op` and writes the three files |
| `angmom_op_mat.{h,cpp}` | `Angmom_op`: `Lx/Ly/Lz` through ladder operators |
| `unittests/` | unit tests, see below |

`velocity_op` is also the kernel used by `module_rt` (TDDFT) and `module_lr`,
so its interface is consumed outside this directory.

## Where new code goes

This directory owns **operator matrices whose elements come from the orbital
tables and projectors** — not from the Hamiltonian or the density matrix.

- Add a new operator here when its matrix is basis-intrinsic (no SCF quantities
  needed). Create `<name>_op_*.{h,cpp}` and register the source in
  `CMakeLists.txt` (`objects` list) and in `source/Makefile.Objects`
  (`OBJS_IO_LCAO`, LCAO-only files must not land in `OBJS_IO`).
- Keep the INPUT switch out of this module: the caller owns the decision
  (`ctrl_scf_lcao.cpp`, `esolver_gets.cpp`) and passes `Input_para` plus every
  dependency explicitly. Do not read `PARAM`, `GlobalV` or `GlobalC` here.
- Output-file naming belongs with the writer of the same operator; reuse
  `source_io/module_hs/hsr_writer.h` (`write_hcontainer_csr`) instead of
  duplicating the CSR header.
- H(R)/S(R)/T(k)/dH(R)/Vxc writers stay in `source_io/module_hs`; they either
  depend on the Hamiltonian or belong to the SCF output family. Today
  `module_hs/output_mat_sparse.cpp` includes `pos_op_mat.h` while this module
  includes `module_hs` CSR helpers, so keep new cross-module dependencies
  one-directional and header-light rather than letting the two directories
  grow into an include cycle.
- Unit tests live in `unittests/` as `test_<file_under_test>.cpp`, with a
  matching target in `unittests/CMakeLists.txt`. Write scratch files to
  `testing::TempDir()`, never to the working directory.

## Conventions and limits

- LCAO only: built inside `if(ENABLE_LCAO)`.
- The NAO basis is real, so `p` and `L` matrix elements are purely imaginary.
  The stored real part is identically zero by construction.
- `nspin = 4` is not supported for `out_mat_p`: the matrix is written as one
  scalar (`npol = 1`) block. The parameter parser rejects the combination.
- Ionic steps in filenames and headers are one-based. `get_s` has no ionic
  step; it passes `istep = -1`, which suppresses the `g{N}` suffix and writes
  the step-less names.
- Units are Hartree atomic units; the `/2` in the nonlocal commutator kernel
  converts the D matrix from Ry to Hartree.

## Build and test

```bash
cmake -S . -B build
cmake --build build -j --target abacus MODULE_LCAO_pos_op_csr_test MODULE_LCAO_angmom_op_mat_test
ctest --test-dir build -R "MODULE_LCAO_pos_op_csr_test|MODULE_LCAO_angmom_op_mat_test" --output-on-failure
```

The integration coverage for the momentum output is the Autotest case
`tests/03_NAO_multik/get_s`, which compares `prx/pry/prz_nao.csr` against the
stored references.
