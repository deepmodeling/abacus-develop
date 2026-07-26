# Issue: Decouple nspin from HContainer template parameter TR

## Background

In the current rt-TDDFT LCAO implementation, the real-space Hamiltonian container `HContainer<TR>` uses template parameter `TR` to distinguish between gamma-point calculations (`TR = double`) and multi-k-point/SOC calculations (`TR = std::complex<double>`). This type-based dispatch implicitly assumes that `TR = complex` is only needed when `nspin = 4` (non-collinear spin with SOC), leading to hard-coded `nspin=4` assumptions throughout the operator chain.

This design prevents using `HContainer<std::complex<double>>` for `nspin = 1` or `nspin = 2` calculations, which is required for:

1. **Hybrid functional TDDFT**: The EXX operator (`OperatorEXX`) needs a unified `complex<double>` HContainer to store all contributions (EXX + TD corrections) before folding to k-space, regardless of nspin.
2. **Consistent TD phase folding**: The TD phase `exp(i * A(t) * tau)` is spin-independent but mathematically complex. Folding `HContainer<double>` through `add_to_matrix(complex*, ...)` works numerically but loses type consistency and blocks future EXX integration.
3. **Unified operator chain**: Currently `TDNonlocal` and `TDEkinetic` hardcode `HContainer<complex<double>>* hR_tmp` internally, creating a type mismatch with the parent `OperatorLCAO<TK, TR>` where `TR` may be `double`.

## Problem Analysis

### Core Issue
`nspin` is conflated with template parameter `TR`. Three instances of `HContainer<T>` exist:
- `HContainer<double>` — gamma-only, nspin=1/2
- `HContainer<complex>` — multi-k, nspin=4
- No `HContainer<complex>` for nspin=1/2 (missing type combination)

### Affected Files (Hardcoded nspin=4)

| File | Line | Pattern | Severity |
|------|------|---------|----------|
| `module_operator_lcao/veff_lcao.cpp` | 125-127 | `vector<double*>(4)`, `for(is=0; is<4)` in `complex,complex` specialization | Critical |
| `module_gint/gint_vl_nspin4.h` | 39 | `const int nspin_ = 4;` hardcoded member | Critical |
| `module_gint/gint_interface.h` | 14-16 | Overload resolution ties `HContainer<complex>*` to nspin=4 path | Critical |
| `module_operator_lcao/nonlocal_force_stress.hpp` | 250, 258 | `nlm1.size() / 4`, `for(is=0; is<4)` in force/stress specializations | High |
| `module_operator_lcao/nonlocal_force_stress.hpp` | 317, 325 | Same pattern in `cal_stress_IJR` | High |
| `module_operator_lcao/operator_lcao.cpp` | 20-56 | Missing `get_hs_pointers()` specialization for `TK=double, TR=complex` | Medium |
| `module_operator_lcao/td_nonlocal_lcao.h` | 59 | `HContainer<complex<double>>* hR_tmp` hardcoded type | Medium |
| `module_operator_lcao/td_ekinetic_lcao.h` | 96 | Same as above | Medium |
| `module_operator_lcao/td_ekinetic_lcao.cpp` | 236 | `static_cast<HContainer<complex<double>>*>` on void* | Medium |

## Proposed Solution

Replace type-based dispatch with runtime `nspin` branching while preserving template flexibility for `HContainer<TR>`:

1. **Merge `Veff` specializations**: Replace three full template specializations with one generic `contributeHR()` that branches on `this->nspin` at runtime.
2. **Unify `cal_gint_vl` interface**: Add `nspin` parameter to allow dynamic selection of computational kernel instead of relying on function overload resolution.
3. **Templatize `hR_tmp`**: Change TD operators to use `HContainer<TR>* hR_tmp` consistently with the parent class template parameter.
4. **Add missing template instantiation**: Provide `get_hs_pointers()` for `OperatorLCAO<double, complex<double>>`.

## Todo List

- [ ] **v1.1** `veff_lcao.cpp`: Merge 3 full template specializations of `contributeHR()` into a single generic implementation with `if (nspin == 4)` / `else if (nspin == 2)` / `else` branching. Remove hardcoded `vector(..., 4)` and `for(is=0; is<4)`.
- [ ] **v1.2** `gint_interface.h` + `gint_interface.cpp`: Add unified `cal_gint_vl(nspin, ...)` interface that accepts either `const double*` or `vector<const double*>` and dispatches to the correct kernel based on runtime `nspin` value.
- [ ] **v1.3** `gint_vl_nspin4.h` + `gint_vl_nspin4.cpp`: Replace hardcoded `const int nspin_ = 4` with constructor parameter. Rename class to `Gint_vl_nc` (non-collinear) to reflect actual purpose.
- [ ] **v1.4** `nonlocal_force_stress.hpp`: Replace `is < 4` and `nlm1.size() / 4` with `is < npol * npol` (or `nspin`) in both `cal_force_IJR` and `cal_stress_IJR` specializations for `complex,complex`.
- [ ] **v1.5** `operator_lcao.cpp`: Add `get_hs_pointers()` specialization for `OperatorLCAO<double, complex<double>>` to prevent link errors.
- [ ] **v1.6** `td_nonlocal_lcao.h` + `.cpp`: Change `HContainer<complex<double>>* hR_tmp` to `HContainer<TR>* hR_tmp`. Update `initialize_HR_tmp`, `set_HR_fixed`, and `cal_HR_IJR` accordingly.
- [ ] **v1.7** `td_ekinetic_lcao.h` + `.cpp`: Same as v1.6 for TDEkinetic.
- [ ] **v1.8** `td_info.h`: Replace hardcoded `HContainer<complex<double>>` in `current_term`, `velocity_HR`, and all accessor methods with template parameter or abstract base.
- [ ] **v1.9** Add explicit template instantiations for all new type combinations in respective `.cpp` files.
- [ ] **v2.0** Write unit tests: Verify `HContainer<complex<double>>` works correctly with `nspin=1` and `nspin=2` configurations.
- [ ] **v2.1** Integration test: Hybrid functional TDDFT (`td_stype=2`) with `nspin=1` produces correct results after refactoring.

## Acceptance Criteria

1. Code compiles with `TR = complex<double>` and `nspin = 1` without errors or runtime crashes.
2. All existing tests pass: gamma-only (nspin=1, TR=double), multi-k (nspin=1/2, TR=complex), SOC (nspin=4, TR=complex).
3. `HContainer<complex<double>>` with `nspin=1/2` produces numerically identical results to `HContainer<double>` (imaginary parts should be zero).
4. Hybrid TDDFT (`td_stype=2`) with `nspin=1` can use unified complex HContainer for EXX + TD contribution folding.
5. No `static_cast<HContainer<complex<double>>*>` on `void*` remains in the TD operator chain.
