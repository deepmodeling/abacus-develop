# DMR pairing conventions in LCAO force/stress (DM/DMK/DMR chain)

- Date: 2026-08-09
- Related code: `source/source_estate/module_dm/` (DMK/DMR),
  `source/source_hamilt/module_hcontainer/func_folding.cpp` (folding_HR),
  `source/source_lcao/module_operator_lcao/` (overlap/ekinetic/nonlocal
  force and stress)
- This document fixes the DMK->DMR Fourier sign convention and the
  contraction mapping between the real-space density matrix and the two-center
  integrals (including the nonlocal R2-R1 pairing) used by the force/stress
  code. It is the derivation referenced by the code comments (the
  "See Sec. 3 ..." comments in `density_matrix.cpp`).

## 1. Index and R-lattice conventions

- Wave-function coefficients `C_{μn}(k)`: μ = orbital (row), n = band
  (column), stored as `psi::Psi`'s `wfc(ib, iw)`.
- Density-matrix elements `D_{μν}`: μ = row (bra) orbital, ν = column (ket)
  orbital.
- Lattice vector R: the matrix block of `HContainer::AtomPair(iat1, iat2, R)`
  is `O_{μν}(R)` with μ ∈ iat1 (home cell) and ν ∈ iat2 (cell R); the vector
  satisfies `dtau = τ(iat2) + R·L − τ(iat1)` (`UnitCell::cal_dtau`).
- k and R are both in **direct coordinates**; the phase `k·R` carries a
  factor of `2π`.

## 2. k-space and real-space Fourier pair

- Forward transform (operators, `folding_HR`, `func_folding.cpp`):

  ```
  O(k) = Σ_R e^{+ikR} O(R)
  ```

- Inverse transform (density matrix, `DensityMatrix_Tools::cal_DMR*`,
  `density_matrix.cpp`):

  ```
  D(R) = Σ_k e^{-ikR} DMK(k)
  ```

- Normalization: the k-point weights `w_k` are embedded in `DMK`
  (`wg = w_k·occ`, `occupy.cpp`), so `cal_DMR` must **not** multiply by
  `1/Nk` again.
- `DMK(μ,ν;k) = Σ_n f_nk C*_{μn}(k) C_{νn}(k) = (C f C†)^T = D_std^T`
  (`cal_dm_psi.cpp`, `zgemm('N','T')` on pre-conjugated coefficients).

## 3. Forward and inverse transforms are mutual inverses (Sec. 3 convention)

For the same set of R vectors, `e^{+ikR}` (folding_HR) and `e^{-ikR}`
(cal_DMR) form a self-consistent forward/inverse Fourier pair. If the
`-sinp` in `density_matrix.cpp` were changed back to `+sinp`, then `D(R)`
would actually store the forward transform `Σ_k e^{+ikR} DMK(k)`, which is no
longer the inverse of `folding_HR`; non-diagonal / non-symmetric consumers
(such as the nspin=4 nonlocal force) then become wrong. This sign is locked
by the unit test `T1_fourier_round_trip` (see below).

## 4. General form of the force/stress

The Feynman–Hellmann-type forces (overlap/ekinetic/nonlocal etc.) read:

```
F = Σ_k Tr[ D(k) · ∂O(k)/∂τ ]
  = Σ_k Σ_{μν} D_{μν}(k) ∂O_{νμ}(k)/∂τ
```

Substituting `∂O(k) = Σ_R e^{+ikR} ∂O(R)`:

```
F = Σ_R Σ_{μν} ∂O_{νμ}(R) D_{μν}(−R)          (*)
   with D(−R) := Σ_k e^{+ikR} DMK(k)
```

The code stores `D(R) = Σ_k e^{-ikR} DMK(k)`; the two are equivalent through
the closed-trace argument of Sec. 5.

## 5. Closed-trace conjugate protection (Sec. 5 pairing)

In the two-center integral path the code contracts same-index pairs for every
`(iat1,iat2,R)`:

```
W_code = Σ_{μνR} D(μ,ν;R) ∂O(μ,ν;R)
```

Substituting `D(R) = Σ_k e^{-ikR} DMK(k)`:

```
W_code = Σ_k Σ_{μν} DMK(μ,ν;k) [Σ_R e^{-ikR} ∂O(μ,ν;R)]
       = Σ_k Tr( DMK(k) · ∂O'(−k) )
```

where `∂O'(−k) := Σ_R e^{-ikR} ∂O(R)` is the folded operator evaluated at
`−k`. If `∂O` is stored in the HContainer with **full-direction pairing**
(for every `(iat1,iat2,R)` there is `(iat2,iat1,−R)` with
`∂O(iat2,iat1,−R) = ∂O(iat1,iat2,R)†`), then `∂O'(k) = Σ_R e^{+ikR} ∂O(R)`
is Hermitian for each k and `∂O'(−k) = ∂O'(k)†`. Hence

```
W_code = Σ_k w_k Tr( DMK(k) · ∂O'(k)† )
       = Σ_k w_k conj( Tr( DMK(k) ∂O'(k) ) )
       = Σ_k w_k Tr( DMK(k) ∂O'(k) )        [trace of Hermitian pair is real]
       = Σ_k w_k Re Tr( DMK(k) ∂O'(k) )
```

That is, the **total is protected by the conjugate** rather than by
"term-by-term equality". The protection holds under three conditions:

1. The contraction is a **closed Frobenius trace over all orbitals**
   (no element-wise consumption);
2. The paired operator is **Hermitian per k** (HContainer full-direction
   paired storage);
3. Only the **real part** is kept.

Any consumer that breaks one of these conditions is exposed to an O(1)
pairing error.

### 5.1 overlap / ekinetic (two-center integral path)

`cal_force_stress_2center` in `operator_fs_utils.hpp` contracts
`D(μ,ν;R) ∂O(μ,ν;R)` with same indices for every `(iat1,iat2,R)`, and
`(iat2,iat1,−R)` is visited once as well (full R-set iteration), hence
`finalize_force_stress(force_factor = 1.0)` (no factor 2 needed).

### 5.2 nonlocal (β projection, R2−R1 pairing)

`nonlocal_fs.cpp::cal_force_stress` enumerates two neighbor atoms relative to
the central atom `iat0`: `iat1` (cell `R1`) and `iat2` (cell `R2`), and
consumes the DMR block

```
dmR->find_matrix(iat1, iat2, R2 − R1)
```

i.e. the pairing structure is `R_vector = R2 − R1` (`R_index2 − R_index1` in
`nonlocal_fs.cpp`), which is not a simple same-R closed trace. The force on
this path reads

```
F(iat0) = Σ_{iat1,iat2} <∂β_{iat1,R1}/∂τ_iat0 | D_{iat1,iat2}(R2−R1) | β_{iat2,R2}>
```

Each term consumes the complex DMR block element-wise by
`(iat1,iat2,R2−R1)` (for nspin=4 the block is provided by `cal_DMR_full`
without dropping the imaginary part). This path **does not** satisfy the
same-R closed-trace protection of Sec. 5; its correctness relies on:

1. `cal_DMR_full` using `e^{-ikR}` (mutually inverse with `folding_HR`'s
   `e^{+ikR}`);
2. the DMR block taken at `(iat1,iat2,R2−R1)` matching the index order of
   the β projection;
3. numerically locked by the SOC finite-difference force integration test
   (`tests/integrate/240_NO_KP_15_SO_FD`).

If a half-set storage optimization is introduced in the future (storing only
`(iat1,iat2,R)` without `(iat2,iat1,−R)`), the Sec. 5 protection argument
and the `force_factor = 1.0` convention both break; the factor must be
restored to 2 and the T2/T3/T8 guards re-run.
