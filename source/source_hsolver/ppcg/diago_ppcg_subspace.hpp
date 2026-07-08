#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace hsolver {

//==============================================================================
// BLOCK_SUBSPACE STRATEGY
//==============================================================================

// ---------------------------------------------------------------------------
// Lock converged eigenpairs: columns with residual below threshold
// ---------------------------------------------------------------------------
template <typename T, typename Device>
void DiagoPPCG<T, Device>::lock_epairs(
    const std::vector<T>& residual,
    const std::vector<double>& ethr_band,
    std::vector<int>& active_cols) const
{
    active_cols.clear();
    active_cols.reserve(n_band_);
    std::vector<double> nrm2_all(n_band_, 0.0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (n_dim_ * n_band_ > 4096)
#endif
    for (int j = 0; j < n_band_; ++j)
    {
        double nrm2 = 0.0;
        for (int ig = 0; ig < n_dim_; ++ig)
            nrm2 += static_cast<double>(std::norm(residual[idx(ig, j, ld_psi_)]));
        nrm2_all[j] = nrm2;
    }
    reduce_pool_if_mpi_ready(nrm2_all.data(), n_band_);
    for (int j = 0; j < n_band_; ++j)
    {
        const Real rnrm = std::sqrt(std::max(static_cast<Real>(nrm2_all[j]),
                                             static_cast<Real>(0)));
        const Real thr = std::max(static_cast<Real>(ethr_band[j]), diag_thr_);
        if (rnrm > thr)
            active_cols.push_back(j);
    }
}

// ---------------------------------------------------------------------------
// Build K = V^H H V and M = V^H S V where V = [psi, w]
// ---------------------------------------------------------------------------
template <typename T, typename Device>
void DiagoPPCG<T, Device>::build_small_subspace(
    const T* psi,
    const std::vector<int>& cols,
    SmallSubspace& subspace) const
{
    const int l = static_cast<int>(cols.size());
    const int dim = 2 * l;
    subspace.k.resize(dim * dim);
    subspace.m.resize(dim * dim);
    subspace.eval.resize(dim);

    copy_cols(psi, cols, subspace.psi_l);
    copy_cols(spsi_.data(), cols, subspace.spsi_l);
    copy_cols(hpsi_.data(), cols, subspace.hpsi_l);
    copy_cols(w_.data(), cols, subspace.w_l);
    copy_cols(sw_.data(), cols, subspace.sw_l);
    copy_cols(hw_.data(), cols, subspace.hw_l);

    // ---------------------------------------------------------------------------
    // Normalize w columns to unit S-norm for numerical stability.
    //
    // The w block of the Gram matrix M has entries O(||w||^2) which become
    // tiny when residuals are small, making M nearly singular and causing
    // sygvd to produce garbage eigenvectors.
    //
    // Scaling to unit S-norm keeps M well-conditioned (diagonal ~1) without
    // changing the subspace. The same scaled basis is reused in update_one_block.
    // ---------------------------------------------------------------------------
    auto scale_to_unit_snorm = [this](std::vector<T>& x,
                                      std::vector<T>& sx,
                                      std::vector<T>& hx,
                                      int lcols) {
        std::vector<double> sn2_all(lcols, 0.0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (n_dim_ * lcols > 4096)
#endif
        for (int j = 0; j < lcols; ++j) {
            double sn2 = 0.0;
            for (int ig = 0; ig < n_dim_; ++ig)
                sn2 += static_cast<double>(std::real(std::conj(x[idx(ig, j, ld_psi_)])
                                                     * sx[idx(ig, j, ld_psi_)]));
            sn2_all[j] = sn2;
        }
        reduce_pool_if_mpi_ready(sn2_all.data(), lcols);
        for (int j = 0; j < lcols; ++j) {
            Real sn = std::sqrt(std::max(static_cast<Real>(sn2_all[j]),
                                         static_cast<Real>(1e-30)));
            // Only scale if the norm is non-negligible; a near-zero
            // column is a converged band whose contribution is harmless.
            if (sn > static_cast<Real>(1e-15)) {
                Real inv = static_cast<Real>(1) / sn;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (n_dim_ > 4096)
#endif
                for (int ig = 0; ig < n_dim_; ++ig) {
                    x[ idx(ig, j, ld_psi_)]  *= inv;
                    sx[idx(ig, j, ld_psi_)] *= inv;
                    hx[idx(ig, j, ld_psi_)] *= inv;
                }
            }
        }
    };
    scale_to_unit_snorm(subspace.w_l,
                        subspace.sw_l,
                        subspace.hw_l,
                        l);

    auto copy_block = [&](const std::vector<T>& src,
                          const int col0,
                          std::vector<T>& dst)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (ld_psi_ * l > 4096)
#endif
        for (int j = 0; j < l; ++j)
            std::copy(src.begin() + j * ld_psi_,
                      src.begin() + (j + 1) * ld_psi_,
                      dst.begin() + (col0 + j) * ld_psi_);
    };

    auto hermitize = [&](std::vector<T>& mat)
    {
        for (int j = 0; j < dim; ++j)
        {
            mat[j + j * dim] = T(std::real(mat[j + j * dim]), 0);
            for (int i = j + 1; i < dim; ++i)
            {
                const T avg = (mat[i + j * dim] + std::conj(mat[j + i * dim]))
                            * static_cast<Real>(0.5);
                mat[i + j * dim] = avg;
                mat[j + i * dim] = std::conj(avg);
            }
        }
    };

    subspace.basis.resize(ld_psi_ * dim);
    subspace.hbasis.resize(ld_psi_ * dim);
    subspace.sbasis.resize(ld_psi_ * dim);
    copy_block(subspace.psi_l, 0, subspace.basis);
    copy_block(subspace.hpsi_l, 0, subspace.hbasis);
    copy_block(subspace.spsi_l, 0, subspace.sbasis);
    copy_block(subspace.w_l, l, subspace.basis);
    copy_block(subspace.hw_l, l, subspace.hbasis);
    copy_block(subspace.sw_l, l, subspace.sbasis);

    gram(subspace.basis.data(), subspace.hbasis.data(), dim, dim, subspace.k, dim);
    gram(subspace.basis.data(), subspace.sbasis.data(), dim, dim, subspace.m, dim);
    hermitize(subspace.k);
    hermitize(subspace.m);
}

// ---------------------------------------------------------------------------
// Solve K v = λ M v (small generalized eigenvalue problem)
// ---------------------------------------------------------------------------
template <typename T, typename Device>
void DiagoPPCG<T, Device>::solve_small_generalized(
    int dim, SmallSubspace& subspace) const
{
    // Try with increasing diagonal shifts; fall back to identity (no update)
    // if the subspace is too ill-conditioned.
    // Save originals; sygvd modifies both matrices in-place before it may
    // fail.
    const std::vector<T> k0 = subspace.k;
    const std::vector<T> m0 = subspace.m;
    const Real shifts[] = {static_cast<Real>(0),
                           static_cast<Real>(1e-10),
                           static_cast<Real>(1e-8),
                           static_cast<Real>(1e-6)};
    for (const Real shift : shifts)
    {
        subspace.k = k0;
        subspace.m = m0;
        for (int i = 0; i < dim; ++i)
            subspace.m[i + i * dim] += T(shift);

        try
        {
            HermitianLapack<T>::sygvd(dim, subspace.k.data(),
                                      subspace.m.data(),
                                      subspace.eval.data());
            return;
        }
        catch (const std::runtime_error&)
        {
            // Try the next diagonal shift.
        }
    }
    // All attempts failed — set eigenvectors to identity (no update).
    std::fill(subspace.k.begin(), subspace.k.end(), T(0));
    for (int i = 0; i < dim; ++i)
    {
        subspace.k[i + i * dim] = T(1);
        subspace.eval[i] = static_cast<Real>(std::real(k0[i + i * dim]))
                         / std::max(static_cast<Real>(std::real(m0[i + i * dim])),
                                    static_cast<Real>(1e-30));
    }
}

// ---------------------------------------------------------------------------
// Update wavefunctions from small subspace eigenvectors
// ---------------------------------------------------------------------------
template <typename T, typename Device>
void DiagoPPCG<T, Device>::update_one_block(
    T* psi,
    const std::vector<int>& cols,
    int l,
    SmallSubspace& subspace)
{
    const int dim = 2 * l;
    const T* eigvec = subspace.k.data();

    subspace.psi_new.assign(ld_psi_ * l, T(0));
    subspace.spsi_new.assign(ld_psi_ * l, T(0));
    subspace.hpsi_new.assign(ld_psi_ * l, T(0));

    subspace.coeff_state.resize(dim * l);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (l * l > 4096)
#endif
    for (int j = 0; j < l; ++j)
    {
        for (int i = 0; i < l; ++i)
        {
            subspace.coeff_state[i + j * dim] = eigvec[i + j * dim];
            subspace.coeff_state[(l + i) + j * dim] = eigvec[(l + i) + j * dim];
        }
    }

    auto fill_basis = [&](const std::vector<T>& a,
                          const std::vector<T>& b,
                          std::vector<T>& basis)
    {
        basis.resize(ld_psi_ * dim);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (ld_psi_ * l > 4096)
#endif
        for (int j = 0; j < l; ++j)
        {
            std::copy(a.begin() + j * ld_psi_,
                      a.begin() + (j + 1) * ld_psi_,
                      basis.begin() + j * ld_psi_);
            std::copy(b.begin() + j * ld_psi_,
                      b.begin() + (j + 1) * ld_psi_,
                      basis.begin() + (l + j) * ld_psi_);
        }
    };

    auto combine = [&](const std::vector<T>& basis,
                       const std::vector<T>& coeff,
                       std::vector<T>& out)
    {
        const T one = T(1);
        const T zero = T(0);
        ModuleBase::gemm_op<T, Device>()('N',
                                         'N',
                                         n_dim_,
                                         l,
                                         dim,
                                         &one,
                                         basis.data(),
                                         ld_psi_,
                                         coeff.data(),
                                         dim,
                                         &zero,
                                         out.data(),
                                         ld_psi_);
    };

    fill_basis(subspace.psi_l, subspace.w_l, subspace.basis);
    fill_basis(subspace.spsi_l, subspace.sw_l, subspace.sbasis);
    fill_basis(subspace.hpsi_l, subspace.hw_l, subspace.hbasis);

    combine(subspace.basis, subspace.coeff_state, subspace.psi_new);
    combine(subspace.sbasis, subspace.coeff_state, subspace.spsi_new);
    combine(subspace.hbasis, subspace.coeff_state, subspace.hpsi_new);

    scatter_cols(psi, cols, subspace.psi_new);
    scatter_cols(spsi_.data(), cols, subspace.spsi_new);
    scatter_cols(hpsi_.data(), cols, subspace.hpsi_new);
}

} // namespace hsolver
