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
// Build K = V^H H V and M = V^H S V where V = [psi, w, p]
// ---------------------------------------------------------------------------
template <typename T, typename Device>
void DiagoPPCG<T, Device>::build_small_subspace(
    const T* psi,
    const std::vector<int>& cols,
    bool use_p,
    SmallSubspace& subspace) const
{
    const int l = static_cast<int>(cols.size());
    const int nblk = use_p ? 3 : 2;
    const int dim = nblk * l;
    subspace.k.assign(dim * dim, T(0));
    subspace.m.assign(dim * dim, T(0));
    subspace.eval.assign(dim, static_cast<Real>(0));
    subspace.w_scale.assign(l, static_cast<Real>(1));
    subspace.p_scale.assign(l, static_cast<Real>(1));

    std::vector<T> psi_l, spsi_l, hpsi_l;
    std::vector<T> w_l, sw_l, hw_l;
    std::vector<T> p_l, sp_l, hp_l;
    copy_cols(psi, cols, psi_l);
    copy_cols(spsi_.data(), cols, spsi_l);
    copy_cols(hpsi_.data(), cols, hpsi_l);
    copy_cols(w_.data(), cols, w_l);
    copy_cols(sw_.data(), cols, sw_l);
    copy_cols(hw_.data(), cols, hw_l);
    if (use_p)
    {
        copy_cols(p_.data(), cols, p_l);
        copy_cols(sp_.data(), cols, sp_l);
        copy_cols(hp_.data(), cols, hp_l);
    }

    // ---------------------------------------------------------------------------
    // Normalize w and p columns to unit S-norm for numerical stability.
    //
    // The [w, p] block of the Gram matrix M has entries O(||w||²) which
    // become tiny when residuals are small, making M nearly singular and
    // causing sygvd to produce garbage eigenvectors.
    //
    // Scaling to unit S-norm keeps M well-conditioned (diagonal ~1) without
    // changing the subspace.  The Ritz values are identical and the Ritz
    // vector coefficients in update_one_block automatically compensate.
    // ---------------------------------------------------------------------------
    auto scale_to_unit_snorm = [this](std::vector<T>& x, std::vector<T>& sx,
                                       std::vector<T>& hx, int lcols,
                                       std::vector<Real>& scale) {
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
                scale[j] = inv;
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
    scale_to_unit_snorm(w_l, sw_l, hw_l, l, subspace.w_scale);
    if (use_p)
        scale_to_unit_snorm(p_l, sp_l, hp_l, l, subspace.p_scale);

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

    std::vector<T> basis(ld_psi_ * dim, T(0));
    std::vector<T> hbasis(ld_psi_ * dim, T(0));
    std::vector<T> sbasis(ld_psi_ * dim, T(0));
    copy_block(psi_l, 0, basis);
    copy_block(hpsi_l, 0, hbasis);
    copy_block(spsi_l, 0, sbasis);
    copy_block(w_l, l, basis);
    copy_block(hw_l, l, hbasis);
    copy_block(sw_l, l, sbasis);
    if (use_p)
    {
        copy_block(p_l, 2 * l, basis);
        copy_block(hp_l, 2 * l, hbasis);
        copy_block(sp_l, 2 * l, sbasis);
    }

    gram(basis.data(), hbasis.data(), dim, dim, subspace.k, dim);
    gram(basis.data(), sbasis.data(), dim, dim, subspace.m, dim);
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
    bool use_p,
    const SmallSubspace& subspace)
{
    const int dim = (use_p ? 3 : 2) * l;
    const T* eigvec = subspace.k.data();

    std::vector<T> psi_l, spsi_l, hpsi_l;
    std::vector<T> w_l, sw_l, hw_l;
    std::vector<T> p_l, sp_l, hp_l;
    copy_cols(psi, cols, psi_l);
    copy_cols(spsi_.data(), cols, spsi_l);
    copy_cols(hpsi_.data(), cols, hpsi_l);
    copy_cols(w_.data(), cols, w_l);
    copy_cols(sw_.data(), cols, sw_l);
    copy_cols(hw_.data(), cols, hw_l);
    if (use_p)
    {
        copy_cols(p_.data(), cols, p_l);
        copy_cols(sp_.data(), cols, sp_l);
        copy_cols(hp_.data(), cols, hp_l);
    }

    std::vector<T> psi_new(ld_psi_ * l, T(0));
    std::vector<T> spsi_new(ld_psi_ * l, T(0));
    std::vector<T> hpsi_new(ld_psi_ * l, T(0));
    std::vector<T> p_new(ld_psi_ * l, T(0));
    std::vector<T> sp_new(ld_psi_ * l, T(0));
    std::vector<T> hp_new(ld_psi_ * l, T(0));

    std::vector<T> coeff_state(dim * l, T(0));
    std::vector<T> coeff_dir(dim * l, T(0));
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (l * l > 4096)
#endif
    for (int j = 0; j < l; ++j)
    {
        for (int i = 0; i < l; ++i)
        {
            coeff_state[i + j * dim] = eigvec[i + j * dim];
            const T cw = eigvec[(l + i) + j * dim] * subspace.w_scale[i];
            coeff_state[(l + i) + j * dim] = cw;
            coeff_dir[(l + i) + j * dim] = cw;
            if (use_p)
            {
                const T cp = eigvec[(2*l + i) + j * dim] * subspace.p_scale[i];
                coeff_state[(2*l + i) + j * dim] = cp;
                coeff_dir[(2*l + i) + j * dim] = cp;
            }
        }
    }

    auto fill_basis = [&](const std::vector<T>& a,
                          const std::vector<T>& b,
                          const std::vector<T>& c,
                          std::vector<T>& basis)
    {
        basis.assign(ld_psi_ * dim, T(0));
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
            if (use_p)
            {
                std::copy(c.begin() + j * ld_psi_,
                          c.begin() + (j + 1) * ld_psi_,
                          basis.begin() + (2 * l + j) * ld_psi_);
            }
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

    std::vector<T> psi_basis;
    std::vector<T> spsi_basis;
    std::vector<T> hpsi_basis;
    fill_basis(psi_l, w_l, p_l, psi_basis);
    fill_basis(spsi_l, sw_l, sp_l, spsi_basis);
    fill_basis(hpsi_l, hw_l, hp_l, hpsi_basis);

    combine(psi_basis, coeff_state, psi_new);
    combine(spsi_basis, coeff_state, spsi_new);
    combine(hpsi_basis, coeff_state, hpsi_new);
    combine(psi_basis, coeff_dir, p_new);
    combine(spsi_basis, coeff_dir, sp_new);
    combine(hpsi_basis, coeff_dir, hp_new);

    scatter_cols(psi, cols, psi_new);
    scatter_cols(spsi_.data(), cols, spsi_new);
    scatter_cols(hpsi_.data(), cols, hpsi_new);
    scatter_cols(p_.data(), cols, p_new);
    scatter_cols(sp_.data(), cols, sp_new);
    scatter_cols(hp_.data(), cols, hp_new);
}

} // namespace hsolver
