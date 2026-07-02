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
    for (int j = 0; j < n_band_; ++j)
    {
        Real nrm2 = 0;
        for (int ig = 0; ig < n_dim_; ++ig)
            nrm2 += static_cast<Real>(std::norm(residual[idx(ig, j, ld_psi_)]));
        const Real rnrm = std::sqrt(std::max(nrm2, static_cast<Real>(0)));
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
        for (int j = 0; j < lcols; ++j) {
            Real sn2 = 0;
            for (int ig = 0; ig < n_dim_; ++ig)
                sn2 += std::real(std::conj(x[idx(ig, j, ld_psi_)])
                                 * sx[idx(ig, j, ld_psi_)]);
            Real sn = std::sqrt(std::max(sn2, static_cast<Real>(1e-30)));
            // Only scale if the norm is non-negligible; a near-zero
            // column is a converged band whose contribution is harmless.
            if (sn > static_cast<Real>(1e-15)) {
                Real inv = static_cast<Real>(1) / sn;
                scale[j] = inv;
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

    auto fill_sym = [&](const std::vector<T>& a, const std::vector<T>& b,
                        int r0, int c0, std::vector<T>& mat)
    {
        std::vector<T> g;
        gram(a.data(), b.data(), l, l, g, l);
        for (int j = 0; j < l; ++j)
            for (int i = 0; i < l; ++i)
            {
                mat[(r0 + i) + (c0 + j) * dim] = g[i + j * l];
                mat[(c0 + j) + (r0 + i) * dim] = std::conj(g[i + j * l]);
            }
    };

    fill_sym(psi_l, hpsi_l, 0,   0,   subspace.k);
    fill_sym(psi_l, spsi_l, 0,   0,   subspace.m);
    fill_sym(w_l,   hw_l,   l,   l,   subspace.k);
    fill_sym(w_l,   sw_l,   l,   l,   subspace.m);
    fill_sym(psi_l, hw_l,   0,   l,   subspace.k);
    fill_sym(psi_l, sw_l,   0,   l,   subspace.m);

    if (use_p)
    {
        fill_sym(p_l, hp_l, 2*l, 2*l, subspace.k);
        fill_sym(p_l, sp_l, 2*l, 2*l, subspace.m);
        fill_sym(psi_l, hp_l, 0,   2*l, subspace.k);
        fill_sym(psi_l, sp_l, 0,   2*l, subspace.m);
        fill_sym(w_l,   hp_l, l,   2*l, subspace.k);
        fill_sym(w_l,   sp_l, l,   2*l, subspace.m);
    }
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
        subspace.k[i + i * dim] = T(1);
    std::fill(subspace.eval.begin(), subspace.eval.end(), static_cast<Real>(0));
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

    for (int j = 0; j < l; ++j)
    {
        for (int i = 0; i < l; ++i)
        {
            const T cpsi = eigvec[i + j * dim];
            const T cw   = eigvec[(l + i) + j * dim] * subspace.w_scale[i];

            for (int ig = 0; ig < n_dim_; ++ig)
            {
                psi_new[idx(ig, j, ld_psi_)]  += psi_l[idx(ig, i, ld_psi_)] * cpsi
                                                + w_l[ idx(ig, i, ld_psi_)] * cw;
                spsi_new[idx(ig, j, ld_psi_)] += spsi_l[idx(ig, i, ld_psi_)] * cpsi
                                                + sw_l[ idx(ig, i, ld_psi_)] * cw;
                hpsi_new[idx(ig, j, ld_psi_)] += hpsi_l[idx(ig, i, ld_psi_)] * cpsi
                                                + hw_l[ idx(ig, i, ld_psi_)] * cw;
                p_new[idx(ig, j, ld_psi_)]    += w_l[ idx(ig, i, ld_psi_)] * cw;
                sp_new[idx(ig, j, ld_psi_)]   += sw_l[ idx(ig, i, ld_psi_)] * cw;
                hp_new[idx(ig, j, ld_psi_)]   += hw_l[ idx(ig, i, ld_psi_)] * cw;
            }

            if (use_p)
            {
                const T cp = eigvec[(2*l + i) + j * dim] * subspace.p_scale[i];
                for (int ig = 0; ig < n_dim_; ++ig)
                {
                    psi_new[idx(ig, j, ld_psi_)]  += p_l[ idx(ig, i, ld_psi_)] * cp;
                    spsi_new[idx(ig, j, ld_psi_)] += sp_l[idx(ig, i, ld_psi_)] * cp;
                    hpsi_new[idx(ig, j, ld_psi_)] += hp_l[idx(ig, i, ld_psi_)] * cp;
                    p_new[idx(ig, j, ld_psi_)]    += p_l[ idx(ig, i, ld_psi_)] * cp;
                    sp_new[idx(ig, j, ld_psi_)]   += sp_l[idx(ig, i, ld_psi_)] * cp;
                    hp_new[idx(ig, j, ld_psi_)]   += hp_l[idx(ig, i, ld_psi_)] * cp;
                }
            }
        }
    }

    scatter_cols(psi, cols, psi_new);
    scatter_cols(spsi_.data(), cols, spsi_new);
    scatter_cols(hpsi_.data(), cols, hpsi_new);
    scatter_cols(p_.data(), cols, p_new);
    scatter_cols(sp_.data(), cols, sp_new);
    scatter_cols(hp_.data(), cols, hp_new);
}

} // namespace hsolver
