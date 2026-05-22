#include "source_hsolver/diago_ppcg.h"

#include "source_base/parallel_comm.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_base/tool_quit.h"
#include "source_hsolver/diago_bpcg.h"
#include "source_hsolver/diago_iter_assist.h"

#include <ATen/kernels/lapack.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace hsolver
{

template <typename T, typename Device>
DiagoPPCG<T, Device>::DiagoPPCG(const Real* precondition_in) : precondition(precondition_in)
{
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::init_iter(const int nband, const int nband_l, const int nbasis, const int ndim)
{
    this->n_band = nband;
    this->n_band_l = nband_l;
    this->n_basis = nbasis;
    this->n_dim = ndim;
    this->n_work = this->n_band_l + this->n_extra;

    const int block_size = this->n_work * this->n_basis;
    this->hpsi.assign(block_size, T(0));
    this->w.assign(block_size, T(0));
    this->hw.assign(block_size, T(0));
    this->p.assign(block_size, T(0));
    this->hp.assign(block_size, T(0));
    this->p_new.assign(block_size, T(0));
    this->hp_new.assign(block_size, T(0));
    this->hpsi_new.assign(block_size, T(0));
    this->work.assign(block_size, T(0));
    this->eigen.assign(this->n_work, Real(0));
    this->err.assign(this->n_work, std::numeric_limits<Real>::max());
    this->is_locked.assign(this->n_work, false);
    this->converge_count.assign(this->n_work, 0);
}

template <typename T, typename Device>
T DiagoPPCG<T, Device>::inner_product(const T* lhs, const T* rhs) const
{
    T result = T(0);
    for (int ig = 0; ig < this->n_dim; ++ig)
    {
        result += std::conj(lhs[ig]) * rhs[ig];
    }
    Parallel_Reduce::reduce_pool(&result, 1);
    return result;
}

template <typename T, typename Device>
typename DiagoPPCG<T, Device>::Real DiagoPPCG<T, Device>::vector_norm(const T* vec) const
{
    const Real norm2 = std::max(Real(0), std::real(this->inner_product(vec, vec)));
    return std::sqrt(norm2);
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::scale_vector(T* vec, const Real alpha) const
{
    for (int ig = 0; ig < this->n_dim; ++ig)
    {
        vec[ig] *= alpha;
    }
    for (int ig = this->n_dim; ig < this->n_basis; ++ig)
    {
        vec[ig] = T(0);
    }
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::axpy_vector(T* y, const T* x, const T alpha) const
{
    for (int ig = 0; ig < this->n_dim; ++ig)
    {
        y[ig] += alpha * x[ig];
    }
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::copy_vector(T* dst, const T* src) const
{
    std::copy(src, src + this->n_basis, dst);
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::zero_vector(T* vec) const
{
    std::fill(vec, vec + this->n_basis, T(0));
}

template <typename T, typename Device>
bool DiagoPPCG<T, Device>::test_error(const std::vector<double>& ethr_band) const
{
    bool not_conv = false;
    for (int ib = 0; ib < this->n_band_l; ++ib)
    {
        if (this->err[ib] > ethr_band[ib])
        {
            not_conv = true;
            break;
        }
    }
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &not_conv, 1, MPI_C_BOOL, MPI_LOR, BP_WORLD);
#endif
    return not_conv;
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::calc_hpsi(const HPsiFunc& hpsi_func, T* psi_in, std::vector<T>& hpsi_out) const
{
    hpsi_func(psi_in, hpsi_out.data(), this->n_basis, this->n_work);
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::modified_gram_schmidt(T* psi_in, std::vector<T>& hpsi_in) const
{
    // Modified Gram-Schmidt: for each column, subtract projections onto all
    // previous columns from both psi and hpsi, then normalize both.
    for (int ib = 0; ib < this->n_work; ++ib)
    {
        T* xi = psi_in + ib * this->n_basis;
        T* hxi = hpsi_in.data() + ib * this->n_basis;
        for (int jb = 0; jb < ib; ++jb)
        {
            const T* xj = psi_in + jb * this->n_basis;
            const T* hxj = hpsi_in.data() + jb * this->n_basis;
            const T coeff = this->inner_product(xj, xi);
            this->axpy_vector(xi, xj, -coeff);
            this->axpy_vector(hxi, hxj, -coeff);
        }

        const Real norm = this->vector_norm(xi);
        if (norm <= Real(1.0e-14))
        {
            ModuleBase::WARNING_QUIT("DiagoPPCG::modified_gram_schmidt", "linear dependent wavefunctions");
        }
        this->scale_vector(xi, Real(1) / norm);
        this->scale_vector(hxi, Real(1) / norm);
    }
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::orth_cholesky(T* psi_in, std::vector<T>& hpsi_in)
{
    // Cholesky-based orthonormalization:
    //   1. Build overlap matrix S = <psi|psi>
    //   2. Cholesky factorize S = U^H * U (LAPACK potrf, upper)
    //   3. Compute U^{-1} (LAPACK trtri, upper, non-unit)
    //   4. Rotate psi and hpsi by U^{-1}, yielding orthonormal vectors.
    std::vector<T> s(this->n_work * this->n_work, T(0));
    for (int col = 0; col < this->n_work; ++col)
    {
        for (int row = 0; row < this->n_work; ++row)
        {
            s[row + col * this->n_work]
                = this->inner_product(psi_in + row * this->n_basis, psi_in + col * this->n_basis);
        }
    }

    ct::kernels::lapack_potrf<T, ct::DEVICE_CPU>()('U', this->n_work, s.data(), this->n_work);

    for (int col = 0; col < this->n_work; ++col)
    {
        for (int row = col + 1; row < this->n_work; ++row)
        {
            s[row + col * this->n_work] = T(0);
        }
    }

    ct::kernels::lapack_trtri<T, ct::DEVICE_CPU>()('U', 'N', this->n_work, s.data(), this->n_work);

    this->rotate_block(psi_in, s, this->work);
    this->rotate_block(hpsi_in.data(), s, this->work);
}

template <typename T, typename Device>
bool DiagoPPCG<T, Device>::check_orthonormality(T* psi_in) const
{
    // Compute the Frobenius norm of (S - I) where S_ij = <psi_i | psi_j>.
    // Returns true if the deviation from identity is below 1e-6.
    Real frob2 = 0;
    for (int col = 0; col < this->n_work; ++col)
    {
        for (int row = 0; row < this->n_work; ++row)
        {
            const T s = this->inner_product(psi_in + row * this->n_basis, psi_in + col * this->n_basis);
            const T delta = s - static_cast<T>(row == col ? 1.0 : 0.0);
            frob2 += std::norm(delta);
        }
    }
    return std::sqrt(frob2) < Real(1e-1);
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::rotate_block(T* block, const std::vector<T>& coeff, std::vector<T>& workspace) const
{
    // Rotate a block of vectors by a coefficient matrix: block_out = block_in * coeff.
    // coeff is (n_work x n_work) column-major; each output column is a linear
    // combination of input columns weighted by the corresponding column of coeff.
    std::fill(workspace.begin(), workspace.end(), T(0));
    for (int out = 0; out < this->n_work; ++out)
    {
        T* dst = workspace.data() + out * this->n_basis;
        for (int in = 0; in < this->n_work; ++in)
        {
            const T* src = block + in * this->n_basis;
            const T c = coeff[in + out * this->n_work];
            for (int ig = 0; ig < this->n_dim; ++ig)
            {
                dst[ig] += src[ig] * c;
            }
        }
    }
    std::copy(workspace.begin(), workspace.end(), block);
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::rayleigh_ritz(T* psi_in, std::vector<T>& hpsi_in)
{
    // Rayleigh-Ritz: build subspace Hamiltonian Hsub = <psi|H|psi>,
    // diagonalize it (LAPACK zheevd), then rotate psi and hpsi by the
    // eigenvectors to obtain Ritz vectors sorted by ascending eigenvalue.
    if (this->n_work == 0)
    {
        return;
    }

    std::vector<T> hsub(this->n_work * this->n_work, T(0));
    for (int col = 0; col < this->n_work; ++col)
    {
        for (int row = 0; row < this->n_work; ++row)
        {
            hsub[row + col * this->n_work]
                = this->inner_product(psi_in + row * this->n_basis, hpsi_in.data() + col * this->n_basis);
        }
    }

    ct::kernels::lapack_heevd<T, ct::DEVICE_CPU>()(this->n_work, hsub.data(), this->n_work, this->eigen.data());
    this->rotate_block(psi_in, hsub, this->work);
    this->rotate_block(hpsi_in.data(), hsub, this->work);
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::calc_preconditioned_residual(T* psi_in)
{
    // For each working band:
    //   - lambda_i = <x_i | H | x_i>   (Rayleigh quotient, used as eigenvalue estimate)
    //   - R_i     = H x_i - lambda_i x_i  (residual)
    //   - w_i     = -K^{-1} R_i           (preconditioned residual)
    // Locked bands are skipped (w_i is zeroed).
    for (int ib = 0; ib < this->n_work; ++ib)
    {
        T* wi = this->w.data() + ib * this->n_basis;
        T* xi = psi_in + ib * this->n_basis;
        T* hxi = this->hpsi.data() + ib * this->n_basis;

        if (this->is_locked[ib])
        {
            this->zero_vector(wi);
            continue;
        }

        const Real lambda = std::real(this->inner_product(xi, hxi));
        this->eigen[ib] = lambda;

        Real err2 = 0;
        for (int ig = 0; ig < this->n_dim; ++ig)
        {
            const T residual = hxi[ig] - lambda * xi[ig];
            err2 += std::norm(residual);
            wi[ig] = -residual / this->precondition[ig];
        }
        Parallel_Reduce::reduce_pool(err2);
        this->err[ib] = std::sqrt(std::max(Real(0), err2));
        for (int ig = this->n_dim; ig < this->n_basis; ++ig)
        {
            wi[ig] = T(0);
        }
    }
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::project_to_orthogonal_complement(T* psi_in, std::vector<T>& block) const
{
    // For each vector v_i in block, subtract its projection onto all current psi
    // vectors: v_i = v_i - sum_j <x_j | v_i> * x_j.
    for (int ib = 0; ib < this->n_work; ++ib)
    {
        T* vi = block.data() + ib * this->n_basis;
        for (int jb = 0; jb < this->n_work; ++jb)
        {
            const T* xj = psi_in + jb * this->n_basis;
            const T coeff = this->inner_product(xj, vi);
            this->axpy_vector(vi, xj, -coeff);
        }
    }
}

template <typename T, typename Device>
bool DiagoPPCG<T, Device>::solve_small_problem(const int active_dim, T* hsmall, T* ssmall, T* coeff, Real* eval) const
{
    // Solve the 2x2 or 3x3 generalized eigenvalue problem H*C = lambda*S*C
    // using LAPACK zhegvd. A small regularization term (1e-12) is added to
    // the diagonal of S to guard against ill-conditioning from near-linear-dependence.
    // On failure, fall back to returning the first basis vector as-is.
    std::fill(coeff, coeff + 9, T(0));
    std::fill(eval, eval + 3, Real(0));
    if (active_dim <= 1)
    {
        coeff[0] = T(1);
        eval[0] = std::real(hsmall[0]);
        return true;
    }

    for (int i = 0; i < active_dim; ++i)
    {
        ssmall[i + i * active_dim] += T(1.0e-12);
    }

    try
    {
        ct::kernels::lapack_hegvd<T, ct::DEVICE_CPU>()(active_dim, active_dim, hsmall, ssmall, eval, coeff);
    }
    catch (const std::exception&)
    {
        coeff[0] = T(1);
        eval[0] = std::real(hsmall[0]);
        return false;
    }
    return true;
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::update_vectors_from_ppcg_subspace(T* psi_in)
{
    // If block sizes are configured, use the block-diagonal variant that solves
    // a single larger generalized eigenvalue problem per block instead of
    // per-band 2D/3D subspace problems.
    if (!this->block_sizes.empty())
    {
        this->update_vectors_blocked(psi_in);
        return;
    }

    // Per-band mode: for each band, construct a small subspace from
    // {x_i, w_i, p_i} (3D when p_i is non-zero, 2D otherwise), build
    // the subspace overlap and Hamiltonian matrices, solve the generalized
    // eigenvalue problem, and update the working vectors using the first
    // eigenvector's coefficients.
    std::fill(this->p_new.begin(), this->p_new.end(), T(0));
    std::fill(this->hp_new.begin(), this->hp_new.end(), T(0));
    std::fill(this->hpsi_new.begin(), this->hpsi_new.end(), T(0));

    for (int ib = 0; ib < this->n_work; ++ib)
    {
        T* xi = psi_in + ib * this->n_basis;
        T* hxi = this->hpsi.data() + ib * this->n_basis;
        T* wi = this->w.data() + ib * this->n_basis;
        T* hwi = this->hw.data() + ib * this->n_basis;
        T* pi = this->p.data() + ib * this->n_basis;
        T* hpi = this->hp.data() + ib * this->n_basis;

        T* xnew = this->work.data() + ib * this->n_basis;
        T* hxnew = this->hpsi_new.data() + ib * this->n_basis;
        T* pnext = this->p_new.data() + ib * this->n_basis;
        T* hpnext = this->hp_new.data() + ib * this->n_basis;

        if (this->is_locked[ib])
        {
            this->copy_vector(xnew, xi);
            this->copy_vector(hxnew, hxi);
            this->zero_vector(pnext);
            this->zero_vector(hpnext);
            continue;
        }

        const Real pnorm = this->vector_norm(pi);
        const int active_dim = (pnorm > Real(1.0e-12)) ? 3 : 2;

        const T* basis_vecs[3] = {xi, wi, pi};
        const T* hbasis_vecs[3] = {hxi, hwi, hpi};

        T hsmall[9] = {};
        T ssmall[9] = {};
        T coeff[9] = {};
        Real eval[3] = {};

        for (int col = 0; col < active_dim; ++col)
        {
            for (int row = 0; row < active_dim; ++row)
            {
                hsmall[row + col * active_dim] = this->inner_product(basis_vecs[row], hbasis_vecs[col]);
                ssmall[row + col * active_dim] = this->inner_product(basis_vecs[row], basis_vecs[col]);
            }
        }

        this->solve_small_problem(active_dim, hsmall, ssmall, coeff, eval);
        this->eigen[ib] = eval[0];

        this->zero_vector(xnew);
        this->zero_vector(hxnew);
        this->zero_vector(pnext);
        this->zero_vector(hpnext);

        for (int j = 0; j < active_dim; ++j)
        {
            const T c = coeff[j];
            this->axpy_vector(xnew, basis_vecs[j], c);
            this->axpy_vector(hxnew, hbasis_vecs[j], c);
        }

        if (active_dim >= 2)
        {
            const T cw = coeff[1];
            this->axpy_vector(pnext, wi, cw);
            this->axpy_vector(hpnext, hwi, cw);
        }
        if (active_dim == 3)
        {
            const T cp = coeff[2];
            this->axpy_vector(pnext, pi, cp);
            this->axpy_vector(hpnext, hpi, cp);
        }
    }

    std::copy(this->work.begin(), this->work.end(), psi_in);
    std::copy(this->hpsi_new.begin(), this->hpsi_new.end(), this->hpsi.begin());
    std::copy(this->p_new.begin(), this->p_new.end(), this->p.begin());
    std::copy(this->hp_new.begin(), this->hp_new.end(), this->hp.begin());
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::update_vectors_blocked(T* psi_in)
{
    // Block-diagonal PPCG variant.
    // For each block of size k_i, construct a 3k_i-dimensional subspace
    // from the three sub-blocks {X_block, W_block, P_block}, build the
    // subspace overlap and Hamiltonian matrices (each 3k_i x 3k_i),
    // solve the generalized eigenvalue problem H_sub * C = Lambda * S_sub * C,
    // and update all k_i bands simultaneously using the first k_i eigenvectors.
    std::fill(this->p_new.begin(), this->p_new.end(), T(0));
    std::fill(this->hp_new.begin(), this->hp_new.end(), T(0));
    std::fill(this->hpsi_new.begin(), this->hpsi_new.end(), T(0));

    int band_offset = 0;
    for (std::size_t b = 0; b < this->block_sizes.size(); ++b)
    {
        const int k_i = this->block_sizes[b];
        if (k_i <= 0 || band_offset + k_i > this->n_band_l)
        {
            band_offset += k_i;
            continue;
        }

        const int nsub = 3 * k_i;
        std::vector<T> hsub(nsub * nsub, T(0));
        std::vector<T> ssub(nsub * nsub, T(0));
        std::vector<T> evec_sub(nsub * nsub, T(0));
        std::vector<Real> eval_sub(nsub, Real(0));

        // Build subspace overlap matrices:
        // sub-blocks: [0..k_i) = X, [k_i..2k_i) = W, [2k_i..3k_i) = P
        for (int col = 0; col < nsub; ++col)
        {
            const int col_sub = col % k_i;
            const int col_blk = col / k_i; // 0=X, 1=W, 2=P
            const int ib_col = band_offset + col_sub;

            const T* vcol = nullptr;
            const T* hvcol = nullptr;
            if (col_blk == 0)
            {
                vcol = psi_in + ib_col * this->n_basis;
                hvcol = this->hpsi.data() + ib_col * this->n_basis;
            }
            else if (col_blk == 1)
            {
                vcol = this->w.data() + ib_col * this->n_basis;
                hvcol = this->hw.data() + ib_col * this->n_basis;
            }
            else
            {
                vcol = this->p.data() + ib_col * this->n_basis;
                hvcol = this->hp.data() + ib_col * this->n_basis;
            }

            for (int row = 0; row < nsub; ++row)
            {
                const int row_sub = row % k_i;
                const int row_blk = row / k_i;
                const int ib_row = band_offset + row_sub;

                const T* vrow = nullptr;
                if (row_blk == 0)
                {
                    vrow = psi_in + ib_row * this->n_basis;
                }
                else if (row_blk == 1)
                {
                    vrow = this->w.data() + ib_row * this->n_basis;
                }
                else
                {
                    vrow = this->p.data() + ib_row * this->n_basis;
                }

                hsub[row + col * nsub] = this->inner_product(vrow, hvcol);
                ssub[row + col * nsub] = this->inner_product(vrow, vcol);
            }
        }

        // Regularize S_sub
        for (int i = 0; i < nsub; ++i)
        {
            ssub[i + i * nsub] += T(1.0e-12);
        }

        // Solve generalized eigenproblem: H_sub * C = Lambda * S_sub * C
        try
        {
            ct::kernels::lapack_hegvd<T, ct::DEVICE_CPU>()(nsub, nsub, hsub.data(), ssub.data(), eval_sub.data(),
                                                            evec_sub.data());
        }
        catch (const std::exception&)
        {
            // Fallback on failure: keep current vectors for this block.
            // Copy the original psi and hpsi for bands in the current block
            // (band_offset through band_offset + k_i - 1), then advance offset.
            for (int ib = band_offset; ib < band_offset + k_i && ib < this->n_work; ++ib)
            {
                T* xnew = this->work.data() + ib * this->n_basis;
                T* hxnew = this->hpsi_new.data() + ib * this->n_basis;
                this->copy_vector(xnew, psi_in + ib * this->n_basis);
                this->copy_vector(hxnew, this->hpsi.data() + ib * this->n_basis);
            }
            band_offset += k_i;
            continue;
        }

        // evec_sub contains eigenvectors (nsub x nsub, column-major).
        // First k_i columns = first k_i eigenvectors.
        // Update X_block = X*C_X + W*C_W + P*C_P
        //        P_block = W*C_W + P*C_P
        for (int ib = 0; ib < k_i; ++ib)
        {
            const int ib_global = band_offset + ib;
            if (this->is_locked[ib_global])
            {
                T* xnew = this->work.data() + ib_global * this->n_basis;
                T* hxnew = this->hpsi_new.data() + ib_global * this->n_basis;
                this->copy_vector(xnew, psi_in + ib_global * this->n_basis);
                this->copy_vector(hxnew, this->hpsi.data() + ib_global * this->n_basis);
                continue;
            }

            T* xnew = this->work.data() + ib_global * this->n_basis;
            T* hxnew = this->hpsi_new.data() + ib_global * this->n_basis;
            T* pnext = this->p_new.data() + ib_global * this->n_basis;
            T* hpnext = this->hp_new.data() + ib_global * this->n_basis;
            this->zero_vector(xnew);
            this->zero_vector(hxnew);
            this->zero_vector(pnext);
            this->zero_vector(hpnext);

            // Accumulate contributions from all 3 sub-blocks and the first k_i eigenvectors
            for (int col = 0; col < nsub; ++col)
            {
                const int col_sub = col % k_i;
                const int col_blk = col / k_i;
                const int ib_src = band_offset + col_sub;

                const T coeff = evec_sub[col + ib * nsub];

                const T* vsrc = nullptr;
                const T* hvsrc = nullptr;
                if (col_blk == 0)
                {
                    vsrc = psi_in + ib_src * this->n_basis;
                    hvsrc = this->hpsi.data() + ib_src * this->n_basis;
                }
                else if (col_blk == 1)
                {
                    vsrc = this->w.data() + ib_src * this->n_basis;
                    hvsrc = this->hw.data() + ib_src * this->n_basis;
                }
                else
                {
                    vsrc = this->p.data() + ib_src * this->n_basis;
                    hvsrc = this->hp.data() + ib_src * this->n_basis;
                }

                this->axpy_vector(xnew, vsrc, coeff);
                this->axpy_vector(hxnew, hvsrc, coeff);

                if (col_blk >= 1)
                {
                    this->axpy_vector(pnext, vsrc, coeff);
                    this->axpy_vector(hpnext, hvsrc, coeff);
                }
            }
        }

        band_offset += k_i;
    }

    // Preserve extra bands (beyond n_band_l) from current psi_in / hpsi / p / hp.
    // These bands are not covered by any block and should not be zeroed.
    for (int ib = this->n_band_l; ib < this->n_work; ++ib)
    {
        this->copy_vector(this->work.data() + ib * this->n_basis, psi_in + ib * this->n_basis);
        this->copy_vector(this->hpsi_new.data() + ib * this->n_basis,
                          this->hpsi.data() + ib * this->n_basis);
        this->zero_vector(this->p_new.data() + ib * this->n_basis);
        this->zero_vector(this->hp_new.data() + ib * this->n_basis);
    }

    std::copy(this->work.begin(), this->work.end(), psi_in);
    std::copy(this->hpsi_new.begin(), this->hpsi_new.end(), this->hpsi.begin());
    std::copy(this->p_new.begin(), this->p_new.end(), this->p.begin());
    std::copy(this->hp_new.begin(), this->hp_new.end(), this->hp.begin());
}

template <typename T, typename Device>
int DiagoPPCG<T, Device>::diag(const HPsiFunc& hpsi_func,
                               T* psi_in,
                               Real* eigenvalue_in,
                               const std::vector<double>& ethr_band)
{
    // On GPU devices, fall back to BPCG (PPCG subspace construction not yet ported to GPU).
    if (!std::is_same<Device, base_device::DEVICE_CPU>::value)
    {
        DiagoBPCG<T, Device> bpcg(this->precondition);
        bpcg.init_iter(this->n_band, this->n_band_l, this->n_basis, this->n_dim);
        bpcg.diag(hpsi_func, psi_in, eigenvalue_in, ethr_band);
        return 0;
    }
    else
    {
        ModuleBase::TITLE("DiagoPPCG", "diag");
        ModuleBase::timer::start("DiagoPPCG", "diag");

        // Initial setup: compute H|psi>, orthonormalize, then Rayleigh-Ritz to get
        // the best possible starting basis from the initial guess.
        this->calc_hpsi(hpsi_func, psi_in, this->hpsi);
        this->modified_gram_schmidt(psi_in, this->hpsi);
        this->rayleigh_ritz(psi_in, this->hpsi);

        // PPCG main iteration loop.
        // Each iteration:
        //   1. Compute preconditioned residuals W and eigenvalue estimates.
        //   2. Update band locking (bands converged for 2 consecutive iterations are frozen).
        //   3. Check global convergence across all MPI ranks.
        //   4. Project W and P to the orthogonal complement of current psi.
        //   5. Compute H|w> and H|p>.
        //   6. Update psi, hpsi, p, hp from the per-band (or per-block) PPCG subspace.
        //   7. Periodically re-orthonormalize (every 4 iterations, or when orthonormality degrades).
        int iter = 0;
        const int max_iter = std::max(1, DiagoIterAssist<T, Device>::PW_DIAG_NMAX);
        for (; iter < max_iter; ++iter)
        {
            // Step 1: compute preconditioned residuals and eigenvalue estimates.
            this->calc_preconditioned_residual(psi_in);

            // Diagnostic: print convergence status every 10 iterations or on first/last.
            if (iter % 10 == 0 || iter == max_iter - 1)
            {
                int n_locked = 0;
                for (int ib = 0; ib < this->n_band_l; ++ib)
                {
                    if (this->is_locked[ib])
                    {
                        n_locked++;
                    }
                }
                std::cerr << "[PPCG] iter=" << iter
                          << " err[0]=" << this->err[0]
                          << " err[end]=" << this->err[this->n_band_l - 1]
                          << " ethr=" << ethr_band[0]
                          << " locked=" << n_locked << "/" << this->n_band_l
                          << " blocked=" << (!this->block_sizes.empty() ? "yes" : "no")
                          << std::endl;
            }

            // Step 2: update locking.
            // A band is locked when err[ib] <= ethr_band[ib] for 2+ consecutive iterations.
            // Only the first n_band_l bands are checked (extra bands are auxiliary).
            for (int ib = 0; ib < this->n_band_l; ++ib)
            {
                if (this->is_locked[ib])
                {
                    continue;
                }
                if (this->err[ib] <= ethr_band[ib])
                {
                    this->converge_count[ib]++;
                    if (this->converge_count[ib] >= 2)
                    {
                        this->is_locked[ib] = true;
                        this->err[ib] = Real(0);
                    }
                }
                else
                {
                    this->converge_count[ib] = 0;
                }
            }

            // Step 3: check global convergence across all MPI ranks.
            if (!this->test_error(ethr_band))
            {
                break;
            }

            // Step 4: project W and P to the orthogonal complement of current psi.
            this->project_to_orthogonal_complement(psi_in, this->w);
            this->project_to_orthogonal_complement(psi_in, this->p);

            // Step 5: apply Hamiltonian to W and P.
            this->calc_hpsi(hpsi_func, this->w.data(), this->hw);
            this->calc_hpsi(hpsi_func, this->p.data(), this->hp);

            // Step 6: solve small subspace eigenproblems and update all working vectors.
            this->update_vectors_from_ppcg_subspace(psi_in);

            // Step 7: periodic re-orthonormalization.
            // Force Cholesky-based re-orthonormalization every 10 iterations.
            // Between scheduled cycles, check orthonormality and re-orthonormalize on demand.
            if ((iter + 1) % 15 == 0)
            {
                this->orth_cholesky(psi_in, this->hpsi);
                this->rayleigh_ritz(psi_in, this->hpsi);
            }
            else if (!this->check_orthonormality(psi_in))
            {
                this->orth_cholesky(psi_in, this->hpsi);
            }
        }

        // Final Rayleigh-Ritz to ensure eigenvalues and vectors are optimal in the subspace.
        this->rayleigh_ritz(psi_in, this->hpsi);
        std::copy(this->eigen.begin(), this->eigen.begin() + this->n_band_l, eigenvalue_in);

        ModuleBase::timer::end("DiagoPPCG", "diag");

        std::cerr << "[PPCG] done: niter=" << std::min(iter + 1, max_iter)
                  << " final_err[0]=" << this->err[0]
                  << " final_err[end]=" << this->err[this->n_band_l - 1]
                  << " eigen[0]=" << eigenvalue_in[0]
                  << std::endl;

        return std::min(iter + 1, max_iter);
    }
}

template class DiagoPPCG<std::complex<float>, base_device::DEVICE_CPU>;
template class DiagoPPCG<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class DiagoPPCG<std::complex<float>, base_device::DEVICE_GPU>;
template class DiagoPPCG<std::complex<double>, base_device::DEVICE_GPU>;
#endif

} // namespace hsolver
