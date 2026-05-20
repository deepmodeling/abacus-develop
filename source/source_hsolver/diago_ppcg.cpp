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
    return std::sqrt(frob2) < Real(1e-6);
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::rotate_block(T* block, const std::vector<T>& coeff, std::vector<T>& workspace) const
{
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
    // Block diagonal mode: solve per-block instead of per-band
    if (!this->block_sizes.empty())
    {
        this->update_vectors_blocked(psi_in);
        return;
    }

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
            // Fallback on failure: keep current vectors for this block
            band_offset += k_i;
            for (int ib = band_offset; ib < band_offset + k_i && ib < this->n_work; ++ib)
            {
                T* xnew = this->work.data() + ib * this->n_basis;
                T* hxnew = this->hpsi_new.data() + ib * this->n_basis;
                this->copy_vector(xnew, psi_in + ib * this->n_basis);
                this->copy_vector(hxnew, this->hpsi.data() + ib * this->n_basis);
            }
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

        this->calc_hpsi(hpsi_func, psi_in, this->hpsi);
        this->modified_gram_schmidt(psi_in, this->hpsi);
        this->rayleigh_ritz(psi_in, this->hpsi);

        int iter = 0;
        const int max_iter = std::max(1, DiagoIterAssist<T, Device>::PW_DIAG_NMAX);
        for (; iter < max_iter; ++iter)
        {
            this->calc_preconditioned_residual(psi_in);

            // Update locking: bands converged for 2+ consecutive iterations are locked
            // Only check the first n_band_l bands (extra bands are auxiliary)
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

            if (!this->test_error(ethr_band))
            {
                break;
            }

            this->project_to_orthogonal_complement(psi_in, this->w);
            this->project_to_orthogonal_complement(psi_in, this->p);

            this->calc_hpsi(hpsi_func, this->w.data(), this->hw);
            this->calc_hpsi(hpsi_func, this->p.data(), this->hp);

            this->update_vectors_from_ppcg_subspace(psi_in);

            if ((iter + 1) % 4 == 0)
            {
                this->orth_cholesky(psi_in, this->hpsi);
                this->rayleigh_ritz(psi_in, this->hpsi);
            }
            else if (!this->check_orthonormality(psi_in))
            {
                this->orth_cholesky(psi_in, this->hpsi);
            }
        }

        this->rayleigh_ritz(psi_in, this->hpsi);
        std::copy(this->eigen.begin(), this->eigen.begin() + this->n_band_l, eigenvalue_in);

        ModuleBase::timer::end("DiagoPPCG", "diag");
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
