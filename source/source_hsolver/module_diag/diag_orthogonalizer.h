#ifndef DIAG_ORTHOGONALIZER_H_
#define DIAG_ORTHOGONALIZER_H_

#include "source_base/kernels/math_kernel_op.h"
#include "source_base/macros.h"
#include "source_base/module_device/device.h"
#include "source_base/module_device/memory_op.h"
#include "source_base/module_device/types.h"
#include "source_base/parallel_reduce.h"
#include "source_base/para_gemm.h"
#include "source_base/tool_quit.h"
#include "source_hsolver/para_linear_transform.h"

#include <ATen/core/tensor_types.h>
#include <ATen/kernels/lapack.h>

#include <cmath>
#include <complex>
#include <vector>

namespace hsolver
{

template <typename T>
struct DiagOrthScalar
{
    static const T* one()
    {
        static const T value = static_cast<T>(1.0);
        return &value;
    }

    static const T* zero()
    {
        static const T value = static_cast<T>(0.0);
        return &value;
    }

    static const T* neg_one()
    {
        static const T value = static_cast<T>(-1.0);
        return &value;
    }
};

/**
 * Shared orthogonalization kernels for iterative diagonalizers.
 *
 * The class intentionally knows only about dense block vectors and BLAS-like
 * operations. Algorithm classes decide when to orthogonalize; this helper owns
 * the repeated mechanics: overlap matrices, projection, normalization checks,
 * modified Gram-Schmidt, and Cholesky orthogonalization.
 */
template <typename T, typename Device>
class DiagOrthogonalizer
{
  private:
    using Real = typename GetTypeReal<T>::type;
    using ct_Device = typename ct::PsiToContainer<Device>::type;

    using resmem_op = base_device::memory::resize_memory_op<T, Device>;
    using delmem_op = base_device::memory::delete_memory_op<T, Device>;
    using setmem_op = base_device::memory::set_memory_op<T, Device>;
    using syncmem_op = base_device::memory::synchronize_memory_op<T, Device, Device>;
    using syncmem_d2h = base_device::memory::synchronize_memory_op<T, base_device::DEVICE_CPU, Device>;
    using syncmem_h2d = base_device::memory::synchronize_memory_op<T, Device, base_device::DEVICE_CPU>;

  public:
    DiagOrthogonalizer(const int dim, const int lda) : dim_(dim), lda_(lda)
    {
    }

    Real vector_norm(const T* vec) const
    {
        Real norm = ModuleBase::dot_real_op<T, Device>()(this->dim_, vec, vec, false);
#ifdef __MPI
        Parallel_Reduce::reduce_pool(norm);
#endif
        return std::sqrt(norm);
    }

    void scale_vector(T* vec, const Real alpha) const
    {
        ModuleBase::vector_mul_real_op<T, Device>()(this->dim_, vec, vec, alpha);
        if (this->lda_ > this->dim_)
        {
            setmem_op()(vec + this->dim_, 0, this->lda_ - this->dim_);
        }
    }

    void rotate_block(T* block, const int ncol, const T* coeff, T* workspace) const
    {
        T* d_coeff = nullptr;
        resmem_op()(d_coeff, ncol * ncol);
        syncmem_h2d()(d_coeff, coeff, ncol * ncol);

        ModuleBase::gemm_op<T, Device>()('N',
                                         'N',
                                         this->dim_,
                                         ncol,
                                         ncol,
                                         DiagOrthScalar<T>::one(),
                                         block,
                                         this->lda_,
                                         d_coeff,
                                         ncol,
                                         DiagOrthScalar<T>::zero(),
                                         workspace,
                                         this->lda_);
        delmem_op()(d_coeff);
        syncmem_op()(block, workspace, this->lda_ * ncol);
    }

    void modified_gram_schmidt(T* block, T* hblock, const int ncol) const
    {
        for (int ib = 0; ib < ncol; ++ib)
        {
            T* xi = block + ib * this->lda_;
            T* hxi = hblock == nullptr ? nullptr : hblock + ib * this->lda_;

            if (ib > 0)
            {
                T* d_lag = nullptr;
                resmem_op()(d_lag, ib);
                setmem_op()(d_lag, 0, ib);
                ModuleBase::gemv_op<T, Device>()('C',
                                                 this->dim_,
                                                 ib,
                                                 DiagOrthScalar<T>::one(),
                                                 block,
                                                 this->lda_,
                                                 xi,
                                                 1,
                                                 DiagOrthScalar<T>::zero(),
                                                 d_lag,
                                                 1);

                std::vector<T> lag(ib);
                syncmem_d2h()(lag.data(), d_lag, ib);
                delmem_op()(d_lag);
#ifdef __MPI
                Parallel_Reduce::reduce_pool(lag.data(), ib);
#endif

                T* d_lag_reduced = nullptr;
                resmem_op()(d_lag_reduced, ib);
                syncmem_h2d()(d_lag_reduced, lag.data(), ib);

                ModuleBase::gemv_op<T, Device>()('N',
                                                 this->dim_,
                                                 ib,
                                                 DiagOrthScalar<T>::neg_one(),
                                                 block,
                                                 this->lda_,
                                                 d_lag_reduced,
                                                 1,
                                                 DiagOrthScalar<T>::one(),
                                                 xi,
                                                 1);
                if (hxi != nullptr)
                {
                    ModuleBase::gemv_op<T, Device>()('N',
                                                     this->dim_,
                                                     ib,
                                                     DiagOrthScalar<T>::neg_one(),
                                                     hblock,
                                                     this->lda_,
                                                     d_lag_reduced,
                                                     1,
                                                     DiagOrthScalar<T>::one(),
                                                     hxi,
                                                     1);
                }
                delmem_op()(d_lag_reduced);
            }

            const Real norm = this->vector_norm(xi);
            if (norm <= Real(1.0e-14))
            {
                ModuleBase::WARNING_QUIT("DiagOrthogonalizer::modified_gram_schmidt",
                                         "linear dependent wavefunctions");
            }
            this->scale_vector(xi, Real(1) / norm);
            if (hxi != nullptr)
            {
                this->scale_vector(hxi, Real(1) / norm);
            }
        }
    }

    void cholesky_orth(T* block, T* hblock, T* workspace, const int ncol) const
    {
        T* d_s = nullptr;
        resmem_op()(d_s, ncol * ncol);
        setmem_op()(d_s, 0, ncol * ncol);
        ModuleBase::gemm_op<T, Device>()('C',
                                         'N',
                                         ncol,
                                         ncol,
                                         this->dim_,
                                         DiagOrthScalar<T>::one(),
                                         block,
                                         this->lda_,
                                         block,
                                         this->lda_,
                                         DiagOrthScalar<T>::zero(),
                                         d_s,
                                         ncol);

        std::vector<T> s(ncol * ncol);
        syncmem_d2h()(s.data(), d_s, ncol * ncol);
        delmem_op()(d_s);
#ifdef __MPI
        Parallel_Reduce::reduce_pool(s.data(), ncol * ncol);
#endif

        ct::kernels::lapack_potrf<T, ct::DEVICE_CPU>()('U', ncol, s.data(), ncol);
        for (int col = 0; col < ncol; ++col)
        {
            for (int row = col + 1; row < ncol; ++row)
            {
                s[row + col * ncol] = T(0);
            }
        }
        ct::kernels::lapack_trtri<T, ct::DEVICE_CPU>()('U', 'N', ncol, s.data(), ncol);

        this->rotate_block(block, ncol, s.data(), workspace);
        if (hblock != nullptr)
        {
            this->rotate_block(hblock, ncol, s.data(), workspace);
        }
    }

    bool check_orthonormality(const T* block, const int ncol, const Real tolerance) const
    {
        T* d_s = nullptr;
        resmem_op()(d_s, ncol * ncol);
        setmem_op()(d_s, 0, ncol * ncol);
        ModuleBase::gemm_op<T, Device>()('C',
                                         'N',
                                         ncol,
                                         ncol,
                                         this->dim_,
                                         DiagOrthScalar<T>::one(),
                                         block,
                                         this->lda_,
                                         block,
                                         this->lda_,
                                         DiagOrthScalar<T>::zero(),
                                         d_s,
                                         ncol);

        std::vector<T> s(ncol * ncol);
        syncmem_d2h()(s.data(), d_s, ncol * ncol);
        delmem_op()(d_s);
#ifdef __MPI
        Parallel_Reduce::reduce_pool(s.data(), ncol * ncol);
#endif

        Real frob2 = 0;
        for (int col = 0; col < ncol; ++col)
        {
            for (int row = 0; row < ncol; ++row)
            {
                const T delta = s[row + col * ncol] - static_cast<T>(row == col ? 1.0 : 0.0);
                frob2 += std::norm(delta);
            }
        }
        return std::sqrt(frob2) < tolerance;
    }

    void project_out(const T* basis, T* block, const int basis_cols, const int block_cols) const
    {
        T* d_coeff = nullptr;
        resmem_op()(d_coeff, basis_cols * block_cols);
        setmem_op()(d_coeff, 0, basis_cols * block_cols);
        ModuleBase::gemm_op<T, Device>()('C',
                                         'N',
                                         basis_cols,
                                         block_cols,
                                         this->dim_,
                                         DiagOrthScalar<T>::one(),
                                         basis,
                                         this->lda_,
                                         block,
                                         this->lda_,
                                         DiagOrthScalar<T>::zero(),
                                         d_coeff,
                                         basis_cols);

        std::vector<T> coeff(basis_cols * block_cols);
        syncmem_d2h()(coeff.data(), d_coeff, basis_cols * block_cols);
        delmem_op()(d_coeff);
#ifdef __MPI
        Parallel_Reduce::reduce_pool(coeff.data(), basis_cols * block_cols);
#endif

        T* d_coeff_reduced = nullptr;
        resmem_op()(d_coeff_reduced, basis_cols * block_cols);
        syncmem_h2d()(d_coeff_reduced, coeff.data(), basis_cols * block_cols);

        ModuleBase::gemm_op<T, Device>()('N',
                                         'N',
                                         this->dim_,
                                         block_cols,
                                         basis_cols,
                                         DiagOrthScalar<T>::neg_one(),
                                         basis,
                                         this->lda_,
                                         d_coeff_reduced,
                                         basis_cols,
                                         DiagOrthScalar<T>::one(),
                                         block,
                                         this->lda_);
        delmem_op()(d_coeff_reduced);
    }

    void overlap_with_metric(const T* basis,
                             const T* metric_block,
                             T* coeff,
                             const int basis_cols,
                             const int block_cols) const
    {
        if (basis_cols <= 0 || block_cols <= 0)
        {
            return;
        }
        ModuleBase::gemm_op<T, Device>()('C',
                                         'N',
                                         basis_cols,
                                         block_cols,
                                         this->dim_,
                                         DiagOrthScalar<T>::one(),
                                         basis,
                                         this->lda_,
                                         metric_block,
                                         this->lda_,
                                         DiagOrthScalar<T>::zero(),
                                         coeff,
                                         basis_cols);
#ifdef __MPI
        Parallel_Reduce::reduce_pool(coeff, basis_cols * block_cols);
#endif
    }

    void project_out_with_coeff(const T* basis,
                                const T* coeff,
                                T* block,
                                const int basis_cols,
                                const int block_cols) const
    {
        if (basis_cols <= 0 || block_cols <= 0)
        {
            return;
        }
        ModuleBase::gemm_op<T, Device>()('N',
                                         'N',
                                         this->dim_,
                                         block_cols,
                                         basis_cols,
                                         DiagOrthScalar<T>::neg_one(),
                                         basis,
                                         this->lda_,
                                         coeff,
                                         basis_cols,
                                         DiagOrthScalar<T>::one(),
                                         block,
                                         this->lda_);
    }

    Real schmidt_orthogonalize_s_metric(const T* basis,
                                        const T* s_target,
                                        T* target,
                                        T* coeff,
                                        const int current_col,
                                        const Real min_norm = Real(0),
                                        const char* warning_source
                                        = "DiagOrthogonalizer::schmidt_orthogonalize_s_metric") const
    {
        this->overlap_with_metric(basis, s_target, coeff, current_col + 1, 1);
        this->project_out_with_coeff(basis, coeff, target, current_col, 1);

        T raw_norm = T(0);
        syncmem_d2h()(&raw_norm, coeff + current_col, 1);
        Real norm2 = static_cast<Real>(std::real(raw_norm))
                     - ModuleBase::dot_real_op<T, Device>()(current_col, coeff, coeff, false);
        if (norm2 <= Real(0))
        {
            ModuleBase::WARNING_QUIT("DiagOrthogonalizer::schmidt_orthogonalize_s_metric",
                                     "psi_norm <= 0.0");
        }

        const Real norm = std::sqrt(norm2);
        if (norm <= min_norm)
        {
            ModuleBase::WARNING_QUIT(warning_source, "psi_norm is below the orthogonalization threshold");
        }
        this->scale_vector(target, Real(1) / norm);
        return norm;
    }

    void project_out_parallel(const T* basis,
                              T* block,
                              T* coeff,
                              ModuleBase::PGemmCN<T, Device>& pmmcn,
                              PLinearTransform<T, Device>& plintrans) const
    {
        pmmcn.multiply(1.0, basis, block, 0.0, coeff);
        plintrans.act(-1.0, basis, coeff, 1.0, block);
    }

    void cholesky_orth_parallel(T* workspace,
                                T* block,
                                T* hblock,
                                T* coeff,
                                const int ncol,
                                ModuleBase::PGemmCN<T, Device>& pmmcn,
                                PLinearTransform<T, Device>& plintrans) const
    {
        pmmcn.multiply(1.0, block, block, 0.0, coeff);

        ct::kernels::set_matrix<T, ct_Device>()('L', coeff, ncol);
        ct::kernels::lapack_potrf<T, ct_Device>()('U', ncol, coeff, ncol);
        ct::kernels::lapack_trtri<T, ct_Device>()('U', 'N', ncol, coeff, ncol);

        this->rotate_parallel(block, coeff, workspace, ncol, plintrans);
        this->rotate_parallel(hblock, coeff, workspace, ncol, plintrans);
    }

  private:
    void rotate_parallel(T* block,
                         T* coeff,
                         T* workspace,
                         const int ncol,
                         PLinearTransform<T, Device>& plintrans) const
    {
        plintrans.act(1.0, block, coeff, 0.0, workspace);
        syncmem_op()(block, workspace, this->lda_ * ncol);
    }

    int dim_ = 0;
    int lda_ = 0;
};

} // namespace hsolver

#endif // DIAG_ORTHOGONALIZER_H_
