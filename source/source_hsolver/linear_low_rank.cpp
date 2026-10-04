#include "source_hsolver/linear_low_rank.h"

#include "source_hsolver/kernels/linear_op.h"
#include "source_hsolver/linear_workspace.h"

namespace hsolver
{
template <typename T, typename Device>
void LinearResponse<T, Device>::update(LinearAlgebra<T, Device>& algebra,
                                       int ld,
                                       int dim,
                                       int bands,
                                       const T* solution,
                                       const T* seed,
                                       const T* seed_residual,
                                       double tolerance,
                                       ct::Tensor* workspace)
{
    rank_ = 0;
    const int64_t size = std::max<int64_t>(1, static_cast<int64_t>(ld) * bands);
    linear_buffer<T, Device>(&directions_, size);
    linear_buffer<T, Device>(&images_, size);
    linear_buffer<T, Device>(workspace, 4 * size);
    T* z = workspace->template data<T>();
    T* w = z + size;
    using Real = typename GetTypeReal<T>::type;
    const double cutoff = std::max(100 * tolerance, 100 * static_cast<double>(std::numeric_limits<Real>::epsilon()));
    linear_op<T, Device>().batch(ld, dim, bands, z, solution, seed, T(1), T(-1), nullptr, nullptr, nullptr);
    if (dim > 0)
    {
        base_device::memory::synchronize_memory_2d_op<T, Device, Device>()(w, ld, seed_residual, ld, dim, bands);
    }
    // Form a rank-revealing block basis, then repair its Gram matrix once.
    // This avoids a device transfer and MPI reduction for every history column.
    std::vector<std::complex<double>> gram = algebra.cross(ld, dim, bands, bands, w, w);
    std::vector<std::complex<double>> c = linear_gram_basis(gram, bands, cutoff, &rank_);
    if (rank_ == 0)
    {
        return;
    }
    T* rz = z + 2 * size;
    T* rw = rz + size;
    algebra.expand(ld, dim, bands, rank_, z, c, rz, T(0));
    algebra.expand(ld, dim, bands, rank_, w, c, rw, T(0));
    gram = algebra.cross(ld, dim, rank_, rank_, rw, rw);
    int refined_rank = 0;
    c = linear_gram_basis(gram, rank_, 1e-6, &refined_rank);
    if (refined_rank > 0)
    {
        algebra.expand(ld, dim, rank_, refined_rank, rz, c, directions_.template data<T>(), T(0));
        algebra.expand(ld, dim, rank_, refined_rank, rw, c, images_.template data<T>(), T(0));
    }
    rank_ = refined_rank;
    if (rank_ > 0)
    {
        gram = algebra.cross(ld, dim, rank_, rank_, images(), images());
        const double limit = std::max(1e-8, 100.0 * std::numeric_limits<Real>::epsilon());
        for (int j = 0; j < rank_; ++j)
        {
            for (int i = 0; i < rank_; ++i)
            {
                if (!std::isfinite(std::abs(gram[i + j * rank_]))
                    || std::abs(gram[i + j * rank_] - (i == j ? std::complex<double>(1) : std::complex<double>(0))) > limit)
                {
                    rank_ = 0;
                    return;
                }
            }
        }
    }
}

template <typename T, typename Device>
void LinearLowRank<T, Device>::prepare_response(int ld, int rank, const T* directions, const T* images)
{
    prepare(ld, rank, directions, images, images, nullptr);
}

template <typename T, typename Device>
void LinearLowRank<T, Device>::prepare_subspace(int ld, int rank, const T* basis, const T* images, const LinearSmallLU& factor)
{
    prepare(ld, rank, basis, images, basis, &factor);
}

template <typename T, typename Device>
void LinearLowRank<T, Device>::prepare(int ld, int rank, const T* z, const T* w, const T* test, const LinearSmallLU* factor)
{
    rank_ = rank;
    test_ = test;
    factor_ = factor;
    if (rank == 0)
    {
        return;
    }
    linear_buffer<T, Device>(&correction_, static_cast<int64_t>(ld) * rank);
    T* correction = correction_.template data<T>();
    if (diagonal_)
    {
        linear_op<T, Device>().diagonal(ld, dim_, rank, diagonal_, w, correction);
    }
    else if (dim_ > 0)
    {
        base_device::memory::synchronize_memory_2d_op<T, Device, Device>()(correction, ld, w, ld, dim_, rank);
    }
    linear_op<T, Device>().batch(ld, dim_, rank, correction, z, correction, T(1), T(-1), nullptr, nullptr, nullptr);
}

template <typename T, typename Device>
void LinearLowRank<T, Device>::apply(const T* x, T* y, int ld, int nvec) const
{
    if (diagonal_)
    {
        linear_op<T, Device>().diagonal(ld, dim_, nvec, diagonal_, x, y);
    }
    else if (dim_ > 0)
    {
        base_device::memory::synchronize_memory_2d_op<T, Device, Device>()(y, ld, x, ld, dim_, nvec);
    }
    if (rank_ == 0)
    {
        return;
    }
    std::vector<std::complex<double>> coefficients = algebra_.projection_cross(ld, dim_, rank_, nvec, test_, x);
    if (factor_ && !factor_->solve(&coefficients, nvec))
    {
        // Changing the preconditioner inside a recurrence invalidates BiCGSTAB and CGS.
        throw LinearPreconditionerError();
    }
    algebra_.expand(ld, dim_, rank_, nvec, correction_.template data<T>(), coefficients, y, T(1));
}

template class LinearResponse<std::complex<float>, base_device::DEVICE_CPU>;
template class LinearResponse<std::complex<double>, base_device::DEVICE_CPU>;
#if defined(__CUDA) || defined(__ROCM)
template class LinearResponse<std::complex<float>, base_device::DEVICE_GPU>;
template class LinearResponse<std::complex<double>, base_device::DEVICE_GPU>;
#endif
template class LinearLowRank<std::complex<float>, base_device::DEVICE_CPU>;
template class LinearLowRank<std::complex<double>, base_device::DEVICE_CPU>;
#if defined(__CUDA) || defined(__ROCM)
template class LinearLowRank<std::complex<float>, base_device::DEVICE_GPU>;
template class LinearLowRank<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace hsolver
