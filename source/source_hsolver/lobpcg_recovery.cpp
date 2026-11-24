#include "diago_lobpcg.h"
#include "lobpcg_internal.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_base/memory_recorder.h"
#include "source_base/module_container/ATen/kernels/lapack.h"
#include "source_base/module_container/base/third_party/blas.h"
#include "source_base/timer.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace hsolver
{
using lobpcg_internal::host_tensor;

namespace
{
bool finite_projection(const ct::Tensor& matrix)
{
    const std::complex<double>* values = matrix.data<std::complex<double>>();
    for (int64_t i = 0; i < matrix.NumElements(); ++i)
    {
        if (!std::isfinite(values[i].real()) || !std::isfinite(values[i].imag()))
        {
            return false;
        }
    }
    return true;
}

// For G = U^H U, ||G||_1 bounds lambda_max and
// 1 / (||U^-1||_1 ||U^-1||_inf) bounds lambda_min from below.
// Use Cholesky for well-conditioned full-rank blocks; otherwise use EVD
// to determine the rank. The caller checks the transformed block.
bool cholesky_rank_transform(const ct::Tensor& gram, const int m, const double relative, ct::Tensor& transform)
{
    using C = std::complex<double>;
    double upper = 0.0;
    for (int j = 0; j < m; ++j)
    {
        double sum = 0.0;
        for (int i = 0; i < m; ++i)
        {
            sum += std::abs(gram.data<std::complex<double>>()[i + j * m]);
        }
        upper = std::max(upper, sum);
    }
    transform = gram;
    int info = 0;
    ct::lapackConnector::potrf('U', m, transform.data<std::complex<double>>(), m, info);
    if (info != 0)
    {
        return false;
    }
    ct::lapackConnector::trtri('U', 'N', m, transform.data<std::complex<double>>(), m, info);
    if (info != 0)
    {
        return false;
    }
    double norm_one = 0.0;
    double norm_inf = 0.0;
    for (int j = 0; j < m; ++j)
    {
        for (int i = j + 1; i < m; ++i)
        {
            transform.data<std::complex<double>>()[i + j * m] = C(0);
        }
    }
    for (int j = 0; j < m; ++j)
    {
        double column = 0.0;
        double row = 0.0;
        for (int i = 0; i < m; ++i)
        {
            column += std::abs(transform.data<std::complex<double>>()[i + j * m]);
            row += std::abs(transform.data<std::complex<double>>()[j + i * m]);
        }
        norm_one = std::max(norm_one, column);
        norm_inf = std::max(norm_inf, row);
    }
    const double inverse_bound = norm_one * norm_inf;
    const double cutoff = relative * upper;
    const double condition_bound = upper * inverse_bound;
    return std::isfinite(condition_bound) && condition_bound <= 1e4 && inverse_bound > 0.0 &&
           1.0 / inverse_bound > cutoff;
}

} // namespace

/// After GV metric failure, refresh true H/S images and retry; if needed, repair only the increments.
template <typename T, typename Device>
bool DiagoLOBPCG<T, Device>::recover_projected(const HSOperator<T, Device>& op,
    const VectorImages& basis,
    int& working_columns,
    const int wanted,
    ct::Tensor& projected,
    ct::Tensor& gram,
    SerialGvWorkspace& workspace)
{
    ModuleBase::timer::Scope timer("DiagoLOBPCG", "recover_projected");
    const int offset = n_max_ - wanted;
    int columns = working_columns - offset;
    // Reapply the physical operators: independently rotated images can lose
    // accuracy when trial directions nearly cancel.
    op.hpsi(basis.x, basis.hx, n_dim_, columns);
    project_padding(basis.hx, columns);
    if (use_overlap_)
    {
        op.spsi(basis.x, basis.sx, n_dim_, columns);
        project_padding(basis.sx, columns);
    }
    form_projected_matrices(basis, columns, nullptr, projected, gram);
    if (solve_projected(projected, gram, columns, wanted, workspace))
    {
        return true;
    }
    const int added = working_columns - n_max_;
    const int kept = repair_basis(op, added);
    if (kept <= 0)
    {
        error_message_ += "; basis recovery failed, retained increments=" + std::to_string(kept);
        return false;
    }
    working_columns = n_max_ + kept;
    columns = working_columns - offset;
    form_projected_matrices(basis, columns, nullptr, projected, gram);
    return solve_projected(projected, gram, columns, wanted, workspace);
}

/// Preserve X, project increments with its actual S metric, then repair their rank and refresh H images.
template <typename T, typename Device>
int DiagoLOBPCG<T, Device>::repair_basis(const HSOperator<T, Device>& op, const int added)
{
    ModuleBase::timer::Scope function_timer("DiagoLOBPCG", "repair_basis");
    const int n = n_dim_;
    const int reference_columns = n_max_;
    const VectorImages reference = basis_images(0);
    const VectorImages increments = basis_images(reference_columns);
    const double epsilon = std::numeric_limits<Real>::epsilon();
    const Real tolerance = std::max(static_cast<Real>(1e-10), static_cast<Real>(100 * epsilon));
    const T minus_one(-1);
    int kept = added;
    if (use_overlap_)
    {
        op.spsi(reference.x, reference.sx, n_dim_, reference_columns);
        project_padding(reference.sx, reference_columns);
    }
    // Use the actual metric M = X^H S X for projection in finite precision.
    ct::Tensor reference_factor =
        reference.sx == nullptr
            ? algebra_.gram_tensor(n, n, reference_columns, reference.x)
            : algebra_.cross_tensor(n, n, reference_columns, reference_columns, reference.x, reference.sx);
    const std::size_t reference_bytes = reference_factor.buffer().GetAllocatedBytes();
    ModuleBase::Memory::record("LOBPCG::recovery_reference", reference_bytes, false);
    int factor_info = 0;
    if (comm_.rank == 0)
    {
        if (finite_projection(reference_factor))
        {
            ct::lapackConnector::potrf(
                'U', reference_columns, reference_factor.data<Wide>(), reference_columns, factor_info);
        }
        else
        {
            factor_info = -1;
        }
    }
    Real shared_info = static_cast<Real>(factor_info);
    broadcast(&shared_info, 1);
    if (shared_info != 0)
    {
        return -1;
    }
    if (use_overlap_ && kept > 0)
    {
        op.spsi(increments.x, increments.sx, n_dim_, kept);
        project_padding(increments.sx, kept);
    }
    for (int pass = 0; pass < 3 && kept > 0; ++pass)
    {
        const T* metric_image = increments.sx != nullptr ? increments.sx : increments.x;
        ct::Tensor cross = algebra_.cross_tensor(n, n, reference_columns, kept, reference.x, metric_image);
        if (comm_.rank == 0)
        {
            // M = R^H R: solve M C = X^H S D without forming an inverse.
            const char left = 'L';
            const char upper = 'U';
            const char adjoint = 'C';
            const char normal = 'N';
            const Wide one(1);
            ct::BlasConnector::trsm(left,
                upper,
                adjoint,
                normal,
                reference_columns,
                kept,
                one,
                reference_factor.data<Wide>(),
                reference_columns,
                cross.template data<Wide>(),
                reference_columns);
            ct::BlasConnector::trsm(left,
                upper,
                normal,
                normal,
                reference_columns,
                kept,
                one,
                reference_factor.data<Wide>(),
                reference_columns,
                cross.template data<Wide>(),
                reference_columns);
        }
        const int cross_size = reference_columns * kept;
        broadcast(cross.template data<Wide>(), cross_size);
        const T* coefficients = prepare_rotation_coefficients(cross);
        algebra_.expand_prepared(n, n, reference_columns, kept, reference.x, coefficients, increments.x, minus_one);
        if (use_overlap_)
        {
            op.spsi(increments.x, increments.sx, n_dim_, kept);
            project_padding(increments.sx, kept);
        }
        kept = rank_orthonormalize(kept, increments.x, increments.sx);
        if (kept <= 0)
        {
            return kept;
        }
        if (use_overlap_)
        {
            op.spsi(increments.x, increments.sx, n_dim_, kept);
            project_padding(increments.sx, kept);
        }
        metric_image = increments.sx != nullptr ? increments.sx : increments.x;
        const ct::Tensor overlap = algebra_.cross_tensor(n, n, reference_columns, kept, reference.x, metric_image);
        double coupling = 0.0;
        for (int64_t i = 0; i < overlap.NumElements(); ++i)
        {
            const Wide& value = overlap.template data<Wide>()[i];
            coupling = std::max(coupling, std::abs(value));
        }
        // Recovery acceptance: both coupling to X and the self metric must pass.
        bool orthonormal = true;
        {
            const ct::Tensor gram = increments.sx == nullptr
                                        ? algebra_.gram_tensor(n, n, kept, increments.x)
                                        : algebra_.cross_tensor(n, n, kept, kept, increments.x, increments.sx);
            for (int j = 0; j < kept && orthonormal; ++j)
            {
                for (int i = 0; i < kept; ++i)
                {
                    const Wide expected = i == j ? Wide(1) : Wide(0);
                    const Real deviation = std::abs(gram.data<Wide>()[i + j * kept] - expected);
                    if (!std::isfinite(deviation) || deviation > tolerance)
                    {
                        orthonormal = false;
                        break;
                    }
                }
            }
        }
        if (coupling <= tolerance && orthonormal)
        {
            op.hpsi(increments.x, increments.hx, n_dim_, kept);
            project_padding(increments.hx, kept);
            return kept;
        }
    }
    return -1;
}

/// Recovery-only metric normalization: use Cholesky for a well-conditioned block, otherwise spectral rank truncation.
template <typename T, typename Device>
int DiagoLOBPCG<T, Device>::rank_orthonormalize(const int m, T* x, T* sx)
{
    ModuleBase::timer::Scope function_timer("DiagoLOBPCG", "rank_orthonormalize");
    ct::Tensor gram =
        sx == nullptr ? algebra_.gram_tensor(n_dim_, n_dim_, m, x) : algebra_.cross_tensor(n_dim_, n_dim_, m, m, x, sx);
    ct::Tensor transform = host_tensor<Wide>(0, Wide(0));
    const double eps = std::numeric_limits<Real>::epsilon();
    // Set the rank cutoff by wavefunction precision before inverse-square-root
    // scaling amplifies roundoff. Explicit S images can introduce O(eps)
    // negative Gram eigenvalues.
    const double relative = 64 * std::max(std::numeric_limits<double>::epsilon(), eps);
    int rank = -1;
    bool triangular = false;
    std::size_t orthogonal_work = 0;
    if (comm_.rank == 0)
    {
        try
        {
            for (int64_t i = 0; i < gram.NumElements(); ++i)
            {
                const Wide& value = gram.data<Wide>()[i];
                if (!std::isfinite(std::abs(value)))
                {
                    throw std::runtime_error("Non-finite LOBPCG Gram matrix.");
                }
            }
            if (cholesky_rank_transform(gram, m, relative, transform))
            {
                rank = m;
                triangular = true;
            }
            else
            {
                ct::Tensor u(gram);
                ct::Tensor values = host_tensor<double>(m, double(0));
                orthogonal_work = gram.buffer().GetAllocatedBytes() + transform.buffer().GetAllocatedBytes() +
                                  2 * u.buffer().GetAllocatedBytes() + values.buffer().GetAllocatedBytes();
                ct::kernels::lapack_heevd<Wide, ct::DEVICE_CPU> solve;
                solve(m, u.data<Wide>(), m, values.data<double>());
                const double cutoff = relative * std::max(0.0, values.data<double>()[values.NumElements() - 1]);
                if (sx != nullptr && values.data<double>()[0] < -cutoff)
                {
                    throw std::runtime_error("LOBPCG requires a positive overlap metric.");
                }
                int first = 0;
                while (first < m && values.data<double>()[first] <= cutoff)
                {
                    ++first;
                }
                rank = m - first;
                const int64_t scaled_size = static_cast<int64_t>(m) * rank;
                ct::Tensor scaled = host_tensor<Wide>(scaled_size, Wide(0));
                for (int j = 0; j < rank; ++j)
                {
                    const int source = first + j;
                    const double scale = 1 / std::sqrt(values.data<double>()[source]);
                    for (int i = 0; i < m; ++i)
                    {
                        scaled.data<Wide>()[i + j * m] = u.data<Wide>()[i + source * m] * scale;
                    }
                }
                if (rank < m)
                {
                    transform = std::move(scaled);
                }
                else
                {
                    // Preserve U Lambda^-1/2 U^H at full rank, U_r Lambda_r^-1/2 otherwise.
                    const int64_t transform_size = static_cast<int64_t>(m) * m;
                    transform = host_tensor<Wide>(transform_size, Wide(0));
                    const Wide one(1);
                    const Wide zero(0);
                    ModuleBase::gemm_op<Wide, base_device::DEVICE_CPU>()('N',
                        'C',
                        m,
                        m,
                        m,
                        &one,
                        scaled.data<Wide>(),
                        m,
                        u.data<Wide>(),
                        m,
                        &zero,
                        transform.data<Wide>(),
                        m);
                }
            }
        }
        catch (const std::exception&)
        {
            rank = -1;
        }
    }
    ModuleBase::Memory::record("LOBPCG::ortho_work", orthogonal_work, false);
    // Broadcast rank and transform type to select TRMM or GEMM on every rank.
    double metadata[2] = {static_cast<double>(rank), static_cast<double>(triangular)};
    broadcast(metadata, 2);
    rank = static_cast<int>(metadata[0]);
    triangular = metadata[1] != 0.0;
    if (rank <= 0)
    {
        return rank;
    }
    const int64_t transform_size = static_cast<int64_t>(m) * rank;
    transform.resize({transform_size});
    const int count = m * rank;
    broadcast(transform.data<Wide>(), count);
    const T* coefficients = prepare_rotation_coefficients(transform);
    ct::Tensor rotated(t_type_, device_type_, {0});
    if (!triangular)
    {
        const int64_t elements = static_cast<int64_t>(n_dim_) * rank;
        rotated.resize({elements});
    }
    const std::size_t scratch_bytes = rotated.buffer().GetAllocatedBytes();
    if (device_type_ == ct::DeviceType::CpuDevice)
    {
        ModuleBase::Memory::record("LOBPCG::basis_scratch", scratch_bytes, false);
    }
#if defined(__CUDA) || defined(__ROCM)
    else
    {
        ModuleBase::Memory::record_gpu("LOBPCG::basis_scratch", scratch_bytes, false);
    }
#endif
    for (T* image : {x, sx})
    {
        if (image == nullptr)
        {
            continue;
        }
        if (triangular)
        {
            const T one(1);
            ModuleBase::trmm_op<T, Device>()('R', 'U', 'N', 'N', n_dim_, m, &one, coefficients, m, image, n_dim_);
        }
        else
        {
            algebra_.expand_prepared(n_dim_, n_dim_, m, rank, image, coefficients, rotated.data<T>(), zero_);
            const int elements = n_dim_ * rank;
            copy_op_(elements, rotated.data<T>(), 1, image, 1);
        }
    }
    const std::size_t ortho_bytes = gram.buffer().GetAllocatedBytes() + transform.buffer().GetAllocatedBytes();
    ModuleBase::Memory::record("LOBPCG::ortho_work", ortho_bytes, false);
    return rank;
}

template int DiagoLOBPCG<std::complex<float>, base_device::DEVICE_CPU>::rank_orthonormalize(
    const int, std::complex<float>*, std::complex<float>*);
template bool DiagoLOBPCG<std::complex<float>, base_device::DEVICE_CPU>::recover_projected(
    const HSOperator<std::complex<float>, base_device::DEVICE_CPU>&,
    const VectorImages&,
    int&,
    const int,
    ct::Tensor&,
    ct::Tensor&,
    SerialGvWorkspace&);
template int DiagoLOBPCG<std::complex<float>, base_device::DEVICE_CPU>::repair_basis(
    const HSOperator<std::complex<float>, base_device::DEVICE_CPU>&, const int);

template int DiagoLOBPCG<std::complex<double>, base_device::DEVICE_CPU>::rank_orthonormalize(
    const int, std::complex<double>*, std::complex<double>*);
template bool DiagoLOBPCG<std::complex<double>, base_device::DEVICE_CPU>::recover_projected(
    const HSOperator<std::complex<double>, base_device::DEVICE_CPU>&,
    const VectorImages&,
    int&,
    const int,
    ct::Tensor&,
    ct::Tensor&,
    SerialGvWorkspace&);
template int DiagoLOBPCG<std::complex<double>, base_device::DEVICE_CPU>::repair_basis(
    const HSOperator<std::complex<double>, base_device::DEVICE_CPU>&, const int);

#if defined(__CUDA) || defined(__ROCM)

template int DiagoLOBPCG<std::complex<float>, base_device::DEVICE_GPU>::rank_orthonormalize(
    const int, std::complex<float>*, std::complex<float>*);
template bool DiagoLOBPCG<std::complex<float>, base_device::DEVICE_GPU>::recover_projected(
    const HSOperator<std::complex<float>, base_device::DEVICE_GPU>&,
    const VectorImages&,
    int&,
    const int,
    ct::Tensor&,
    ct::Tensor&,
    SerialGvWorkspace&);
template int DiagoLOBPCG<std::complex<float>, base_device::DEVICE_GPU>::repair_basis(
    const HSOperator<std::complex<float>, base_device::DEVICE_GPU>&, const int);

template int DiagoLOBPCG<std::complex<double>, base_device::DEVICE_GPU>::rank_orthonormalize(
    const int, std::complex<double>*, std::complex<double>*);
template bool DiagoLOBPCG<std::complex<double>, base_device::DEVICE_GPU>::recover_projected(
    const HSOperator<std::complex<double>, base_device::DEVICE_GPU>&,
    const VectorImages&,
    int&,
    const int,
    ct::Tensor&,
    ct::Tensor&,
    SerialGvWorkspace&);
template int DiagoLOBPCG<std::complex<double>, base_device::DEVICE_GPU>::repair_basis(
    const HSOperator<std::complex<double>, base_device::DEVICE_GPU>&, const int);
#endif

} // namespace hsolver
