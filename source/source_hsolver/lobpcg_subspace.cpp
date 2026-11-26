#include "diago_lobpcg.h"
#include "lobpcg_internal.h"
#include "dense_eigensolver.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_base/memory_recorder.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#ifdef __MPI
#include "source_base/module_parallel/para_mpi_func.h"
#endif

#include <algorithm>
#include <exception>
#include <type_traits>

namespace hsolver
{
using lobpcg_internal::host_tensor;
using lobpcg_internal::host_tensor_cast_or_view;

// Projected K/G, generalized RR, and the retained X/P projection cache.

/// Form K = B^H H B and G = B^H S B; with a retained [X,P] block, compute only the new W panel.
template <typename T, typename Device>
void DiagoLOBPCG<T, Device>::form_projected_matrices(
    const VectorImages& basis, const int columns, const ProjectionCache* retained, ct::Tensor& k, ct::Tensor& g)
{
    ModuleBase::timer::Scope timer("DiagoLOBPCG", "form_projected_matrices");
    const int keep = retained == nullptr ? 0 : retained->columns;
    const ct::Tensor* saved[] = {
        retained == nullptr ? nullptr : &retained->k, retained == nullptr ? nullptr : &retained->g};
    ct::Tensor* outputs[] = {&k, &g};
    const T* images[] = {basis.hx, basis.sx == nullptr ? basis.x : basis.sx};
    const int added = columns - keep;
    std::size_t local_bytes = 0;
    for (int matrix = 0; matrix < 2; ++matrix)
    {
        ct::Tensor& result = *outputs[matrix];
        const bool nc_gram = matrix == 1 && basis.sx == nullptr && std::is_same<Device, base_device::DEVICE_CPU>::value;
        if (keep == 0)
        {
            if (nc_gram)
            {
                result = algebra_.gram_root_tensor(n_dim_, n_dim_, columns, basis.x);
            }
            else
            {
                result = algebra_.cross_root_tensor(n_dim_, n_dim_, columns, columns, basis.x, images[matrix]);
            }
            continue;
        }
        const int64_t elements = comm_.rank == 0 ? static_cast<int64_t>(columns) * columns : 0;
        result = host_tensor<Wide>(elements, Wide(0));
        if (comm_.rank == 0)
        {
            for (int j = 0; j < keep; ++j)
            {
                std::copy_n(saved[matrix]->template data<Wide>() + j * keep, keep, result.data<Wide>() + j * columns);
            }
        }
        if (added == 0)
        {
            continue;
        }
        ct::Tensor panel = host_tensor<Wide>(0, Wide(0));
        if (nc_gram)
        {
            // Compute W^H W with HERK, then reduce the entire new panel once.
            const int n = n_dim_;
            const int64_t basis_elements = static_cast<int64_t>(n) * columns;
            const ct::Tensor basis_wide = host_tensor_cast_or_view<Wide>(basis.x, basis_elements);
            const Wide* wide_basis = basis_wide.data<Wide>();
            const Wide* w = wide_basis + keep * n;
            panel = host_tensor<Wide>(static_cast<std::size_t>(columns) * added, Wide(0));
            ct::Tensor self = host_tensor<Wide>(static_cast<std::size_t>(added) * added, Wide(0));
            const std::size_t panel_bytes = panel.buffer().GetAllocatedBytes() + self.buffer().GetAllocatedBytes() +
                                            basis_wide.buffer().GetAllocatedBytes();
            local_bytes = std::max(local_bytes, panel_bytes);
            const Wide one(1);
            const Wide zero(0);
            ModuleBase::gemm_op<Wide, base_device::DEVICE_CPU>()(
                'C', 'N', keep, added, n, &one, wide_basis, n, w, n, &zero, panel.data<Wide>(), columns);
            ModuleBase::gram_herk_cpu(n, n, added, w, self.data<Wide>());
            linear_complete_hermitian_upper(added, self.data<Wide>());
            for (int j = 0; j < added; ++j)
            {
                std::copy_n(self.data<Wide>() + j * added, added, panel.data<Wide>() + keep + j * columns);
            }
#ifdef __MPI
            const int panel_size = columns * added;
            Parallel_Common::reduce_data_to_root(panel.data<Wide>(), panel_size, comm_.comm, 0);
#endif
        }
        else
        {
            const T* new_image = images[matrix] + keep * n_dim_;
            panel = algebra_.cross_root_tensor(n_dim_, n_dim_, columns, added, basis.x, new_image);
            if (comm_.rank == 0)
            {
                const std::size_t panel_bytes = panel.buffer().GetAllocatedBytes();
                local_bytes = std::max(local_bytes, panel_bytes);
            }
        }
        if (comm_.rank == 0)
        {
            // Retained blocks are already global; only the new panel was reduced.
            Wide* output = result.data<Wide>();
            std::copy_n(panel.data<Wide>(), panel.NumElements(), output + keep * columns);
            for (int j = 0; j < keep; ++j)
            {
                for (int i = keep; i < columns; ++i)
                {
                    output[i + j * columns] = std::conj(output[j + i * columns]);
                }
            }
        }
    }
    ModuleBase::Memory::record("LOBPCG::projection_local", local_bytes, false);
}

/// Solve K C = G C Theta for the wanted Ritz pairs and share the result; return false for metric failure.
template <typename T, typename Device>
bool DiagoLOBPCG<T, Device>::solve_projected(const ct::Tensor& k,
    const ct::Tensor& g,
    const int dimension,
    const int wanted,
    SerialGvWorkspace& serial_workspace)
{
    ModuleBase::timer::Scope function_timer("DiagoLOBPCG", "solve_projected");
    error_message_.clear();
    ct::Tensor values = host_tensor<double>(dimension, double(0));
    const int64_t coefficient_count = static_cast<int64_t>(dimension) * wanted;
    ct::Tensor coefficients(ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::CpuDevice, {0});
    if (std::is_same<T, Wide>::value)
    {
        // Reuse ritz_coefficients_ as the solve and gather/broadcast output buffer.
        coefficients = host_tensor_cast_or_view<Wide>(ritz_coefficients_.data<T>(), coefficient_count);
    }
    else
    {
        coefficients = host_tensor<Wide>(coefficient_count, Wide(0));
    }
    const std::size_t output_bytes = values.buffer().GetAllocatedBytes() + coefficients.buffer().GetAllocatedBytes();
    ModuleBase::Memory::record("LOBPCG::rr_output", output_bytes, false);
    // Root-only backends publish one outcome before any rank returns/exits:
    // 1 = success, 0 = non-positive metric, -1 = fatal backend failure.
    double result = 1;
    const bool gpu_rr = std::is_same<Device, base_device::DEVICE_GPU>::value && dense::GpuDenseEigensolver::available();
    try
    {
        if (para_rr_method_ == 0)
        {
            if (comm_.rank == 0)
            {
                if (gpu_rr)
                {
                    const bool positive = gpu_rr_solver_.solve_generalized(k.data<Wide>(),
                        g.data<Wide>(),
                        dimension,
                        wanted,
                        values.data<double>(),
                        coefficients.data<Wide>());
                    result = positive ? 1 : 0;
                    if (!positive)
                    {
                        error_message_ = "CUDA generalized RR: overlap is not positive definite.";
                    }
                }
                else
                {
                    const int info = serial_workspace.solve(k.data<Wide>(),
                        g.data<Wide>(),
                        dimension,
                        wanted,
                        values.data<double>(),
                        coefficients.data<Wide>());
                    if (info != 0)
                    {
                        error_message_ = "Serial HEGVX failed, info=" + std::to_string(info);
                        result = info > dimension ? 0 : -1;
                    }
                }
            }
        }
        else if (para_rr_method_ == 2)
        {
#ifdef __MPI
            rr_solver_.solve_generalized(k.data<Wide>(),
                g.data<Wide>(),
                dimension,
                wanted,
                values.data<double>(),
                coefficients.data<Wide>(),
                comm_.comm,
                para_rr_bs_);
#else
            ModuleBase::WARNING_QUIT("DiagoLOBPCG::solve_projected", "Parallel LOBPCG RR requires an MPI build.");
#endif
        }
        else
        {
            ModuleBase::WARNING_QUIT("DiagoLOBPCG::solve_projected", "LOBPCG supports diag_subspace=0 or 2. Please select 0 or 2.");
        }
    }
    catch (const std::exception& error)
    {
        error_message_ = error.what();
        if (para_rr_method_ == 0)
        {
            result = -1; // All ranks must receive the root failure before exiting.
        }
        else
        {
            const std::string reason(error.what());
            const bool recoverable = reason.find("info=16") != std::string::npos ||
                                     reason.find("not positive definite") != std::string::npos;
            result = recoverable ? 0 : -1;
        }
    }
    if (para_rr_method_ == 0)
    {
        // Only root solves the local RR; every rank records the same selected
        // backend label, including failures and zero capacity on nonowners.
        if (gpu_rr)
        {
            ModuleBase::Memory::record("RR::gpu_work", 0, false);
#if defined(__CUDA) || defined(__ROCM)
            const std::size_t bytes = comm_.rank == 0 ? gpu_rr_solver_.storage_bytes() : 0;
            ModuleBase::Memory::record_gpu("RR::gpu_work", bytes, false);
#endif
        }
        else
        {
            const std::size_t bytes = serial_workspace.storage_bytes();
            ModuleBase::Memory::record("RR::serial_gv", bytes, false);
        }
        broadcast(&result, 1);
    }
    if (result <= 0)
    {
#ifdef __MPI
        if (para_rr_method_ == 0 && comm_.nproc > 1)
        {
            const Parallel::ParaWorld world = Parallel::ParaWorld::make_mpi("lobpcg_rr", comm_.comm);
            Parallel::bcast_string(error_message_, world, 0);
        }
#endif
        error_message_ += " (dimension=" + std::to_string(dimension) + ", wanted=" + std::to_string(wanted) + ")";
        if (result < 0)
        {
            ModuleBase::WARNING_QUIT("DiagoLOBPCG::solve_projected", error_message_);
        }
        return false;
    }
    if (comm_.rank == 0)
    {
        for (int i = 0; i < wanted; ++i)
        {
            e_red_.data<Real>()[i] = static_cast<Real>(values.data<double>()[i]);
        }
    }
    broadcast(e_red_.data<Real>(), wanted);
    broadcast(coefficients.data<Wide>(), dimension * wanted);
    if (!std::is_same<T, Wide>::value)
    {
        for (int i = 0; i < dimension * wanted; ++i)
        {
            ritz_coefficients_.data<T>()[i] = static_cast<T>(coefficients.data<Wide>()[i]);
        }
    }
    return true;
}

/// Cache projections for the new [X,P] using coefficient-space products and the Ritz XX identities.
template <typename T, typename Device>
void DiagoLOBPCG<T, Device>::retain_pair(const ct::Tensor& k,
    const ct::Tensor& g,
    const int rows,
    const ct::Tensor& x_coefficients,
    const int first_x,
    const int x_columns,
    const int p_columns,
    const Real* ritz_values,
    const int zero_prefix,
    ProjectionCache& cache)
{
    ModuleBase::timer::Scope function_timer("DiagoLOBPCG", "retain_pair");
    cache.columns = x_columns + p_columns;
    if (comm_.rank == 0)
    {
        const Wide* cx = x_coefficients.data<Wide>() + first_x * rows;
        const T* cp_t = u_p_.data<T>();
        const int keep = x_columns + p_columns;
        const int64_t packed_size = static_cast<int64_t>(rows) * keep;
        cache.structural_coefficients.resize({packed_size});
        Wide* packed = cache.structural_coefficients.template data<Wide>();
        std::copy_n(cx, rows * x_columns, packed);
        for (int i = 0; i < rows * p_columns; ++i)
        {
            packed[rows * x_columns + i] = static_cast<Wide>(cp_t[i]);
        }
        if (!std::is_same<Device, base_device::DEVICE_CPU>::value && p_columns > 0)
        {
            linear_buffer<Wide, Device>(&cache.device_coefficients, packed_size);
            using Upload = base_device::memory::synchronize_memory_op<Wide, Device, base_device::DEVICE_CPU>;
            Upload()(cache.device_coefficients.template data<Wide>(), packed, packed_size);
        }
        const ct::Tensor units = host_tensor<Real>(x_columns, static_cast<Real>(1));
        const ct::Tensor* matrices[] = {&g, &k};
        ct::Tensor* outputs[] = {&cache.g, &cache.k};
        const Real* diagonals[] = {units.template data<Real>(), ritz_values};
        for (int index = 0; index < 2; ++index)
        {
            const ct::Tensor& matrix = *matrices[index];
            ct::Tensor& retained_matrix = *outputs[index];
            const Real* diagonal = diagonals[index];
            const int nonzero_rows = rows - zero_prefix;
            const int64_t retained_elements = static_cast<int64_t>(keep) * keep;
            retained_matrix.resize({retained_elements});
            retained_matrix.zero();
            for (int i = 0; i < x_columns; ++i)
            {
                retained_matrix.template data<Wide>()[i + i * keep] = diagonal[i];
            }
            if (p_columns == 0)
            {
                continue;
            }
            const int64_t product_size = static_cast<int64_t>(rows) * p_columns;
            const int64_t result_size = static_cast<int64_t>(keep) * p_columns;
            Wide* xp = retained_matrix.template data<Wide>() + x_columns * keep;
            const Wide* matrix_data = matrix.data<Wide>();
            const Wide* coefficients = cache.structural_coefficients.template data<Wide>();
            Wide* product = nullptr;
            Wide* result = xp;
            if (std::is_same<Device, base_device::DEVICE_CPU>::value)
            {
                cache.h_product.resize({product_size});
                product = cache.h_product.template data<Wide>();
            }
            else
            {
                // Host projections stay authoritative; only the GEMMs use Device.
                const int64_t matrix_size = matrix.NumElements();
                linear_buffer<Wide, Device>(&cache.device_matrix, matrix_size);
                linear_buffer<Wide, Device>(&cache.device_product, product_size);
                linear_buffer<Wide, Device>(&cache.device_result, result_size);
                using Upload = base_device::memory::synchronize_memory_op<Wide, Device, base_device::DEVICE_CPU>;
                Wide* device_data = cache.device_matrix.template data<Wide>();
                Upload()(device_data, matrix_data, matrix_size);
                matrix_data = device_data;
                coefficients = cache.device_coefficients.template data<Wide>();
                product = cache.device_product.template data<Wide>();
                result = cache.device_result.template data<Wide>();
            }
            const Wide one(1);
            const Wide zero(0);
            const Wide* panel = matrix_data + zero_prefix * rows;
            const Wide* cp = coefficients + x_columns * rows + zero_prefix;
            // Form Q C_P, then C^H (Q C_P).
            ModuleBase::gemm_op<Wide, Device>()(
                'N', 'N', rows, p_columns, nonzero_rows, &one, panel, rows, cp, rows, &zero, product, rows);
            ModuleBase::gemm_op<Wide, Device>()(
                'C', 'N', keep, p_columns, rows, &one, coefficients, rows, product, rows, &zero, result, keep);
            if (!std::is_same<Device, base_device::DEVICE_CPU>::value)
            {
                using Download = base_device::memory::synchronize_memory_op<Wide, base_device::DEVICE_CPU, Device>;
                Download()(xp, result, result_size);
            }
            for (int j = x_columns; j < keep; ++j)
            {
                for (int i = 0; i < x_columns; ++i)
                {
                    retained_matrix.template data<Wide>()[j + i * keep] =
                        std::conj(retained_matrix.template data<Wide>()[i + j * keep]);
                }
            }
        }
    }
    const std::size_t retained_bytes = cache.k.buffer().GetAllocatedBytes() + cache.g.buffer().GetAllocatedBytes();
    ModuleBase::Memory::record("LOBPCG::retained", retained_bytes, false);
}

template <typename T, typename Device>
DiagoLOBPCG<T, Device>::ProjectionCache::~ProjectionCache()
{
    std::size_t host = h_product.buffer().GetAllocatedBytes() + structural_coefficients.buffer().GetAllocatedBytes();
    std::size_t device = 0;
    const ct::Tensor* buffers[] = {&device_matrix, &device_coefficients, &device_product, &device_result};
    for (const ct::Tensor* buffer : buffers)
    {
        if (buffer->device_type() == ct::DeviceType::CpuDevice)
        {
            host += buffer->buffer().GetAllocatedBytes();
        }
        else
        {
            device += buffer->buffer().GetAllocatedBytes();
        }
    }
    ModuleBase::Memory::record("LOBPCG::cache_work", host, false);
#if defined(__CUDA) || defined(__ROCM)
    ModuleBase::Memory::record_gpu("LOBPCG::cache_work", device, false);
#endif
}

template DiagoLOBPCG<std::complex<float>, base_device::DEVICE_CPU>::ProjectionCache::~ProjectionCache();
template void DiagoLOBPCG<std::complex<float>, base_device::DEVICE_CPU>::form_projected_matrices(
    const VectorImages&, const int, const ProjectionCache*, ct::Tensor&, ct::Tensor&);
template bool DiagoLOBPCG<std::complex<float>, base_device::DEVICE_CPU>::solve_projected(
    const ct::Tensor&, const ct::Tensor&, const int, const int, SerialGvWorkspace&);
template void DiagoLOBPCG<std::complex<float>, base_device::DEVICE_CPU>::retain_pair(const ct::Tensor&,
    const ct::Tensor&,
    const int,
    const ct::Tensor&,
    const int,
    const int,
    const int,
    const Real*,
    const int,
    ProjectionCache&);

template DiagoLOBPCG<std::complex<double>, base_device::DEVICE_CPU>::ProjectionCache::~ProjectionCache();
template void DiagoLOBPCG<std::complex<double>, base_device::DEVICE_CPU>::form_projected_matrices(
    const VectorImages&, const int, const ProjectionCache*, ct::Tensor&, ct::Tensor&);
template bool DiagoLOBPCG<std::complex<double>, base_device::DEVICE_CPU>::solve_projected(
    const ct::Tensor&, const ct::Tensor&, const int, const int, SerialGvWorkspace&);
template void DiagoLOBPCG<std::complex<double>, base_device::DEVICE_CPU>::retain_pair(const ct::Tensor&,
    const ct::Tensor&,
    const int,
    const ct::Tensor&,
    const int,
    const int,
    const int,
    const Real*,
    const int,
    ProjectionCache&);

#if defined(__CUDA) || defined(__ROCM)

template DiagoLOBPCG<std::complex<float>, base_device::DEVICE_GPU>::ProjectionCache::~ProjectionCache();
template void DiagoLOBPCG<std::complex<float>, base_device::DEVICE_GPU>::form_projected_matrices(
    const VectorImages&, const int, const ProjectionCache*, ct::Tensor&, ct::Tensor&);
template bool DiagoLOBPCG<std::complex<float>, base_device::DEVICE_GPU>::solve_projected(
    const ct::Tensor&, const ct::Tensor&, const int, const int, SerialGvWorkspace&);
template void DiagoLOBPCG<std::complex<float>, base_device::DEVICE_GPU>::retain_pair(const ct::Tensor&,
    const ct::Tensor&,
    const int,
    const ct::Tensor&,
    const int,
    const int,
    const int,
    const Real*,
    const int,
    ProjectionCache&);

template DiagoLOBPCG<std::complex<double>, base_device::DEVICE_GPU>::ProjectionCache::~ProjectionCache();
template void DiagoLOBPCG<std::complex<double>, base_device::DEVICE_GPU>::form_projected_matrices(
    const VectorImages&, const int, const ProjectionCache*, ct::Tensor&, ct::Tensor&);
template bool DiagoLOBPCG<std::complex<double>, base_device::DEVICE_GPU>::solve_projected(
    const ct::Tensor&, const ct::Tensor&, const int, const int, SerialGvWorkspace&);
template void DiagoLOBPCG<std::complex<double>, base_device::DEVICE_GPU>::retain_pair(const ct::Tensor&,
    const ct::Tensor&,
    const int,
    const ct::Tensor&,
    const int,
    const int,
    const int,
    const Real*,
    const int,
    ProjectionCache&);
#endif

} // namespace hsolver
