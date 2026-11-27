#include "diago_lobpcg.h"
#include "lobpcg_internal.h"
#include "dense_eigensolver.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_base/memory_recorder.h"
#include "source_base/module_container/ATen/kernels/memory.h"
#include "source_base/parallel_device.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <type_traits>
#include <utility>

namespace hsolver
{
using lobpcg_internal::host_tensor;
using lobpcg_internal::host_tensor_cast_or_view;

/// Store the caller configuration and allocate Device wavefunction and host RR workspaces.
template <typename T, typename Device>
DiagoLOBPCG<T, Device>::DiagoLOBPCG(const Real* precondition,
    const int nband,
    const int ndim,
    const int nmax,
    const LobpcgBasisLayout& layout,
    const diag_comm_info& comm,
    const int para_rr_method,
    const int para_rr_block_size,
    const bool use_overlap,
    const LobpcgLocking locking,
    dense::ParallelDenseEigensolver& rr_solver,
    dense::GpuDenseEigensolver& gpu_rr_solver)
    : n_band_(nband), n_dim_(ndim), n_max_(nmax), layout_(layout), comm_(comm), para_rr_method_(para_rr_method),
      para_rr_bs_(para_rr_block_size), use_overlap_(use_overlap), locking_(locking), rr_solver_(rr_solver),
      gpu_rr_solver_(gpu_rr_solver), algebra_(comm)
{
    if (n_max_ < n_band_)
    {
        ModuleBase::WARNING_QUIT("DiagoLOBPCG", "LOBPCG block must contain all requested bands.");
    }
    r_type_ = ct::DataTypeToEnum<Real>::value;
    t_type_ = ct::DataTypeToEnum<T>::value;
    device_type_ = ct::DeviceTypeToEnum<Device>::value;

    prec_ = ct::Tensor(r_type_, device_type_, {n_dim_});
    syncmem_var_h2d_op()(prec_.data<Real>(), precondition, n_dim_);
    allocate_workspace();
}

/// Record the final owned algebra workspace capacity before its buffers are released.
template <typename T, typename Device>
DiagoLOBPCG<T, Device>::~DiagoLOBPCG()
{
    const std::size_t host_bytes = algebra_.storage_bytes(ct::DeviceType::CpuDevice);
    ModuleBase::Memory::record("LOBPCG::algebra", host_bytes, false);
#if defined(__CUDA) || defined(__ROCM)
    const std::size_t device_bytes = algebra_.storage_bytes(ct::DeviceType::GpuDevice);
    ModuleBase::Memory::record_gpu("LOBPCG::algebra", device_bytes, false);
#endif
}

// Entry and generalized iteration.
/// Solve H X = S X Lambda by generalized RR in [X,P,W], with prefix locking and failure-only repair.
template <typename T, typename Device>
bool DiagoLOBPCG<T, Device>::diag(const HSOperator<T, Device>& op,
    T* psi_in,
    const int ld_psi,
    Real* eigenvalue_in,
    const double tolerance,
    const int max_iter,
    const std::vector<double>& ethr_band)
{
    ModuleBase::timer::Scope timer("DiagoLOBPCG", "diag");
    iterations_ = 0;
    error_message_.clear();

    double global_dim = static_cast<double>(layout_.npw) * layout_.npol;
#ifdef __MPI
    if (comm_.nproc > 1)
    {
        Parallel_Common::reduce_data(&global_dim, 1, comm_.comm);
    }
#endif
    // A full [X,P,W] block needs at least 3*b physical rows; trim only optional buffer bands.
    if (3 * n_max_ > global_dim)
    {
        if (n_band_ > static_cast<int>(global_dim) / 3)
        {
            ModuleBase::WARNING_QUIT("DiagoLOBPCG::diag", "The requested bands cannot form a full-rank LOBPCG space.");
        }
        n_max_ = static_cast<int>(global_dim) / 3;
        allocate_workspace();
    }

    // Register optional labels in a common order across pools with different iteration paths.
    const char* optional_work[] = {"LOBPCG::Cp",
        "LOBPCG::retained",
        "LOBPCG::cache_work",
        "LOBPCG::p_work",
        "LOBPCG::recovery_reference",
        "LOBPCG::ortho_work",
        "LOBPCG::basis_scratch"};
    for (const char* name : optional_work)
    {
        ModuleBase::Memory::record(name, 0, false);
    }
#if defined(__CUDA) || defined(__ROCM)
    ModuleBase::Memory::record_gpu("LOBPCG::cache_work", 0, false);
    ModuleBase::Memory::record_gpu("LOBPCG::basis_scratch", 0, false);
#endif

    // 1. Prepare X, then obtain its initial Ritz pairs and residual.
    initialize(psi_in, ld_psi);
    SerialGvWorkspace serial_workspace;
    bool initially_exact = true;
    { // Initial RR: release temporary K/G, coefficient views and norms before the main loop.
        const VectorImages initial = basis_images(0);
        if (use_overlap_)
        {
            op.spsi(initial.x, initial.sx, n_dim_, n_max_);
            project_padding(initial.sx, n_max_);
        }
        op.hpsi(initial.x, initial.hx, n_dim_, n_max_);
        project_padding(initial.hx, n_max_);
        ct::Tensor k = host_tensor<Wide>(0, Wide(0));
        ct::Tensor g = host_tensor<Wide>(0, Wide(0));
        form_projected_matrices(initial, n_max_, nullptr, k, g);
        ModuleBase::Memory::record("LOBPCG::K", k.buffer().GetAllocatedBytes(), false);
        ModuleBase::Memory::record("LOBPCG::G", g.buffer().GetAllocatedBytes(), false);
        if (!solve_projected(k, g, n_max_, n_max_, serial_workspace))
        {
            ModuleBase::WARNING_QUIT("DiagoLOBPCG::diag", error_message_);
        }
        std::copy_n(e_red_.data<Real>(), n_max_, eig_.data<Real>());
        const int coefficient_count = n_max_ * n_max_;
        const ct::Tensor coefficients = host_tensor_cast_or_view<Wide>(ritz_coefficients_.data<T>(), coefficient_count);
        ModuleBase::Memory::record("LOBPCG::Cx", coefficients.buffer().GetAllocatedBytes(), false);
        const T* prepared = prepare_rotation_coefficients(coefficients);
        rotate_images(initial, n_max_, n_max_, prepared, initial, evec_.data<T>());
        // R = HX - SX*Lambda. Exclude roundoff-only directions before normalizing W;
        // if every requested residual vanishes, these Ritz pairs already solve H/S.
        form_residual(initial, n_max_, eig_.data<Real>(), residual_.data<T>());
        ct::Tensor residual_norm_squared = host_tensor<Real>(n_max_, Real(0));
        for (int i = 0; i < n_max_; ++i)
        {
            const T* column = residual_.data<T>() + i * n_dim_;
            const Real norm = nrm2_op_(n_dim_, column, 1);
            residual_norm_squared.data<Real>()[i] = norm * norm;
        }
#ifdef __MPI
        if (comm_.nproc > 1)
        {
            Parallel_Common::reduce_data(residual_norm_squared.data<Real>(), n_max_, comm_.comm);
        }
#endif
        for (int i = 0; i < n_max_; ++i)
        {
            const Real norm = std::sqrt(residual_norm_squared.data<Real>()[i]);
            const Real scale = std::max(static_cast<Real>(1), std::abs(eig_.data<Real>()[i]));
            const Real roundoff = static_cast<Real>(32) * std::numeric_limits<Real>::epsilon() * scale;
            if (i < n_band_)
            {
                initially_exact = initially_exact && norm <= roundoff;
            }
            if (norm <= roundoff)
            {
                T* column = residual_.data<T>() + i * n_dim_;
                setmem_complex_op()(column, zero_, n_dim_);
            }
        }
    } // End initial RR and roundoff-residual filtering.
    // All requested residuals are roundoff-sized: return these Ritz pairs without constructing W.
    if (initially_exact)
    {
        copy_op_(n_dim_ * n_max_, space_.data<T>(), 1, x_new_.data<T>(), 1);
        write_results(psi_in, ld_psi, eigenvalue_in);
        return true;
    }
    const bool soft = locking_ == LobpcgLocking::Soft;
    int first_active = 0;
    int first_w = n_max_;
    int n_w = build_search(op, 0, n_max_, n_max_);
    if (n_w <= 0)
    {
        error_message_ = "No nonzero initial search directions.";
        ModuleBase::WARNING_QUIT("DiagoLOBPCG::diag", error_message_);
    }
    int working_columns = n_max_ + n_w;

    ct::Tensor eig_prev(eig_);
    ct::Tensor tol_band = host_tensor<Real>(n_max_, static_cast<Real>(tolerance));
    for (int i = 0; i < n_band_ && i < static_cast<int>(ethr_band.size()); ++i)
    {
        tol_band.data<Real>()[i] = static_cast<Real>(ethr_band[i]);
    }

    // H is fixed only within this call. Start with a full projection;
    // normal updates reuse the retained X/P blocks of K and G.
    ProjectionCache cache;

    for (int iter = 0; iter < max_iter; ++iter)
    {
        iterations_ = iter + 1;
        // 2. Apply H to new W and assemble the current K/G pair.
        if (n_w > 0)
        {
            const VectorImages search = basis_images(first_w);
            op.hpsi(search.x, search.hx, n_dim_, n_w);
            project_padding(search.hx, n_w);
        }
        const int rr_off = soft ? 0 : first_active;
        int rr_dim = working_columns - rr_off;
        const int x_columns = n_max_ - rr_off;
        const VectorImages rr = basis_images(rr_off);
        ct::Tensor projected = host_tensor<Wide>(0, Wide(0));
        ct::Tensor gram = host_tensor<Wide>(0, Wide(0));
        form_projected_matrices(rr, rr_dim, &cache, projected, gram);
        ModuleBase::Memory::record("LOBPCG::K", projected.buffer().GetAllocatedBytes(), false);
        ModuleBase::Memory::record("LOBPCG::G", gram.buffer().GetAllocatedBytes(), false);
        // 3. Solve generalized RR; only actual failure enters recovery.
        if (!solve_projected(projected, gram, rr_dim, x_columns, serial_workspace))
        {
            // Release stale projections while retaining reusable scratch capacity.
            cache.columns = 0;
            cache.k = host_tensor<Wide>(0, Wide(0));
            cache.g = host_tensor<Wide>(0, Wide(0));
            if (!recover_projected(op, rr, working_columns, x_columns, projected, gram, serial_workspace))
            {
                ModuleBase::WARNING_QUIT("DiagoLOBPCG::diag", error_message_);
            }
            rr_dim = working_columns - rr_off;
        }
        // 4. Rotate X/SX and decide convergence before any next-step work.
        Real* active_eigenvalues = eig_.data<Real>() + rr_off;
        std::copy_n(e_red_.data<Real>(), x_columns, active_eigenvalues);

        const int active_offset = rr_off * n_dim_;
        T* active_sx = use_overlap_ ? sx_new_.data<T>() + active_offset : nullptr;
        const VectorImages active = {x_new_.data<T>() + active_offset, hx_new_.data<T>() + active_offset, active_sx};
        const int coefficient_count = rr_dim * x_columns;
        const ct::Tensor coefficients = host_tensor_cast_or_view<Wide>(ritz_coefficients_.data<T>(), coefficient_count);
        ModuleBase::Memory::record("LOBPCG::Cx", coefficients.buffer().GetAllocatedBytes(), false);
        // Rotate X/SX now; only the needed HX columns are rotated below.
        const VectorImages ritz_source = {rr.x, nullptr, rr.sx};
        // Share the prepared coefficients across X/SX and HX rotations.
        const T* prepared = prepare_rotation_coefficients(coefficients);
        rotate_images(ritz_source, rr_dim, x_columns, prepared, active, nullptr);

        if (lock_converged(eig_prev, tol_band, first_active))
        {
            write_results(psi_in, ld_psi, eigenvalue_in);
            return true;
        }

        const int active_columns = n_max_ - first_active;
        if (iter + 1 == max_iter)
        {
            // Return the current Ritz pairs when the iteration budget is exhausted.
            break;
        }
        // 5. Form the active residual, rotating only the required HX.
        // Newly locked columns never enter the next residual/search block.
        const int search_offset = first_active * n_dim_;
        T* search_sx = use_overlap_ ? sx_new_.data<T>() + search_offset : nullptr;
        const VectorImages search_ritz = {
            x_new_.data<T>() + search_offset, hx_new_.data<T>() + search_offset, search_sx};
        const int hx_columns = soft ? x_columns : active_columns;
        const int first_hx_column = soft ? 0 : x_columns - active_columns;
        const VectorImages hx_target = soft ? active : search_ritz;
        const T* hx_coefficients = prepared + first_hx_column * rr_dim;
        algebra_.expand_prepared(n_dim_, n_dim_, rr_dim, hx_columns, rr.hx, hx_coefficients, hx_target.hx, zero_);
        const Real* search_eigenvalues = eig_.data<Real>() + first_active;
        T* active_residual = residual_.data<T>() + first_active * n_dim_;
        form_residual(search_ritz, active_columns, search_eigenvalues, active_residual);

        // 6. Build P/HP/SP and their consistently scaled coefficients.
        const int n_p = build_direction(rr, rr_dim, x_columns, active_columns, coefficients, gram);
        // 7. Retain the new X/P projections before overwriting the old basis.
        const Real* retained_values = active_eigenvalues + first_hx_column;
        // Expansion zeroes these rows exactly; normalization preserves them.
        // New locks remove this structure.
        const int p_zero_prefix = x_columns == active_columns ? x_columns : 0;
        retain_pair(projected,
            gram,
            rr_dim,
            coefficients,
            first_hx_column,
            hx_columns,
            n_p,
            retained_values,
            p_zero_prefix,
            cache);
        first_w = n_max_ + n_p;

        const int x_elements = n_dim_ * x_columns;
        copy_op_(x_elements, active.x, 1, rr.x, 1);
        // Soft locking retains every X in RR; their complete H images must
        // remain valid even while only the active suffix generates residuals.
        const int hx_elements = n_dim_ * hx_columns;
        const int next_hx_first = soft ? 0 : first_active;
        const VectorImages next_active = basis_images(next_hx_first);
        copy_op_(hx_elements, hx_target.hx, 1, next_active.hx, 1);
        if (use_overlap_)
        {
            copy_op_(x_elements, active.sx, 1, rr.sx, 1);
        }
        // 8. Prepare W for the next iteration. Its H image is deferred.
        n_w = build_search(op, first_active, active_columns, first_w);
        if (n_w == 0 && n_p == 0)
        {
            error_message_ = "No nonzero search or history directions remain.";
            write_results(psi_in, ld_psi, eigenvalue_in);
            ModuleBase::WARNING_QUIT("DiagoLOBPCG::diag", error_message_);
        }
        working_columns = n_max_ + n_p + n_w;
    }
    iterations_ = max_iter;
    for (int i = 0; i < n_band_; ++i)
    {
        if (!std::isfinite(eig_.data<Real>()[i]))
        {
            error_message_ = "Non-finite target Ritz eigenvalue.";
            ModuleBase::WARNING_QUIT("DiagoLOBPCG::diag", error_message_);
        }
    }
    write_results(psi_in, ld_psi, eigenvalue_in);
    return false;
}

// Main algorithm operations, in execution order.
/// Copy the input X, seed optional buffer bands, and zero inactive rows in each spinor component.
template <typename T, typename Device>
void DiagoLOBPCG<T, Device>::initialize(T* psi_in, const int ld_psi)
{
    ModuleBase::timer::Scope function_timer("DiagoLOBPCG", "initialize");
    T* x = space_.data<T>();
    syncmem_complex_2d_op()(x, n_dim_, psi_in, ld_psi, n_dim_, n_band_);
    const int stride = n_dim_ / layout_.npol;
    for (int j = n_band_; j < n_max_; ++j)
    {
        ct::Tensor seed = host_tensor<T>(n_dim_, zero_);
        std::seed_seq seed_key{comm_.rank, j};
        std::mt19937 generator(seed_key);
        std::uniform_real_distribution<double> uniform(-0.5, 0.5);
        for (int component = 0; component < layout_.npol; ++component)
        {
            for (int ig = 0; ig < layout_.npw; ++ig)
            {
                seed.data<T>()[component * stride + ig] = static_cast<T>(uniform(generator));
            }
        }
        syncmem_complex_h2d_op()(x + j * n_dim_, seed.data<T>(), n_dim_);
    }
    // Seeded buffer columns already have zero padding; only the caller's columns need clearing.
    project_padding(x, n_band_);
}

/// Form R = HX - SX*Lambda for the supplied block; a null SX means S = I.
template <typename T, typename Device>
void DiagoLOBPCG<T, Device>::form_residual(
    const VectorImages& vectors, const int count, const Real* eigenvalues, T* residual)
{
    ModuleBase::timer::Scope function_timer("DiagoLOBPCG", "form_residual");
    const int elements = n_dim_ * count;
    copy_op_(elements, vectors.hx, 1, residual, 1);
    const T* metric_x = vectors.sx == nullptr ? vectors.x : vectors.sx;
    for (int i = 0; i < count; ++i)
    {
        const T alpha = static_cast<T>(-eigenvalues[i]);
        const T* column = metric_x + i * n_dim_;
        T* r = residual + i * n_dim_;
        axpy_op_(n_dim_, &alpha, column, 1, r, 1);
    }
}

/// Precondition the active residual, project out the locked prefix, then normalize and compact W.
template <typename T, typename Device>
int DiagoLOBPCG<T, Device>::build_search(
    const HSOperator<T, Device>& op, const int locked, const int count, const int first_w)
{
    ModuleBase::timer::Scope function_timer("DiagoLOBPCG", "build_search");
    const int n = n_dim_;
    const VectorImages block = basis_images(first_w);
    const T* residual = residual_.data<T>() + locked * n;
    const Real* values = eig_.data<Real>() + locked;
    const int elements = n * count;
    copy_op_(elements, residual, 1, block.x, 1);
    syncmem_var_h2d_op()(eig_device_.data<Real>(), values, count);
    // Residual padding is zero; division by the positive preconditioner preserves it.
    precondition_op_(n_dim_, space_.data<T>(), first_w, count, prec_.data<Real>(), eig_device_.data<Real>());
    if (use_overlap_)
    {
        op.spsi(block.x, block.sx, n_dim_, count);
        project_padding(block.sx, count);
    }
    if (locked > 0)
    {
        const VectorImages reference = basis_images(0);
        const T* metric = block.sx == nullptr ? block.x : block.sx;
        const ct::Tensor cross = algebra_.cross_tensor(n, n, locked, count, reference.x, metric);
        const T* prepared = prepare_rotation_coefficients(cross);
        const T minus_one(-1);
        algebra_.expand_prepared(n, n, locked, count, reference.x, prepared, block.x, minus_one);
        if (block.sx != nullptr)
        {
            algebra_.expand_prepared(n, n, locked, count, reference.sx, prepared, block.sx, minus_one);
        }
    }
    const VectorImages search = {block.x, nullptr, block.sx};
    return normalize_columns(search, count, nullptr, 0);
}

/// Lock the leading bands satisfying the eigenvalue-change tolerance; soft locking can reactivate them.
template <typename T, typename Device>
bool DiagoLOBPCG<T, Device>::lock_converged(ct::Tensor& eig_prev, const ct::Tensor& tol_band, int& locked)
{
    ModuleBase::timer::Scope function_timer("DiagoLOBPCG", "lock_converged");
    const int previous_locked = locked;
    locked = 0;
    const bool soft = locking_ == LobpcgLocking::Soft;
    bool prefix = true;
    for (int i = 0; i < n_max_; ++i)
    {
        const double delta = std::abs(eig_.data<Real>()[i] - eig_prev.data<Real>()[i]);
        const double band_tolerance = tol_band.data<Real>()[i];
        const bool band_ok = delta < band_tolerance;
        // Hard retains locked vectors; soft lets their prefix reactivate.
        const bool previous = !soft && i < previous_locked;
        prefix = prefix && (previous || band_ok);
        locked += prefix;
        eig_prev.data<Real>()[i] = eig_.data<Real>()[i];
    }
    return locked >= n_band_;
}

/// Form P from the non-X Ritz coefficients, removing new locked components before normalization.
template <typename T, typename Device>
int DiagoLOBPCG<T, Device>::build_direction(const VectorImages& basis,
    const int dimension,
    const int q,
    const int active,
    const ct::Tensor& x_coefficients,
    const ct::Tensor& g)
{
    ModuleBase::timer::Scope function_timer("DiagoLOBPCG", "build_direction");
    T* direction = u_p_.data<T>();
    const int coefficient_count = dimension * active;
    std::size_t p_bytes = 0;
    {
        const int locked_columns = q - active;
        const Wide* coefficients = x_coefficients.data<Wide>();
        Wide* projected_direction = reinterpret_cast<Wide*>(direction);
        ct::Tensor direction_wide = host_tensor<Wide>(0, Wide(0));
        if (!std::is_same<T, Wide>::value && comm_.rank == 0)
        {
            const int64_t direction_size = static_cast<int64_t>(dimension) * active;
            direction_wide.resize({direction_size});
            projected_direction = direction_wide.data<Wide>();
            p_bytes = direction_wide.buffer().GetAllocatedBytes();
        }
        if (comm_.rank == 0)
        {
            const Wide* active_coefficients = coefficients + locked_columns * dimension;
            std::copy_n(active_coefficients, dimension * active, projected_direction);
            for (int j = 0; j < active; ++j)
            {
                std::fill_n(projected_direction + j * dimension, q, Wide(0));
            }
            if (locked_columns > 0)
            {
                ct::Tensor left = host_tensor<Wide>(static_cast<std::size_t>(locked_columns) * dimension, Wide(0));
                ct::Tensor cross = host_tensor<Wide>(static_cast<std::size_t>(locked_columns) * active, Wide(0));
                p_bytes = left.buffer().GetAllocatedBytes() + cross.buffer().GetAllocatedBytes() +
                          direction_wide.buffer().GetAllocatedBytes();
                const Wide one(1);
                const Wide zero(0);
                const Wide minus_one(-1);
                ModuleBase::gemm_op<Wide, base_device::DEVICE_CPU>()('C',
                    'N',
                    locked_columns,
                    dimension,
                    dimension,
                    &one,
                    coefficients,
                    dimension,
                    g.data<Wide>(),
                    dimension,
                    &zero,
                    left.data<Wide>(),
                    locked_columns);
                ModuleBase::gemm_op<Wide, base_device::DEVICE_CPU>()('N',
                    'N',
                    locked_columns,
                    active,
                    dimension,
                    &one,
                    left.data<Wide>(),
                    locked_columns,
                    projected_direction,
                    dimension,
                    &zero,
                    cross.template data<Wide>(),
                    locked_columns);
                ModuleBase::gemm_op<Wide, base_device::DEVICE_CPU>()('N',
                    'N',
                    dimension,
                    active,
                    locked_columns,
                    &minus_one,
                    coefficients,
                    dimension,
                    cross.template data<Wide>(),
                    locked_columns,
                    &one,
                    projected_direction,
                    dimension);
            }
        }
        if (!std::is_same<T, Wide>::value && comm_.rank == 0)
        {
            for (int i = 0; i < dimension * active; ++i)
            {
                direction[i] = static_cast<T>(projected_direction[i]);
            }
        }
        ModuleBase::Memory::record("LOBPCG::p_work", p_bytes, false);
        broadcast(direction, coefficient_count);
    }
    const ct::Tensor p_coefficients = host_tensor_cast_or_view<Wide>(direction, coefficient_count);
    ModuleBase::Memory::record("LOBPCG::Cp", p_coefficients.buffer().GetAllocatedBytes(), false);
    const VectorImages history = basis_images(n_max_);
    const T* prepared = prepare_rotation_coefficients(p_coefficients);
    rotate_images(basis, dimension, active, prepared, history, evec_.data<T>());
    return normalize_columns(history, active, direction, dimension);
}

/// Copy only the requested Ritz pairs back to the caller, preserving its leading dimension.
template <typename T, typename Device>
void DiagoLOBPCG<T, Device>::write_results(T* psi_in, const int ld_psi, Real* eigenvalue_in)
{
    syncmem_complex_2d_op()(psi_in, ld_psi, x_new_.data<T>(), n_dim_, n_dim_, n_band_);
    std::copy_n(eig_.data<Real>(), n_band_, eigenvalue_in);
}

// Shared vector transformations and storage.
/// Apply the same right coefficient matrix to each supplied X/HX/SX image; scratch permits overlap.
template <typename T, typename Device>
void DiagoLOBPCG<T, Device>::rotate_images(const VectorImages& source,
    const int source_columns,
    const int target_columns,
    const T* coefficients,
    const VectorImages& target,
    T* scratch)
{
    ModuleBase::timer::Scope function_timer("DiagoLOBPCG", "rotate_images");
    const T* inputs[] = {source.x, source.hx, source.sx};
    T* outputs[] = {target.x, target.hx, target.sx};
    const int elements = n_dim_ * target_columns;
    for (int image = 0; image < 3; ++image)
    {
        if (inputs[image] == nullptr)
        {
            continue;
        }
        T* destination = scratch == nullptr ? outputs[image] : scratch;
        algebra_.expand_prepared(
            n_dim_, n_dim_, source_columns, target_columns, inputs[image], coefficients, destination, zero_);
        if (scratch != nullptr)
        {
            copy_op_(elements, scratch, 1, outputs[image], 1);
        }
    }
}

/// Scale nonzero columns by their global Euclidean norm and compact images and coefficients together.
template <typename T, typename Device>
int DiagoLOBPCG<T, Device>::normalize_columns(
    const VectorImages& block, const int columns, T* coefficients, const int coefficient_rows)
{
    ModuleBase::timer::Scope function_timer("DiagoLOBPCG", "normalize_columns");
    const int n = n_dim_;
    ct::Tensor norms = host_tensor<double>(columns, 0.0);
    if (std::is_same<Device, base_device::DEVICE_CPU>::value)
    {
        for (int j = 0; j < columns; ++j)
        {
            for (int i = 0; i < n; ++i)
            {
                const Wide value = static_cast<Wide>(block.x[i + j * n]);
                norms.data<double>()[j] += std::norm(value);
            }
        }
#ifdef __MPI
        Parallel_Common::reduce_data(norms.data<double>(), columns, comm_.comm);
#endif
    }
    else
    {
        // dots_tensor returns globally reduced products for corresponding columns.
        const int stride = n * columns;
        const ct::Tensor dots = algebra_.dots_tensor(n, n, columns, 1, stride, block.x, block.x);
        for (int j = 0; j < columns; ++j)
        {
            norms.data<double>()[j] = dots.template data<Wide>()[j].real();
        }
    }
    int kept = 0;
    for (int j = 0; j < columns; ++j)
    {
        if (!std::isfinite(norms.data<double>()[j]))
        {
            ModuleBase::WARNING_QUIT("DiagoLOBPCG::normalize_columns", "Nonfinite LOBPCG search direction.");
        }
        if (norms.data<double>()[j] <= std::numeric_limits<Real>::min())
        {
            continue;
        }
        const Real scale = static_cast<Real>(1.0 / std::sqrt(norms.data<double>()[j]));
        // Compaction is monotone (kept <= j); no unread column is overwritten.
        for (T* image : {block.x, block.hx, block.sx})
        {
            if (image == nullptr)
            {
                continue;
            }
            T* target = image + kept * n;
            const T* source = image + j * n;
            if (std::is_same<Device, base_device::DEVICE_CPU>::value)
            {
                // A serial loop avoids the vector kernel's per-column OpenMP overhead.
                for (int i = 0; i < n; ++i)
                {
                    target[i] = source[i] * scale;
                }
            }
            else
            {
                ModuleBase::vector_mul_real_op<T, Device>()(n, target, source, scale);
            }
        }
        if (coefficients != nullptr)
        {
            for (int i = 0; i < coefficient_rows; ++i)
            {
                coefficients[i + kept * coefficient_rows] = coefficients[i + j * coefficient_rows] * scale;
            }
        }
        ++kept;
    }
    return kept;
}

/// Borrow host double coefficients on CPU, or convert/upload them into reusable algebra scratch.
template <typename T, typename Device>
const T* DiagoLOBPCG<T, Device>::prepare_rotation_coefficients(const ct::Tensor& coefficients)
{
    // CPU double BLAS can consume the caller-owned coefficients directly.
    // The caller retains them through every synchronous local rotation.
    const Wide* input = coefficients.data<Wide>();
    if (std::is_same<Device, base_device::DEVICE_CPU>::value && std::is_same<T, Wide>::value)
    {
        return reinterpret_cast<const T*>(input);
    }
    const int64_t count = coefficients.NumElements();
    return algebra_.prepare_coefficients(input, count);
}

/// Return non-owning X/HX/SX pointers at one column offset in the shared trial-basis storage.
template <typename T, typename Device>
typename DiagoLOBPCG<T, Device>::VectorImages DiagoLOBPCG<T, Device>::basis_images(const int first)
{
    const int offset = first * n_dim_;
    T* sx = use_overlap_ ? sspace_.data<T>() + offset : nullptr;
    return {space_.data<T>() + offset, hspace_.data<T>() + offset, sx};
}

/// Zero unused rows of every spinor component so full-stride BLAS products exclude nonphysical entries.
template <typename T, typename Device>
void DiagoLOBPCG<T, Device>::project_padding(T* x, const int columns) const
{
    const int stride = n_dim_ / layout_.npol;
    for (int j = 0; j < columns; ++j)
    {
        for (int component = 0; component < layout_.npol; ++component)
        {
            T* first = x + j * n_dim_ + component * stride;
            const int padding = stride - layout_.npw;
            if (padding > 0)
            {
                setmem_complex_op()(first + layout_.npw, zero_, padding);
            }
        }
    }
}

/// Allocate wavefunction blocks and small RR buffers for the current internal width, and record capacity.
template <typename T, typename Device>
void DiagoLOBPCG<T, Device>::allocate_workspace()
{
    ModuleBase::timer::Scope function_timer("DiagoLOBPCG", "allocate_workspace");
    const int len_space = 3 * n_max_;
    const ct::DeviceType host = ct::DeviceType::CpuDevice;
    eig_ = ct::Tensor(r_type_, host, {n_max_});
    eig_device_ = ct::Tensor(r_type_, device_type_, {n_max_});
    evec_ = ct::Tensor(t_type_, device_type_, {n_dim_, n_max_});
    space_ = ct::Tensor(t_type_, device_type_, {n_dim_, len_space});
    hspace_ = ct::Tensor(t_type_, device_type_, {n_dim_, len_space});
    // C has at most (3*b) rows and b requested Ritz columns.
    ritz_coefficients_ = ct::Tensor(t_type_, host, {len_space, n_max_});
    e_red_ = ct::Tensor(r_type_, host, {len_space});
    x_new_ = ct::Tensor(t_type_, device_type_, {n_dim_, n_max_});
    hx_new_ = ct::Tensor(t_type_, device_type_, {n_dim_, n_max_});
    if (use_overlap_)
    {
        sspace_ = ct::Tensor(t_type_, device_type_, {n_dim_, len_space});
        sx_new_ = ct::Tensor(t_type_, device_type_, {n_dim_, n_max_});
        sspace_.zero();
        sx_new_.zero();
    }
    residual_ = ct::Tensor(t_type_, device_type_, {n_dim_, n_max_});
    u_p_ = ct::Tensor(t_type_, host, {len_space, n_max_});

    for (ct::Tensor* buffer :
        {&eig_, &space_, &hspace_, &ritz_coefficients_, &e_red_, &x_new_, &hx_new_, &residual_, &u_p_})
    {
        buffer->zero();
    }

    const std::pair<const char*, const ct::Tensor*> allocations[] = {{"LOBPCG::prec", &prec_},
        {"LOBPCG::eig", &eig_},
        {"LOBPCG::eig_device", &eig_device_},
        {"LOBPCG::evec", &evec_},
        {"LOBPCG::space", &space_},
        {"LOBPCG::hspace", &hspace_},
        {"LOBPCG::sspace", &sspace_},
        {"LOBPCG::ritz_coefficients", &ritz_coefficients_},
        {"LOBPCG::e_red", &e_red_},
        {"LOBPCG::x_new", &x_new_},
        {"LOBPCG::hx_new", &hx_new_},
        {"LOBPCG::sx_new", &sx_new_},
        {"LOBPCG::residual", &residual_},
        {"LOBPCG::u_p", &u_p_}};
    for (const std::pair<const char*, const ct::Tensor*>& allocation : allocations)
    {
        const std::size_t bytes = allocation.second->buffer().GetAllocatedBytes();
        if (allocation.second->device_type() == ct::DeviceType::CpuDevice)
        {
            ModuleBase::Memory::record(allocation.first, bytes, false);
        }
#if defined(__CUDA) || defined(__ROCM)
        else
        {
            ModuleBase::Memory::record_gpu(allocation.first, bytes, false);
        }
#endif
    }
    const std::size_t host_bytes = algebra_.storage_bytes(ct::DeviceType::CpuDevice);
    ModuleBase::Memory::record("LOBPCG::algebra", host_bytes, false);
#if defined(__CUDA) || defined(__ROCM)
    const std::size_t device_bytes = algebra_.storage_bytes(ct::DeviceType::GpuDevice);
    ModuleBase::Memory::record_gpu("LOBPCG::algebra", device_bytes, false);
#endif
}

template class DiagoLOBPCG<std::complex<float>, base_device::DEVICE_CPU>;
template class DiagoLOBPCG<std::complex<double>, base_device::DEVICE_CPU>;
#if defined(__CUDA) || defined(__ROCM)
template class DiagoLOBPCG<std::complex<float>, base_device::DEVICE_GPU>;
template class DiagoLOBPCG<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace hsolver
