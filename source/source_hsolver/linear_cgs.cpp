#include "source_hsolver/linear_cgs.h"

#include <numeric>

namespace hsolver
{
template <typename T, typename Device>
LinearCGS<T, Device>::LinearCGS(const double tolerance, const diag_comm_info& comm) : tolerance_(tolerance), work_(comm, 9)
{
}

template <typename T, typename Device>
void LinearCGS<T, Device>::retire_converged()
{
    work_.dot_pair(ld_,
                   dim_,
                   active_,
                   work_.data(residual_slot),
                   work_.data(residual_slot),
                   work_.data(shadow_slot),
                   work_.data(residual_slot),
                   norm_.data(),
                   rho_.data());
    swaps_.clear();
    for (int band = active_ - 1; band >= 0; --band)
    {
        const double error = linear_norm(norm_[band]);
        if (!std::isfinite(error) || error > threshold_[band])
        {
            continue;
        }
        const int last = --active_;
        if (band != last)
        {
            swaps_.push_back(band);
            swaps_.push_back(last);
            std::swap(rho_[band], rho_[last]);
            std::swap(rho_prev_[band], rho_prev_[last]);
            std::swap(threshold_[band], threshold_[last]);
            std::swap(original_[band], original_[last]);
        }
    }
    work_.swap_columns(ld_, dim_, swaps_);
}

template <typename T, typename Device>
bool LinearCGS<T, Device>::iterate(const LinearOperator<T, Device>& op,
                                   const LinearOperator<T, Device>& preconditioner,
                                   const bool first_iteration,
                                   LinearSolveResult* result)
{
    T* r = work_.data(residual_slot);
    T* shadow = work_.data(shadow_slot);
    T* p = work_.data(direction_slot);
    T* q = work_.data(q_slot);
    T* u = work_.data(u_slot);
    T* v = work_.data(image_slot);
    T* scratch = work_.data(scratch_slot);
    T* direction = work_.data(update_slot);
    T* solution = work_.data(solution_slot);
    const bool identity = preconditioner.is_identity();
    for (int band = 0; band < active_; ++band)
    {
        if (linear_bad_divisor(rho_[band]) || linear_bad_divisor(rho_prev_[band]))
        {
            result->status = LinearSolveStatus::breakdown;
            result->failed_band = original_[band];
            return false;
        }
        beta_[band] = first_iteration ? T(0) : rho_[band] / rho_prev_[band];
    }
    work_.cgs_direction(ld_, dim_, active_, p, u, r, q, beta_.data());
    if (!identity)
    {
        preconditioner.apply(p, scratch, ld_, active_);
    }
    const T* preconditioned_direction = identity ? p : scratch;
    work_.apply(op, preconditioned_direction, v, ld_, active_);
    work_.dot(ld_, dim_, active_, shadow, v, denominator_.data());
    for (int band = 0; band < active_; ++band)
    {
        if (linear_bad_divisor(denominator_[band]))
        {
            result->status = LinearSolveStatus::breakdown;
            result->failed_band = original_[band];
            return false;
        }
        alpha_[band] = rho_[band] / denominator_[band];
    }
    work_.cgs_alpha(ld_, dim_, active_, q, direction, u, v, alpha_.data());
    if (!identity)
    {
        preconditioner.apply(direction, scratch, ld_, active_);
    }
    const T* update = identity ? direction : scratch;
    work_.apply(op, update, v, ld_, active_);
    work_.cgs_finish(ld_, dim_, active_, solution, r, update, v);
    rho_prev_ = rho_;
    return true;
}

template <typename T, typename Device>
LinearSolveResult LinearCGS<T, Device>::solve(const LinearOperator<T, Device>& op,
                                              const LinearOperator<T, Device>& preconditioner,
                                              int ld,
                                              int nband,
                                              int dim,
                                              T* x,
                                              const T* b,
                                              const T* initial_residual,
                                              const int max_iterations)
{
    const LinearSolveTimer timer("LinearCGS");
    work_.prepare(ld, dim, nband, x, b, tolerance_, max_iterations);
    work_.reset_statistics();
    ld_ = ld;
    dim_ = dim;
    beta_.resize(nband);
    rho_.resize(nband);
    rho_prev_.resize(nband);
    alpha_.resize(nband);
    denominator_.resize(nband);
    norm_.resize(nband);
    threshold_.resize(nband);
    original_.resize(nband);
    work_.dot(ld, dim, nband, b, b, norm_.data());
    for (int band = 0; band < nband; ++band)
    {
        const double rhs_norm = linear_norm(norm_[band]);
        threshold_[band] = tolerance_ * (rhs_norm == 0.0 ? 1.0 : rhs_norm);
    }
    const std::vector<double> original_threshold = threshold_;
    LinearSolveResult result;
    while (true)
    {
        work_.clear();
        active_ = nband;
        threshold_ = original_threshold;
        std::iota(original_.begin(), original_.end(), 0);
        rho_prev_.assign(nband, T(1));
        work_.copy(ld, dim, nband, x, work_.data(solution_slot));
        if (initial_residual && result.restarts == 0)
        {
            work_.copy(ld, dim, nband, initial_residual, work_.data(residual_slot));
        }
        else
        {
            work_.residual(op, ld, dim, nband, work_.data(solution_slot), b, work_.data(residual_slot));
        }
        work_.copy(ld, dim, nband, work_.data(residual_slot), work_.data(shadow_slot));
        retire_converged();
        result.status = LinearSolveStatus::max_iterations;
        const int cycle_start = result.iterations;
        try
        {
            while (result.iterations < max_iterations && active_ > 0)
            {
                const bool first_iteration = result.iterations == cycle_start;
                ++result.iterations;
                if (!iterate(op, preconditioner, first_iteration, &result))
                {
                    break;
                }
                retire_converged();
            }
        }
        catch (const LinearPreconditionerError&)
        {
            result.status = LinearSolveStatus::preconditioner_failure;
        }
        work_.restore(ld, dim, nband, original_, work_.data(solution_slot), x);
        if (active_ == 0)
        {
            result.status = LinearSolveStatus::residual_mismatch;
        }
        work_.verify(op, ld, dim, nband, x, b, original_threshold, work_.data(residual_slot), &result);
        if (result.status != LinearSolveStatus::residual_mismatch || result.iterations >= max_iterations || result.iterations == cycle_start
            || !std::isfinite(result.max_residual))
        {
            return result;
        }
        ++result.restarts;
    }
}

template class LinearCGS<std::complex<float>, base_device::DEVICE_CPU>;
template class LinearCGS<std::complex<double>, base_device::DEVICE_CPU>;
#if defined(__CUDA) || defined(__ROCM)
template class LinearCGS<std::complex<float>, base_device::DEVICE_GPU>;
template class LinearCGS<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace hsolver
