#include "source_hsolver/linear_bicgstab.h"

#include "source_hsolver/kernels/linear_op.h"

#include <numeric>

namespace hsolver
{
template <typename T, typename Device>
LinearBiCGSTAB<T, Device>::LinearBiCGSTAB(const double tolerance, const int max_iter, const diag_comm_info& comm)
    : tolerance_(tolerance), max_iter_(max_iter), work_(comm)
{
}

template <typename T, typename Device>
bool LinearBiCGSTAB<T, Device>::breakdown(const T divisor, const int band, LinearSolveResult* result) const
{
    if (!linear_bad_divisor(divisor))
    {
        return false;
    }
    result->status = LinearSolveStatus::breakdown;
    result->failed_band = original_[band];
    return true;
}

template <typename T, typename Device>
void LinearBiCGSTAB<T, Device>::retire(const int band)
{
    const int last = active_ - 1;
    if (band != last)
    {
        swaps_.push_back(band);
        swaps_.push_back(last);
        std::swap(rho_[band], rho_[last]);
        std::swap(rho_prev_[band], rho_prev_[last]);
        std::swap(alpha_[band], alpha_[last]);
        std::swap(omega_[band], omega_[last]);
        std::swap(threshold_[band], threshold_[last]);
        std::swap(original_[band], original_[last]);
    }
    --active_;
}

template <typename T, typename Device>
void LinearBiCGSTAB<T, Device>::retire_converged(const bool alpha_step, const bool identity)
{
    const T* residual = work_.data(alpha_step ? intermediate_slot : residual_slot);
    if (alpha_step)
    {
        work_.dot(ld_, dim_, active_, residual, residual, norm_.data());
    }
    else
    {
        work_.dot_pair(ld_, dim_, active_, residual, residual, work_.data(shadow_slot), residual, norm_.data(), rho_.data());
    }
    swaps_.clear();
    skip_.assign(active_, 1);
    bool any_converged = false;
    for (int band = 0; band < active_; ++band)
    {
        const double error = linear_norm(norm_[band]);
        if (std::isfinite(error) && error <= threshold_[band])
        {
            skip_[band] = 0;
            any_converged = true;
        }
    }
    if (!any_converged)
    {
        return;
    }
    if (alpha_step)
    {
        work_.batch(ld_,
                    dim_,
                    active_,
                    work_.data(solution_slot),
                    work_.data(solution_slot),
                    work_.data(identity ? direction_slot : precond_direction_slot),
                    T(1),
                    T(1),
                    nullptr,
                    alpha_.data(),
                    skip_.data());
    }
    for (int band = active_ - 1; band >= 0; --band)
    {
        const double error = linear_norm(norm_[band]);
        if (std::isfinite(error) && error <= threshold_[band])
        {
            retire(band);
        }
    }
    work_.swap_columns(ld_, dim_, swaps_);
}

template <typename T, typename Device>
bool LinearBiCGSTAB<T, Device>::iterate(const LinearOperator<T, Device>& op,
                                        const LinearOperator<T, Device>& preconditioner,
                                        LinearSolveResult* result)
{
    T* r = work_.data(residual_slot);
    const bool identity = preconditioner.is_identity();
    T* shadow = work_.data(shadow_slot);
    T* p = work_.data(direction_slot);
    T* v = work_.data(image_slot);
    T* s = work_.data(intermediate_slot);
    T* y = work_.data(identity ? direction_slot : precond_direction_slot);
    T* z = work_.data(identity ? intermediate_slot : precond_intermediate_slot);
    T* t = work_.data(intermediate_image_slot);
    T* solution = work_.data(solution_slot);

    for (int band = 0; band < active_; ++band)
    {
        if (breakdown(rho_[band], band, result) || breakdown(rho_prev_[band], band, result) || breakdown(omega_[band], band, result))
        {
            return false;
        }
        beta_[band] = (rho_[band] / rho_prev_[band]) * (alpha_[band] / omega_[band]);
    }
    work_.bicg_update(ld_, dim_, active_, true, p, r, v, beta_.data(), omega_.data());
    if (!identity)
    {
        preconditioner.apply(p, y, ld_, active_);
    }
    work_.apply(op, y, v, ld_, active_);
    work_.dot(ld_, dim_, active_, shadow, v, denominator_.data());
    for (int band = 0; band < active_; ++band)
    {
        if (breakdown(denominator_[band], band, result))
        {
            return false;
        }
        alpha_[band] = rho_[band] / denominator_[band];
    }
    work_.batch(ld_, dim_, active_, s, r, v, T(1), T(-1), nullptr, alpha_.data(), nullptr);
    // Retire exact alpha steps before evaluating omega = <t,s>/<t,t>.
    retire_converged(true, identity);
    if (active_ == 0)
    {
        return true;
    }
    if (!identity)
    {
        preconditioner.apply(s, z, ld_, active_);
    }
    work_.apply(op, z, t, ld_, active_);
    work_.dot_pair(ld_, dim_, active_, t, s, t, t, denominator_.data(), tt_.data());
    for (int band = 0; band < active_; ++band)
    {
        if (breakdown(tt_[band], band, result))
        {
            return false;
        }
        omega_[band] = denominator_[band] / tt_[band];
        if (breakdown(omega_[band], band, result))
        {
            return false;
        }

        rho_prev_[band] = rho_[band];
    }
    work_.bicg_update(ld_, dim_, active_, false, solution, y, z, alpha_.data(), omega_.data());
    work_.batch(ld_, dim_, active_, r, s, t, T(1), T(-1), nullptr, omega_.data(), nullptr);
    retire_converged(false, identity);
    return true;
}

template <typename T, typename Device>
LinearSolveResult LinearBiCGSTAB<T, Device>::solve(const LinearOperator<T, Device>& op,
                                                   const LinearOperator<T, Device>& preconditioner,
                                                   int ld,
                                                   int nband,
                                                   int dim,
                                                   T* x,
                                                   const T* b,
                                                   const T* initial_residual)
{
    const LinearSolveTimer timer("LinearBiCGSTAB");
    work_.prepare(ld, dim, nband, x, b, tolerance_, max_iter_);
    work_.reset_statistics();
    ld_ = ld;
    dim_ = dim;
    active_ = nband;
    beta_.resize(nband);
    rho_.assign(nband, T(0));
    rho_prev_.assign(nband, T(1));
    alpha_.assign(nband, T(1));
    omega_.assign(nband, T(1));
    denominator_.resize(nband);
    tt_.resize(nband);
    norm_.resize(nband);
    threshold_.resize(nband);
    original_.resize(nband);
    std::iota(original_.begin(), original_.end(), 0);
    work_.dot(ld, dim, nband, b, b, norm_.data());
    for (int band = 0; band < nband; ++band)
    {
        threshold_[band] = tolerance_ * std::max(1.0, linear_norm(norm_[band]));
    }
    const std::vector<double> original_threshold = threshold_;
    LinearSolveResult result;
    while (true)
    {
        // A restart keeps x, but discards the recursively accumulated residual
        // and search directions. All cycles share one iteration budget.
        work_.clear();
        active_ = nband;
        threshold_ = original_threshold;
        std::iota(original_.begin(), original_.end(), 0);
        rho_prev_.assign(nband, T(1));
        alpha_.assign(nband, T(1));
        omega_.assign(nband, T(1));
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
        retire_converged(false, preconditioner.is_identity());
        result.status = LinearSolveStatus::max_iterations;
        const int cycle_start = result.iterations;
        while (result.iterations < max_iter_ && active_ > 0)
        {
            ++result.iterations;
            if (!iterate(op, preconditioner, &result))
            {
                break;
            }
        }
        work_.restore(ld, dim, nband, original_, work_.data(solution_slot), x);
        if (active_ == 0)
        {
            result.status = LinearSolveStatus::residual_mismatch;
        }
        work_.verify(op, ld, dim, nband, x, b, original_threshold, work_.data(residual_slot), &result);
        if (result.status != LinearSolveStatus::residual_mismatch || result.iterations >= max_iter_ || result.iterations == cycle_start
            || !std::isfinite(result.max_residual))
        {
            return result;
        }
        ++result.restarts;
    }
}

template class LinearBiCGSTAB<std::complex<float>, base_device::DEVICE_CPU>;
template class LinearBiCGSTAB<std::complex<double>, base_device::DEVICE_CPU>;
#if defined(__CUDA) || defined(__ROCM)
template class LinearBiCGSTAB<std::complex<float>, base_device::DEVICE_GPU>;
template class LinearBiCGSTAB<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace hsolver
