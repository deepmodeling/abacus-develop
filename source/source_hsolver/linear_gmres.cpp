#include "source_hsolver/linear_gmres.h"

#include <numeric>
#include <stdexcept>

namespace hsolver
{
template <typename T, typename Device>
LinearGMRES<T, Device>::LinearGMRES(double tolerance, const LinearSolveOptions& options, const diag_comm_info& comm)
    : tolerance_(tolerance), max_iter_(options.max_iterations), restart_(0), reconstruct_(options.reconstruct), work_(comm, 2),
      algebra_(comm)
{
    restart_ = std::min(options.restart, std::max(1, options.max_iterations));
    if (restart_ <= 0 || restart_ > (std::numeric_limits<int>::max() - 1) / 3)
    {
        throw std::invalid_argument("Invalid GMRES restart dimension.");
    }
}

template <typename T, typename Device>
T* LinearGMRES<T, Device>::vector(int slot)
{
    return krylov_.template data<T>() + static_cast<int64_t>(slot) * stride_;
}

template <typename T, typename Device>
void LinearGMRES<T, Device>::retire(int band)
{
    const int last = --active_;
    if (band == last)
    {
        return;
    }
    swaps_.push_back(band);
    swaps_.push_back(last);
    std::swap(order_[band], order_[last]);
    std::swap(h_[band], h_[last]);
    std::swap(g_[band], g_[last]);
    std::swap(sine_[band], sine_[last]);
    std::swap(cosine_[band], cosine_[last]);
}

template <typename T, typename Device>
void LinearGMRES<T, Device>::apply_swaps(int last)
{
    work_.swap_columns(ld_, dim_, swaps_);
    work_.swap_vectors(ld_, dim_, last + 2, stride_, basis(0), swaps_);
    if (last >= 0)
    {
        work_.swap_vectors(ld_, dim_, last + 1, stride_, direction(0), swaps_);
        work_.swap_vectors(ld_, dim_, last + 1, stride_, image(0), swaps_);
    }
    swaps_.clear();
}

template <typename T, typename Device>
T* LinearGMRES<T, Device>::residual()
{
    return work_.data(0);
}

template <typename T, typename Device>
T* LinearGMRES<T, Device>::solution()
{
    return work_.data(1);
}

template <typename T, typename Device>
T* LinearGMRES<T, Device>::basis(int index)
{
    return vector(index);
}

template <typename T, typename Device>
T* LinearGMRES<T, Device>::direction(int index)
{
    return vector(restart_ + 1 + index);
}

template <typename T, typename Device>
T* LinearGMRES<T, Device>::image(int index)
{
    return vector(2 * restart_ + 1 + index);
}

template <typename T, typename Device>
void LinearGMRES<T, Device>::orthogonalize(int j, T* next, std::vector<T>* coefficients)
{
    std::vector<T>& coeff = *coefficients;
    // Two-pass classical Gram-Schmidt batches the global reductions.
    for (int pass = 0; pass < 2; ++pass)
    {
        const std::vector<Wide> dots = algebra_.arnoldi_dots(ld_, dim_, active_, j + 1, stride_, basis(0), next);
        for (int i = 0; i <= j; ++i)
        {
            for (int b = 0; b < active_; ++b)
            {
                h_[b][i + j * (restart_ + 1)] += dots[i * active_ + b];
                coeff[b] = T(dots[i * active_ + b]);
            }
            work_.batch(ld_, dim_, active_, next, next, basis(i), T(1), T(-1), nullptr, coeff.data(), nullptr);
        }
    }
}

template <typename T, typename Device>
bool LinearGMRES<T, Device>::update_qr(int b, int j, double norm)
{
    Wide* column = h_[b].data() + j * (restart_ + 1);
    column[j + 1] = norm;
    for (int i = 0; i < j; ++i)
    {
        const Wide upper = cosine_[b][i] * column[i] + sine_[b][i] * column[i + 1];
        column[i + 1] = -std::conj(sine_[b][i]) * column[i] + cosine_[b][i] * column[i + 1];
        column[i] = upper;
    }
    const double magnitude = std::hypot(std::abs(column[j]), norm);
    if (magnitude == 0 || !std::isfinite(magnitude))
    {
        return false;
    }
    const Wide phase = std::abs(column[j]) == 0 ? Wide(1) : column[j] / std::abs(column[j]);
    cosine_[b][j] = std::abs(column[j]) / magnitude;
    sine_[b][j] = phase * norm / magnitude;
    column[j] = phase * magnitude;
    column[j + 1] = 0;
    g_[b][j + 1] = -std::conj(sine_[b][j]) * g_[b][j];
    g_[b][j] *= cosine_[b][j];
    return true;
}

template <typename T, typename Device>
bool LinearGMRES<T, Device>::back_substitute(int b, int j, std::vector<Wide>* weights)
{
    weights->assign(g_[b].begin(), g_[b].begin() + j + 1);
    for (int i = j; i >= 0; --i)
    {
        (*weights)[i] /= h_[b][i + i * (restart_ + 1)];
        if (!std::isfinite(std::abs((*weights)[i])))
        {
            return false;
        }
        for (int k = 0; k < i; ++k)
        {
            (*weights)[k] -= h_[b][k + i * (restart_ + 1)] * (*weights)[i];
        }
    }
    return true;
}

template <typename T, typename Device>
void LinearGMRES<T, Device>::update_solution(int j,
                                             const std::vector<int>& done,
                                             const std::vector<std::vector<Wide>>& weights,
                                             bool reconstruct,
                                             std::vector<T>* coefficients)
{
    std::vector<T>& coeff = *coefficients;
    const bool any_done = std::find(done.begin(), done.end(), 1) != done.end();
    for (int i = 0; any_done && i <= j; ++i)
    {
        for (int b = 0; b < active_; ++b)
        {
            coeff[b] = done[b] ? T(weights[b][i]) : T(0);
        }
        work_.batch(ld_, dim_, active_, solution(), solution(), direction(i), T(1), T(1), nullptr, coeff.data(), nullptr);
        if (reconstruct)
        {
            work_.batch(ld_, dim_, active_, residual(), residual(), image(i), T(1), T(-1), nullptr, coeff.data(), nullptr);
        }
    }
}

template <typename T, typename Device>
bool LinearGMRES<T, Device>::start_cycle(const std::vector<double>& threshold, std::vector<T>* coefficients)
{
    h_.assign(bands_, std::vector<Wide>((restart_ + 1) * restart_, 0));
    g_.assign(bands_, std::vector<Wide>(restart_ + 1, 0));
    sine_.assign(bands_, std::vector<Wide>(restart_, 0));
    cosine_.assign(bands_, std::vector<double>(restart_, 0));
    std::vector<Wide> norms = algebra_.dots(ld_, dim_, bands_, 1, stride_, residual(), residual());
    std::vector<T>& coeff = *coefficients;
    active_ = bands_;
    swaps_.clear();
    std::iota(order_.begin(), order_.end(), 0);
    for (int b = 0; b < bands_; ++b)
    {
        const double norm = linear_norm(norms[b]);
        if (!std::isfinite(norm))
        {
            return false;
        }
        g_[b][0] = norm;
        coeff[b] = norm == 0 ? T(0) : T(1.0 / norm);
    }
    work_.batch(ld_, dim_, bands_, basis(0), residual(), residual(), T(1), T(0), coeff.data(), nullptr, nullptr);
    for (int b = bands_ - 1; b >= 0; --b)
    {
        if (std::abs(g_[b][0]) <= 0.8 * threshold[order_[b]])
        {
            retire(b);
        }
    }
    apply_swaps(-1);
    return true;
}

template <typename T, typename Device>
bool LinearGMRES<T, Device>::cycle(const LinearOperator<T, Device>& op,
                                   const LinearOperator<T, Device>& preconditioner,
                                   const std::vector<double>& threshold,
                                   bool reconstruct,
                                   LinearSolveResult* result)
{
    std::vector<T> coeff(bands_);
    if (!start_cycle(threshold, &coeff))
    {
        return false;
    }
    for (int j = 0; j < restart_ && active_ > 0 && result->iterations < max_iter_; ++j)
    {
        T* z = direction(j);
        T* raw = image(j);
        T* next = basis(j + 1);
        preconditioner.apply(basis(j), z, ld_, active_);
        work_.apply(op, z, raw, ld_, active_);
        work_.copy(ld_, dim_, active_, raw, next);
        orthogonalize(j, next, &coeff);
        const std::vector<Wide> norms = algebra_.dots(ld_, dim_, active_, 1, stride_, next, next);
        ++result->iterations;
        std::vector<int> done(active_, 0);
        std::vector<std::vector<Wide>> weights(active_);
        for (int b = 0; b < active_; ++b)
        {
            const double norm = linear_norm(norms[b]);
            if (!std::isfinite(norm))
            {
                return false;
            }
            if (!update_qr(b, j, norm))
            {
                return false;
            }
            done[b] = std::abs(g_[b][j + 1]) <= 0.8 * threshold[order_[b]] || norm <= std::numeric_limits<Real>::min() || j + 1 == restart_
                      || result->iterations == max_iter_;
            coeff[b] = done[b] ? T(0) : T(1.0 / norm);
            if (!done[b])
            {
                continue;
            }
            if (!back_substitute(b, j, &weights[b]))
            {
                return false;
            }
        }
        work_.batch(ld_, dim_, active_, next, next, next, T(1), T(0), coeff.data(), nullptr, nullptr);
        update_solution(j, done, weights, reconstruct, &coeff);
        for (int b = active_ - 1; b >= 0; --b)
        {
            if (done[b])
            {
                retire(b);
            }
        }
        apply_swaps(j);
    }
    return true;
}

template <typename T, typename Device>
LinearSolveResult LinearGMRES<T, Device>::solve(const LinearOperator<T, Device>& op,
                                                const LinearOperator<T, Device>& preconditioner,
                                                int ld,
                                                int nvec,
                                                int dim,
                                                T* x,
                                                const T* b,
                                                const T* initial_residual,
                                                bool force_check)
{
    const LinearSolveTimer timer("LinearGMRES");
    work_.prepare(ld, dim, nvec, x, b, tolerance_, max_iter_);
    work_.reset_statistics();
    ld_ = ld;
    dim_ = dim;
    bands_ = nvec;
    stride_ = std::max(1, ld * nvec);
    linear_buffer<T, Device>(&krylov_, static_cast<int64_t>(3 * restart_ + 1) * stride_);
    order_.resize(nvec);
    const std::vector<Wide> rhs_norms = algebra_.dots(ld, dim, nvec, 1, stride_, b, b);
    std::vector<double> threshold(nvec);
    for (int i = 0; i < nvec; ++i)
    {
        threshold[i] = tolerance_ * std::max(1.0, linear_norm(rhs_norms[i]));
    }
    bool reconstruct = reconstruct_;
    force_check = force_check || tolerance_ < 100 * std::numeric_limits<Real>::epsilon();
    LinearSolveResult result;
    while (true)
    {
        work_.clear();
        work_.copy(ld, dim, nvec, x, solution());
        if (initial_residual && result.restarts == 0)
        {
            work_.copy(ld, dim, nvec, initial_residual, residual());
        }
        else
        {
            work_.residual(op, ld, dim, nvec, x, b, residual());
        }
        const int start = result.iterations;
        bool regular = false;
        try
        {
            regular = cycle(op, preconditioner, threshold, reconstruct, &result);
            result.status = regular ? LinearSolveStatus::max_iterations : LinearSolveStatus::breakdown;
        }
        catch (const LinearPreconditionerError&)
        {
            // Unlike completed Arnoldi steps, this failed attempt has not yet been counted.
            ++result.iterations;
            result.status = LinearSolveStatus::preconditioner_failure;
        }
        work_.restore(ld, dim, nvec, order_, solution(), x);
        bool accepted = regular && reconstruct;
        bool invalid_reconstruction = false;
        if (accepted)
        {
            const std::vector<Wide> norms = algebra_.dots(ld, dim, nvec, 1, stride_, residual(), residual());
            result.max_residual = 0;
            for (int i = 0; i < nvec; ++i)
            {
                const double error = linear_norm(norms[i]);
                result.max_residual = std::max(result.max_residual, error);
                invalid_reconstruction = invalid_reconstruction || !std::isfinite(error);
                if (!std::isfinite(error) || error > threshold[order_[i]])
                {
                    accepted = false;
                }
            }
            if (accepted && !force_check)
            {
                result.status = LinearSolveStatus::converged;
                result.failed_band = -1;
                result.reconstructed = true;
                work_.statistics(&result);
                return result;
            }
        }
        work_.residual(op, ld, dim, nvec, x, b, residual());
        ++result.true_checks;
        work_.statistics(&result);
        const std::vector<Wide> true_norms = algebra_.dots(ld, dim, nvec, 1, stride_, residual(), residual());
        result.max_residual = 0;
        result.failed_band = -1;
        for (int i = 0; i < nvec; ++i)
        {
            const double error = linear_norm(true_norms[i]);
            result.max_residual = std::max(result.max_residual, error);
            if ((!std::isfinite(error) || error > threshold[i]) && result.failed_band < 0)
            {
                result.failed_band = i;
            }
        }
        if (result.failed_band < 0)
        {
            result.status = LinearSolveStatus::converged;
        }
        if (result.status == LinearSolveStatus::converged)
        {
            return result;
        }
        if (result.status == LinearSolveStatus::preconditioner_failure)
        {
            // Let the caller retry with a different preconditioner and the remaining budget.
            return result;
        }
        const bool recover_anomaly = reconstruct && !regular;
        // A full restart cycle is not a reconstruction failure.
        if (reconstruct && (accepted || !regular || invalid_reconstruction))
        {
            ++result.reconstruction_fallbacks;
            reconstruct = false;
        }
        if (!regular && !recover_anomaly)
        {
            return result;
        }
        const bool retry_initial = initial_residual && result.restarts == 0;
        if (result.iterations >= max_iter_ || (result.iterations == start && !retry_initial) || !std::isfinite(result.max_residual))
        {
            return result;
        }
        ++result.restarts;
    }
}
template class LinearGMRES<std::complex<float>, base_device::DEVICE_CPU>;
template class LinearGMRES<std::complex<double>, base_device::DEVICE_CPU>;
#if defined(__CUDA) || defined(__ROCM)
template class LinearGMRES<std::complex<float>, base_device::DEVICE_GPU>;
template class LinearGMRES<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace hsolver
