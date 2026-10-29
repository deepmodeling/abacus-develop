#include "source_hsolver/module_linalg/exp_lanczos.h"

#include "source_base/module_container/base/third_party/lapack.h"
#include "source_base/parallel_device.h"
#include "source_base/timer.h"
#include "source_hsolver/module_linalg/kernels/linear_op.h"

#include <limits>

namespace hsolver
{
namespace
{
using Wide = std::complex<double>;

void add_work(ExpResult& total, const ExpResult& part)
{
    total.hpsi_calls += part.hpsi_calls;
    total.hpsi_cols += part.hpsi_cols;
    total.bisect_count += part.bisect_count;
    total.reorth_count += part.reorth_count;
    total.order_or_dim = std::max(total.order_or_dim, part.order_or_dim);
}

ExpResult failed_band(ExpResult result, int band, bool retry)
{
    result.failed_band = band;
    result.state = retry ? ExpState::retry : ExpState::failed;
    result.has_err_est = false;
    return result;
}
} // namespace

template <typename T, typename Device>
ExpLanczos<T, Device>::ExpLanczos(const ExpControls& exp_ctrl, const diag_comm_info& comm)
    : ExpAction<T, Device>(exp_ctrl, comm), comm_(comm), work_(comm, 1)
{
}

template <typename T, typename Device>
int ExpLanczos<T, Device>::batch_size(const OrbLayout& layout) const
{
    // The layout can differ between ranks, but operator column counts must agree.
    double rows = std::max(1, layout.ld);
#ifdef __MPI
    if (comm_.nproc > 1)
    {
        std::vector<double> sizes(comm_.nproc, 0.0);
        sizes[comm_.rank] = rows;
        Parallel_Common::reduce_data(sizes.data(), comm_.nproc, comm_.comm);
        rows = *std::max_element(sizes.begin(), sizes.end());
    }
#endif
    const double vectors = static_cast<double>(this->exp_ctrl_.krylov_max) + 6.0;
    const double bytes_per_band = sizeof(T) * rows * vectors;
    const double budget = 256.0 * 1024.0 * 1024.0;
    const int memory_limit = static_cast<int>(std::min(64.0, budget / bytes_per_band));
    return std::max(1, std::min(layout.bands, memory_limit));
}

template <typename T, typename Device>
void ExpLanczos<T, Device>::combine(const OrbLayout& layout, int count, const T* basis, const std::vector<T>& coefficients, T* output)
{
    ModuleBase::timer::start("ExpLanczos", "combine");
    linear_buffer<T, Device>(&coefficients_, coefficients.size());
    T* data = coefficients_.template data<T>();
    base_device::memory::synchronize_memory_op<T, Device, base_device::DEVICE_CPU>()(data, coefficients.data(), coefficients.size());
    const int dim = layout.npol == 2 ? layout.ld : layout.npw;
    const int stride = layout.ld * layout.bands;
    linear_op<T, Device>().combine(layout.ld, dim, layout.bands, count, stride, basis, data, output);
    ModuleBase::timer::end("ExpLanczos", "combine");
}

template <typename T, typename Device>
std::vector<Wide> ExpLanczos<T, Device>::basis_gram(const OrbLayout& layout, int count, const T* basis)
{
    ModuleBase::timer::start("ExpLanczos", "basis_gram");
    const int size = layout.bands * count * count;
    linear_buffer<Wide, Device>(&products_, size);
    Wide* data = products_.template data<Wide>();
    const int dim = layout.npol == 2 ? layout.ld : layout.npw;
    const int stride = layout.ld * layout.bands;
    linear_op<T, Device>().block_gram(layout.ld, dim, layout.bands, count, stride, basis, data);
    std::vector<Wide> gram(size);
    base_device::memory::synchronize_memory_op<Wide, base_device::DEVICE_CPU, Device>()(gram.data(), data, size);
#ifdef __MPI
    if (comm_.nproc > 1)
    {
        Parallel_Common::reduce_data(gram.data(), size, comm_.comm);
    }
#endif
    ModuleBase::timer::end("ExpLanczos", "basis_gram");
    return gram;
}

template <typename T, typename Device>
bool ExpLanczos<T, Device>::diagonalize(const std::vector<double>& diagonal,
                                        const std::vector<double>& off_diagonal,
                                        int band,
                                        int bands,
                                        int count)
{
    ModuleBase::timer::start("ExpLanczos", "diagonalize");
    small_vectors_.assign(count * count, 0.0);
    small_values_.resize(count);
    const int work_size = std::max(1, 3 * count);
    small_work_.resize(work_size);
    for (int i = 0; i < count; ++i)
    {
        small_vectors_[i * count + i] = diagonal[i * bands + band];
        if (i + 1 < count)
        {
            small_vectors_[(i + 1) * count + i] = off_diagonal[i * bands + band];
        }
    }
    // Only the real symmetric tridiagonal projection is diagonalized.
    const char job = 'V';
    const char upper = 'U';
    int info = 0;
    ct::lapackConnector::heev(job, upper, count, small_vectors_.data(), count, small_values_.data(), small_work_.data(), work_size, info);
    ModuleBase::timer::end("ExpLanczos", "diagonalize");
    return info == 0;
}

template <typename T, typename Device>
ExpResult ExpLanczos<T, Device>::attempt(const HSOperator<T, Device>& op,
                                         const OrbLayout& layout,
                                         double alpha,
                                         double rel_tol,
                                         const T* input,
                                         T* output,
                                         T* full,
                                         bool reorthogonalize,
                                         bool* lost_orthogonality)
{
    using Real = typename T::value_type;
    const int ld = layout.ld;
    const int bands = layout.bands;
    const int dim = layout.npol == 2 ? ld : layout.npw;
    const int stride = ld * bands;
    const int limit = this->exp_ctrl_.krylov_max;
    const int64_t size = std::max(1, stride);
    const double time = full ? 2.0 * std::abs(alpha) : std::abs(alpha);
    *lost_orthogonality = false;
    ExpResult result;
    result.has_err_est = true;
    result.reorth_count = reorthogonalize ? 1 : 0;
    const int64_t basis_size = size * static_cast<int64_t>(limit);
    if (basis_.NumElements() < basis_size)
    {
        // Avoid holding two large Krylov allocations while growing the workspace.
        basis_ = ct::Tensor();
    }
    linear_buffer<T, Device>(&basis_, basis_size);
    linear_buffer<T, Device>(&image_, size);
    linear_buffer<T, Device>(&packed_, size);
    linear_buffer<T, Device>(&hpacked_, size);
    linear_buffer<int, Device>(&permutation_, bands);
    T* q = basis_.template data<T>();
    T* w = image_.template data<T>();
    T* packed = packed_.template data<T>();
    T* hpacked = hpacked_.template data<T>();
    int* map = permutation_.template data<int>();
    base_device::memory::set_memory_op<T, Device>()(output, 0, stride);
    if (full)
    {
        base_device::memory::set_memory_op<T, Device>()(full, 0, stride);
    }
    base_device::memory::synchronize_memory_op<T, Device, Device>()(q, input, stride);
    zero_pad<T, Device>(q, layout);
    const std::vector<std::complex<double>> norms = this->algebra_.dots(ld, dim, bands, 1, stride, q, q);
    std::vector<double> scale(bands);
    std::vector<T> factors(bands);
    std::vector<bool> done(bands, false);
    std::vector<double> diagonal(static_cast<size_t>(bands) * limit);
    std::vector<double> beta(static_cast<size_t>(bands) * limit);
    std::vector<double> log_bound(bands, 0.0);
    std::vector<double> correction_sq(bands, 0.0);
    std::vector<double> errors(bands, 0.0);
    std::vector<int> active;
    active.reserve(bands);
    for (int b = 0; b < bands; ++b)
    {
        if (!std::isfinite(norms[b].real()) || norms[b].real() < 0.0)
        {
            return failed_band(result, b, false);
        }
        scale[b] = std::sqrt(norms[b].real());
        done[b] = scale[b] == 0.0;
        factors[b] = done[b] ? T(0) : T(1.0 / scale[b]);
    }
    this->work_.batch(ld, dim, bands, q, q, q, T(1), T(0), factors.data(), nullptr, nullptr);
    for (int j = 0; j < limit; ++j)
    {
        active.clear();
        for (int b = 0; b < bands; ++b)
        {
            if (!done[b])
            {
                active.push_back(b);
            }
        }
        if (active.empty())
        {
            return result;
        }
        const int count = active.size();
        T* column = q + static_cast<int64_t>(j) * stride;
        if (count == bands)
        {
            op.hpsi(column, w, ld, bands);
        }
        else
        {
            base_device::memory::synchronize_memory_op<int, Device, base_device::DEVICE_CPU>()(map, active.data(), count);
            linear_op<T, Device>().gather(ld, ld, count, 1, stride, column, packed, map);
            op.hpsi(packed, hpacked, ld, count);
            base_device::memory::set_memory_op<T, Device>()(w, 0, stride);
            linear_op<T, Device>().scatter(ld, ld, count, hpacked, w, map);
        }
        zero_pad<T, Device>(w, layout);
        ++result.hpsi_calls;
        result.hpsi_cols += count;
        const std::vector<std::complex<double>> projection = this->algebra_.dots(ld, dim, bands, 1, stride, column, w);
        const int m = j + 1;
        result.order_or_dim = m;
        const int first = std::max(0, j - 1);
        const int columns = j - first + 1;
        std::vector<T> recurrence(columns * bands, T(0));
        for (int b: active)
        {
            if (!std::isfinite(std::abs(projection[b])))
            {
                return failed_band(result, b, false);
            }
            diagonal[j * bands + b] = projection[b].real();
            recurrence[(columns - 1) * bands + b] = T(-projection[b].real());
            if (j > 0)
            {
                recurrence[b] = T(-beta[(j - 1) * bands + b]);
            }
        }
        const T* previous = q + static_cast<int64_t>(first) * stride;
        combine(layout, columns, previous, recurrence, w);
        if (reorthogonalize)
        {
            std::vector<double> removed(bands, 0.0);
            for (int pass = 0; pass < 2; ++pass)
            {
                const std::vector<std::complex<double>> projections = this->algebra_.dots(ld, dim, bands, m, stride, q, w);
                std::vector<T> correction(m * bands, T(0));
                for (int i = 0; i < m; ++i)
                {
                    for (int b: active)
                    {
                        const Wide value = projections[i * bands + b];
                        if (!std::isfinite(std::abs(value)))
                        {
                            return failed_band(result, b, false);
                        }
                        removed[b] += std::abs(value);
                        correction[i * bands + b] = T(-value);
                    }
                }
                combine(layout, m, q, correction, w);
            }
            for (int b: active)
            {
                // ||Q*C|| <= sum_i |C_i| for normalized columns. Accumulate
                // a Frobenius bound on the discarded recurrence corrections.
                correction_sq[b] += removed[b] * removed[b];
            }
        }
        const std::vector<std::complex<double>> residuals = this->algebra_.dots(ld, dim, bands, 1, stride, w, w);
        std::vector<int> candidates;
        for (int b: active)
        {
            if (!std::isfinite(residuals[b].real()) || residuals[b].real() < 0.0)
            {
                return failed_band(result, b, false);
            }
            const double residual = std::sqrt(residuals[b].real());
            beta[j * bands + b] = residual;
            if (residual == 0.0 || time == 0.0)
            {
                log_bound[b] = -std::numeric_limits<double>::infinity();
            }
            else
            {
                log_bound[b] += std::log(time) + std::log(residual) - std::log(static_cast<double>(m));
            }
            // Jawecki et al., Theorem 1: |t|^m product(beta_1,...,beta_m)/m!.
            // Duhamel adds |t|*||F|| for the discarded reorthogonalization terms.
            errors[b] = std::exp(log_bound[b]) + time * std::sqrt(correction_sq[b]);
            if (errors[b] <= rel_tol)
            {
                candidates.push_back(b);
            }
            else if (residual <= std::numeric_limits<Real>::min())
            {
                return failed_band(result, b, true);
            }
        }
        if (!candidates.empty())
        {
            const std::vector<std::complex<double>> gram = basis_gram(layout, m, q);
            // The defect bound follows the three-term relation, not Q-adjoint*Q=I.
            // Rebuilding for insignificant loss of orthogonality wastes H applications.
            // Check the actual output norm separately at the requested precision.
            const double roundoff = std::numeric_limits<Real>::epsilon();
            const double orth_tol = std::max(std::sqrt(roundoff), 0.1 * rel_tol);
            const double norm_tol = std::max(100.0 * m * std::numeric_limits<Real>::epsilon(), 0.1 * rel_tol);
            for (int b: candidates)
            {
                double defect = 0.0;
                for (int row = 0; row < m; ++row)
                {
                    double row_sum = 0.0;
                    for (int col = 0; col < m; ++col)
                    {
                        const Wide identity = row == col ? Wide(1) : Wide(0);
                        const double value = std::abs(gram[(col * m + row) * bands + b] - identity);
                        if (!std::isfinite(value))
                        {
                            return failed_band(result, b, false);
                        }
                        row_sum += value;
                    }
                    defect = std::max(defect, row_sum);
                }
                if (defect > orth_tol)
                {
                    *lost_orthogonality = true;
                    return failed_band(result, b, true);
                }
            }
            std::vector<T> coefficients(m * bands, T(0));
            std::vector<T> full_coefficients;
            if (full)
            {
                full_coefficients.assign(m * bands, T(0));
            }
            for (int b: candidates)
            {
                if (!diagonalize(diagonal, beta, b, bands, m))
                {
                    return failed_band(result, b, false);
                }
                std::vector<Wide> half_values(m, Wide(0));
                std::vector<Wide> full_values(m, Wide(0));
                for (int k = 0; k < m; ++k)
                {
                    const double weight = scale[b] * small_vectors_[k * m];
                    const Wide half_phase = weight * std::exp(Wide(0, -alpha * small_values_[k]));
                    Wide full_phase(0);
                    if (full)
                    {
                        full_phase = weight * std::exp(Wide(0, -2.0 * alpha * small_values_[k]));
                    }
                    for (int i = 0; i < m; ++i)
                    {
                        const double component = small_vectors_[k * m + i];
                        half_values[i] += component * half_phase;
                        if (full)
                        {
                            full_values[i] += component * full_phase;
                        }
                    }
                }
                for (int i = 0; i < m; ++i)
                {
                    coefficients[i * bands + b] = T(half_values[i]);
                    if (full)
                    {
                        full_coefficients[i * bands + b] = T(full_values[i]);
                    }
                }
            }
            combine(layout, m, q, coefficients, output);
            if (full)
            {
                combine(layout, m, q, full_coefficients, full);
            }
            const std::vector<std::complex<double>> output_norms = this->algebra_.dots(ld, dim, bands, 1, stride, output, output);
            std::vector<Wide> full_norms;
            if (full)
            {
                full_norms = this->algebra_.dots(ld, dim, bands, 1, stride, full, full);
            }
            for (int b: candidates)
            {
                const double norm_error = std::abs(std::sqrt(output_norms[b].real()) / scale[b] - 1.0);
                double full_error = 0.0;
                if (full)
                {
                    full_error = std::abs(std::sqrt(full_norms[b].real()) / scale[b] - 1.0);
                }
                if (!std::isfinite(norm_error) || !std::isfinite(full_error))
                {
                    return failed_band(result, b, false);
                }
                if (norm_error > norm_tol || full_error > norm_tol)
                {
                    *lost_orthogonality = true;
                    return failed_band(result, b, true);
                }
                done[b] = true;
                result.rel_err_est = std::max(result.rel_err_est, errors[b]);
            }
        }
        bool finished = true;
        for (int b = 0; b < bands; ++b)
        {
            factors[b] = T(0);
            if (!done[b])
            {
                finished = false;
                factors[b] = T(1.0 / beta[j * bands + b]);
                result.failed_band = b;
            }
        }
        if (finished)
        {
            result.failed_band = -1;
            return result;
        }
        if (m < limit)
        {
            T* next = q + static_cast<int64_t>(m) * stride;
            base_device::memory::set_memory_op<T, Device>()(next, 0, stride);
            this->work_.batch(ld, dim, bands, next, w, w, T(1), T(0), factors.data(), nullptr, nullptr);
        }
    }
    result.state = ExpState::retry;
    result.rel_err_est = *std::max_element(errors.begin(), errors.end());
    return result;
}

template <typename T, typename Device>
ExpResult ExpLanczos<T, Device>::interval(const HSOperator<T, Device>& op,
                                          const OrbLayout& layout,
                                          double alpha,
                                          double rel_tol,
                                          const T* input,
                                          T* output,
                                          T* full,
                                          int depth)
{
    bool lost_orthogonality = false;
    ExpResult trial = attempt(op, layout, alpha, rel_tol, input, output, full, false, &lost_orthogonality);
    if (lost_orthogonality)
    {
        const ExpResult previous = trial;
        trial = attempt(op, layout, alpha, rel_tol, input, output, full, true, &lost_orthogonality);
        add_work(trial, previous);
    }
    if (trial.state != ExpState::retry)
    {
        return trial;
    }
    if (depth == 10)
    {
        trial.state = ExpState::failed;
        return trial;
    }
    ++trial.bisect_count;
    ct::Tensor midpoint;
    const int64_t size = static_cast<int64_t>(layout.ld) * layout.bands;
    T* middle = output;
    T* endpoint = full;
    double sub_alpha = alpha;
    if (!full)
    {
        linear_buffer<T, Device>(&midpoint, size);
        middle = midpoint.template data<T>();
        endpoint = output;
        sub_alpha = alpha / 2.0;
    }
    const double sub_tol = rel_tol / 2.0;
    const int next_depth = depth + 1;
    ExpResult left = interval(op, layout, sub_alpha, sub_tol, input, middle, nullptr, next_depth);
    add_work(left, trial);
    if (left.state != ExpState::success)
    {
        return left;
    }
    ExpResult right = interval(op, layout, sub_alpha, sub_tol, middle, endpoint, nullptr, next_depth);
    add_work(right, left);
    if (right.state == ExpState::success)
    {
        right.rel_err_est += left.rel_err_est;
    }
    return right;
}

template <typename T, typename Device>
ExpResult ExpLanczos<T, Device>::evaluate(const HSOperator<T, Device>& op,
                                          const OrbLayout& layout,
                                          double alpha,
                                          const T* input,
                                          T* output,
                                          T* full)
{
    ModuleBase::timer::start("ExpLanczos", "evaluate");
    ExpResult result;
    result.has_err_est = true;
    const int batch_limit = batch_size(layout);
    for (int first = 0; first < layout.bands; first += batch_limit)
    {
        const int bands = std::min(batch_limit, layout.bands - first);
        const OrbLayout batch{layout.ld, layout.npw, layout.npol, bands};
        const int64_t offset = static_cast<int64_t>(first) * layout.ld;
        const T* src = input + offset;
        T* dst = output + offset;
        T* full_dst = full ? full + offset : nullptr;
        const ExpResult local = interval(op, batch, alpha, this->exp_ctrl_.rel_tol, src, dst, full_dst, 0);
        result.has_err_est = result.has_err_est && local.has_err_est;
        add_work(result, local);
        result.rel_err_est = std::max(result.rel_err_est, local.rel_err_est);
        if (local.state != ExpState::success)
        {
            result.state = local.state;
            result.failed_band = first + local.failed_band;
            break;
        }
    }
    ModuleBase::timer::end("ExpLanczos", "evaluate");
    return result;
}

template <typename T, typename Device>
ExpResult ExpLanczos<T, Device>::apply(const HSOperator<T, Device>& op, const OrbLayout& layout, double alpha, const T* input, T* output)
{
    return evaluate(op, layout, alpha, input, output, nullptr);
}

template <typename T, typename Device>
ExpResult ExpLanczos<T, Device>::apply_pair(const HSOperator<T, Device>& op,
                                            const OrbLayout& layout,
                                            double alpha,
                                            const T* input,
                                            T* half,
                                            T* full)
{
    return evaluate(op, layout, alpha, input, half, full);
}

template class ExpLanczos<std::complex<float>, base_device::DEVICE_CPU>;
template class ExpLanczos<std::complex<double>, base_device::DEVICE_CPU>;
#if defined(__CUDA) || defined(__ROCM)
template class ExpLanczos<std::complex<float>, base_device::DEVICE_GPU>;
template class ExpLanczos<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace hsolver
