#include "source_hsolver/module_rt_pw/exp_lanczos.h"

#include "source_base/module_container/base/third_party/lapack.h"
#include "source_base/timer.h"

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
    total.order_or_dim = std::max(total.order_or_dim, part.order_or_dim);
}

bool diagonalize(std::vector<Wide>* vectors, int n, std::vector<double>* values)
{
    const char job = 'V';
    const char upper = 'U';
    const int lwork = std::max(1, n * n + 2 * n);
    const int lrwork = std::max(1, 1 + 5 * n + 2 * n * n);
    const int liwork = std::max(1, 3 + 5 * n);
    std::vector<Wide> work(lwork);
    std::vector<double> rwork(lrwork);
    std::vector<int> iwork(liwork);
    values->resize(n);
    int info = 0;
    zheevd_(&job,
            &upper,
            &n,
            vectors->data(),
            &n,
            values->data(),
            work.data(),
            &lwork,
            rwork.data(),
            &lrwork,
            iwork.data(),
            &liwork,
            &info);
    return info == 0;
}

std::vector<double> gauss_rule(int order)
{
    std::vector<double> rule(2 * order);
    const double pi = std::acos(-1.0);
    for (int i = 0; i < order; ++i)
    {
        double x = std::cos(pi * (i + 0.75) / (order + 0.5));
        double derivative = 0.0;
        for (int iteration = 0; iteration < 20; ++iteration)
        {
            double p = 1.0;
            double previous = 0.0;
            for (int k = 1; k <= order; ++k)
            {
                const double next = ((2 * k - 1) * x * p - (k - 1) * previous) / k;
                previous = p;
                p = next;
            }
            derivative = order * (x * p - previous) / (x * x - 1);
            const double change = p / derivative;
            x -= change;
            if (std::abs(change) < 1e-15)
            {
                break;
            }
        }
        rule[i] = x;
        rule[order + i] = 2.0 / ((1.0 - x * x) * derivative * derivative);
    }
    return rule;
}

double defect_integral(const std::vector<double>& values,
                       const std::vector<Wide>& vectors,
                       double time,
                       int panels,
                       const std::vector<double>& rule,
                       double sign)
{
    const int n = values.size();
    const int order = rule.size() / 2;
    const double width = time / panels;
    double integral = 0.0;
    for (int panel = 0; panel < panels; ++panel)
    {
        for (int node = 0; node < order; ++node)
        {
            const double t = sign * (panel + 0.5 * (rule[node] + 1.0)) * width;
            Wide amplitude(0);
            for (int k = 0; k < n; ++k)
            {
                // Removing a scalar phase improves evaluation for shifted spectra.
                amplitude += vectors[n - 1 + k * n] * std::conj(vectors[k * n]) * std::exp(Wide(0, -t * (values[k] - values[0])));
            }
            integral += rule[order + node] * std::abs(amplitude);
        }
    }
    return 0.5 * width * integral;
}

double estimate_error(const std::vector<double>& values, const std::vector<Wide>& vectors, double alpha, double beta, double rel_tol)
{
    static const std::vector<double> rule8 = gauss_rule(8);
    static const std::vector<double> rule16 = gauss_rule(16);
    const double time = std::abs(alpha);
    const double sign = alpha < 0 ? -1.0 : 1.0;
    const double phase = time * (values.back() - values.front());
    const double requested = std::ceil(phase / (std::acos(-1.0) / 4.0));
    if (!std::isfinite(requested) || requested > 4096)
    {
        return std::numeric_limits<double>::infinity();
    }
    int panels = std::max(1, static_cast<int>(requested));
    for (; panels <= 4096; panels *= 2)
    {
        const double low = beta * defect_integral(values, vectors, time, panels, rule8, sign);
        const double high = beta * defect_integral(values, vectors, time, panels, rule16, sign);
        if (std::abs(low - high) <= rel_tol / 8.0)
        {
            return 2.0 * std::max(low, high);
        }
    }
    return std::numeric_limits<double>::infinity();
}
} // namespace

template <typename T, typename Device>
ExpLanczos<T, Device>::ExpLanczos(const ExpOpts& exp_opts, const diag_comm_info& comm) : ExpAction<T, Device>(exp_opts, comm)
{
}

template <typename T, typename Device>
ExpResult ExpLanczos<T, Device>::attempt(const HSOperator<T, Device>& op,
                                    const OrbLayout& layout,
                                    double alpha,
                                    double rel_tol,
                                    const T* input,
                                    T* output)
{
    ExpResult result;
    using Real = typename T::value_type;
    const int ld = layout.ld;
    const int bands = layout.bands;
    const int dim = layout.npol == 2 ? ld : layout.npw;
    const int stride = ld * bands;
    const int limit = this->exp_opts_.krylov_max;
    const int64_t size = std::max(1, stride);
    linear_buffer<T, Device>(&basis_, size * (limit + 1));
    linear_buffer<T, Device>(&image_, size);
    linear_buffer<T, Device>(&packed_, size);
    linear_buffer<T, Device>(&hpacked_, size);
    T* q = basis_.template data<T>();
    T* w = image_.template data<T>();
    T* packed = packed_.template data<T>();
    T* hpacked = hpacked_.template data<T>();
    base_device::memory::set_memory_op<T, Device>()(q, 0, size * (limit + 1));
    base_device::memory::set_memory_op<T, Device>()(output, 0, stride);
    base_device::memory::synchronize_memory_op<T, Device, Device>()(q, input, stride);
    zero_pad<T, Device>(q, layout);
    const auto norms = this->algebra_.dots(ld, dim, bands, 1, stride, q, q);
    std::vector<double> scale(bands);
    std::vector<T> factors(bands);
    std::vector<bool> done(bands, false);
    std::vector<Wide> matrix(static_cast<size_t>(bands) * limit * limit, Wide(0));
    result.has_err_est = true;
    for (int b = 0; b < bands; ++b)
    {
        scale[b] = std::sqrt(std::max(0.0, norms[b].real()));
        if (!std::isfinite(norms[b].real()))
        {
            result.state = ExpState::failed;
            result.has_err_est = false;
            result.failed_band = b;
            return result;
        }
        done[b] = scale[b] == 0.0;
        factors[b] = done[b] ? T(0) : T(1.0 / scale[b]);
    }
    this->work_.batch(ld, dim, bands, q, q, q, T(1), T(0), factors.data(), nullptr, nullptr);
    for (int j = 0; j < limit; ++j)
    {
        std::vector<int> active;
        T* column = q + static_cast<int64_t>(j) * stride;
        for (int b = 0; b < bands; ++b)
        {
            if (done[b])
            {
                continue;
            }
            T* destination = packed + active.size() * ld;
            base_device::memory::synchronize_memory_op<T, Device, Device>()(destination, column + b * ld, ld);
            active.push_back(b);
        }
        if (active.empty())
        {
            result.failed_band = -1;
            return result;
        }
        const int count = active.size();
        op.hpsi(packed, hpacked, ld, count);
        const OrbLayout batch_layout{ld, layout.npw, layout.npol, count};
        zero_pad<T, Device>(hpacked, batch_layout);
        ++result.hpsi_calls;
        result.hpsi_cols += count;
        base_device::memory::set_memory_op<T, Device>()(w, 0, stride);
        for (int a = 0; a < count; ++a)
        {
            base_device::memory::synchronize_memory_op<T, Device, Device>()(w + active[a] * ld, hpacked + a * ld, ld);
        }
        // Two complete projections retain small off-tridiagonal corrections.
        for (int pass = 0; pass < 2; ++pass)
        {
            const auto products = this->algebra_.dots(ld, dim, bands, j + 1, stride, q, w);
            for (int i = 0; i <= j; ++i)
            {
                for (int b = 0; b < bands; ++b)
                {
                    const Wide value = products[i * bands + b];
                    matrix[(static_cast<size_t>(b) * limit + j) * limit + i] += value;
                    factors[b] = T(-value);
                }
                const T* qi = q + static_cast<int64_t>(i) * stride;
                this->work_.batch(ld, dim, bands, w, w, qi, T(1), T(1), nullptr, factors.data(), nullptr);
            }
        }
        const auto residuals = this->algebra_.dots(ld, dim, bands, 1, stride, w, w);
        const int m = j + 1;
        result.order_or_dim = std::max(result.order_or_dim, m);
        double max_error = 0.0;
        for (int b: active)
        {
            const double beta = std::sqrt(std::max(0.0, residuals[b].real()));
            if (!std::isfinite(residuals[b].real()))
            {
                result.state = ExpState::failed;
                result.has_err_est = false;
                result.failed_band = b;
                return result;
            }
            const bool breakdown = beta <= 100.0 * std::numeric_limits<Real>::epsilon();
            const bool evaluate = breakdown || m == limit || (m >= 4 && m % 2 == 0);
            if (evaluate)
            {
                std::vector<Wide> vectors(m * m);
                const Wide* block = matrix.data() + static_cast<size_t>(b) * limit * limit;
                for (int col = 0; col < m; ++col)
                {
                    for (int row = 0; row <= col; ++row)
                    {
                        // The upper triangle contains explicitly computed projections.
                        Wide value = block[row + col * limit];
                        if (row == col)
                        {
                            value = Wide(value.real(), 0);
                        }
                        vectors[row + col * m] = value;
                        vectors[col + row * m] = std::conj(value);
                    }
                }
                std::vector<double> values;
                if (!diagonalize(&vectors, m, &values))
                {
                    result.state = ExpState::failed;
                    result.has_err_est = false;
                    result.failed_band = b;
                    return result;
                }
                const double error = estimate_error(values, vectors, alpha, beta, rel_tol);
                max_error = std::max(max_error, error);
                if (error <= rel_tol)
                {
                    std::vector<Wide> coefficients(m, Wide(0));
                    for (int k = 0; k < m; ++k)
                    {
                        const Wide phase = scale[b] * std::exp(Wide(0, -alpha * values[k])) * std::conj(vectors[k * m]);
                        for (int i = 0; i < m; ++i)
                        {
                            coefficients[i] += vectors[i + k * m] * phase;
                        }
                    }
                    this->algebra_.expand(stride, dim, m, 1, q + b * ld, coefficients, output + b * ld, T(0));
                    result.rel_err_est = std::max(result.rel_err_est, error);
                    done[b] = true;
                }
            }
            factors[b] = done[b] || beta == 0.0 ? T(0) : T(1.0 / beta);
        }
        for (int b = 0; b < bands; ++b)
        {
            if (done[b])
            {
                factors[b] = T(0);
            }
        }
        T* next = q + static_cast<int64_t>(m) * stride;
        this->work_.batch(ld, dim, bands, next, w, w, T(1), T(0), factors.data(), nullptr, nullptr);
        if (m == limit)
        {
            result.rel_err_est = std::max(result.rel_err_est, max_error);
        }
    }
    for (int b = 0; b < bands; ++b)
    {
        if (!done[b])
        {
            result.state = ExpState::retry;
            result.failed_band = b;
            return result;
        }
    }
    result.failed_band = -1;
    return result;
}

template <typename T, typename Device>
ExpResult ExpLanczos<T, Device>::interval(const HSOperator<T, Device>& op,
                                     const OrbLayout& layout,
                                     double alpha,
                                     double rel_tol,
                                     const T* input,
                                     T* output,
                                     int depth)
{
    ExpResult trial = attempt(op, layout, alpha, rel_tol, input, output);
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
    linear_buffer<T, Device>(&midpoint, size);
    T* middle = midpoint.template data<T>();
    const double half_alpha = alpha / 2.0;
    const double half_tol = rel_tol / 2.0;
    const int next_depth = depth + 1;
    ExpResult left = interval(op, layout, half_alpha, half_tol, input, middle, next_depth);
    add_work(left, trial);
    if (left.state != ExpState::success)
    {
        return left;
    }
    ExpResult right = interval(op, layout, half_alpha, half_tol, middle, output, next_depth);
    add_work(right, left);
    if (right.state == ExpState::success)
    {
        right.rel_err_est += left.rel_err_est;
    }
    return right;
}

template <typename T, typename Device>
ExpResult ExpLanczos<T, Device>::apply(const HSOperator<T, Device>& op, const OrbLayout& layout, double alpha, const T* input, T* output)
{
    ModuleBase::timer::start("ExpLanczos", "apply");
    ExpResult result;
    result.has_err_est = true;
    // Limit basis storage independently of the number of occupied bands.
    for (int first = 0; first < layout.bands; first += 8)
    {
        const int bands = std::min(8, layout.bands - first);
        const OrbLayout batch{layout.ld, layout.npw, layout.npol, bands};
        const int64_t offset = static_cast<int64_t>(first) * layout.ld;
        const T* src = input + offset;
        T* dst = output + offset;
        const ExpResult local = interval(op, batch, alpha, this->exp_opts_.rel_tol, src, dst, 0);
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
    ModuleBase::timer::end("ExpLanczos", "apply");
    return result;
}

template class ExpLanczos<std::complex<float>, base_device::DEVICE_CPU>;
template class ExpLanczos<std::complex<double>, base_device::DEVICE_CPU>;
#if defined(__CUDA) || defined(__ROCM)
template class ExpLanczos<std::complex<float>, base_device::DEVICE_GPU>;
template class ExpLanczos<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace hsolver
