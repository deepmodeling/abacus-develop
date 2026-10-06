#include "source_base/kernels/math_kernel_op.h"
#include "source_base/module_container/base/third_party/lapack.h"
#include "source_base/timer.h"
#include "source_hsolver/orthonormal.h"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace hsolver
{
namespace
{
using Wide = std::complex<double>;

std::vector<Wide> identity(int n)
{
    std::vector<Wide> result(n * n, Wide(0));
    for (int i = 0; i < n; ++i)
    {
        result[i + i * n] = Wide(1);
    }
    return result;
}

std::vector<Wide> multiply(const std::vector<Wide>& a, const std::vector<Wide>& b, int n)
{
    std::vector<Wide> result(n * n);
    const Wide one(1);
    const Wide zero(0);
    ModuleBase::gemm_op<Wide, base_device::DEVICE_CPU>()('N', 'N', n, n, n, &one, a.data(), n, b.data(), n, &zero, result.data(), n);
    return result;
}

bool cholesky(const std::vector<Wide>& g, int n, std::vector<Wide>* c)
{
    ModuleBase::timer::start("Orthonormal", "cholesky");
    *c = g;
    const char upper = 'U';
    const char nonunit = 'N';
    int info = 0;
    zpotrf_(&upper, &n, c->data(), &n, &info);
    if (info == 0)
    {
        ztrtri_(&upper, &nonunit, &n, c->data(), &n, &info);
    }
    for (int j = 0; j < n; ++j)
    {
        for (int i = j + 1; i < n; ++i)
        {
            (*c)[i + j * n] = Wide(0);
        }
    }
    ModuleBase::timer::end("Orthonormal", "cholesky");
    return info == 0;
}

bool lowdin(const std::vector<Wide>& g, int n, std::vector<Wide>* c)
{
    ModuleBase::timer::start("Orthonormal", "lowdin");
    std::vector<Wide> u(g);
    std::vector<double> values(n);
    const char vectors = 'V';
    const char upper = 'U';
    const int lwork = std::max(1, n * n + 2 * n);
    const int lrwork = std::max(1, 1 + 5 * n + 2 * n * n);
    const int liwork = std::max(1, 3 + 5 * n);
    std::vector<Wide> work(lwork);
    std::vector<double> rwork(lrwork);
    std::vector<int> iwork(liwork);
    int info = 0;
    zheevd_(&vectors, &upper, &n, u.data(), &n, values.data(), work.data(), &lwork, rwork.data(), &lrwork, iwork.data(), &liwork, &info);
    bool valid = info == 0;
    for (double value: values)
    {
        valid = valid && std::isfinite(value) && value > 0.0;
    }
    if (valid)
    {
        std::vector<Wide> scaled(u);
        for (int j = 0; j < n; ++j)
        {
            for (int i = 0; i < n; ++i)
            {
                scaled[i + j * n] /= std::sqrt(values[j]);
            }
        }
        c->resize(n * n);
        const Wide one(1);
        const Wide zero(0);
        ModuleBase::gemm_op<Wide, base_device::DEVICE_CPU>()('N', 'C', n, n, n, &one, scaled.data(), n, u.data(), n, &zero, c->data(), n);
    }
    ModuleBase::timer::end("Orthonormal", "lowdin");
    return valid;
}

bool newton_schulz(const std::vector<Wide>& g, int n, std::vector<Wide>* c)
{
    ModuleBase::timer::start("Orthonormal", "newton_schulz");
    constexpr double iteration_tol = 1e-12;
    constexpr double max_deviation = 0.1;
    constexpr int max_updates = 5;
    const std::vector<Wide> unit = identity(n);
    double bound = 0.0;
    for (int i = 0; i < n; ++i)
    {
        double sum = 0.0;
        for (int j = 0; j < n; ++j)
        {
            sum += std::abs(g[i + j * n] - unit[i + j * n]);
        }
        bound = std::max(bound, sum);
    }
    bool valid = false;
    *c = unit;
    if (bound <= max_deviation)
    {
        for (int iteration = 0; iteration <= max_updates; ++iteration)
        {
            const std::vector<Wide> square = multiply(*c, *c, n);
            std::vector<Wide> residual = multiply(g, square, n);
            if (orth_error(residual, n) <= iteration_tol)
            {
                valid = true;
                break;
            }
            if (iteration == max_updates)
            {
                break;
            }
            for (int i = 0; i < n * n; ++i)
            {
                residual[i] = 0.5 * (3.0 * unit[i] - residual[i]);
            }
            *c = multiply(*c, residual, n);
        }
    }
    ModuleBase::timer::end("Orthonormal", "newton_schulz");
    return valid;
}
} // namespace

const char* orth_method_name(OrthMethod method)
{
    switch (method)
    {
    case OrthMethod::cholesky:
        return "cholesky";
    case OrthMethod::lowdin:
        return "lowdin";
    case OrthMethod::newton_schulz:
        return "newton_schulz";
    default:
        return "none";
    }
}

OrthMethod parse_orth_method(const std::string& name)
{
    for (OrthMethod method: {OrthMethod::none, OrthMethod::cholesky, OrthMethod::lowdin, OrthMethod::newton_schulz})
    {
        if (name == orth_method_name(method))
        {
            return method;
        }
    }
    throw std::invalid_argument("Unsupported orthonormalization method: " + name);
}

double orth_error(const std::vector<Wide>& g, int n)
{
    double error = 0.0;
    for (int j = 0; j < n; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            const Wide target = i == j ? Wide(1) : Wide(0);
            const double value = std::abs(g[i + j * n] - target);
            if (!std::isfinite(value))
            {
                return std::numeric_limits<double>::infinity();
            }
            error = std::max(error, value);
        }
    }
    return error;
}

const char* orth_failure_name(OrthFailure failure)
{
    switch (failure)
    {
    case OrthFailure::nonfinite_gram:
        return "nonfinite Gram matrix";
    case OrthFailure::nonfinite_norm:
        return "nonfinite orbital norm";
    case OrthFailure::nonpositive_norm:
        return "nonpositive orbital norm";
    case OrthFailure::factorization_failed:
        return "factorization/iteration failed";
    case OrthFailure::invalid_candidate:
        return "nonfinite or zero-norm candidate";
    case OrthFailure::no_improvement:
        return "candidate did not improve orthogonality";
    case OrthFailure::tolerance_not_met:
        return "orthogonality tolerance not met";
    default:
        return "none";
    }
}

bool orth_transform(const std::vector<Wide>& g, int n, OrthMethod method, std::vector<Wide>* c)
{
    ModuleBase::timer::start("Orthonormal", "orth_transform");
    bool valid = false;
    switch (method)
    {
    case OrthMethod::cholesky:
        valid = cholesky(g, n, c);
        break;
    case OrthMethod::lowdin:
        valid = lowdin(g, n, c);
        break;
    case OrthMethod::newton_schulz:
        valid = newton_schulz(g, n, c);
        break;
    default:
        valid = false;
    }
    for (const Wide& value: *c)
    {
        valid = valid && std::isfinite(std::abs(value));
    }
    ModuleBase::timer::end("Orthonormal", "orth_transform");
    return valid;
}
} // namespace hsolver
