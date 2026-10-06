#include "source_hsolver/orthonormal.h"

#include "source_base/parallel_device.h"
#include "source_base/timer.h"

#include <cmath>

namespace hsolver
{
namespace
{
bool positive_norms(const std::vector<std::complex<double>>& gram, int bands)
{
    for (int band = 0; band < bands; ++band)
    {
        if (!(gram[band + band * bands].real() > 0.0))
        {
            return false;
        }
    }
    return true;
}
} // namespace

template <typename T, typename Device>
Orthonormal<T, Device>::Orthonormal(const diag_comm_info& comm) : comm_(comm), algebra_(comm)
{
}

template <typename T, typename Device>
std::vector<std::complex<double>> Orthonormal<T, Device>::gram(const T* input, int ld, int dim, int bands)
{
    ModuleBase::timer::start("Orthonormal", "gram");
    std::vector<std::complex<double>> result = algebra_.gram(ld, dim, bands, input);
    ModuleBase::timer::end("Orthonormal", "gram");
    return result;
}

template <typename T, typename Device>
bool Orthonormal<T, Device>::factor(const std::vector<std::complex<double>>& g,
                                    int bands,
                                    OrthMethod method,
                                    std::vector<std::complex<double>>* transform)
{
    ModuleBase::timer::start("Orthonormal", "factor");
    bool valid = false;
    if (comm_.rank == 0)
    {
        valid = orth_transform(g, bands, method, transform);
    }
    double status = static_cast<double>(valid);
#ifdef __MPI
    if (comm_.nproc > 1)
    {
        Parallel_Common::bcast_data(&status, 1, comm_.comm, 0);
    }
#endif
    valid = status != 0.0;
    if (valid)
    {
        transform->resize(static_cast<std::size_t>(bands) * bands);
#ifdef __MPI
        if (comm_.nproc > 1)
        {
            // Every partition must apply the same coefficients and take the same recovery branch.
            const int elements = bands * bands;
            Parallel_Common::bcast_data(transform->data(), elements, comm_.comm, 0);
        }
#endif
    }
    ModuleBase::timer::end("Orthonormal", "factor");
    return valid;
}

template <typename T, typename Device>
void Orthonormal<T, Device>::rotate(const T* input,
                                    T* output,
                                    int ld,
                                    int dim,
                                    int bands,
                                    const std::vector<std::complex<double>>& transform)
{
    ModuleBase::timer::start("Orthonormal", "rotate");
    std::vector<std::complex<double>> delta(transform);
    for (int i = 0; i < bands; ++i)
    {
        delta[i + i * bands] -= 1.0;
    }
    if (dim > 0)
    {
        base_device::memory::synchronize_memory_2d_op<T, Device, Device>()(output, ld, input, ld, dim, bands);
    }
    algebra_.expand(ld, dim, bands, bands, input, delta, output, T(1));
    ModuleBase::timer::end("Orthonormal", "rotate");
}

template <typename T, typename Device>
bool Orthonormal<T, Device>::try_candidate(T* input,
                                           int ld,
                                           int dim,
                                           int bands,
                                           const std::vector<std::complex<double>>& transform,
                                           std::vector<std::complex<double>>* g,
                                           OrthResult* result)
{
    ModuleBase::timer::start("Orthonormal", "try_candidate");
    T* candidate = candidate_.template data<T>();
    rotate(input, candidate, ld, dim, bands, transform);
    ++result->passes;
    const std::vector<std::complex<double>> check = gram(candidate, ld, dim, bands);
    const double error = orth_error(check, bands);
    const bool valid = std::isfinite(error) && positive_norms(check, bands);
    const bool improved = valid && error < result->after;
    if (improved)
    {
        if (dim > 0)
        {
            base_device::memory::synchronize_memory_2d_op<T, Device, Device>()(input, ld, candidate, ld, dim, bands);
        }
        result->after = error;
        *g = check;
        result->failure = OrthFailure::tolerance_not_met;
    }
    else if (!valid || result->after > orth_tolerance<T>())
    {
        result->failure = valid ? OrthFailure::no_improvement : OrthFailure::invalid_candidate;
        result->reason += std::string(orth_method_name(result->actual)) + ": " + orth_failure_name(result->failure) + "; ";
        ++result->rejected;
    }
    ModuleBase::timer::end("Orthonormal", "try_candidate");
    return improved;
}

template <typename T, typename Device>
void Orthonormal<T, Device>::correct(T* input,
                                     int ld,
                                     int dim,
                                     int bands,
                                     OrthMethod method,
                                     std::vector<std::complex<double>>* g,
                                     OrthResult* result)
{
    ModuleBase::timer::start("Orthonormal", "correct");
    std::vector<OrthMethod> methods{method};
    if (method != OrthMethod::cholesky)
    {
        methods.push_back(OrthMethod::cholesky);
    }
    if (method != OrthMethod::lowdin)
    {
        methods.push_back(OrthMethod::lowdin);
    }
    const int64_t elements = static_cast<int64_t>(ld) * bands;
    linear_buffer<T, Device>(&candidate_, elements);
    bool changed = false;
    std::vector<std::complex<double>> transform;
    for (int pass = 0; pass < 2; ++pass)
    {
        bool improved = false;
        for (std::size_t attempt = 0; attempt < methods.size(); ++attempt)
        {
            result->actual = methods[attempt];
            if (attempt > 0)
            {
                ++result->fallbacks;
            }
            const bool factored = factor(*g, bands, result->actual, &transform);
            if (factored)
            {
                improved = try_candidate(input, ld, dim, bands, transform, g, result);
            }
            else
            {
                result->failure = OrthFailure::factorization_failed;
                result->reason += std::string(orth_method_name(result->actual)) + ": " + orth_failure_name(result->failure) + "; ";
            }
            if (improved || result->after <= orth_tolerance<T>())
            {
                break;
            }
        }
        changed = changed || improved;
        if (result->after <= orth_tolerance<T>())
        {
            result->status = changed ? OrthStatus::accepted : OrthStatus::unchanged;
            result->failure = OrthFailure::none;
            break;
        }
        // Another pass is useful only after an accepted update changed the Gram matrix.
        if (!improved)
        {
            break;
        }
    }
    ModuleBase::timer::end("Orthonormal", "correct");
}

template <typename T, typename Device>
OrthResult Orthonormal<T, Device>::inspect(const T* input, int ld, int dim, int bands, bool full_gram, bool collect_norms)
{
    ModuleBase::timer::start("Orthonormal", "inspect");
    OrthResult result;
    result.gram_checked = full_gram;
    if (full_gram)
    {
        const std::vector<std::complex<double>> g = gram(input, ld, dim, bands);
        result.before = orth_error(g, bands);
        result.after = result.before;
        for (int band = 0; collect_norms && band < bands; ++band)
        {
            result.norms.push_back(g[band + band * bands].real());
        }
        if (!std::isfinite(result.after))
        {
            result.failure = OrthFailure::nonfinite_gram;
        }
        else if (!positive_norms(g, bands))
        {
            result.failure = OrthFailure::nonpositive_norm;
        }
        else
        {
            result.status = OrthStatus::inspected;
        }
    }
    else
    {
        // Corresponding-band products promote operands before multiplication, without a full FP64 copy.
        std::vector<std::complex<double>> norms;
        if (bands > 0)
        {
            norms = algebra_.dots(ld, dim, bands, 1, 0, input, input);
        }
        result.status = OrthStatus::disabled;
        for (int band = 0; band < bands; ++band)
        {
            const std::complex<double> value = norms[band];
            if (collect_norms)
            {
                result.norms.push_back(value.real());
            }
            if (result.status == OrthStatus::failed)
            {
                continue;
            }
            if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
            {
                result.failure = OrthFailure::nonfinite_norm;
            }
            else if (value.real() <= 0.0)
            {
                result.failure = OrthFailure::nonpositive_norm;
            }
            if (result.failure != OrthFailure::none)
            {
                result.status = OrthStatus::failed;
                result.reason = "band=" + std::to_string(band);
            }
        }
    }
    ModuleBase::timer::end("Orthonormal", "inspect");
    return result;
}

template <typename T, typename Device>
OrthResult Orthonormal<T, Device>::apply(T* input, int ld, int dim, int bands, OrthMethod method, bool collect_norms)
{
    ModuleBase::timer::start("Orthonormal", "apply");
    if (method == OrthMethod::none)
    {
        const OrthResult result = inspect(input, ld, dim, bands, false, collect_norms);
        ModuleBase::timer::end("Orthonormal", "apply");
        return result;
    }
    OrthResult result;
    result.gram_checked = true;
    std::vector<std::complex<double>> g = gram(input, ld, dim, bands);
    result.before = orth_error(g, bands);
    result.after = result.before;
    if (!std::isfinite(result.before))
    {
        result.failure = OrthFailure::nonfinite_gram;
    }
    else if (!positive_norms(g, bands))
    {
        result.failure = OrthFailure::nonpositive_norm;
    }
    else if (bands == 0)
    {
        result.status = OrthStatus::unchanged;
    }
    else
    {
        correct(input, ld, dim, bands, method, &g, &result);
    }
    for (int i = 0; collect_norms && i < bands; ++i)
    {
        result.norms.push_back(g[i + i * bands].real());
    }
    ModuleBase::timer::end("Orthonormal", "apply");
    return result;
}

template class Orthonormal<std::complex<float>, base_device::DEVICE_CPU>;
template class Orthonormal<std::complex<double>, base_device::DEVICE_CPU>;
#if defined(__CUDA) || defined(__ROCM)
template class Orthonormal<std::complex<float>, base_device::DEVICE_GPU>;
template class Orthonormal<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace hsolver
