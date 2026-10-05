#include "source_hsolver/orthonormal.h"

#include "source_base/parallel_device.h"
#include "source_base/timer.h"
#include "source_hsolver/kernels/linear_op.h"

#include <chrono>
#include <cmath>
#include <type_traits>

namespace hsolver
{
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
                                    std::vector<std::complex<double>>* transform,
                                    OrthResult* result)
{
    ModuleBase::timer::start("Orthonormal", "factor");
    OrthResult local;
    bool valid = false;
    if (comm_.rank == 0)
    {
        valid = orth_transform(g, bands, method, transform, &local);
    }
    double status[] = {static_cast<double>(valid),
                       static_cast<double>(local.actual),
                       static_cast<double>(local.fallbacks),
                       static_cast<double>(local.events[0]),
                       static_cast<double>(local.events[1]),
                       static_cast<double>(local.events[2])};
#ifdef __MPI
    if (comm_.nproc > 1)
    {
        Parallel_Common::bcast_data(status, 6, comm_.comm, 0);
    }
#endif
    valid = status[0] != 0.0;
    result->actual = static_cast<OrthMethod>(static_cast<int>(status[1]));
    result->fallbacks += static_cast<int>(status[2]);
    for (int i = 0; i < 3; ++i)
    {
        result->events[i] += static_cast<int>(status[i + 3]);
    }
    if (status[2] > 0.0)
    {
        const int index = 4 + static_cast<int>(result->actual);
        ++result->events[index];
    }
    result->reason += local.reason;
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
OrthResult Orthonormal<T, Device>::apply(T* input, int ld, int dim, int bands, OrthMethod method)
{
    ModuleBase::timer::start("Orthonormal", "apply");
    const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    OrthResult result;
    std::vector<std::complex<double>> g = gram(input, ld, dim, bands);
    result.before = orth_error(g, bands);
    result.after = result.before;
    const bool finite_gram = std::isfinite(result.before);
    if (!finite_gram)
    {
        // Distinguish nonfinite coefficients from overflow of otherwise finite products.
        double invalid = 0.0;
        std::vector<T> column(dim);
        for (int band = 0; band < bands && dim > 0; ++band)
        {
            const T* source = input + static_cast<int64_t>(band) * ld;
            base_device::memory::synchronize_memory_op<T, base_device::DEVICE_CPU, Device>()(column.data(), source, dim);
            for (const T& value: column)
            {
                if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
                {
                    invalid = 1.0;
                }
            }
        }
#ifdef __MPI
        Parallel_Common::reduce_data(&invalid, 1, comm_.comm);
#endif
        result.invalid_input = invalid > 0.0;
        result.skipped = method != OrthMethod::none;
        result.reason = "nonfinite Gram matrix; ";
        ++result.events[3];
    }
    const double threshold = std::is_same<T, std::complex<double>>::value ? 1e-12 : 1e-6;
    if (method != OrthMethod::none && finite_gram && bands > 0)
    {
        const int64_t elements = static_cast<int64_t>(ld) * bands;
        linear_buffer<T, Device>(&candidate_, elements);
        T* candidate = candidate_.template data<T>();
        for (int pass = 0; pass < 2; ++pass)
        {
            std::vector<std::complex<double>> c;
            if (!factor(g, bands, method, &c, &result))
            {
                result.skipped = true;
                break;
            }
            rotate(input, candidate, ld, dim, bands, c);
            ++result.passes;
            const std::vector<std::complex<double>> check = gram(candidate, ld, dim, bands);
            const double error = orth_error(check, bands);
            bool nonzero_columns = true;
            for (int band = 0; band < bands; ++band)
            {
                nonzero_columns = nonzero_columns && check[band + band * bands].real() > 0.0;
            }
            if (std::isfinite(error) && nonzero_columns && error < result.after)
            {
                if (dim > 0)
                {
                    base_device::memory::synchronize_memory_2d_op<T, Device, Device>()(input, ld, candidate, ld, dim, bands);
                }
                result.after = error;
                g = check;
            }
            else if (result.after > threshold)
            {
                result.skipped = true;
                result.reason += "candidate did not improve finite input; ";
                ++result.events[4];
            }
            if (result.after <= threshold)
            {
                break;
            }
        }
    }
    for (int i = 0; i < bands; ++i)
    {
        result.norms.push_back(g[i + i * bands].real());
    }
    linear_op<T, Device>().synchronize();
    const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::vector<double> times(comm_.nproc, 0.0);
    times[comm_.rank] = elapsed;
#ifdef __MPI
    Parallel_Common::reduce_data(times.data(), times.size(), comm_.comm);
#endif
    result.seconds = *std::max_element(times.begin(), times.end());
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
