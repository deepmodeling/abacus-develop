#ifndef STERNHEIMER_CHANNEL_PARALLEL_H
#define STERNHEIMER_CHANNEL_PARALLEL_H

#include <exception>
#include <stdexcept>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ModuleRI
{

template <typename Result, typename Worker>
std::vector<Result> run_sternheimer_channel_tasks(const int num_channels,
                                                  Worker worker,
                                                  const int requested_threads = 0)
{
    if (num_channels < 0)
    {
        throw std::invalid_argument("Sternheimer channel task count must be non-negative.");
    }
    if (requested_threads < 0)
    {
        throw std::invalid_argument("Sternheimer channel thread count must be non-negative.");
    }

    std::vector<Result> results(static_cast<std::size_t>(num_channels));
    std::vector<std::exception_ptr> errors(static_cast<std::size_t>(num_channels));

#ifdef _OPENMP
    const int worker_threads = requested_threads > 0 ? requested_threads : omp_get_max_threads();
#pragma omp parallel for schedule(dynamic) num_threads(worker_threads)
#endif
    for (int channel_index = 0; channel_index < num_channels; ++channel_index)
    {
        try
        {
            results[static_cast<std::size_t>(channel_index)] = worker(channel_index);
        }
        catch (...)
        {
            errors[static_cast<std::size_t>(channel_index)] = std::current_exception();
        }
    }

    for (const std::exception_ptr& error: errors)
    {
        if (error)
        {
            std::rethrow_exception(error);
        }
    }
    return results;
}

} // namespace ModuleRI

#endif
