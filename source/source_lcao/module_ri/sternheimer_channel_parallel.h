#ifndef STERNHEIMER_CHANNEL_PARALLEL_H
#define STERNHEIMER_CHANNEL_PARALLEL_H

#include <algorithm>
#include <exception>
#include <stdexcept>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ModuleRI
{

template <typename Result, typename Worker>
std::vector<Result> run_sternheimer_channel_tasks(const int num_channels, Worker worker, const int max_workers = 0)
{
    if (num_channels < 0)
    {
        throw std::invalid_argument("Sternheimer channel task count must be non-negative.");
    }
    if (max_workers < 0)
    {
        throw std::invalid_argument("Sternheimer channel maximum worker count must be non-negative.");
    }

    std::vector<Result> results(static_cast<std::size_t>(num_channels));
    std::vector<std::exception_ptr> errors(static_cast<std::size_t>(num_channels));
    if (num_channels == 0)
    {
        return results;
    }

#ifdef _OPENMP
    const int worker_count
        = max_workers > 0 ? std::min(num_channels, max_workers) : std::min(num_channels, omp_get_max_threads());
#pragma omp parallel for schedule(dynamic) num_threads(worker_count)
#endif
    for (int channel_index = 0; channel_index < num_channels; ++channel_index)
    {
        // A worker may run concurrently; shared captures must be read-only or indexed by channel.
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
