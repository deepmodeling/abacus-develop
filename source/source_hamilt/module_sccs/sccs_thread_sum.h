#ifndef SCCS_THREAD_SUM_H
#define SCCS_THREAD_SUM_H

#include <array>
#include <cstddef>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ModuleSccs
{
// Sums over the grid points [0, size): each thread runs accumulate(i, partial)
// over a fixed static share of the points in index order, starting from zero,
// and the partials are added in thread order. Unlike an OpenMP reduction, whose
// partials meet in run-dependent order, a given thread count gives the same
// sums on every run, and one thread gives the serial loop. The sqrt-CG
// iteration count depends on the last bits of its inner products.
template <std::size_t Count, typename Accumulate>
std::array<double, Count> thread_sums(std::size_t size, const Accumulate& accumulate)
{
#ifdef _OPENMP
    const int thread_count = omp_get_max_threads();
#else
    const int thread_count = 1;
#endif
    std::array<double, Count> zero;
    zero.fill(0.0);
    std::vector<std::array<double, Count>> partials(thread_count, zero);
#ifdef _OPENMP
#pragma omp parallel num_threads(thread_count)
#endif
    {
#ifdef _OPENMP
        const int thread = omp_get_thread_num();
#else
        const int thread = 0;
#endif
        std::array<double, Count> partial = zero;
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
        for (std::size_t i = 0; i < size; ++i)
        {
            accumulate(i, partial);
        }
        partials[thread] = partial;
    }
    std::array<double, Count> sums = zero;
    for (int thread = 0; thread < thread_count; ++thread)
    {
        for (std::size_t k = 0; k < Count; ++k)
        {
            sums[k] += partials[thread][k];
        }
    }
    return sums;
}
} // namespace ModuleSccs

#endif
