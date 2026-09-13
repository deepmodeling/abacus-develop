#ifndef STERNHEIMER_WEAK_Q_DISPATCH_H
#define STERNHEIMER_WEAK_Q_DISPATCH_H

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <limits>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef __MKL
#include <mkl_service.h>
#endif

namespace ModuleRI
{

struct SternheimerWeakQLayout
{
    int workers;
    int inner_threads;
    int allocated_threads;

    void validate() const
    {
        if (workers <= 0 || inner_threads <= 0 || allocated_threads <= 0
            || workers > allocated_threads / inner_threads)
            throw std::invalid_argument("Weak q dispatch requires positive workers/inner/allocation and workers*inner <= allocation.");
    }
};

namespace sternheimer_weak_q_detail
{
#ifdef _OPENMP
class DynamicThreadsGuard
{
  public:
    DynamicThreadsGuard() : previous_dynamic_(omp_get_dynamic()) { omp_set_dynamic(0); }
    ~DynamicThreadsGuard() { omp_set_dynamic(previous_dynamic_); }
    DynamicThreadsGuard(const DynamicThreadsGuard&) = delete;
    DynamicThreadsGuard& operator=(const DynamicThreadsGuard&) = delete;

  private:
    int previous_dynamic_;
};
#endif

// Settings belong to the calling OpenMP task / MKL thread, never the process.
class InnerThreadsGuard
{
  public:
    explicit InnerThreadsGuard(int threads)
    {
#ifdef _OPENMP
        previous_threads_ = omp_get_max_threads();
        previous_levels_ = omp_get_max_active_levels();
        omp_set_num_threads(threads);
        const int inner_level = omp_get_active_level() + 1;
        if (threads > 1 && previous_levels_ < inner_level)
            omp_set_max_active_levels(inner_level);
#endif
#ifdef __MKL
        previous_mkl_ = mkl_set_num_threads_local(threads);
#endif
        (void)threads;
    }

    ~InnerThreadsGuard()
    {
#ifdef __MKL
        mkl_set_num_threads_local(previous_mkl_);
#endif
#ifdef _OPENMP
        omp_set_max_active_levels(previous_levels_);
        omp_set_num_threads(previous_threads_);
#endif
    }

    InnerThreadsGuard(const InnerThreadsGuard&) = delete;
    InnerThreadsGuard& operator=(const InnerThreadsGuard&) = delete;

  private:
#ifdef _OPENMP
    int previous_threads_;
    int previous_levels_;
#endif
#ifdef __MKL
    int previous_mkl_;
#endif
};
} // namespace sternheimer_weak_q_detail

// factory(index) returns a worker value (including move-only handles); callback
// receives (zero-based task, worker&, index). Both callables must support concurrent
// invocation; mutable solve state belongs in the worker, not shared captures.
// Every task is attempted once, even after callback failure. Factory failure starts
// no tasks. Like run_sternheimer_channel_tasks, exceptions leave only after the team.
// Counts have layout.workers entries. A runtime-limited team is rejected before
// any factory or callback runs; an empty task range creates no team or workers.
// Dynamic team adjustment is disabled for the dispatch and restored on exit.
template <typename Factory, typename Callback>
std::vector<std::size_t> run_sternheimer_weak_q_tasks(const std::int64_t task_count,
                                                    const SternheimerWeakQLayout& layout,
                                                    Factory&& factory,
                                                    Callback&& callback)
{
    layout.validate();
    if (task_count < 0
        || static_cast<std::uint64_t>(task_count) > std::numeric_limits<std::size_t>::max())
        throw std::invalid_argument("Weak q dispatch task count must be non-negative and fit size_t.");
#ifndef _OPENMP
    if (layout.workers != 1)
        throw std::invalid_argument("Weak q dispatch with multiple workers requires OpenMP.");
#endif
    std::vector<std::size_t> counts(static_cast<std::size_t>(layout.workers), 0);
    if (task_count == 0) return counts;

    using FactoryResult = decltype(factory(0));
    static_assert(!std::is_reference<FactoryResult>::value, "Weak q factory must return a worker value, not a reference.");
    using Worker = typename std::decay<FactoryResult>::type;
    std::vector<std::exception_ptr> errors(static_cast<std::size_t>(layout.workers));
    std::atomic<bool> factory_failed{false};

#ifdef _OPENMP
    sternheimer_weak_q_detail::DynamicThreadsGuard dynamic_threads;
#pragma omp parallel num_threads(layout.workers)
#endif
    {
        const int index =
#ifdef _OPENMP
            omp_get_thread_num();
#else
            0;
#endif
        sternheimer_weak_q_detail::InnerThreadsGuard inner_threads(layout.inner_threads);
        std::unique_ptr<Worker> worker;
        try
        {
#ifdef _OPENMP
            if (omp_get_num_threads() != layout.workers)
                throw std::runtime_error("Weak q dispatch actual OpenMP team does not match requested workers.");
#endif
            worker.reset(new Worker(factory(index)));
        }
        catch (...)
        {
            errors[index] = std::current_exception();
            factory_failed.store(true, std::memory_order_relaxed);
        }
        // All threads must agree to enter the worksharing loop, including when a
        // factory throws. Keep this flag separate from callback exception slots.
#ifdef _OPENMP
#pragma omp barrier
#endif
        if (!factory_failed.load(std::memory_order_relaxed))
        {
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
            for (std::int64_t task = 0; task < task_count; ++task)
            {
                try
                {
                    callback(task, *worker, index);
                }
                catch (...)
                {
                    if (!errors[index]) errors[index] = std::current_exception();
                }
                ++counts[index];
            }
        }
    }
    for (const auto& error : errors)
        if (error) std::rethrow_exception(error);
    return counts;
}

} // namespace ModuleRI
#endif
