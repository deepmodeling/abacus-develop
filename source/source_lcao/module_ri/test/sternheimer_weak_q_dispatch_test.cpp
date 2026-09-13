#include <gtest/gtest.h>

#include "source_lcao/module_ri/sternheimer_weak_q_dispatch.h"

#include <array>
#include <atomic>
#include <limits>
#include <memory>
#include <numeric>
#include <string>

using ModuleRI::SternheimerWeakQLayout;
using ModuleRI::run_sternheimer_weak_q_tasks;

namespace
{
int thread_index()
{
#ifdef _OPENMP
    return omp_get_thread_num();
#else
    return 0;
#endif
}

struct Worker
{
    explicit Worker(int index) : index(index) {}
    Worker(const Worker&) = delete;
    Worker& operator=(const Worker&) = delete;
    int index;
    std::size_t calls = 0;
};

void check_coverage(const SternheimerWeakQLayout& layout, int task_count)
{
    std::vector<std::atomic<int>> visits(task_count);
    for (auto& count : visits) count.store(0);
    std::vector<int> factories(layout.workers, 0);
    std::vector<int> team_sizes(layout.workers, 0);
    std::vector<int> owners(task_count, -1);
    std::vector<std::size_t> observed(layout.workers, 0);
    const auto counts = run_sternheimer_weak_q_tasks(task_count, layout,
        [&](int index) {
            EXPECT_EQ(thread_index(), index);
            ++factories[index];
#ifdef _OPENMP
            team_sizes[index] = omp_get_num_threads();
            EXPECT_EQ(team_sizes[index], layout.workers);
            EXPECT_FALSE(omp_get_dynamic());
            EXPECT_EQ(omp_get_max_threads(), layout.inner_threads);
#else
            team_sizes[index] = 1;
#endif
#ifdef __MKL
            EXPECT_EQ(mkl_get_max_threads(), layout.inner_threads);
            EXPECT_EQ(mkl_set_num_threads_local(layout.inner_threads), layout.inner_threads);
#endif
            return std::unique_ptr<Worker>(new Worker(index));
        },
        [&](std::int64_t task, std::unique_ptr<Worker>& worker, int index) {
            EXPECT_EQ(thread_index(), index);
            EXPECT_EQ(worker->index, index);
#ifdef _OPENMP
            EXPECT_FALSE(omp_get_dynamic());
#endif
#ifdef __MKL
            EXPECT_EQ(mkl_get_max_threads(), layout.inner_threads);
            EXPECT_EQ(mkl_set_num_threads_local(layout.inner_threads), layout.inner_threads);
#endif
            ++worker->calls;
            observed[index] = worker->calls;
            visits[task].fetch_add(1);
            owners[task] = index;
        });
    ASSERT_EQ(counts.size(), static_cast<std::size_t>(layout.workers));
    ASSERT_EQ(team_sizes[0], layout.workers);
    EXPECT_EQ(std::accumulate(factories.begin(), factories.end(), 0), team_sizes[0]);
    EXPECT_EQ(std::accumulate(counts.begin(), counts.end(), std::size_t(0)),
              static_cast<std::size_t>(task_count));
    for (int i = 0; i < task_count; ++i)
    {
        EXPECT_EQ(visits[i].load(), 1);
        EXPECT_GE(owners[i], 0);
        EXPECT_LT(owners[i], layout.workers);
    }
    for (int i = 0; i < layout.workers; ++i)
    {
        EXPECT_EQ(factories[i], 1);
        EXPECT_EQ(counts[i], observed[i]);
        if (counts[i] != 0) EXPECT_EQ(factories[i], 1);
    }
}
} // namespace

TEST(SternheimerWeakQDispatch, ValidatesAll48CpuCandidatesWithoutOverflow)
{
    for (int workers : {1, 4, 8, 16, 24, 48})
        EXPECT_NO_THROW((SternheimerWeakQLayout{workers, 48 / workers, 48}.validate()));
    EXPECT_NO_THROW((SternheimerWeakQLayout{3, 5, 48}.validate()));
    for (const auto& layout : std::array<SternheimerWeakQLayout, 9>{{
             {0, 1, 48}, {-1, 1, 48}, {1, 0, 48}, {1, -1, 48},
             {1, 1, 0}, {1, 1, -1}, {49, 1, 48}, {8, 7, 48},
             {std::numeric_limits<int>::max(), 2, std::numeric_limits<int>::max()}}})
        EXPECT_THROW(layout.validate(), std::invalid_argument);
}

TEST(SternheimerWeakQDispatch, RejectsInvalidInputBeforeCreatingWorkers)
{
    int created = 0;
    auto factory = [&](int) { ++created; return 0; };
    auto callback = [](std::int64_t, int&, int) {};
    EXPECT_THROW(run_sternheimer_weak_q_tasks(-1, {1, 1, 1}, factory, callback), std::invalid_argument);
    EXPECT_THROW(run_sternheimer_weak_q_tasks(0, {0, 1, 1}, factory, callback), std::invalid_argument);
    EXPECT_THROW(run_sternheimer_weak_q_tasks(10, {2, 2, 3}, factory, callback), std::invalid_argument);
    EXPECT_EQ(created, 0);
}

TEST(SternheimerWeakQDispatch, EmptyTasksDoNotCreateWorkers)
{
    int created = 0;
    const auto counts = run_sternheimer_weak_q_tasks(0, {1, 1, 1},
        [&](int) { ++created; return 0; }, [](std::int64_t, int&, int) { FAIL(); });
    EXPECT_EQ(created, 0);
    EXPECT_EQ(counts, std::vector<std::size_t>{0});
}

TEST(SternheimerWeakQDispatch, ExactCoverageForMoveOnlyThreadOwnedWorkers)
{
#ifdef _OPENMP
    for (int workers : {1, 4, 8, 16, 24, 48})
    {
        SCOPED_TRACE(workers);
        check_coverage({workers, 48 / workers, 48}, 197);
    }
    check_coverage({8, 1, 8}, 3);
#else
    check_coverage({1, 1, 1}, 197);
    EXPECT_THROW(run_sternheimer_weak_q_tasks(1, {2, 1, 2},
        [](int) { return 0; }, [](std::int64_t, int&, int) {}), std::invalid_argument);
#endif
}

TEST(SternheimerWeakQDispatch, CallbackFailureStillAttemptsEveryTaskOnceAndRethrowsOutsideTeam)
{
    std::array<std::atomic<int>, 37> visits;
    for (auto& count : visits) count.store(0);
    std::atomic<int> alive{0};
    const int workers =
#ifdef _OPENMP
        4;
#else
        1;
#endif
    try
    {
        run_sternheimer_weak_q_tasks(37, {workers, 1, workers},
            [&](int index) {
                ++alive;
                return std::shared_ptr<Worker>(new Worker(index), [&](Worker* worker) {
                    --alive;
                    delete worker;
                });
            },
            [&](std::int64_t task, std::shared_ptr<Worker>&, int) {
                ++visits[task];
                if (task == 3 || task == 17) throw std::runtime_error("callback failure");
            });
        FAIL() << "Callback exception was lost";
    }
    catch (const std::runtime_error& error)
    {
        EXPECT_EQ(std::string(error.what()), "callback failure");
#ifdef _OPENMP
        EXPECT_FALSE(omp_in_parallel());
#endif
    }
    EXPECT_EQ(alive.load(), 0);
    for (const auto& count : visits) EXPECT_EQ(count.load(), 1);
}

TEST(SternheimerWeakQDispatch, FactoryFailurePreventsAllCallbacksAndDoesNotDeadlock)
{
    std::atomic<int> callbacks{0};
    std::atomic<int> alive{0};
    const int workers =
#ifdef _OPENMP
        4;
#else
        1;
#endif
    EXPECT_THROW(run_sternheimer_weak_q_tasks(37, {workers, 1, workers},
        [&](int index) {
            if (index == 0) throw std::logic_error("factory failure");
            ++alive;
            return std::shared_ptr<Worker>(new Worker(index), [&](Worker* worker) {
                --alive;
                delete worker;
            });
        }, [&](std::int64_t, std::shared_ptr<Worker>&, int) { ++callbacks; }), std::logic_error);
    EXPECT_EQ(callbacks.load(), 0);
    EXPECT_EQ(alive.load(), 0);
}

TEST(SternheimerWeakQDispatch, PreservesNonStandardExceptionType)
{
    EXPECT_THROW(run_sternheimer_weak_q_tasks(2, {1, 1, 1},
        [](int) -> int { throw 17; }, [](std::int64_t, int&, int) {}), int);
    EXPECT_THROW(run_sternheimer_weak_q_tasks(2, {1, 1, 1},
        [](int) { return 0; }, [](std::int64_t, int&, int) { throw 23; }), int);
}

TEST(SternheimerWeakQDispatch, RestoresThreadSettingsOnSuccessAndBothExceptionPaths)
{
#ifdef _OPENMP
    const int previous_threads = omp_get_max_threads();
    const int previous_levels = omp_get_max_active_levels();
    const int previous_dynamic = omp_get_dynamic();
    omp_set_dynamic(1);
#endif
#ifdef __MKL
    const int previous_mkl = mkl_set_num_threads_local(3);
#endif
    for (int failure : {0, 1, 2})
    {
        try
        {
            run_sternheimer_weak_q_tasks(3, {1, 2, 2},
                [&](int) {
#ifdef _OPENMP
                    EXPECT_EQ(omp_get_max_threads(), 2);
                    EXPECT_FALSE(omp_get_dynamic());
#endif
#ifdef __MKL
                    EXPECT_EQ(mkl_get_max_threads(), 2);
                    EXPECT_EQ(mkl_set_num_threads_local(2), 2);
#endif
                    if (failure == 1) throw std::runtime_error("factory");
                    return 0;
                }, [&](std::int64_t, int&, int) {
#ifdef _OPENMP
                    EXPECT_FALSE(omp_get_dynamic());
#endif
#ifdef __MKL
                    EXPECT_EQ(mkl_get_max_threads(), 2);
                    EXPECT_EQ(mkl_set_num_threads_local(2), 2);
#endif
                    if (failure == 2) throw std::runtime_error("callback");
                });
            EXPECT_EQ(failure, 0);
        }
        catch (const std::runtime_error&) { EXPECT_NE(failure, 0); }
#ifdef _OPENMP
        EXPECT_EQ(omp_get_max_threads(), previous_threads);
        EXPECT_EQ(omp_get_max_active_levels(), previous_levels);
        EXPECT_TRUE(omp_get_dynamic());
#endif
#ifdef __MKL
        EXPECT_EQ(mkl_get_max_threads(), 3);
        EXPECT_EQ(mkl_set_num_threads_local(3), 3);
#endif
    }
#ifdef __MKL
    EXPECT_EQ(mkl_set_num_threads_local(previous_mkl), 3);
#endif
#ifdef _OPENMP
    omp_set_dynamic(previous_dynamic);
#endif
}

#ifdef _OPENMP
TEST(SternheimerWeakQDispatch, RuntimeSerializedTeamFailsBeforeCreatingWorkers)
{
    const int previous_levels = omp_get_max_active_levels();
    const int previous_dynamic = omp_get_dynamic();
    omp_set_max_active_levels(0);
    omp_set_dynamic(1);
    std::atomic<int> factories{0};
    std::atomic<int> callbacks{0};
    EXPECT_THROW(run_sternheimer_weak_q_tasks(37, {4, 1, 4},
        [&](int) { ++factories; return 0; },
        [&](std::int64_t, int&, int) { ++callbacks; }), std::runtime_error);
    EXPECT_EQ(factories.load(), 0);
    EXPECT_EQ(callbacks.load(), 0);
    EXPECT_EQ(omp_get_max_active_levels(), 0);
    EXPECT_TRUE(omp_get_dynamic());
    omp_set_dynamic(previous_dynamic);
    omp_set_max_active_levels(previous_levels);
}

TEST(SternheimerWeakQDispatch, InnerOpenMPTeamUsesRequestedThreadCount)
{
    const int previous_dynamic = omp_get_dynamic();
    omp_set_dynamic(1);
    run_sternheimer_weak_q_tasks(8, {4, 2, 8}, [](int) { return 0; },
        [](std::int64_t, int&, int) {
            EXPECT_EQ(omp_get_num_threads(), 4);
            EXPECT_FALSE(omp_get_dynamic());
            int inner_size = 0;
#pragma omp parallel shared(inner_size)
            {
#pragma omp single
                inner_size = omp_get_num_threads();
            }
            EXPECT_EQ(inner_size, 2);
        });
    EXPECT_TRUE(omp_get_dynamic());
    omp_set_dynamic(previous_dynamic);
}
#endif

#ifdef __MKL
TEST(SternheimerWeakQDispatch, RestoresMklGlobalDefaultSentinel)
{
    const int previous_mkl = mkl_set_num_threads_local(0);
    for (int failure : {0, 1, 2})
    {
        try
        {
            run_sternheimer_weak_q_tasks(1, {1, 2, 2},
                [&](int) {
                    EXPECT_EQ(mkl_get_max_threads(), 2);
                    EXPECT_EQ(mkl_set_num_threads_local(2), 2);
                    if (failure == 1) throw std::runtime_error("factory");
                    return 0;
                }, [&](std::int64_t, int&, int) {
                    EXPECT_EQ(mkl_get_max_threads(), 2);
                    EXPECT_EQ(mkl_set_num_threads_local(2), 2);
                    if (failure == 2) throw std::runtime_error("callback");
                });
            EXPECT_EQ(failure, 0);
        }
        catch (const std::runtime_error&) { EXPECT_NE(failure, 0); }
        EXPECT_EQ(mkl_set_num_threads_local(0), 0);
    }
    EXPECT_EQ(mkl_set_num_threads_local(previous_mkl), 0);
}
#endif
