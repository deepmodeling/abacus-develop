#include "../sccs_thread_sum.h"

#include <array>
#include <cstddef>
#include <set>
#include <vector>

#include <gtest/gtest.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace
{
// Sets the OpenMP team size for one test and turns off dynamic teams, so that
// the runtime cannot shrink the team. Both settings are restored afterwards, so
// the other tests of the binary keep OMP_NUM_THREADS and OMP_DYNAMIC.
class ThreadCount
{
public:
    explicit ThreadCount(int count)
    {
#ifdef _OPENMP
        previous_ = omp_get_max_threads();
        previous_dynamic_ = omp_get_dynamic();
        omp_set_dynamic(0);
        omp_set_num_threads(count);
#else
        (void)count;
#endif
    }

    ~ThreadCount()
    {
#ifdef _OPENMP
        omp_set_num_threads(previous_);
        omp_set_dynamic(previous_dynamic_);
#endif
    }

private:
    int previous_ = 1;
    int previous_dynamic_ = 0;
};
} // namespace

// Three threads over 1001 points: the static shares are uneven. Integer terms
// sum exactly, 1/(i + 1) matches the serial loop to rounding, and a second run
// with the same team gives the same bits.
TEST(SccsThreadSum, UnevenSharesMatchTheSerialSumAndRepeat)
{
    const ThreadCount team(3);
    const std::size_t size = 1001;
    std::vector<int> owner(size, 0);
    const auto add_terms = [&owner](std::size_t i, std::array<double, 2>& sum) {
#ifdef _OPENMP
        owner[i] = omp_get_thread_num();
#endif
        const double index = static_cast<double>(i);
        sum[0] += index;
        sum[1] += 1.0 / (index + 1.0);
    };
    const std::array<double, 2> first = ModuleSccs::thread_sums<2>(size, add_terms);
    const std::array<double, 2> second = ModuleSccs::thread_sums<2>(size, add_terms);
    double serial = 0.0;
    for (std::size_t i = 0; i < size; ++i)
    {
        const double index = static_cast<double>(i);
        serial += 1.0 / (index + 1.0);
    }
    const double point_count = static_cast<double>(size);
    const double integer_sum = 0.5 * point_count * (point_count - 1.0);
    EXPECT_EQ(first[0], integer_sum);
    EXPECT_NEAR(first[1], serial, 1e-13);
    EXPECT_EQ(first[0], second[0]);
    EXPECT_EQ(first[1], second[1]);
#ifdef _OPENMP
    // The partial sums of more than one thread were combined.
    const std::set<int> threads(owner.begin(), owner.end());
    EXPECT_GT(threads.size(), 1u);
#endif
}

// Four threads over two points leave two threads without points, and an empty
// range gives zero.
TEST(SccsThreadSum, ThreadsWithoutPointsAddZero)
{
    const ThreadCount team(4);
    const auto add_point = [](std::size_t i, std::array<double, 1>& sum) {
        const double index = static_cast<double>(i);
        sum[0] += 0.25 + index;
    };
    const std::array<double, 1> pair = ModuleSccs::thread_sums<1>(2, add_point);
    const std::array<double, 1> empty = ModuleSccs::thread_sums<1>(0, add_point);
    EXPECT_EQ(pair[0], 1.5);
    EXPECT_EQ(empty[0], 0.0);
}
