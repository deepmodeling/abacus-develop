#ifdef NDEBUG
#undef NDEBUG
#endif

#include "source_base/module_container/NDArray/dist_ndarray.h"

#include "gtest/gtest.h"

#include <cstddef>
#include <memory>
#include <type_traits>

static_assert(!std::is_copy_constructible<DistNDArray<double>>::value,
              "DistNDArray must not be copy constructible");
static_assert(!std::is_copy_assignable<DistNDArray<double>>::value,
              "DistNDArray must not be copy assignable");
static_assert(std::is_nothrow_move_constructible<DistNDArray<double>>::value,
              "DistNDArray move construction must be noexcept");
static_assert(std::is_nothrow_move_assignable<DistNDArray<double>>::value,
              "DistNDArray move assignment must be noexcept");

namespace
{
std::shared_ptr<BlockCyclicDistribution> make_serial_distribution(
    const int nrglb,
    const int ncglb)
{
    std::shared_ptr<BlockCyclicDistribution> distribution
        = std::make_shared<BlockCyclicDistribution>();
    distribution->set_serial(nrglb, ncglb);
    return distribution;
}

void check_linearization(const int nrglb,
                         const int ncglb,
                         const StorageMajorOrder major_order)
{
    std::shared_ptr<BlockCyclicDistribution> distribution
        = make_serial_distribution(nrglb, ncglb);
    const int nloc
        = major_order == StorageMajorOrder::ColumnMajor ? nrglb : ncglb;
    const int ld = nloc > 1 ? nloc : 1;

    DistNDArray<double> array(distribution, major_order, ld);

    EXPECT_EQ(array.distribution().get(), distribution.get());
    EXPECT_EQ(array.major_order(), major_order);
    EXPECT_EQ(array.ld(), ld);
    EXPECT_EQ(array.nrloc(), nrglb);
    EXPECT_EQ(array.ncloc(), ncglb);

    const std::size_t ld_size = static_cast<std::size_t>(array.ld());

    for (int i = 0; i < nrglb; ++i)
    {
        for (int j = 0; j < ncglb; ++j)
        {
            const double value = static_cast<double>(i * ncglb + j);
            array.at_local(i, j) = value;
        }
    }

    for (int i = 0; i < nrglb; ++i)
    {
        for (int j = 0; j < ncglb; ++j)
        {
            std::size_t iflat = 0;
            if (major_order == StorageMajorOrder::ColumnMajor)
            {
                iflat = static_cast<std::size_t>(j) * ld_size
                    + static_cast<std::size_t>(i);
            }
            else
            {
                iflat = static_cast<std::size_t>(i) * ld_size
                    + static_cast<std::size_t>(j);
            }

            const double* const expected_address = array.data() + iflat;
            EXPECT_EQ(&array.at_local(i, j), expected_address);

            const double value = static_cast<double>(i * ncglb + j);
            EXPECT_EQ(array.at_local(i, j), value);

            const int il = distribution->irloc(i);
            const int jl = distribution->icloc(j);
            ASSERT_NE(il, -1);
            ASSERT_NE(jl, -1);
            EXPECT_EQ(&array.at_global(i, j), &array.at_local(il, jl));
            EXPECT_EQ(array.at_global(i, j), array.at_local(il, jl));

            const bool owned = distribution->in_this_processor(i, j);
            EXPECT_EQ(array.owns_global_index(i, j), owned);
        }
    }
}
} // namespace

TEST(DistNDArrayTest, ColumnMajorLinearization)
{
    const int nrglb = 8;
    const int ncglb = 6;
    check_linearization(nrglb, ncglb, StorageMajorOrder::ColumnMajor);
}

TEST(DistNDArrayTest, RowMajorLinearization)
{
    const int nrglb = 2;
    const int ncglb = 3;
    check_linearization(nrglb, ncglb, StorageMajorOrder::RowMajor);
}

TEST(DistNDArrayTest, AtGlobalAssertsOutsideGlobalBounds)
{
    const int nrglb = 4;
    const int ncglb = 3;
    std::shared_ptr<BlockCyclicDistribution> distribution
        = make_serial_distribution(nrglb, ncglb);
    const int ld = nrglb > 1 ? nrglb : 1;
    DistNDArray<double> array(distribution, StorageMajorOrder::ColumnMajor, ld);

    EXPECT_DEATH(array.at_global(nrglb, 0), "");
}

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
