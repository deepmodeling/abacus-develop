#include "sccs_test.h"
#include "../sccs_pw_coulomb.h"

#include "source_base/parallel_global.h"

#include <limits>

namespace SccsTest
{
int pool_size = 1;
int pool_rank = 0;
}

using SccsPwCoulombTest = SccsTest::PwTest;

TEST_F(SccsPwCoulombTest, FourierModeAndConstantBackground)
{
    const std::vector<double> mode = cosine_mode(0);
    std::vector<double> charge = mode;
    for (double& value : charge)
    {
        value += 0.3;
    }
    ModuleSccs::PeriodicCoulombOperator coulomb(basis, tpiba);
    std::vector<double> potential;
    ASSERT_TRUE(coulomb.apply_potential(charge, potential, error)) << error;
    const double kernel = ModuleBase::FOUR_PI / (tpiba * tpiba);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const double expected = kernel * mode[ir];
        EXPECT_NEAR(potential[ir], expected, 1e-11);
    }
    charge.assign(basis.nrxx, 1.0);
    ASSERT_TRUE(coulomb.apply_potential(charge, potential, error)) << error;
    for (double value : potential)
    {
        EXPECT_NEAR(value, 0.0, 1e-12);
    }
}

TEST_F(SccsPwCoulombTest, RankLocalInvalidInputReturnsCollectively)
{
    std::vector<double> charge(basis.nrxx, 0.0);
    if (basis.poolrank == 0)
    {
        charge[0] = std::numeric_limits<double>::quiet_NaN();
    }
    ModuleSccs::PeriodicCoulombOperator coulomb(basis, tpiba);
    std::vector<double> potential(1, 12.0);
    EXPECT_FALSE(coulomb.apply_potential(charge, potential, error));
    EXPECT_FALSE(error.empty());
    ASSERT_EQ(potential.size(), 1u);
    EXPECT_DOUBLE_EQ(potential[0], 12.0);
    charge.assign(basis.nrxx, 0.0);
    if (basis.poolrank == 0)
    {
        charge.pop_back();
    }
    EXPECT_FALSE(coulomb.apply_potential(charge, potential, error));
    EXPECT_DOUBLE_EQ(potential[0], 12.0);
}

int main(int argc, char** argv)
{
    int threads = 1;
    Parallel_Global::read_pal_param(argc, argv, SccsTest::pool_size, threads, SccsTest::pool_rank);
    int band_size = 1;
    int band_rank = 0;
    int band_group = 0;
    int pool = 0;
#ifdef __MPI
    // These tests create only pool/band groups, not diagonalization/grid groups.
    GRID_WORLD = MPI_COMM_NULL;
    DIAG_WORLD = MPI_COMM_NULL;
    Parallel_Global::init_pools(SccsTest::pool_size,
                               SccsTest::pool_rank,
                               1,
                               1,
                               band_size,
                               band_rank,
                               band_group,
                               SccsTest::pool_size,
                               SccsTest::pool_rank,
                               pool);
#endif
    testing::InitGoogleTest(&argc, argv);
    const int status = RUN_ALL_TESTS();
#ifdef __MPI
    Parallel_Global::finalize_mpi();
#endif
    return status;
}
