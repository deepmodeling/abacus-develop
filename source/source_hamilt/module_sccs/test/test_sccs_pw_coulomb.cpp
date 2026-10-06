#include "sccs_test.h"
#include "../sccs_pw_coulomb.h"

#include "source_base/parallel_global.h"

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
    coulomb.apply_potential(charge, potential);
    const double kernel = ModuleBase::FOUR_PI / (tpiba * tpiba);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const double expected = kernel * mode[ir];
        EXPECT_NEAR(potential[ir], expected, 1e-11);
    }
    charge.assign(basis.nrxx, 1.0);
    coulomb.apply_potential(charge, potential);
    for (double value : potential)
    {
        EXPECT_NEAR(value, 0.0, 1e-12);
    }
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
