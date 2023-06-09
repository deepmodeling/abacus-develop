#include "module_basis/module_nao/atomic_radials.h"

#include "gtest/gtest.h"

#ifdef __MPI
#include <mpi.h>
#endif

#include "module_base/constants.h"

/***********************************************************
 *      Unit test of class "AtomicRadials"
 ***********************************************************/
/*!
 *  Tested functions:
 *
 *  - read
 *      - parse an orbital file and initialize the NumericalRadial objects
 *
 *  - all "getters"
 *      - Get access to private members.
 *
 *  - all "batch setters"
 *      - Set a property for all NumericalRadial objects at once
 *                                                                      */
class AtomicRadialsTest : public ::testing::Test
{
  protected:
    void SetUp(){};
    void TearDown(){};

    AtomicRadials Ti_radials;                                         //!< object under test
    std::string file = "../../../../../tests/PP_ORB/Ti_gga_10au_100Ry_4s2p2d1f.orb"; //!< orbital file to read from
    std::string log_file = "./test_files/atomic_orbital.log";         //!< file for logging

    double tol = 1e-12; //!< numerical tolerance for grid & values
};

TEST_F(AtomicRadialsTest, ReadAndGet)
{

    Ti_radials.build(file);

    EXPECT_EQ(Ti_radials.lmax(), 3);
    EXPECT_EQ(Ti_radials.nzeta(0), 4);
    EXPECT_EQ(Ti_radials.nzeta(1), 2);
    EXPECT_EQ(Ti_radials.nzeta(2), 2);
    EXPECT_EQ(Ti_radials.nzeta(3), 1);
    EXPECT_EQ(Ti_radials.nzeta_max(), 4);
    EXPECT_EQ(Ti_radials.nchi(), 9);
    EXPECT_DOUBLE_EQ(Ti_radials.rcut_max(), 10.0);
    EXPECT_DOUBLE_EQ(Ti_radials.orb_ecut(), 100.0);

    EXPECT_NEAR(Ti_radials.chi(0, 0).ptr_rvalue()[0], -1.581711853170e-01, tol);
    EXPECT_NEAR(Ti_radials.chi(0, 0).ptr_rvalue()[4], -1.583907030513e-01, tol);
    EXPECT_NEAR(Ti_radials.chi(0, 0).ptr_rvalue()[996], -4.183526380009e-05, tol);
    EXPECT_NEAR(Ti_radials.chi(0, 0).ptr_rvalue()[1000], 0, tol);

    EXPECT_NEAR(Ti_radials.chi(0, 3).ptr_rvalue()[0], -1.166292682541e+00, tol);
    EXPECT_NEAR(Ti_radials.chi(0, 3).ptr_rvalue()[4], -1.164223359672e+00, tol);
    EXPECT_NEAR(Ti_radials.chi(0, 3).ptr_rvalue()[996], -3.183325576529e-04, tol);
    EXPECT_NEAR(Ti_radials.chi(0, 3).ptr_rvalue()[1000], 0, tol);

    EXPECT_NEAR(Ti_radials.chi(3, 0).ptr_rvalue()[0], 0, tol);
    EXPECT_NEAR(Ti_radials.chi(3, 0).ptr_rvalue()[4], 3.744878535962e-05, tol);
    EXPECT_NEAR(Ti_radials.chi(3, 0).ptr_rvalue()[996], 7.495357740660e-05, tol);
    EXPECT_NEAR(Ti_radials.chi(3, 0).ptr_rvalue()[1000], 0, tol);
}

TEST_F(AtomicRadialsTest, BatchSet)
{

    int itype = 5;
    Ti_radials.build(file, itype);

    EXPECT_EQ(Ti_radials.itype(), 5);
    EXPECT_EQ(Ti_radials.chi(0, 0).itype(), 5);
    EXPECT_EQ(Ti_radials.chi(0, 3).itype(), 5);
    EXPECT_EQ(Ti_radials.chi(3, 0).itype(), 5);

    ModuleBase::SphericalBesselTransformer sbt;
    Ti_radials.set_transformer(&sbt);
    EXPECT_EQ(Ti_radials.chi(0, 0).ptr_sbt(), &sbt);
    EXPECT_EQ(Ti_radials.chi(0, 3).ptr_sbt(), &sbt);
    EXPECT_EQ(Ti_radials.chi(3, 0).ptr_sbt(), &sbt);

    Ti_radials.set_uniform_grid(true, 2001, 20.0);
    EXPECT_EQ(Ti_radials.chi(0, 0).nr(), 2001);
    EXPECT_EQ(Ti_radials.chi(0, 3).nr(), 2001);
    EXPECT_EQ(Ti_radials.chi(3, 0).nr(), 2001);
    EXPECT_NEAR(Ti_radials.chi(0, 0).ptr_rgrid()[2000], 20, tol);
    EXPECT_NEAR(Ti_radials.chi(0, 3).ptr_rgrid()[1500], 15, tol);
    EXPECT_NEAR(Ti_radials.chi(3, 0).ptr_rgrid()[1200], 12, tol);
    EXPECT_NEAR(Ti_radials.chi(0, 0).ptr_rvalue()[2000], 0, tol);
    EXPECT_NEAR(Ti_radials.chi(0, 3).ptr_rvalue()[1500], 0, tol);
    EXPECT_NEAR(Ti_radials.chi(3, 0).ptr_rvalue()[1200], 0, tol);

    double grid[5] = {0.0, 1.1, 2.2, 3.3, 4.4};
    Ti_radials.set_grid(true, 5, grid);
    EXPECT_EQ(Ti_radials.chi(0, 0).ptr_rgrid()[1], 1.1);
    EXPECT_EQ(Ti_radials.chi(0, 3).ptr_rgrid()[2], 2.2);
    EXPECT_EQ(Ti_radials.chi(3, 0).ptr_rgrid()[3], 3.3);
}

int main(int argc, char** argv)
{

#ifdef __MPI
    MPI_Init(&argc, &argv);
#endif

    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

#ifdef __MPI
    MPI_Finalize();
#endif

    return result;
}
