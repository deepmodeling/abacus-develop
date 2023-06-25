#include "module_basis/module_nao/beta_radials.h"

#include "gtest/gtest.h"

#ifdef __MPI
#include <mpi.h>
#endif

#include "module_base/constants.h"
#include "module_base/global_variable.h"

/***********************************************************
 *      Unit test of class "BetaRadials"
 ***********************************************************/
/*!
 *  Tested functions:
 *
 *  - build
 *      - parse a pseudopotential file and initialize individual NumericalRadial objects
 *
 *  - all "getters"
 *      - Get access to private members.
 *
 *  - all "batch setters"
 *      - Set a property for all NumericalRadial objects at once
 *                                                                      */
class BetaRadialsTest : public ::testing::Test
{
  protected:
    void SetUp();
    void TearDown() {}

    BetaRadials beta; //!< object under test

    std::string dir = "../../../../../tests/PP_ORB/"; //!< directory with test files
    std::string log_file = "./test_files/atomic_orbital.log"; //!< file for logging
};

void BetaRadialsTest::SetUp() {
#ifdef __MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &GlobalV::MY_RANK);
#endif
}

TEST_F(BetaRadialsTest, ReadAndGet100)
{
    /*
     * Read beta projectors from a UPF file of old format
     *                                                                  */
    std::string file100 = "Zn.LDA.UPF"; //!< a UPF 1.0.0 file to read from
    beta.build(dir+file100, 314, nullptr, GlobalV::MY_RANK);

    EXPECT_EQ(beta.itype(), 314);
    EXPECT_EQ(beta.symbol(), "Zn");
    EXPECT_EQ(beta.lmax(), 2);
    EXPECT_EQ(beta.nzeta(0), 0);
    EXPECT_EQ(beta.nzeta(1), 1);
    EXPECT_EQ(beta.nzeta(2), 1);
    EXPECT_EQ(beta.nzeta_max(), 1);
    EXPECT_EQ(beta.nchi(), 2);

    EXPECT_DOUBLE_EQ(beta.rcut_max(),  3.36331676102);

    EXPECT_DOUBLE_EQ(beta.chi(1, 0).ptr_rgrid()[0]  , 3.21829939971e-05);
    EXPECT_DOUBLE_EQ(beta.chi(1, 0).ptr_rgrid()[4]  , 3.39007851980e-05);
    EXPECT_DOUBLE_EQ(beta.chi(1, 0).ptr_rgrid()[888], 3.31987661585e+00);
    EXPECT_DOUBLE_EQ(beta.chi(1, 0).ptr_rgrid()[889], 3.36331676102e+00);

    EXPECT_DOUBLE_EQ(beta.chi(2, 0).ptr_rgrid()[0]  , 3.21829939971e-05);
    EXPECT_DOUBLE_EQ(beta.chi(2, 0).ptr_rgrid()[4]  , 3.39007851980e-05);
    EXPECT_DOUBLE_EQ(beta.chi(2, 0).ptr_rgrid()[887], 3.27699753773e+00);
    EXPECT_DOUBLE_EQ(beta.chi(2, 0).ptr_rgrid()[888], 3.31987661585e+00);

    EXPECT_DOUBLE_EQ(beta.chi(1, 0).ptr_rvalue()[0]  , -9.73759791529e-09);
    EXPECT_DOUBLE_EQ(beta.chi(1, 0).ptr_rvalue()[4]  , -1.08048430719e-08);
    EXPECT_DOUBLE_EQ(beta.chi(1, 0).ptr_rvalue()[600], -5.66090228695e-02);
    EXPECT_DOUBLE_EQ(beta.chi(1, 0).ptr_rvalue()[888], -1.89079314582e-11);
    EXPECT_DOUBLE_EQ(beta.chi(1, 0).ptr_rvalue()[889],  1.99483148995e-12);

    EXPECT_DOUBLE_EQ(beta.chi(2, 0).ptr_rvalue()[0]  , -2.84316719619e-11);
    EXPECT_DOUBLE_EQ(beta.chi(2, 0).ptr_rvalue()[4]  , -3.32316831459e-11);
    EXPECT_DOUBLE_EQ(beta.chi(2, 0).ptr_rvalue()[600], -3.87224488590e-01);
    EXPECT_DOUBLE_EQ(beta.chi(2, 0).ptr_rvalue()[887],  3.31894955664e-11);
    EXPECT_DOUBLE_EQ(beta.chi(2, 0).ptr_rvalue()[888], -3.22751710952e-12);
}

TEST_F(BetaRadialsTest, ReadAndGetVoid) {
    /*
     * This test reads a UPF file with no beta projectors.
     * The test is to check that the code does not crash.
     *                                                                      */
    std::string file0 = "H.pz-vbc.UPF"; //!< a UPF file with no beta projectors
    beta.build(dir+file0, 5, nullptr, GlobalV::MY_RANK);

    EXPECT_EQ(beta.nchi(), 0);
    EXPECT_EQ(beta.itype(), 5);
    EXPECT_EQ(beta.lmax(), -1);
    EXPECT_EQ(beta.symbol(), "H");
}

TEST_F(BetaRadialsTest, ReadAndGet201)
{
    /*
     * This test read beta projectors from a UPF file of 2.0.1 format
     *                                                                      */
    std::string file201 = "Pb_ONCV_PBE-1.0.upf"; //!< a UPF 2.0.1 file to read from
    beta.build(dir+file201, 999, nullptr, GlobalV::MY_RANK);

    EXPECT_EQ(beta.itype(), 999);
    EXPECT_EQ(beta.symbol(), "Pb");
    EXPECT_EQ(beta.lmax(), 3);
    EXPECT_EQ(beta.nzeta(0), 2);
    EXPECT_EQ(beta.nzeta(1), 2);
    EXPECT_EQ(beta.nzeta(2), 2);
    EXPECT_EQ(beta.nzeta(3), 2);
    EXPECT_EQ(beta.nzeta_max(), 2);
    EXPECT_EQ(beta.nchi(), 8);

    // NOTE: neither "cutoff_radius_index" nor "cutoff_radius" is reliable!
    // the code reads all the values first and then reverse scan to determine the grid size
    EXPECT_DOUBLE_EQ(beta.rcut_max(), 3.68);

    EXPECT_EQ(beta.chi(0,0).rcut(), 3.64);
    EXPECT_EQ(beta.chi(0,0).nr(), 365);
    EXPECT_EQ(beta.chi(0,0).izeta(), 0);
    EXPECT_DOUBLE_EQ(beta.chi(0, 0).ptr_rgrid()[0]  , 0.0000);
    EXPECT_DOUBLE_EQ(beta.chi(0, 0).ptr_rgrid()[8]  , 0.0800);
    EXPECT_DOUBLE_EQ(beta.chi(0, 0).ptr_rgrid()[364], 3.6400);
    EXPECT_DOUBLE_EQ(beta.chi(0, 0).ptr_rvalue()[0]  , 0.0000000000e+00);
    EXPECT_DOUBLE_EQ(beta.chi(0, 0).ptr_rvalue()[4]  , 5.9689893417e-02);
    EXPECT_DOUBLE_EQ(beta.chi(0, 0).ptr_rvalue()[364], -4.5888625103e-07);

    EXPECT_EQ(beta.chi(3,1).rcut(), 3.68);
    EXPECT_EQ(beta.chi(3,1).nr(), 369);
    EXPECT_EQ(beta.chi(3,1).izeta(), 1);
    EXPECT_DOUBLE_EQ(beta.chi(3, 1).ptr_rgrid()[0]  , 0.0000);
    EXPECT_DOUBLE_EQ(beta.chi(3, 1).ptr_rgrid()[8]  , 0.0800);
    EXPECT_DOUBLE_EQ(beta.chi(3, 1).ptr_rgrid()[368], 3.6800);
    EXPECT_DOUBLE_EQ(beta.chi(3, 1).ptr_rvalue()[0]  , 0.0000000000e+00);
    EXPECT_DOUBLE_EQ(beta.chi(3, 1).ptr_rvalue()[4]  , 1.7908487484e-06);
    EXPECT_DOUBLE_EQ(beta.chi(3, 1).ptr_rvalue()[368], -7.0309158570e-06);
}

TEST_F(BetaRadialsTest, BatchSet)
{
    std::string file201 = "Pb_ONCV_PBE-1.0.upf"; //!< a UPF 2.0.1 file to read from
    beta.build(dir+file201, 999, nullptr, GlobalV::MY_RANK);

    ModuleBase::SphericalBesselTransformer sbt;
    beta.set_transformer(&sbt);
    for (int l = 0; l != beta.lmax(); ++l) {
        for (int izeta = 0; izeta != beta.nzeta(l); ++izeta) {
            EXPECT_EQ(beta.chi(l, izeta).ptr_sbt(), &sbt);
        }
    }

    beta.set_uniform_grid(true, 2001, 20.0);
    for (int l = 0; l != beta.lmax(); ++l) {
        for (int izeta = 0; izeta != beta.nzeta(l); ++izeta) {
            EXPECT_EQ(beta.chi(l, izeta).nr(), 2001);
            EXPECT_EQ(beta.chi(l, izeta).rcut(), 20.0);
            EXPECT_EQ(beta.chi(l, izeta).ptr_rgrid()[1500], 15.0);
            EXPECT_EQ(beta.chi(l, izeta).ptr_rvalue()[1500], 0.0);
        }
    }

    double grid[5] = {0.0, 1.1, 2.2, 3.3, 4.4};
    beta.set_grid(true, 5, grid);
    for (int l = 0; l != beta.lmax(); ++l) {
        for (int izeta = 0; izeta != beta.nzeta(l); ++izeta) {
            EXPECT_EQ(beta.chi(l, izeta).ptr_rgrid()[0], 0.0);
            EXPECT_EQ(beta.chi(l, izeta).ptr_rgrid()[1], 1.1);
            EXPECT_EQ(beta.chi(l, izeta).ptr_rgrid()[2], 2.2);
            EXPECT_EQ(beta.chi(l, izeta).ptr_rgrid()[3], 3.3);
            EXPECT_EQ(beta.chi(l, izeta).ptr_rgrid()[4], 4.4);
        }
    }
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
