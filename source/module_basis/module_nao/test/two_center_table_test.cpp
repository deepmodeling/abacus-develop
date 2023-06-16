#include "module_basis/module_nao/two_center_table.h"

#include "gtest/gtest.h"
#include "module_base/spherical_bessel_transformer.h"

#include <chrono>
using iclock = std::chrono::high_resolution_clock;

#ifdef __MPI
#include <mpi.h>
#endif

/***********************************************************
 *      Unit test of class "TwoCenterTable"
 ***********************************************************/
/*!
 *  Tested functions:
 *
 *  - build
 *      - builds a two-center integral radial table from two RadialCollection objects
 *
 *  - 
 *
 *                                                                      */
class TwoCenterTableTest : public ::testing::Test
{
  protected:
    void SetUp();
    void TearDown();

    TwoCenterTable table;

    RadialCollection orb;
    int nfile = 0;                                              // number of orbital files
    std::string* file = nullptr;                                //!< orbital files to read from
    std::string log_file = "./test_files/two_center_table.log"; //!< file for logging
};

void TwoCenterTableTest::SetUp()
{
    std::string dir = "../../../../../tests/PP_ORB/";
    nfile = 8;
    file = new std::string[nfile];
    file[0] = dir + "C_gga_8au_100Ry_2s2p1d.orb";
    file[1] = dir + "H_gga_8au_60Ry_2s1p.orb";
    file[2] = dir + "O_gga_10au_100Ry_2s2p1d.orb";
    file[3] = dir + "Fe_gga_9au_100Ry_4s2p2d1f.orb";
    file[4] = dir + "Cs_gga_10au_100Ry_4s2p1d.orb";
    file[5] = dir + "Pb_gga_7au_100Ry_2s2p2d1f.orb";
    file[6] = dir + "F_gga_7au_100Ry_2s2p1d.orb";
    file[7] = dir + "I_gga_7au_100Ry_2s2p2d1f.orb";
}

void TwoCenterTableTest::TearDown()
{
    delete[] file;
}

TEST_F(TwoCenterTableTest, BuildOverlapAndKinetic)
{
    orb.build(nfile, file, 'o');
    ModuleBase::SphericalBesselTransformer sbt;
    double rmax = orb.rcut_max() * 2.0;
    double dr = 0.01;
    int nr = static_cast<int>(rmax / dr) + 1;
    
    iclock::time_point start = iclock::now();

    orb.set_transformer(&sbt, 0);
    orb.set_uniform_grid(true, nr, rmax, 'i', true);

    table.build(orb, orb, 'S', false);
    table.build(orb, orb, 'S', true);
    table.build(orb, orb, 'T', false);
    table.build(orb, orb, 'T', true);

    std::chrono::duration<double> dur = iclock::now() - start;
    std::cout << "time elapsed = " << dur.count() << " seconds" << std::endl;
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
