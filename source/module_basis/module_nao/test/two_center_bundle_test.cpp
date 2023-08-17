#include "module_basis/module_nao/two_center_bundle.h"

#include "gtest/gtest.h"
#include "module_base/ylm.h"
#include "module_base/global_variable.h"

#ifdef __MPI
#include <mpi.h>
#endif

/***********************************************************
 *      Unit test of class "TwoCenterBundle"
 ***********************************************************/
/*!
 *  Tested functions:
 *
 *  - build
 *      - builds the bundle of TwoCenterIntegrator objects
 *                                                                      */
class TwoCenterBundleTest : public ::testing::Test
{
  protected:
    void SetUp();
    void TearDown();

    TwoCenterBundle bundle;
};

void TwoCenterBundleTest::SetUp()
{
#ifdef __MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &GlobalV::MY_RANK);
#endif


    std::string dir = "../../../../../tests/PP_ORB/";

    int nfile_orb = 3;
    std::string* file_orb = new std::string[nfile_orb];
    file_orb[0] = dir + "C_gga_8au_100Ry_2s2p1d.orb";
    file_orb[1] = dir + "Fe_gga_9au_100Ry_4s2p2d1f.orb";
    file_orb[2] = dir + "O_gga_10au_100Ry_2s2p1d.orb";

    int nfile_pp = 3;
    std::string* file_pp = new std::string[nfile_pp];
    file_pp[0] = dir + "C_ONCV_PBE-1.0.upf";
    file_pp[1] = dir + "Fe_ONCV_PBE-1.0.upf";
    file_pp[2] = dir + "O_ONCV_PBE-1.0.upf";

    int nfile_desc = 0;

    bundle.build(nfile_orb, file_orb, nfile_pp, file_pp, nfile_desc, nullptr);

    delete[] file_orb;
    delete[] file_pp;
}

void TwoCenterBundleTest::TearDown()
{
}

TEST_F(TwoCenterBundleTest, Build)
{
    // transfer ownership to ovl
    TwoCenterIntegrator* ovl = bundle.overlap_orb.release();
    TwoCenterIntegrator* psibeta = bundle.overlap_orb_beta.release();
    TwoCenterIntegrator* kin = bundle.kinetic_orb.release();

    ModuleBase::Vector3<double> vR0 = {0.0, 0.0, 0.0};
    double out;
    ovl->calculate(0, 0, 0, 0, 0, 0, 0, 0, vR0, false, &out);
    std::cout << "out = " << out << std::endl;

    psibeta->calculate(0, 0, 0, 0, 0, 0, 0, 0, vR0, false, &out);
    std::cout << "out = " << out << std::endl;

    kin->calculate(0, 0, 0, 0, 0, 0, 0, 0, vR0, false, &out);
    std::cout << "out = " << out << std::endl;
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
