#include <fstream>
#include "gtest/gtest.h"
#include "module_orbital/ORB_atomic.h"
#include "module_orbital/ORB_atomic_lm.h"

#define private public
#include "module_orbital/ORB_read.h"
#undef private

#ifdef __MPI
#include <mpi.h>
#endif

/***********************************************************
 *          unit test of class "LCAO_Orbitals"
 ***********************************************************/

/** 
 * Tested functions:
 *
 * - read_orb_file
 *   read orbital file from given ifstream and pour the data into
 *   Numerical_Orbital & its member Numerical_Orbital_Lm objects.
 *   This function is the core of the class.
 *
 * - Read_PAO
 *   read
 *
 * - Read_Descriptor
 *   read
 *
 * - Read_Orbitals
 *   read
 *
 * - bcast_files
 *   broadcast
 *
 * - all "getters"
 *   get access to class members
 *
 ***********************************************************/

class LcaoOrbitalsTest : public ::testing::Test
{
protected:

    void SetUp();
    void TearDown();

    // object under unit test
    LCAO_Orbitals lcao_;

    // helper functions
    void init();
    size_t calc_nk(double const& ecutwfc, double const& dk);
};


size_t LcaoOrbitalsTest::calc_nk(double const& ecutwfc, double const& dk) {

    // current formula for calculating nk from ecutwfc & dk
    // see module_orbital/ORB_read.cpp, function "Read_Orbitals"

    size_t nk = 0;

    if(ecutwfc < 20) {
        nk = static_cast<int>( 2 * sqrt(ecutwfc) / dk )  + 4;
    } else {
        nk = static_cast<int>( sqrt(ecutwfc) / dk )  + 4;
    }

    if (nk%2 == 0) {
        ++nk;
    }

    return nk;
}


/*
void LcaoOrbitalsTest::SetUp() {

    ///////////////////////////////////////////////////
    //                  Parameters
    ///////////////////////////////////////////////////

    // mock GlobalV::ofs_running
    std::ofstream ofs_log("ORB_read_test.log");

    // test case directory
    // contains the abacus input files & numerical atomic orbital files
    std::string case_dir = "./lcao_H2O/";

    // these numbers are not read by LCAO_Orbitals; just pass in for simplicity
    int ntype = 2;
    int lmax = 2;
    double lcao_ecut = 100.0;
    double lcao_dk = 0.01;
    double lcao_dr = 0.01;
    double lcao_rmax = 20;

    bool deepks_setorb = true;
    //int out_mat_r = 0; // unused variable?
    bool force_flag = true;
    int my_rank = 0;


    // append numerical atomic orbital files to lcao.orbital_file
    lcao_.read_in_flag = true;
    std::string stru_file = case_dir + "STRU";
    std::ifstream ifs(stru_file);

    ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "NUMERICAL_ORBITAL");
    std::string orb_file;
    for(int it = 0; it < this->ntype; ++it) {
        ifs >> orb_file;
        orb_file = case_dir + orb_file;
        lcao_.orbital_file.push_back(orb_file);
    }

    // set descriptor file
    ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "NUMERICAL_DESCRIPTOR");
    std::string desc_file;
    ifs >> desc_file;
    desc_file = case_dir + desc_file;
    lcao_.descriptor_file = desc_file; 

    // below we mimic ORB_control::read_orb_first
    lcao_.ecutwfc = lcao_ecut;
    lcao_.dk = lcao_dk;
    lcao_.dR = lcao_dr;
    lcao_.Rmax = lcao_rmax;

    lcao_.Read_Orbitals(ofs_log, ntype, lmax, deepks_setorb, 0, 
            force_flag, my_rank);

    ofs_log << "ecutwfc = " << lcao_.ecutwfc << std::endl
        << "kmesh = " << lcao_.get_kmesh() << std::endl
        << "dk = " << lcao_.get_dk() << std::endl
        << "dR = " << lcao_.get_dR() << std::endl
        << "Rmax = " << lcao_.get_Rmax() << std::endl
        << "lmax = " << lcao_.get_lmax() << std::endl
        << "lmax_d = " << lcao_.get_lmax_d() << std::endl
        << "nchimax = " << lcao_.get_nchimax() << std::endl
        << "nchimax_d = " << lcao_.get_nchimax_d() << std::endl
        << "ntype = " << lcao_.get_ntype() << std::endl
        << "dr_uniform = " << lcao_.get_dr_uniform() << std::endl
        << "rcutmax_Phi = " << lcao_.get_rcutmax_Phi() << std::endl;

}
*/

void LcaoOrbitalsTest::SetUp() {

}


void LcaoOrbitalsTest::TearDown() {

}


TEST_F(LcaoOrbitalsTest, ReadOrbFile) {

    // This test checks whether read_orb_file behaves as expected.
    //
    // read_orb_file is supposed to read an orbital file, which contains
    // all the numerical atomic orbital information of a single element.
    //
    // The information includes the element name, the radial mesh, and
    // certain number of radial functions classified by their angular
    // momentum. 

    // mock GlobalV::ofs_running
    std::ofstream ofs_log("ORB_read_test.log");

    // orbital file to read data from
    std::string orb_file = "./lcao_H2O/O_gga_7au_60Ry_2s2p1d.orb";


    std::ifstream ifs(orb_file);

    // index of element type, not meaningful here
    int iat = 7;

    // maximum angular momentum & number of chi for 
    int lmax = 0;

    // maximum number of chi of all angular momemtum
    int nchimax = 0;

    bool force_flag = true;

    // In order to test read_orb_file, we have to reproduce the necessary
    // steps done by Read_Orbitals and Read_PAO/Read_Descriptor before 
    // read_orb_file is called within them.

	delete[] lcao_.Phi;

    lcao_.dk = 0.01;
    lcao_.ecutwfc = 100.0;
    lcao_.kmesh = calc_nk(lcao_.ecutwfc, lcao_.dk);
    lcao_.dr_uniform = 0.001;
    lcao_.Phi = new Numerical_Orbital[iat+1];

    int my_rank = 0;

    lcao_.read_orb_file(ofs_log, ifs, iat, lmax, nchimax, lcao_.Phi,
            force_flag, my_rank);

    // alias
    Numerical_Orbital& ao = lcao_.Phi[iat];

    ////////////////////////////////////////////////////
    //      check data in Numerical_Orbital object
    ////////////////////////////////////////////////////
    EXPECT_EQ(ao.getType(), iat);
    EXPECT_EQ(ao.getLabel(), "O");
    EXPECT_EQ(ao.getLmax(), 2);
    EXPECT_EQ(ao.getNchi(0), 2);
    EXPECT_EQ(ao.getNchi(1), 2);
    EXPECT_EQ(ao.getNchi(2), 1);
    ASSERT_EQ(ao.getTotal_nchi(), 5);


    ////////////////////////////////////////////////////
    //  check data in Numerical_Orbital_Lm object
    ////////////////////////////////////////////////////

    std::vector<int> L_list{0,0,1,1,2};
    std::vector<int> N_list{0,1,0,1,0};

    for (size_t i = 0; i != 5; ++i) {
        int L = L_list[i], N = N_list[i];
        EXPECT_EQ(ao.PhiLN(L,N).getLabel(), "O");
        EXPECT_EQ(ao.PhiLN(L,N).getType(), iat);
        EXPECT_EQ(ao.PhiLN(L,N).getL(), L);
        EXPECT_EQ(ao.PhiLN(L,N).getChi(), N);
        EXPECT_EQ(ao.PhiLN(L,N).getNr(), 701);
        EXPECT_EQ(ao.PhiLN(L,N).getNk(), lcao_.kmesh);
        EXPECT_EQ(ao.PhiLN(L,N).getDk(), lcao_.dk);
        EXPECT_EQ(ao.PhiLN(L,N).getDruniform(), lcao_.dr_uniform);

        for (int ir = 0; ir != 701; ++ir) {
            EXPECT_DOUBLE_EQ(ao.PhiLN(L,N).getRab(ir), 0.01);
            EXPECT_DOUBLE_EQ(ao.PhiLN(L,N).getRadial(ir), 0.01*ir);
        }
    }

    double max_tol = 1e-12;
    EXPECT_NEAR(ao.PhiLN(0,0).getPsi(0), 1.208504975904e+00, max_tol);
    EXPECT_NEAR(ao.PhiLN(0,0).getPsi(1), 1.208605373194e+00, max_tol);
    EXPECT_NEAR(ao.PhiLN(0,0).getPsi(4), 1.210103935461e+00, max_tol);
    EXPECT_NEAR(ao.PhiLN(0,0).getPsi(699), 4.465396560257e-08, max_tol);
    EXPECT_NEAR(ao.PhiLN(0,0).getPsi(700), 0.0, max_tol);

    EXPECT_NEAR(ao.PhiLN(0,1).getPsi(0), 7.254873428942e-01, max_tol);
    EXPECT_NEAR(ao.PhiLN(0,1).getPsi(1), 7.256666701836e-01, max_tol);
    EXPECT_NEAR(ao.PhiLN(0,1).getPsi(4), 7.283448557011e-01, max_tol);
    EXPECT_NEAR(ao.PhiLN(0,1).getPsi(699), -1.916246212603e-06, max_tol);
    EXPECT_NEAR(ao.PhiLN(0,1).getPsi(700), 0.0, max_tol);

    EXPECT_NEAR(ao.PhiLN(1,0).getPsi(0), 0.0, max_tol);
    EXPECT_NEAR(ao.PhiLN(1,0).getPsi(1), 4.626669306440e-02, max_tol);
    EXPECT_NEAR(ao.PhiLN(1,0).getPsi(4), 1.845014292772e-01, max_tol);
    EXPECT_NEAR(ao.PhiLN(1,0).getPsi(699), 2.870401658966e-07, max_tol);
    EXPECT_NEAR(ao.PhiLN(1,0).getPsi(700), 0.0, max_tol);

    EXPECT_NEAR(ao.PhiLN(1,1).getPsi(0), 0.0, max_tol);
    EXPECT_NEAR(ao.PhiLN(1,1).getPsi(1), 3.375340101333e-02, max_tol);
    EXPECT_NEAR(ao.PhiLN(1,1).getPsi(4), 1.346256082234e-01, max_tol);
    EXPECT_NEAR(ao.PhiLN(1,1).getPsi(699), -2.771091616120e-06, max_tol);
    EXPECT_NEAR(ao.PhiLN(1,1).getPsi(700), 0.0, max_tol);

    EXPECT_NEAR(ao.PhiLN(2,0).getPsi(0), 0.0, max_tol);
    EXPECT_NEAR(ao.PhiLN(2,0).getPsi(1), -3.343626342662e-04, max_tol);
    EXPECT_NEAR(ao.PhiLN(2,0).getPsi(4), -5.337546547975e-03, max_tol);
    EXPECT_NEAR(ao.PhiLN(2,0).getPsi(699), 1.396308876444e-06, max_tol);
    EXPECT_NEAR(ao.PhiLN(2,0).getPsi(700), 0.0, max_tol);
}

/*
TEST_F(LcaoOrbitalsTest, Getters) {
    EXPECT_DOUBLE_EQ(lcao_.get_dr_uniform(), 0.001);
    EXPECT_DOUBLE_EQ(lcao_.get_rcutmax_Phi(), 8.0);

    EXPECT_EQ(lcao_.get_lmax_d(), 2);
    EXPECT_EQ(lcao_.get_nchimax_d(), 2);
}
*/


int main(int argc, char **argv)
{

#ifdef __MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&GlobalV::NPROC);
    MPI_Comm_rank(MPI_COMM_WORLD,&GlobalV::MY_RANK);
#endif

    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

#ifdef __MPI
    MPI_Finalize();
#endif

    return result;
}



