#include <fstream>
#include "gtest/gtest.h"
#include "module_base/global_variable.h"
#include "module_basis/module_ao/ORB_atomic.h"
#include "module_basis/module_ao/ORB_atomic_lm.h"

#define private public
#include "module_basis/module_ao/ORB_read.h"
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
 * - Read_Orbitals
 *   read orbital & (optionally) descriptor files and pour the data
 *   into Numerical_Orbital (and its member Numerical_Orbital_Lm) objects
 *
 * - bcast_files
 *   broadcast the orbital_file from rank-0 to all other processes
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

    // initialize lcao_ with parameters below & call Read_Orbitals
    void lcao_read();

    // parameters to pass to lcao_
    std::ofstream ofs_log_;
    int ntype_;
    int lmax_;
    int out_mat_r_; // unused variable
    bool force_flag_;
    int my_rank_;
    bool deepks_setorb_;
    bool read_in_flag_;
    std::string descriptor_file_;
    std::vector<std::string> orbital_file_;
    double ecutwfc_;
    double dk_;
    double dR_;
    double Rmax_;
};


void LcaoOrbitalsTest::SetUp() {
    ofs_log_.open("ORB_read_test.log");
    ntype_ = 2;
    lmax_ = 2;
    out_mat_r_ = 0; // unused variable
    force_flag_ = true;
    my_rank_ = 0;
    deepks_setorb_ = true;

    read_in_flag_ = true;
    descriptor_file_ = "./lcao_H2O/jle.orb";
    orbital_file_.push_back("./lcao_H2O/H_gga_8au_60Ry_2s1p.orb");
    orbital_file_.push_back("./lcao_H2O/O_gga_7au_60Ry_2s2p1d.orb");

    // below we mimic ORB_control::read_orb_first
    ecutwfc_ = 123.0;
    dk_ = 0.01;
    dR_ = 0.01;
    Rmax_ = 20;
}


void LcaoOrbitalsTest::lcao_read() {

    // see UnitCell::read_atom_species in module_cell/read_atoms.cpp
    lcao_.read_in_flag = read_in_flag_;
    lcao_.descriptor_file = descriptor_file_;
    lcao_.orbital_file = orbital_file_;
#ifdef __MPI
    lcao_.bcast_files(ntype_, GlobalV::MY_RANK);
#endif

    // see ORB_control::read_orb_first
    lcao_.ecutwfc = ecutwfc_;
    lcao_.dk = dk_;
    lcao_.dR = dR_;
    lcao_.Rmax = Rmax_;

    lcao_.Read_Orbitals(ofs_log_, ntype_, lmax_, deepks_setorb_, out_mat_r_, 
            force_flag_, my_rank_);
}


void LcaoOrbitalsTest::TearDown() {
}


TEST_F(LcaoOrbitalsTest, ReadInFlag) {
    read_in_flag_ = false;
    EXPECT_EXIT(this->lcao_read(), testing::ExitedWithCode(0), "");
}


TEST_F(LcaoOrbitalsTest, WrongOrbFile) {
    orbital_file_[0] = "./lcao_H2O/H_gga_8au_60Ry_2s1.orb";
    EXPECT_EXIT(this->lcao_read(), testing::ExitedWithCode(0), "");
}


TEST_F(LcaoOrbitalsTest, WrongDescFile) {
    descriptor_file_ = "./lcao_H2O/jl.orb";
    EXPECT_EXIT(this->lcao_read(), testing::ExitedWithCode(0), "");
}


TEST_F(LcaoOrbitalsTest, BcastFiles) {
#ifdef __MPI
    if ( GlobalV::MY_RANK == 0 ) {
        lcao_.orbital_file = orbital_file_;
    }

    if ( GlobalV::MY_RANK != 0) {
        EXPECT_EQ(lcao_.orbital_file, std::vector<std::string>{});
    }

    lcao_.read_in_flag = read_in_flag_;
    lcao_.bcast_files(2, GlobalV::MY_RANK);

    EXPECT_EQ(lcao_.orbital_file, orbital_file_);
#endif
}


TEST_F(LcaoOrbitalsTest, ReadOrbitals) {
    
    this->lcao_read();

    // This test checks whether Read_Orbitals behaves as expected.

    // alias
    Numerical_Orbital& ao0 = lcao_.Phi[0];
    Numerical_Orbital& ao1 = lcao_.Phi[1];
    Numerical_Orbital& aod = lcao_.Alpha[0];

    // we now check whether the contents of the above
    // Numerical_Orbital objects are correct

    // maximum tolerance for double-type comparison
    double max_tol = 1e-12;

    // H
    EXPECT_EQ(ao0.getType(), 0);
    EXPECT_EQ(ao0.getLabel(), "H");
    EXPECT_EQ(ao0.getLmax(), 1);
    EXPECT_EQ(ao0.getNchi(0), 2);
    EXPECT_EQ(ao0.getNchi(1), 1);
    ASSERT_EQ(ao0.getTotal_nchi(), 3);

    std::vector<int> L0_list{0,0,1};
    std::vector<int> N0_list{0,1,0};

    for (size_t i = 0; i != 3; ++i) {
        int L = L0_list[i], N = N0_list[i];
        EXPECT_EQ(ao0.PhiLN(L,N).getLabel(), "H");
        EXPECT_EQ(ao0.PhiLN(L,N).getType(), 0);
        EXPECT_EQ(ao0.PhiLN(L,N).getL(), L);
        EXPECT_EQ(ao0.PhiLN(L,N).getChi(), N);
        EXPECT_EQ(ao0.PhiLN(L,N).getNr(), 801);
        EXPECT_EQ(ao0.PhiLN(L,N).getNk(), lcao_.kmesh);
        EXPECT_EQ(ao0.PhiLN(L,N).getDk(), lcao_.dk);
        EXPECT_EQ(ao0.PhiLN(L,N).getDruniform(), lcao_.dr_uniform);

        for (int ir = 0; ir != 801; ++ir) {
            EXPECT_DOUBLE_EQ(ao0.PhiLN(L,N).getRab(ir), 0.01);
            EXPECT_DOUBLE_EQ(ao0.PhiLN(L,N).getRadial(ir), 0.01*ir);
        }
    }

    EXPECT_NEAR(ao0.PhiLN(0,0).getPsi(0  ), 1.837183001954e+00, max_tol);
    EXPECT_NEAR(ao0.PhiLN(0,0).getPsi(1  ), 1.836944589913e+00, max_tol);
    EXPECT_NEAR(ao0.PhiLN(0,0).getPsi(4  ), 1.833374417163e+00, max_tol);
    EXPECT_NEAR(ao0.PhiLN(0,0).getPsi(799), 3.037233152557e-07, max_tol);
    EXPECT_NEAR(ao0.PhiLN(0,0).getPsi(800), 0.000000000000e+00, max_tol);

    EXPECT_NEAR(ao0.PhiLN(0,1).getPsi(0  ), -2.482045090982e+00, max_tol);
    EXPECT_NEAR(ao0.PhiLN(0,1).getPsi(1  ), -2.481575045574e+00, max_tol);
    EXPECT_NEAR(ao0.PhiLN(0,1).getPsi(4  ), -2.474535579529e+00, max_tol);
    EXPECT_NEAR(ao0.PhiLN(0,1).getPsi(799), 1.115867959482e-06, max_tol);
    EXPECT_NEAR(ao0.PhiLN(0,1).getPsi(800), 0.000000000000e+00, max_tol);

    EXPECT_NEAR(ao0.PhiLN(1,0).getPsi(0  ), 0.000000000000e+00, max_tol);
    EXPECT_NEAR(ao0.PhiLN(1,0).getPsi(1  ), -2.619148756396e-02, max_tol);
    EXPECT_NEAR(ao0.PhiLN(1,0).getPsi(4  ), -1.045849793771e-01, max_tol);
    EXPECT_NEAR(ao0.PhiLN(1,0).getPsi(799), 3.217573100688e-06, max_tol);
    EXPECT_NEAR(ao0.PhiLN(1,0).getPsi(800), 0.000000000000e+00, max_tol);


    // O
    EXPECT_EQ(ao1.getType(), 1);
    EXPECT_EQ(ao1.getLabel(), "O");
    EXPECT_EQ(ao1.getLmax(), 2);
    EXPECT_EQ(ao1.getNchi(0), 2);
    EXPECT_EQ(ao1.getNchi(1), 2);
    EXPECT_EQ(ao1.getNchi(2), 1);
    ASSERT_EQ(ao1.getTotal_nchi(), 5);

    std::vector<int> L1_list{0,0,1,1,2};
    std::vector<int> N1_list{0,1,0,1,0};

    for (size_t i = 0; i != 5; ++i) {
        int L = L1_list[i], N = N1_list[i];
        EXPECT_EQ(ao1.PhiLN(L,N).getLabel(), "O");
        EXPECT_EQ(ao1.PhiLN(L,N).getType(), 1);
        EXPECT_EQ(ao1.PhiLN(L,N).getL(), L);
        EXPECT_EQ(ao1.PhiLN(L,N).getChi(), N);
        EXPECT_EQ(ao1.PhiLN(L,N).getNr(), 701);
        EXPECT_EQ(ao1.PhiLN(L,N).getNk(), lcao_.kmesh);
        EXPECT_EQ(ao1.PhiLN(L,N).getDk(), lcao_.dk);
        EXPECT_EQ(ao1.PhiLN(L,N).getDruniform(), lcao_.dr_uniform);

        for (int ir = 0; ir != 701; ++ir) {
            EXPECT_DOUBLE_EQ(ao1.PhiLN(L,N).getRab(ir), 0.01);
            EXPECT_DOUBLE_EQ(ao1.PhiLN(L,N).getRadial(ir), 0.01*ir);
        }
    }

    EXPECT_NEAR(ao1.PhiLN(0,0).getPsi(0), 1.208504975904e+00, max_tol);
    EXPECT_NEAR(ao1.PhiLN(0,0).getPsi(1), 1.208605373194e+00, max_tol);
    EXPECT_NEAR(ao1.PhiLN(0,0).getPsi(4), 1.210103935461e+00, max_tol);
    EXPECT_NEAR(ao1.PhiLN(0,0).getPsi(699), 4.465396560257e-08, max_tol);
    EXPECT_NEAR(ao1.PhiLN(0,0).getPsi(700), 0.0, max_tol);

    EXPECT_NEAR(ao1.PhiLN(0,1).getPsi(0), 7.254873428942e-01, max_tol);
    EXPECT_NEAR(ao1.PhiLN(0,1).getPsi(1), 7.256666701836e-01, max_tol);
    EXPECT_NEAR(ao1.PhiLN(0,1).getPsi(4), 7.283448557011e-01, max_tol);
    EXPECT_NEAR(ao1.PhiLN(0,1).getPsi(699), -1.916246212603e-06, max_tol);
    EXPECT_NEAR(ao1.PhiLN(0,1).getPsi(700), 0.0, max_tol);

    EXPECT_NEAR(ao1.PhiLN(1,0).getPsi(0), 0.0, max_tol);
    EXPECT_NEAR(ao1.PhiLN(1,0).getPsi(1), 4.626669306440e-02, max_tol);
    EXPECT_NEAR(ao1.PhiLN(1,0).getPsi(4), 1.845014292772e-01, max_tol);
    EXPECT_NEAR(ao1.PhiLN(1,0).getPsi(699), 2.870401658966e-07, max_tol);
    EXPECT_NEAR(ao1.PhiLN(1,0).getPsi(700), 0.0, max_tol);

    EXPECT_NEAR(ao1.PhiLN(1,1).getPsi(0), 0.0, max_tol);
    EXPECT_NEAR(ao1.PhiLN(1,1).getPsi(1), 3.375340101333e-02, max_tol);
    EXPECT_NEAR(ao1.PhiLN(1,1).getPsi(4), 1.346256082234e-01, max_tol);
    EXPECT_NEAR(ao1.PhiLN(1,1).getPsi(699), -2.771091616120e-06, max_tol);
    EXPECT_NEAR(ao1.PhiLN(1,1).getPsi(700), 0.0, max_tol);

    EXPECT_NEAR(ao1.PhiLN(2,0).getPsi(0), 0.0, max_tol);
    EXPECT_NEAR(ao1.PhiLN(2,0).getPsi(1), -3.343626342662e-04, max_tol);
    EXPECT_NEAR(ao1.PhiLN(2,0).getPsi(4), -5.337546547975e-03, max_tol);
    EXPECT_NEAR(ao1.PhiLN(2,0).getPsi(699), 1.396308876444e-06, max_tol);
    EXPECT_NEAR(ao1.PhiLN(2,0).getPsi(700), 0.0, max_tol);


    // Descriptor

    EXPECT_EQ(aod.getType(), 0);
    EXPECT_EQ(aod.getLabel(), "");
    EXPECT_EQ(aod.getLmax(), 2);
    EXPECT_EQ(aod.getNchi(0), 2);
    EXPECT_EQ(aod.getNchi(1), 2);
    EXPECT_EQ(aod.getNchi(2), 2);
    ASSERT_EQ(aod.getTotal_nchi(), 6);

    std::vector<int> Ld_list{0,0,1,1,2,2};
    std::vector<int> Nd_list{0,1,0,1,0,1};

    for (size_t i = 0; i != 6; ++i) {
        int L = Ld_list[i], N = Nd_list[i];
        EXPECT_EQ(aod.PhiLN(L,N).getLabel(), "");
        EXPECT_EQ(aod.PhiLN(L,N).getType(), 0);
        EXPECT_EQ(aod.PhiLN(L,N).getL(), L);
        EXPECT_EQ(aod.PhiLN(L,N).getChi(), N);
        EXPECT_EQ(aod.PhiLN(L,N).getNr(), 205);
        EXPECT_EQ(aod.PhiLN(L,N).getNk(), lcao_.kmesh);
        EXPECT_EQ(aod.PhiLN(L,N).getDk(), lcao_.dk);
        EXPECT_EQ(aod.PhiLN(L,N).getDruniform(), lcao_.dr_uniform);

        for (int ir = 0; ir != 205; ++ir) {
            EXPECT_DOUBLE_EQ(aod.PhiLN(L,N).getRab(ir), 0.01);
            EXPECT_DOUBLE_EQ(aod.PhiLN(L,N).getRadial(ir), 0.01*ir);
        }
    }

    // TODO chi value check is skipped for now
    // orbitals in jle.orb are not normalized
    // getPsi() does not gives the numbers in jle.orb
}


TEST_F(LcaoOrbitalsTest, Getters) {

    this->lcao_read();

    EXPECT_EQ(lcao_.get_ecutwfc(), lcao_.ecutwfc);
    EXPECT_EQ(lcao_.get_kmesh(), lcao_.kmesh);
    EXPECT_EQ(lcao_.get_dk(), lcao_.dk);
    EXPECT_EQ(lcao_.get_dR(), lcao_.dR);
    EXPECT_EQ(lcao_.get_Rmax(), lcao_.Rmax);
    EXPECT_EQ(lcao_.get_lmax(), lcao_.lmax);
    EXPECT_EQ(lcao_.get_lmax_d(), lcao_.lmax_d);
    EXPECT_EQ(lcao_.get_nchimax(), lcao_.nchimax);
    EXPECT_EQ(lcao_.get_nchimax_d(), lcao_.nchimax_d);
    EXPECT_EQ(lcao_.get_ntype(), lcao_.ntype);
    EXPECT_EQ(lcao_.get_dr_uniform(), lcao_.dr_uniform);
    EXPECT_EQ(lcao_.get_rcutmax_Phi(), lcao_.rcutmax_Phi);
}


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


