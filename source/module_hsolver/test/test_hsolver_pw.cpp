#include <gtest/gtest.h>
#include <iostream>
#include <vector>

#define private public
#define protected public

#include "module_hsolver/hsolver_pw.h"
#include "hsolver_supplementary_mock.h"
#include "hsolver_pw_sup.h"

#include "module_base/global_variable.h"

/************************************************
 *  unit test of HSolverPW class
 ***********************************************/

/**
 * Tested function:
 *  - test for template float and double respectively:
 *  - 1. solve()
 *  - 2. initDiagh()
 *  - 3. endDiagh()
 *  - 4. hamiltSolvePsiK()
 *  - 5. updatePsiK()
 *  - 6. update_precondition()
 *  - 7. hsolver::HSolver::diagethr (for cases below)
 * 		- set_diagethr, for setting diagethr;
 *  	- reset_diagethr, for updating diagethr;
 * 		- cal_hsolerror, for calculate actually diagethr;
 */
class TestHSolverPW : public ::testing::Test
{
	public:
    ModulePW::PW_Basis_K pwbk;
    hsolver::HSolverPW<float, psi::DEVICE_CPU> hs_f = hsolver::HSolverPW<float, psi::DEVICE_CPU>(&pwbk, nullptr);
    hsolver::HSolverPW<double, psi::DEVICE_CPU> hs_d = hsolver::HSolverPW<double, psi::DEVICE_CPU>(&pwbk, nullptr);

    hamilt::Hamilt<double> hamilt_test_d;
    hamilt::Hamilt<float> hamilt_test_f;

    psi::Psi<std::complex<double>> psi_test_cd;
    psi::Psi<std::complex<float>> psi_test_cf;

    elecstate::ElecState elecstate_test;

    std::string method_test = "cg";

    std::vector<float> ekb_f;

    std::ofstream temp_ofs;
};

TEST_F(TestHSolverPW, solve)
{
    // initial memory and data
    elecstate_test.ekb.create(1, 2);
    this->ekb_f.resize(2);
    psi_test_cf.resize(1, 2, 3);
    psi_test_cd.resize(1, 2, 3);
    GlobalV::nelec = 1.0;
    // check constructor
    EXPECT_EQ(this->hs_f.classname, "HSolverPW");
    EXPECT_NEAR(this->hs_f.diag_ethr, 0.01, 1.0e-7);
    EXPECT_EQ(this->hs_d.classname, "HSolverPW");
    EXPECT_NEAR(this->hs_d.diag_ethr, 0.01, 1.0e-7);
    // check solve()
    EXPECT_EQ(this->hs_f.initialed_psi, false);
    EXPECT_EQ(this->hs_d.initialed_psi, false);
    this->hs_f.solve(&hamilt_test_f, psi_test_cf, &elecstate_test, method_test, true);
    EXPECT_EQ(this->hs_f.initialed_psi, true);
    for (int i = 0; i < psi_test_cf.size(); i++)
    {
        EXPECT_DOUBLE_EQ(psi_test_cf.get_pointer()[i].real(), i + 3);
    }
    EXPECT_DOUBLE_EQ(elecstate_test.ekb.c[0], 4.0);
    EXPECT_DOUBLE_EQ(elecstate_test.ekb.c[1], 7.0);
    EXPECT_DOUBLE_EQ(hsolver::DiagoIterAssist<float>::avg_iter, 0.0);

    this->hs_d.solve(&hamilt_test_d, psi_test_cd, &elecstate_test, method_test, true);
    EXPECT_EQ(this->hs_d.initialed_psi, true);
    EXPECT_DOUBLE_EQ(hsolver::DiagoIterAssist<double>::avg_iter, 0.0);
    for (int i = 0; i < psi_test_cd.size(); i++)
    {
        EXPECT_DOUBLE_EQ(psi_test_cd.get_pointer()[i].real(), i + 3);
    }
    EXPECT_DOUBLE_EQ(elecstate_test.ekb.c[0], 4.0);
    EXPECT_DOUBLE_EQ(elecstate_test.ekb.c[1], 7.0);

    // check initDiagh()
    this->hs_f.method = "dav";
    this->hs_d.method = "dav";
    this->hs_f.initialed_psi = false;
    this->hs_d.initialed_psi = false;
    this->hs_f.initDiagh();
    this->hs_d.initDiagh();
    // will not change state of initialed_psi in initDiagh
    EXPECT_EQ(this->hs_f.initialed_psi, false);
    EXPECT_EQ(this->hs_d.initialed_psi, false);

    // check hamiltSolvePsiK()
    this->hs_f.hamiltSolvePsiK(&hamilt_test_f, psi_test_cf, ekb_f.data());
    this->hs_d.hamiltSolvePsiK(&hamilt_test_d, psi_test_cd, elecstate_test.ekb.c);
    for (int i = 0; i < psi_test_cf.size(); i++)
    {
        EXPECT_DOUBLE_EQ(psi_test_cf.get_pointer()[i].real(), i + 4);
    }
    for (int i = 0; i < psi_test_cd.size(); i++)
    {
        EXPECT_DOUBLE_EQ(psi_test_cf.get_pointer()[i].real(), i + 4);
    }
    EXPECT_DOUBLE_EQ(ekb_f[0], 5.0);
    EXPECT_DOUBLE_EQ(elecstate_test.ekb.c[0], 5.0);
    EXPECT_DOUBLE_EQ(ekb_f[1], 8.0);
    EXPECT_DOUBLE_EQ(elecstate_test.ekb.c[1], 8.0);

    // check endDiagH()
    this->hs_f.initialed_psi = true;
    this->hs_d.initialed_psi = true;
    this->hs_f.endDiagh();
    this->hs_d.endDiagh();
    // will change state of initialed_psi in endDiagh
    EXPECT_EQ(this->hs_f.initialed_psi, true);
    EXPECT_EQ(this->hs_d.initialed_psi, true);

    // check updatePsiK()
    // skip initializing Psi, Psi will not change
    this->hs_f.updatePsiK(&hamilt_test_f, psi_test_cf, 0);
    this->hs_d.updatePsiK(&hamilt_test_d, psi_test_cd, 0);
    for (int i = 0; i < psi_test_cf.size(); i++)
    {
        EXPECT_DOUBLE_EQ(psi_test_cf.get_pointer()[i].real(), i + 4);
    }
    for (int i = 0; i < psi_test_cd.size(); i++)
    {
        EXPECT_DOUBLE_EQ(psi_test_cd.get_pointer()[i].real(), i + 4);
    }
    // check update_precondition()
    this->hs_f.update_precondition(this->hs_f.precondition, 0, psi_test_cf.get_nbasis());
    this->hs_d.update_precondition(this->hs_d.precondition, 0, psi_test_cd.get_nbasis());
    EXPECT_NEAR(this->hs_f.precondition[0], 2.414213657, 1e-8);
    EXPECT_NEAR(this->hs_f.precondition[1], 3.618033886, 1e-8);
    EXPECT_NEAR(this->hs_f.precondition[2], 6.236067772, 1e-8);
    EXPECT_NEAR(this->hs_d.precondition[0], 2.414213562, 1e-8);
    EXPECT_NEAR(this->hs_d.precondition[1], 3.618033989, 1e-8);
    EXPECT_NEAR(this->hs_d.precondition[2], 6.236067977, 1e-8);

    // check diago_ethr
    GlobalV::init_chg = "atomic";
    GlobalV::PW_DIAG_THR = 1e-7;
    GlobalV::CALCULATION = "scf";
    float test_diagethr = hs_f.set_diagethr(0, 1, 1.0);
    EXPECT_NEAR(hs_f.diag_ethr, 0.01, 1.0e-7);
    EXPECT_NEAR(test_diagethr, 0.01, 1.0e-7);
    GlobalV::CALCULATION = "md";
    GlobalV::init_chg = "file";
    test_diagethr = hs_f.set_diagethr(0, 1, 1.0);
    EXPECT_NEAR(test_diagethr, 1e-5, 1.0e-7);
    test_diagethr = hs_f.set_diagethr(0, 2, 1.0);
    EXPECT_NEAR(test_diagethr, 0.01, 1.0e-7);
    test_diagethr = hs_f.set_diagethr(0, 3, 1.0e-3);
    EXPECT_NEAR(test_diagethr, 0.0001, 1.0e-7);
    test_diagethr = hs_f.reset_diagethr(temp_ofs, 2.0, 1.0);
    EXPECT_NEAR(test_diagethr, 0.1, 1.0e-7);
    test_diagethr = hs_f.reset_diagethr(temp_ofs, 0.5, 1.0);
    EXPECT_NEAR(test_diagethr, 0.1, 1.0e-7);
    test_diagethr = hs_f.cal_hsolerror();
    EXPECT_NEAR(test_diagethr, 0.1, 1.0e-7);

    GlobalV::init_chg = "atomic";
    GlobalV::PW_DIAG_THR = 1e-7;
    GlobalV::CALCULATION = "scf";
    double test_diagethr_d = hs_d.set_diagethr(0, 1, 1.0);
    EXPECT_EQ(hs_d.diag_ethr, 0.01);
    EXPECT_EQ(test_diagethr_d, 0.01);
    GlobalV::CALCULATION = "md";
    GlobalV::init_chg = "file";
    test_diagethr_d = hs_d.set_diagethr(0, 1, 1.0);
    EXPECT_EQ(test_diagethr_d, 1e-5);
    test_diagethr_d = hs_d.set_diagethr(0, 2, 1.0);
    EXPECT_EQ(test_diagethr_d, 0.01);
    test_diagethr_d = hs_d.set_diagethr(0, 3, 1.0e-3);
    EXPECT_EQ(test_diagethr_d, 0.0001);
    test_diagethr_d = hs_d.reset_diagethr(temp_ofs, 2.0, 1.0);
    EXPECT_EQ(test_diagethr_d, 0.1);
    test_diagethr_d = hs_d.reset_diagethr(temp_ofs, 0.5, 1.0);
    EXPECT_EQ(test_diagethr_d, 0.1);
    test_diagethr_d = hs_d.cal_hsolerror();
    EXPECT_EQ(test_diagethr_d, 0.1);
}

/*#include "mpi.h"
#include "module_base/timer.h"
int main(int argc, char **argv)
{
	ModuleBase::timer::disable();
	MPI_Init(&argc, &argv);
	testing::InitGoogleTest(&argc, argv);

	MPI_Comm_size(MPI_COMM_WORLD,&GlobalV::NPROC);
	MPI_Comm_rank(MPI_COMM_WORLD,&GlobalV::MY_RANK);
	int result = RUN_ALL_TESTS();
	
	MPI_Finalize();
	
	return result;
}*/

int main(int argc, char **argv)
{
        MPI_Init(&argc, &argv);
        testing::InitGoogleTest(&argc, argv);

        int nproc;
        int my_rank;

        MPI_Comm_size(MPI_COMM_WORLD,&nproc);
        MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

        int result = RUN_ALL_TESTS();

        MPI_Finalize();

        return result;
}