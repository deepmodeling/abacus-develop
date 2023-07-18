#include <gtest/gtest.h>
#include <iostream>
#include <vector>

#define private public
#define protected public

#include "module_hsolver/hsolver_pw.h"
#include "hsolver_supplementary_mock.h"
#include "hsolver_pw_sup.h"
#include "module_hsolver/hsolver_pw_sdft.h"

#include "module_base/global_variable.h"

//mock for module_sdft
template<typename REAL>
Sto_Func<REAL>::Sto_Func(){}

template class Sto_Func<double>;

Stochastic_hchi::Stochastic_hchi(){};
Stochastic_hchi::~Stochastic_hchi(){};

Stochastic_Iter::Stochastic_Iter()
{
    change = false;
    mu0 = 0;
    method = 2;
}

Stochastic_Iter::~Stochastic_Iter()
{
    delete p_che;
    delete[] spolyv;
    delete[] chiallorder;
}

void Stochastic_Iter::init(int *nchip_in, const int method_in, K_Vectors* pkv, ModulePW::PW_Basis_K *wfc_basis, Stochastic_WF &stowf)
{
    this->nchip = nchip_in;
    this->targetne = 1;
    this->method = method_in;

}

void Stochastic_Iter::orthog(
    const int &ik, 
    psi::Psi<ModuleBase::ComplexMatrix::type, psi::DEVICE_CPU> &psi, 
    Stochastic_WF &stowf
)
{
    //do something to verify this function has been called
    for(int i=0;i<psi.size();i++)
    {
        psi.get_pointer()[i] += 1.1;
    }
    return;
}

void Stochastic_Iter::checkemm(
    const int &ik, 
    int istep, 
    int iter, 
    Stochastic_WF &stowf
)
{
    //do something to verify this function has been called
    stowf.nchi++;
    return;
}

void Stochastic_Iter::calPn(
    const int &ik, 
    Stochastic_WF &stowf
)
{
    //do something to verify this function has been called
    stowf.nbands_diag ++;
    return;
}

void Stochastic_Iter::itermu(
    int iter, 
    elecstate::ElecState *pes
)
{
    //do something to verify this function has been called
    pes->f_en.eband += 1.2;
    return;
}

void Stochastic_Iter::calHsqrtchi(Stochastic_WF &stowf)
{
    //do something to verify this function has been called
    stowf.nchip_max++;
    return;
}

void Stochastic_Iter::sum_stoband(
    Stochastic_WF &stowf, 
    elecstate::ElecState *pes, 
    hamilt::Hamilt<double, psi::DEVICE_CPU> *pHamilt,
    ModulePW::PW_Basis_K* wfc_basis
)
{
    //do something to verify this function has been called
    stowf.nbands_total ++;
    return;
}

Charge::Charge(){};
Charge::~Charge(){};

/************************************************
 *  unit test of HSolverPW_SDFT class
 ***********************************************/

/**
 * Tested function:
 *  - 1. solve()
 *      - with psi;
 *      - without psi;
 *      - skip charge;
 *  - 2. hsolver::HSolverPW_SDFT::diagethr (for cases below)
 * 		- set_diagethr, for setting diagethr;
 * 		- cal_hsolerror, for calculate actually diagethr;
 */
class TestHSolverPW_SDFT : public ::testing::Test
{
	public:
    ModulePW::PW_Basis_K pwbk;
    Stochastic_WF stowf;
    K_Vectors kv;
    wavefunc wf;
    hsolver::HSolverPW_SDFT hs_d = hsolver::HSolverPW_SDFT(&kv, &pwbk, &wf, stowf, 0);

    hamilt::Hamilt<double> hamilt_test_d;

	psi::Psi<std::complex<double>> psi_test_cd;
    psi::Psi<std::complex<double>> psi_test_no;

	elecstate::ElecState elecstate_test;

	std::string method_test = "cg";

	std::ofstream temp_ofs;
};

TEST_F(TestHSolverPW_SDFT, solve)
{
	//initial memory and data
	elecstate_test.ekb.create(1,2);
    elecstate_test.f_en.eband = 0.0;
    stowf.nbands_diag = 0;
    stowf.nbands_total = 0;
    stowf.nchi = 0;
    stowf.nchip_max = 0;
	psi_test_cd.resize(1, 2, 3);
	GlobalV::nelec = 1.0;
    GlobalV::MY_STOGROUP = 0.0;
    int istep = 0;
    int iter = 0;
	//check constructor
	EXPECT_EQ(this->hs_d.classname, "HSolverPW_SDFT");
	EXPECT_NEAR(this->hs_d.diag_ethr, 0.01, 1.0e-7);
	//check solve()
	EXPECT_EQ(this->hs_d.initialed_psi, false);

	this->hs_d.solve(
        &hamilt_test_d, 
        psi_test_cd, 
        &elecstate_test,
        &pwbk, 
        stowf, 
        istep, 
        iter, 
        method_test, 
        false
    );
	EXPECT_EQ(this->hs_d.initialed_psi, true);
	EXPECT_DOUBLE_EQ(hsolver::DiagoIterAssist<double>::avg_iter, 0.0);
	EXPECT_DOUBLE_EQ(elecstate_test.ekb.c[0], 4.0);
	EXPECT_DOUBLE_EQ(elecstate_test.ekb.c[1], 7.0);
    for(int i=0;i<psi_test_cd.size();i++)
	{
        //std::cout<<__FILE__<<__LINE__<<" "<<psi_test_cd.get_pointer()[i]<<std::endl;
		EXPECT_DOUBLE_EQ(psi_test_cd.get_pointer()[i].real(), double(i)+4.1);
	}
    EXPECT_EQ(stowf.nbands_diag, 1);
    EXPECT_EQ(stowf.nbands_total, 1);
    EXPECT_EQ(stowf.nchi, 1);
    EXPECT_EQ(stowf.nchip_max, 1);
    EXPECT_DOUBLE_EQ(elecstate_test.f_en.eband, 1.2);
    /*std::cout<<__FILE__<<__LINE__<<" "<<stowf.nbands_diag<<std::endl;
    std::cout<<__FILE__<<__LINE__<<" "<<stowf.nbands_total<<std::endl;
    std::cout<<__FILE__<<__LINE__<<" "<<stowf.nchi<<std::endl;
    std::cout<<__FILE__<<__LINE__<<" "<<stowf.nchip_max<<std::endl;
    std::cout<<__FILE__<<__LINE__<<" "<<elecstate_test.f_en.eband<<std::endl;*/

    //check diago_ethr
	GlobalV::init_chg = "atomic";
	GlobalV::PW_DIAG_THR = 1e-7;
	GlobalV::CALCULATION = "scf";
	double test_diagethr_d = hs_d.set_diagethr(0, 1, 1.0);
	//std::cout<<__FILE__<<__LINE__<<" "<<test_diagethr_d<<std::endl;
	EXPECT_EQ(hs_d.diag_ethr, 0.01);
	EXPECT_EQ(test_diagethr_d, 0.01);
	GlobalV::CALCULATION = "md";
	GlobalV::init_chg = "file";
	test_diagethr_d = hs_d.set_diagethr(0, 1, 1.0);
	//std::cout<<__FILE__<<__LINE__<<" "<<test_diagethr_d<<std::endl;
    EXPECT_EQ(test_diagethr_d, 1e-5);
	test_diagethr_d = hs_d.set_diagethr(0, 2, 1.0);
	//std::cout<<__FILE__<<__LINE__<<" "<<test_diagethr_d<<std::endl;
	EXPECT_EQ(test_diagethr_d, 0);
	test_diagethr_d = hs_d.set_diagethr(0, 3, 1.0e-3);
	//std::cout<<__FILE__<<__LINE__<<" "<<test_diagethr_d<<std::endl;
	EXPECT_EQ(test_diagethr_d, 0);
    test_diagethr_d = hs_d.cal_hsolerror();
	EXPECT_EQ(test_diagethr_d, 0.0);


}

TEST_F(TestHSolverPW_SDFT, solve_noband_skipcharge)
{
	//initial memory and data
	elecstate_test.ekb.create(1,2);
    elecstate_test.f_en.eband = 0.0;
    stowf.nbands_diag = 0;
    stowf.nbands_total = 0;
    stowf.nchi = 0;
    stowf.nchip_max = 0;
    psi_test_no.nk = 2;
    psi_test_no.nbands = 0;
    psi_test_no.nbasis = 0;
	GlobalV::nelec = 1.0;
    GlobalV::MY_STOGROUP = 0.0;
    GlobalV::NSPIN = 1;
    elecstate_test.charge = new Charge;
    elecstate_test.charge->rho = new double*[1];
    elecstate_test.charge->rho[0] = new double[10];
    elecstate_test.charge->nrxx = 10;
    int istep = 0;
    int iter = 0;
	//check constructor
	EXPECT_EQ(this->hs_d.classname, "HSolverPW_SDFT");
	EXPECT_NEAR(this->hs_d.diag_ethr, 1e-7, 1.0e-10);
	//check solve()
    hs_d.initialed_psi = true;

	this->hs_d.solve(
        &hamilt_test_d, 
        psi_test_no, 
        &elecstate_test, 
        &pwbk,
        stowf, 
        istep, 
        iter, 
        method_test, 
        false
    );
	EXPECT_DOUBLE_EQ(hsolver::DiagoIterAssist<double>::avg_iter, 0.0);
    EXPECT_EQ(stowf.nbands_diag, 2);
    EXPECT_EQ(stowf.nbands_total, 1);
    EXPECT_EQ(stowf.nchi, 2);
    EXPECT_EQ(stowf.nchip_max, 1);
    EXPECT_DOUBLE_EQ(elecstate_test.f_en.eband, 1.2);
    /*std::cout<<__FILE__<<__LINE__<<" "<<stowf.nbands_diag<<std::endl;
    std::cout<<__FILE__<<__LINE__<<" "<<stowf.nbands_total<<std::endl;
    std::cout<<__FILE__<<__LINE__<<" "<<stowf.nchi<<std::endl;
    std::cout<<__FILE__<<__LINE__<<" "<<stowf.nchip_max<<std::endl;
    std::cout<<__FILE__<<__LINE__<<" "<<elecstate_test.f_en.eband<<std::endl;*/

    //test for skip charge
    this->hs_d.solve(
        &hamilt_test_d, 
        psi_test_no, 
        &elecstate_test, 
        &pwbk,
        stowf, 
        istep, 
        iter, 
        method_test, 
        true
    );
    EXPECT_EQ(stowf.nbands_diag, 4);
    EXPECT_EQ(stowf.nbands_total, 1);
    EXPECT_EQ(stowf.nchi, 4);
    EXPECT_EQ(stowf.nchip_max, 2);
    EXPECT_DOUBLE_EQ(elecstate_test.f_en.eband, 2.4);

    delete[] elecstate_test.charge->rho[0];
    delete[] elecstate_test.charge->rho;
    delete elecstate_test.charge;

}

#ifdef __MPI
#include "mpi.h"
#include "module_base/timer.h"
int main(int argc, char **argv)
{
	ModuleBase::timer::disable();
	MPI_Init(&argc, &argv);
	testing::InitGoogleTest(&argc, argv);

	MPI_Comm_size(MPI_COMM_WORLD,&GlobalV::NPROC);
	MPI_Comm_rank(MPI_COMM_WORLD,&GlobalV::MY_RANK);
    MPI_Comm_split(MPI_COMM_WORLD, 0, 1, &PARAPW_WORLD);
	int result = RUN_ALL_TESTS();
	
    MPI_Comm_free(&PARAPW_WORLD);
	MPI_Finalize();
	
	return result;
}

#endif