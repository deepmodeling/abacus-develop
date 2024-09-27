#include <gtest/gtest.h>
#include <iostream>
#define private public
#include "module_parameter/parameter.h"
#undef private
#include <vector>

#define private public
#define protected public
#include "module_hsolver/hsolver_pw.h"
#include "hsolver_supplementary_mock.h"
#include "hsolver_pw_sup.h"
#include "module_hsolver/hsolver_pw_sdft.h"
#include "module_base/global_variable.h"
#undef private
#undef protected

//mock for module_sdft
template<typename REAL>
Sto_Func<REAL>::Sto_Func(){}

template class Sto_Func<double>;

template<typename REAL>
StoChe<REAL>::StoChe(const int& nche, const int& method, const REAL& emax_sto, const REAL& emin_sto)
{
    this->nche = nche;
}
template<typename REAL>
StoChe<REAL>::~StoChe(){}

template class StoChe<double>;

Stochastic_hchi::Stochastic_hchi(){};
Stochastic_hchi::~Stochastic_hchi(){};

Stochastic_Iter::Stochastic_Iter()
{
    change = false;
    mu0 = 0;
    method = 2;
}

Stochastic_Iter::~Stochastic_Iter(){};

void Stochastic_Iter::init(K_Vectors* pkv, ModulePW::PW_Basis_K *wfc_basis, Stochastic_WF &stowf, StoChe<double> &stoche)
{
    this->nchip = stowf.nchip;;
    this->targetne = 1;
    this->method = stoche.method_sto;
}

void Stochastic_Iter::orthog(const int& ik,
                             psi::Psi<ModuleBase::ComplexMatrix::type, base_device::DEVICE_CPU>& psi,
                             Stochastic_WF& stowf)
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

void Stochastic_Iter::sum_stoband(Stochastic_WF& stowf,
                                  elecstate::ElecState* pes,
                                  hamilt::Hamilt<std::complex<double>, base_device::DEVICE_CPU>* pHamilt,
                                  ModulePW::PW_Basis_K* wfc_basis)
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
 */
class TestHSolverPW_SDFT : public ::testing::Test
{
	public:
    TestHSolverPW_SDFT():stoche(8,1,0,0){}
    ModulePW::PW_Basis_K pwbk;
    Stochastic_WF stowf;
    K_Vectors kv;
    wavefunc wf;
    StoChe<double> stoche;
    hsolver::HSolverPW_SDFT hs_d = hsolver::HSolverPW_SDFT(&kv, 
                                                           &pwbk, 
                                                           &wf, 
                                                           stowf, 
                                                           stoche, 
                                                           
                                                           "scf",
                                                           "pw",
                                                           "cg",
                                                           false,
                                                           PARAM.sys.use_uspp,
                                                           PARAM.input.nspin,
                     hsolver::DiagoIterAssist<std::complex<double>>::SCF_ITER,
                     hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_NMAX,
                     hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_THR,
                     hsolver::DiagoIterAssist<std::complex<double>>::need_subspace,
                     false);

    hamilt::Hamilt<std::complex<double>> hamilt_test_d;

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
	PARAM.input.nelec = 1.0;
    GlobalV::MY_STOGROUP = 0.0;
    int istep = 0;
    int iter = 0;

    this->hs_d.solve(&hamilt_test_d,
                     psi_test_cd,
                     &elecstate_test,
                     &pwbk,
                     stowf,
                     istep,
                     iter,
                     false);
	EXPECT_DOUBLE_EQ(hsolver::DiagoIterAssist<std::complex<double>>::avg_iter, 0.0);
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
	PARAM.input.nelec = 1.0;
    GlobalV::MY_STOGROUP = 0.0;
    PARAM.input.nspin = 1;
    elecstate_test.charge = new Charge;
    elecstate_test.charge->rho = new double*[1];
    elecstate_test.charge->rho[0] = new double[10];
    elecstate_test.charge->nrxx = 10;
    int istep = 0;
    int iter = 0;

    this->hs_d.solve(&hamilt_test_d,
                     psi_test_no,
                     &elecstate_test,
                     &pwbk,
                     stowf,
                     istep,
                     iter,
                     false
    );
	EXPECT_DOUBLE_EQ(hsolver::DiagoIterAssist<std::complex<double>>::avg_iter, 0.0);
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
    this->hs_d.solve(&hamilt_test_d,
                     psi_test_no,
                     &elecstate_test,
                     &pwbk,
                     stowf,
                     istep,
                     iter,
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