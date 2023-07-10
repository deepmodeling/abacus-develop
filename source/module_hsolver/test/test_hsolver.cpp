#include <gtest/gtest.h>
#include <iostream>
#include <vector>
#define protected public

#include "module_hsolver/hsolver.h"
#include "hsolver_supplementary_mock.h"

template class hsolver::HSolver<double, psi::DEVICE_CPU>;
template class hsolver::HSolver<float, psi::DEVICE_CPU>;

/************************************************
 *  unit test of HSolver base class
 ***********************************************/

/**
 * Tested function:
 *  - hsolver::HSolver::solve (for cases below)
 * 		- float/double;
 * 		- Psi<double/float/complex<double>/complex<float>>;
 * 		- skip charge or not;
 * 		- stodft case;
 *  - hsolver::HSolver::diagethr (for cases below)
 * 		- set_diagethr, for setting diagethr;
 *  	- reset_diagethr, for updating diagethr;
 * 		- cal_hsolerror, for calculate actually diagethr;
 *  - hsolver::DiagH (for cases below)
 *      - diag() for Psi(FPTYPE) case
 *      - destructor of DiagH and HSolver
 *
 * the definition of supplementary functions is added in hsolver_supplementary_mock.h 
 */

class TestHSolver : public ::testing::Test
{
	public:
	hsolver::HSolver<float, psi::DEVICE_CPU> hs_f;
	hsolver::HSolver<double, psi::DEVICE_CPU> hs_d;

	hamilt::Hamilt<double> hamilt_test_d;
	hamilt::Hamilt<float> hamilt_test_f;

	psi::Psi<std::complex<double>> psi_test_cd;
	psi::Psi<std::complex<float>> psi_test_cf;
	psi::Psi<double> psi_test_d;
	psi::Psi<float> psi_test_f;

	Stochastic_WF stowf_test;
	
	elecstate::ElecState elecstate_test;

	ModulePW::PW_Basis_K* wfcpw;

	std::string method_test = "none";

	std::ofstream temp_ofs;
};

TEST_F(TestHSolver, solve)
{
	hs_f.solve(&hamilt_test_f, psi_test_cf, &elecstate_test, method_test, false);
	hs_f.solve(&hamilt_test_f, psi_test_f, &elecstate_test, method_test, false);
	hs_d.solve(&hamilt_test_d, psi_test_cd, &elecstate_test, method_test, false);
	hs_d.solve(&hamilt_test_d, psi_test_d, &elecstate_test, method_test, false);
	hs_f.solve(&hamilt_test_f, psi_test_cf, &elecstate_test, method_test, true);
	hs_f.solve(&hamilt_test_f, psi_test_f, &elecstate_test, method_test, true);
	hs_d.solve(&hamilt_test_d, psi_test_cd, &elecstate_test, method_test, true);
	hs_d.solve(&hamilt_test_d, psi_test_d, &elecstate_test, method_test, true);
	hs_f.solve(&hamilt_test_f, psi_test_cf, &elecstate_test, wfcpw, stowf_test, 0, 0, method_test, true);
	hs_d.solve(&hamilt_test_d, psi_test_cd, &elecstate_test, wfcpw, stowf_test, 0, 0, method_test, true);
	EXPECT_EQ(hs_f.classname, "none");
	EXPECT_EQ(hs_d.classname, "none");
	EXPECT_EQ(hs_f.method, "none");
	EXPECT_EQ(hs_d.method, "none");
}

TEST_F(TestHSolver, diagethr)
{
    float test_diagethr = hs_f.set_diagethr(0, 0, 0.0);
	EXPECT_EQ(hs_f.diag_ethr, 0.0);
	EXPECT_EQ(test_diagethr, 0.0);
    test_diagethr = hs_f.reset_diagethr(temp_ofs, 0.0, 0.0);
	EXPECT_EQ(hs_f.diag_ethr, 0.0);
	EXPECT_EQ(test_diagethr, 0.0);
    test_diagethr = hs_f.cal_hsolerror();
	EXPECT_EQ(hs_f.diag_ethr, 0.0);
	EXPECT_EQ(test_diagethr, 0.0);

	double test_diagethr_d = hs_d.set_diagethr(0, 0, 0.0);
	EXPECT_EQ(hs_d.diag_ethr, 0.0);
	EXPECT_EQ(test_diagethr_d, 0.0);
    test_diagethr_d = hs_d.reset_diagethr(temp_ofs, 0.0, 0.0);
	EXPECT_EQ(hs_d.diag_ethr, 0.0);
	EXPECT_EQ(test_diagethr_d, 0.0);
    test_diagethr_d = hs_d.cal_hsolerror();
	EXPECT_EQ(hs_d.diag_ethr, 0.0);
	EXPECT_EQ(test_diagethr_d, 0.0);
}
namespace hsolver
{
	template <typename FPTYPE, typename Device = psi::DEVICE_CPU>
	class DiagH_mock : public DiagH<FPTYPE, Device>
	{
		public:
		DiagH_mock(){}
		~DiagH_mock(){}

		void diag(hamilt::Hamilt<FPTYPE, Device> *phm_in, psi::Psi<std::complex<FPTYPE>, Device> &psi, FPTYPE *eigenvalue_in)
		{
			return;
		}
	};
	template class DiagH_mock<double>;
	template class DiagH_mock<float>;
}

TEST_F(TestHSolver, diagh)
{
	this->hs_f.pdiagh = new hsolver::DiagH_mock<float>;
	//test DiagH::diag
	this->hs_f.pdiagh->diag(nullptr, this->psi_test_f, nullptr);
	EXPECT_EQ(this->hs_f.pdiagh->method, "none");
	this->hs_d.pdiagh = new hsolver::DiagH_mock<double>;
	//test DiagH::diag
	this->hs_d.pdiagh->diag(nullptr, this->psi_test_d, nullptr);
	EXPECT_EQ(this->hs_d.pdiagh->method, "none");
	//test HSolver::~HSolver() it will delete pdiagh
}