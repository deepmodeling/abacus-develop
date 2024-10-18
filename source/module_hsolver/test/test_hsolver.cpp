#include <gtest/gtest.h>
#include <iostream>
#include <vector>
#define protected public

#include "module_hsolver/hsolver.h"
#include "hsolver_supplementary_mock.h"

#include <module_base/macros.h>

// template class hsolver::HSolver<std::complex<float>, base_device::DEVICE_CPU>;
// template class hsolver::HSolver<std::complex<double>, base_device::DEVICE_CPU>;

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
 *  - hsolver::DiagH (for cases below)
 *      - diag() for Psi(FPTYPE) case
 *      - destructor of DiagH and HSolver
 *
 * the definition of supplementary functions is added in hsolver_supplementary_mock.h 
 */

class TestHSolver : public ::testing::Test
{
public:
  // hsolver::HSolver<std::complex<float>, base_device::DEVICE_CPU> hs_cf;
  // hsolver::HSolver<std::complex<double>, base_device::DEVICE_CPU> hs_cd;
  // hsolver::HSolver<float, base_device::DEVICE_CPU> hs_f;
  // hsolver::HSolver<double, base_device::DEVICE_CPU> hs_d;

  hamilt::Hamilt<std::complex<double>> hamilt_test_cd;
  hamilt::Hamilt<std::complex<float>> hamilt_test_cf;
  psi::Psi<std::complex<double>> psi_test_cd;
  psi::Psi<std::complex<float>> psi_test_cf;

  hamilt::Hamilt<double> hamilt_test_d;
  hamilt::Hamilt<float> hamilt_test_f;
  psi::Psi<double> psi_test_d;
  psi::Psi<float> psi_test_f;

  Stochastic_WF stowf_test;

  elecstate::ElecState elecstate_test;

  ModulePW::PW_Basis_K* wfcpw;

  std::string method_test = "none";

  std::ofstream temp_ofs;
};

// TEST_F(TestHSolver, solve)
// {
//     hs_cf.solve(&hamilt_test_cf, psi_test_cf, &elecstate_test, method_test, false);
// 	hs_f.solve(&hamilt_test_f, psi_test_f, &elecstate_test, method_test, false);
//     hs_cd.solve(&hamilt_test_cd, psi_test_cd, &elecstate_test, method_test, false);
// 	hs_d.solve(&hamilt_test_d, psi_test_d, &elecstate_test, method_test, false);
//     hs_cf.solve(&hamilt_test_cf, psi_test_cf, &elecstate_test, method_test, true);
// 	hs_f.solve(&hamilt_test_f, psi_test_f, &elecstate_test, method_test, true);
//     hs_cd.solve(&hamilt_test_cd, psi_test_cd, &elecstate_test, method_test, true);
// 	hs_d.solve(&hamilt_test_d, psi_test_d, &elecstate_test, method_test, true);
//     // hs_cf.solve(&hamilt_test_cf, psi_test_cf, &elecstate_test, wfcpw, stowf_test, 0, 0, method_test, true);
//     // hs_cd.solve(&hamilt_test_cd, psi_test_cd, &elecstate_test, wfcpw, stowf_test, 0, 0, method_test, true);
// 	EXPECT_EQ(hs_f.method, "none");
// 	EXPECT_EQ(hs_d.method, "none");
// }

// TEST_F(TestHSolver, diagethr)
// {
//     float test_diagethr = hs_f.set_diagethr(0.0, 0, 0, 0.0);
// 	EXPECT_EQ(test_diagethr, 0.0);

// 	double test_diagethr_d = hs_d.set_diagethr(0.0, 0, 0, 0.0);
// 	EXPECT_EQ(test_diagethr_d, 0.0);
// }