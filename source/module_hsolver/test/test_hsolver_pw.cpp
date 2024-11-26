#include <gtest/gtest.h>
#include <iostream>
#include <vector>

#define private public
#define protected public
#include "module_parameter/parameter.h"
#include "module_hsolver/hsolver_pw.h"
#include "module_hsolver/hsolver_lcaopw.h"
#include "hsolver_supplementary_mock.h"
#include "hsolver_pw_sup.h"
#include "hsolver_supplementary_mock.h"
#include "module_base/global_variable.h"
#include "module_hsolver/hsolver_pw.h"
#undef private
#undef protected
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
 *  - 8. solve()
 *      - lcao_in_pw specific implementation
 */
class TestHSolverPW : public ::testing::Test {
  public:
    ModulePW::PW_Basis_K pwbk;
    hsolver::HSolverPW<std::complex<float>, base_device::DEVICE_CPU> hs_f
        = hsolver::HSolverPW<std::complex<float>, base_device::DEVICE_CPU>(&pwbk,
                                                                           "scf",
                                                                           "pw",
                                                                           "cg",
                                                                           false,
                                                                           PARAM.sys.use_uspp,
                                                                           PARAM.input.nspin,
                     hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_CPU>::SCF_ITER,
                     hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_CPU>::PW_DIAG_NMAX,
                     hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_CPU>::PW_DIAG_THR,
                     hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_CPU>::need_subspace);
    hsolver::HSolverPW<std::complex<double>, base_device::DEVICE_CPU> hs_d
        = hsolver::HSolverPW<std::complex<double>, base_device::DEVICE_CPU>(&pwbk,
                                                                            "scf",
                                                                            "pw",
                                                                            "cg",
                                                                            false,
                                                                            PARAM.sys.use_uspp,
                                                                            PARAM.input.nspin,
                     hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_CPU>::SCF_ITER,
                     hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_CPU>::PW_DIAG_NMAX,
                     hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_CPU>::PW_DIAG_THR,
                     hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_CPU>::need_subspace);

    hamilt::Hamilt<std::complex<double>> hamilt_test_d;
    hamilt::Hamilt<std::complex<float>> hamilt_test_f;

    psi::Psi<std::complex<double>> psi_test_cd;
    psi::Psi<std::complex<float>> psi_test_cf;

    elecstate::ElecState elecstate_test;

    std::string method_test = "cg";

    std::vector<float> ekb_f;

    std::ofstream temp_ofs;
};

// TEST_F(TestHSolverPW, solve) {
//     // initial memory and data
//     elecstate_test.ekb.create(1, 2);
//     elecstate_test.pot = new elecstate::Potential;
//     this->ekb_f.resize(2);
//     psi_test_cf.resize(1, 2, 3);
//     psi_test_cd.resize(1, 2, 3);
//     PARAM.input.nelec = 1.0;

//     // check solve()
//     EXPECT_EQ(this->hs_f.initialed_psi, false);
//     EXPECT_EQ(this->hs_d.initialed_psi, false);

//     this->hs_f.solve(&hamilt_test_f,
//                      psi_test_cf,
//                      &elecstate_test,
//                      elecstate_test.ekb.c,

//                      GlobalV::RANK_IN_POOL,
//                      GlobalV::NPROC_IN_POOL,

//                      true);
//     // EXPECT_EQ(this->hs_f.initialed_psi, true);
//     for (int i = 0; i < psi_test_cf.size(); i++) {
//         EXPECT_DOUBLE_EQ(psi_test_cf.get_pointer()[i].real(), i + 3);
//     }
//     EXPECT_DOUBLE_EQ(elecstate_test.ekb.c[0], 4.0);
//     EXPECT_DOUBLE_EQ(elecstate_test.ekb.c[1], 7.0);
//     EXPECT_DOUBLE_EQ(hsolver::DiagoIterAssist<std::complex<float>>::avg_iter,
//                      0.0);

//     this->hs_d.solve(&hamilt_test_d,
//                      psi_test_cd,
//                      &elecstate_test,
//                      elecstate_test.ekb.c,
                     
//                      GlobalV::RANK_IN_POOL,
//                      GlobalV::NPROC_IN_POOL,

//                      true);
  
//     // EXPECT_EQ(this->hs_d.initialed_psi, true);
//     EXPECT_DOUBLE_EQ(hsolver::DiagoIterAssist<std::complex<double>>::avg_iter,
//                      0.0);
//     for (int i = 0; i < psi_test_cd.size(); i++) {
//         EXPECT_DOUBLE_EQ(psi_test_cd.get_pointer()[i].real(), i + 3);
//     }
//     EXPECT_DOUBLE_EQ(elecstate_test.ekb.c[0], 4.0);
//     EXPECT_DOUBLE_EQ(elecstate_test.ekb.c[1], 7.0);

//     // // check hamiltSolvePsiK()
//     // this->hs_f.hamiltSolvePsiK(&hamilt_test_f, psi_test_cf, this->hs_f.precondition, ekb_f.data());
//     // this->hs_d.hamiltSolvePsiK(&hamilt_test_d,
//     //                            psi_test_cd,
//     //                            this->hs_f.precondition,
//     //                            elecstate_test.ekb.c);
//     // for (int i = 0; i < psi_test_cf.size(); i++) {
//     //     EXPECT_DOUBLE_EQ(psi_test_cf.get_pointer()[i].real(), i + 4);
//     // }
//     // for (int i = 0; i < psi_test_cd.size(); i++) {
//     //     EXPECT_DOUBLE_EQ(psi_test_cf.get_pointer()[i].real(), i + 4);
//     // }
//     // EXPECT_DOUBLE_EQ(ekb_f[0], 5.0);
//     // EXPECT_DOUBLE_EQ(elecstate_test.ekb.c[0], 5.0);
//     // EXPECT_DOUBLE_EQ(ekb_f[1], 8.0);
//     // EXPECT_DOUBLE_EQ(elecstate_test.ekb.c[1], 8.0);

//     // // check endDiagH()
//     // this->hs_f.initialed_psi = true;
//     // this->hs_d.initialed_psi = true;
//     // this->hs_f.endDiagh();
//     // this->hs_d.endDiagh();
//     // // will change state of initialed_psi in endDiagh
//     // EXPECT_EQ(this->hs_f.initialed_psi, true);
//     // EXPECT_EQ(this->hs_d.initialed_psi, true);

//     // // check updatePsiK()
//     // // skip initializing Psi, Psi will not change
//     // this->hs_f.updatePsiK(&hamilt_test_f, psi_test_cf, 0);
//     // this->hs_d.updatePsiK(&hamilt_test_d, psi_test_cd, 0);
//     // for (int i = 0; i < psi_test_cf.size(); i++) {
//     //     EXPECT_DOUBLE_EQ(psi_test_cf.get_pointer()[i].real(), i + 4);
//     // }
//     // for (int i = 0; i < psi_test_cd.size(); i++) {
//     //     EXPECT_DOUBLE_EQ(psi_test_cd.get_pointer()[i].real(), i + 4);
//     // }
//     // // check update_precondition()
//     // this->hs_f.update_precondition(this->hs_f.precondition,
//     //                                0,
//     //                                psi_test_cf.get_nbasis());
//     // this->hs_d.update_precondition(this->hs_d.precondition,
//     //                                0,
//     //                                psi_test_cd.get_nbasis());
//     // EXPECT_NEAR(this->hs_f.precondition[0], 2.414213657, 1e-8);
//     // EXPECT_NEAR(this->hs_f.precondition[1], 3.618033886, 1e-8);
//     // EXPECT_NEAR(this->hs_f.precondition[2], 6.236067772, 1e-8);
//     // EXPECT_NEAR(this->hs_d.precondition[0], 2.414213562, 1e-8);
//     // EXPECT_NEAR(this->hs_d.precondition[1], 3.618033989, 1e-8);
//     // EXPECT_NEAR(this->hs_d.precondition[2], 6.236067977, 1e-8);

//     // // check diago_ethr
//     // PARAM.input.init_chg = "atomic";
//     // GlobalV::PW_DIAG_THR = 1e-7;
//     // PARAM.input.calculation = "scf";
//     // float test_diagethr = hs_f.set_diagethr(hs_f.diag_ethr, 0, 1, 1.0);
//     // EXPECT_NEAR(hs_f.diag_ethr, 0.01, 1.0e-7);
//     // EXPECT_NEAR(test_diagethr, 0.01, 1.0e-7);
//     // PARAM.input.calculation = "md";
//     // PARAM.input.init_chg = "file";
//     // test_diagethr = hs_f.set_diagethr(hs_f.diag_ethr, 0, 1, 1.0);
//     // EXPECT_NEAR(test_diagethr, 1e-5, 1.0e-7);
//     // test_diagethr = hs_f.set_diagethr(hs_f.diag_ethr, 0, 2, 1.0);
//     // EXPECT_NEAR(test_diagethr, 0.01, 1.0e-7);
//     // test_diagethr = hs_f.set_diagethr(hs_f.diag_ethr, 0, 3, 1.0e-3);
//     // EXPECT_NEAR(test_diagethr, 0.0001, 1.0e-7);

//     // PARAM.input.init_chg = "atomic";
//     // GlobalV::PW_DIAG_THR = 1e-7;
//     // PARAM.input.calculation = "scf";
//     // double test_diagethr_d = hs_d.set_diagethr(hs_d.diag_ethr, 0, 1, 1.0);
//     // EXPECT_EQ(hs_d.diag_ethr, 0.01);
//     // EXPECT_EQ(test_diagethr_d, 0.01);
//     // PARAM.input.calculation = "md";
//     // PARAM.input.init_chg = "file";
//     // test_diagethr_d = hs_d.set_diagethr(hs_d.diag_ethr, 0, 1, 1.0);
//     // EXPECT_EQ(test_diagethr_d, 1e-5);
//     // test_diagethr_d = hs_d.set_diagethr(hs_d.diag_ethr, 0, 2, 1.0);
//     // EXPECT_EQ(test_diagethr_d, 0.01);
//     // test_diagethr_d = hs_d.set_diagethr(hs_d.diag_ethr, 0, 3, 1.0e-3);
//     // EXPECT_EQ(test_diagethr_d, 0.0001);
// }

TEST_F(TestHSolverPW, SolveLcaoInPW) {
    pwbk.nks = 1;
    // initial memory and data
    elecstate_test.ekb.create(1, 2);
    elecstate_test.pot = new elecstate::Potential;
    // 1 kpt, 2 bands, 3 basis
    psi_test_cf.resize(1, 2, 3);
    psi_test_cd.resize(1, 2, 3);

    psi::Psi<std::complex<double>> transform_test_cd;
    psi::Psi<std::complex<float>> transform_test_cf;
    // transform psi, the old wanf2, has 2 local basis and 3 pw basis.
    // so in principle the hcc has dimension 3*3 to diagonalize
    // 2 lowest eigenvalues will be selected and save to psi
    transform_test_cd.resize(1, 3, 3);
    transform_test_cf.resize(1, 3, 3);

    std::complex<double> psi_value_d = {0.0, 0.0};
    std::complex<float> psi_value_f = {0.0, 0.0};
    for (int iband = 0; iband < transform_test_cd.get_nbands(); iband++) {
        for (int ibasis = 0; ibasis < transform_test_cd.get_nbasis();
             ibasis++) {
            transform_test_cd
                .get_pointer()[iband * transform_test_cd.get_nbasis() + ibasis]
                = psi_value_d;
            transform_test_cf
                .get_pointer()[iband * transform_test_cf.get_nbasis() + ibasis]
                = psi_value_f;
            psi_value_d += std::complex<double>(1.0, 0.0);
            psi_value_f += std::complex<float>(1.0, 0.0);
        }
    }
    PARAM.input.nelec = 1.0;

    // check solve()
    elecstate_test.ekb.c[0] = 1.0;
    elecstate_test.ekb.c[1] = 2.0;

    hsolver::HSolverLIP<std::complex<float>> hs_f_lip
        = hsolver::HSolverLIP<std::complex<float>>(&pwbk);
    hsolver::HSolverLIP<std::complex<double>> hs_d_lip
        = hsolver::HSolverLIP<std::complex<double>>(&pwbk);
    hs_f_lip.solve(&hamilt_test_f, psi_test_cf, &elecstate_test,
        transform_test_cf, true);
    EXPECT_DOUBLE_EQ(hsolver::DiagoIterAssist<std::complex<float>>::avg_iter, 0.0);
    for (int i = 0; i < psi_test_cf.size(); i++)
    {
        EXPECT_DOUBLE_EQ(psi_test_cf.get_pointer()[i].real(), i);
    }
    EXPECT_DOUBLE_EQ(elecstate_test.ekb.c[0], 0.0);
    EXPECT_DOUBLE_EQ(elecstate_test.ekb.c[1], 0.0);

    elecstate_test.ekb.c[0] = 1.0;
    elecstate_test.ekb.c[1] = 2.0;
    hs_d_lip.solve(&hamilt_test_d, psi_test_cd, &elecstate_test, transform_test_cd, true);
    EXPECT_DOUBLE_EQ(hsolver::DiagoIterAssist<std::complex<double>>::avg_iter, 0.0);
    for (int i = 0; i < psi_test_cd.size(); i++)
    {
        EXPECT_DOUBLE_EQ(psi_test_cd.get_pointer()[i].real(), i);
    }
    EXPECT_DOUBLE_EQ(elecstate_test.ekb.c[0], 0.0);
    EXPECT_DOUBLE_EQ(elecstate_test.ekb.c[1], 0.0);
}
