#include "../diago_iter_assist.h"
#include "../diago_params.h"

#include <complex>

#include <gtest/gtest.h>

TEST(DiagoParamsPWTest, AppliesPwDiagNmaxForScf)
{
    Input_para input;
    input.calculation = "scf";
    input.ks_solver = "ppcg";
    input.pw_diag_nmax = 73;

    typedef hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_CPU> Assist;
    const int saved_nmax = Assist::PW_DIAG_NMAX;
    Assist::PW_DIAG_NMAX = 19;

    hsolver::setup_diago_params_pw<std::complex<double>, base_device::DEVICE_CPU>(
        0,
        1,
        1.0e-8,
        input);

    EXPECT_EQ(Assist::PW_DIAG_NMAX, 73);
    Assist::PW_DIAG_NMAX = saved_nmax;
}

TEST(DiagoParamsPWTest, PreservesSharedNmaxForNscf)
{
    Input_para input;
    input.calculation = "nscf";
    input.pw_diag_nmax = 73;

    typedef hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_CPU> Assist;
    const int saved_nmax = Assist::PW_DIAG_NMAX;
    Assist::PW_DIAG_NMAX = 19;

    const char* solvers[] = {"cg", "dav", "dav_subspace", "bpcg", "ppcg"};
    for (const char* solver : solvers)
    {
        input.ks_solver = solver;
        hsolver::setup_diago_params_pw<std::complex<double>, base_device::DEVICE_CPU>(
            0,
            1,
            1.0e-8,
            input);
        EXPECT_EQ(Assist::PW_DIAG_NMAX, 19) << solver;
    }
    Assist::PW_DIAG_NMAX = saved_nmax;
}
