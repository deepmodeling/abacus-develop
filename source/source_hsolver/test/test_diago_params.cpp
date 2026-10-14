#include "../diago_iter_assist.h"
#include "../diago_params.h"

#include <complex>

#include <gtest/gtest.h>

template <typename T>
class DiagoParamsPWTest : public ::testing::Test
{
  protected:
    using Assist = hsolver::DiagoIterAssist<T, base_device::DEVICE_CPU>;
    using Real = typename GetTypeReal<T>::type;

    void SetUp() override
    {
        saved_nmax = Assist::PW_DIAG_NMAX;
        saved_thr = Assist::PW_DIAG_THR;
        saved_scf_iter = Assist::SCF_ITER;
        saved_need_subspace = Assist::need_subspace;
    }

    void TearDown() override
    {
        Assist::PW_DIAG_NMAX = saved_nmax;
        Assist::PW_DIAG_THR = saved_thr;
        Assist::SCF_ITER = saved_scf_iter;
        Assist::need_subspace = saved_need_subspace;
    }

    void setup(const Input_para& input)
    {
        hsolver::setup_diago_params_pw<T, base_device::DEVICE_CPU>(0, 1, 1.0e-8, input);
    }

    int saved_nmax;
    Real saved_thr;
    int saved_scf_iter;
    bool saved_need_subspace;
};

using DiagoParamsTypes = ::testing::Types<std::complex<float>, std::complex<double>, double>;
TYPED_TEST_SUITE(DiagoParamsPWTest, DiagoParamsTypes);

TYPED_TEST(DiagoParamsPWTest, AppliesPwDiagNmaxForScf)
{
    using Assist = typename TestFixture::Assist;
    Input_para input;
    input.calculation = "scf";
    input.pw_diag_nmax = 73;

    const char* solvers[] = {"cg", "dav", "dav_subspace", "ppcg"};
    for (const char* solver : solvers)
    {
        input.ks_solver = solver;
        Assist::PW_DIAG_NMAX = 19;
        this->setup(input);
        EXPECT_EQ(Assist::PW_DIAG_NMAX, 73) << solver;
    }
}

TYPED_TEST(DiagoParamsPWTest, AppliesPwDiagNmaxForNscf)
{
    using Assist = typename TestFixture::Assist;
    Input_para input;
    input.calculation = "nscf";
    input.pw_diag_nmax = 73;

    // These solvers consume the shared limit; BPCG has its own update budget.
    const char* solvers[] = {"cg", "dav", "dav_subspace", "ppcg"};
    for (const char* solver : solvers)
    {
        input.ks_solver = solver;
        Assist::PW_DIAG_NMAX = 19;
        this->setup(input);
        EXPECT_EQ(Assist::PW_DIAG_NMAX, 73) << solver;
    }
}

TYPED_TEST(DiagoParamsPWTest, NscfBudgetDoesNotDependOnPreviousCalculation)
{
    using Assist = typename TestFixture::Assist;
    Input_para input;

    const char* solvers[] = {"cg", "dav", "dav_subspace", "ppcg"};
    for (const char* solver : solvers)
    {
        input.ks_solver = solver;
        input.calculation = "scf";
        input.pw_diag_nmax = 91;
        this->setup(input);
        EXPECT_EQ(Assist::PW_DIAG_NMAX, 91) << solver;

        input.calculation = "nscf";
        input.pw_diag_nmax = 7;
        this->setup(input);
        EXPECT_EQ(Assist::PW_DIAG_NMAX, 7) << solver;

        input.pw_diag_nmax = 200;
        this->setup(input);
        EXPECT_EQ(Assist::PW_DIAG_NMAX, 200) << solver;

        Input_para default_input;
        default_input.calculation = "nscf";
        default_input.ks_solver = solver;
        this->setup(default_input);
        EXPECT_EQ(Assist::PW_DIAG_NMAX, 50) << solver;
    }
}
