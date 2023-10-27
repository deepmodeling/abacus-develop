#include <algorithm>
#include <string>

#include "../spin_constrain.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

/************************************************
 *  unit test of functions in template_helpers.cpp
 ***********************************************/

/**
 * - Tested functions:
 *  - Functions in template_helpers.cpp are defined by template specialization.
 *    but they are not used in the code.
 *    So, we just test if they can be called without error.
 */

K_Vectors::K_Vectors()
{
}
K_Vectors::~K_Vectors()
{
}

class SpinConstrainTest : public testing::Test
{
  protected:
    SpinConstrain<double, psi::DEVICE_CPU>& sc = SpinConstrain<double, psi::DEVICE_CPU>::getScInstance();
};

TEST_F(SpinConstrainTest, TemplatHelpers)
{
    // this is a trivial test as the double version is not used
    std::vector<std::complex<double>> Sloc2;
    EXPECT_NO_THROW(sc.cal_h_lambda(nullptr, Sloc2, true));
    EXPECT_NO_THROW(sc.cal_mw_from_lambda(0));
    EXPECT_NO_THROW(sc.cal_MW_k(nullptr, std::vector<std::vector<std::complex<double>>>(0)));
    EXPECT_NO_THROW(sc.cal_MW(0, nullptr, false));
    ModuleBase::matrix orbMulP;
    EXPECT_NO_THROW(sc.convert(orbMulP));
    EXPECT_NO_THROW(sc.run_lambda_loop(0));
    std::vector<std::vector<std::vector<double>>> AorbMulP;
    EXPECT_NO_THROW(sc.calculate_MW(AorbMulP));
    ModuleBase::ComplexMatrix mud;
    ModuleBase::matrix MecMulP;
    EXPECT_NO_THROW(sc.collect_MW(MecMulP, mud, 0));
    EXPECT_FALSE(sc.check_rms_stop(0, 0, 0.0));
    EXPECT_NO_THROW(sc.print_termination());
    std::vector<ModuleBase::Vector3<double>> new_spin, old_spin, new_delta_lambda, old_delta_lambda;
    EXPECT_FALSE(sc.check_gradient_decay(new_spin, old_spin, new_delta_lambda, old_delta_lambda));
}