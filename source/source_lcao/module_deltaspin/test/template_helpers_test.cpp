#include <algorithm>
#include <string>

#include "source_lcao/module_deltaspin/lambda_loop_helper.h"
#include "source_lcao/module_deltaspin/spin_constrain.h"
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
#include "source_cell/klist.h"

class SpinConstrainTest : public testing::Test
{
  protected:
    spinconstrain::SpinConstrain<double>& sc
        = spinconstrain::SpinConstrain<double>::getScInstance();
};

TEST_F(SpinConstrainTest, TemplatHelpers)
{
    // this is a trivial test as the double version is not used
    std::vector<std::complex<double>> Sloc2;
    EXPECT_NO_THROW(sc.cal_mw_from_lambda(0));
    EXPECT_NO_THROW(sc.cal_mi_lcao(0,false));
    EXPECT_NO_THROW(sc.run_lambda_loop(0, true, std::cout));
    EXPECT_FALSE(spinconstrain::check_rms_stop(sc, 0, 0, 0.0, 0.0, 0.0, std::cout));
    EXPECT_NO_THROW(spinconstrain::print_termination(sc, std::cout));
    EXPECT_NO_THROW(spinconstrain::print_header(sc, std::cout));
    std::vector<ModuleBase::Vector3<double>> new_spin, old_spin, new_delta_lambda, old_delta_lambda;
    EXPECT_FALSE(spinconstrain::check_gradient_decay(sc, new_spin, old_spin, new_delta_lambda, old_delta_lambda, true, std::cout));
    double alpha = 0.0;
    EXPECT_NO_THROW(spinconstrain::check_restriction(sc, new_spin, alpha, std::cout));
    EXPECT_EQ(spinconstrain::cal_alpha_opt(sc, new_spin, old_spin, alpha), 0.0);
}