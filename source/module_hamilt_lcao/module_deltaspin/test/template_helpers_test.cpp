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
UnitCell::UnitCell()
{
}
UnitCell::~UnitCell()
{
}
InfoNonlocal::InfoNonlocal()
{
}
InfoNonlocal::~InfoNonlocal()
{
}
Magnetism::Magnetism()
{
}
Magnetism::~Magnetism()
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
    EXPECT_NO_THROW(sc.cal_h_lambda(nullptr, Sloc2));
    EXPECT_NO_THROW(sc.cal_mw_from_lambda(0));
    EXPECT_NO_THROW(sc.cal_MW_k(nullptr, std::vector<std::vector<std::complex<double>>>(0)));
    UnitCell ucell;
    EXPECT_NO_THROW(sc.cal_MW(0, nullptr, ucell, false));
    ModuleBase::matrix orbMulP;
    EXPECT_NO_THROW(sc.convert(orbMulP));
    EXPECT_NO_THROW(sc.run_lambda_loop(0));
    K_Vectors kv;
    EXPECT_NO_THROW(sc.init_sc(0.0,
                               0,
                               0,
                               0.0,
                               0.0,
                               false,
                               ucell,
                               "",
                               0,
                               nullptr,
                               0,
                               kv,
                               "",
                               nullptr,
                               nullptr,
                               nullptr,
                               nullptr,
                               nullptr));
    EXPECT_NO_THROW(sc.set_input_parameters(0.0, 0, 0, 0.0, 0.0, false));
    EXPECT_NO_THROW(sc.set_ParaV(nullptr));
    EXPECT_NO_THROW(sc.set_solver_parameters(0, kv, nullptr, nullptr, nullptr, nullptr, "", nullptr));
    EXPECT_NO_THROW(sc.bcast_ScData("", 0, 0));
    std::map<int, int> atomCounts;
    std::map<int, int> orbitalCounts;
    EXPECT_NO_THROW(sc.set_orb_counts(atomCounts, orbitalCounts));
}