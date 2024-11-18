#include <algorithm>
#include <random>
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#define private public
#include "module_parameter/parameter.h"
#undef private

#include "../sc_lambda_lcao.h"
#include "module_cell/klist.h"
#include "module_hamilt_lcao/module_deltaspin/spin_constrain.h"

// mockes
K_Vectors::K_Vectors()
{
}
K_Vectors::~K_Vectors()
{
}

template <>
void hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>::init(const int ik_in)
{
}

template <>
void hamilt::OperatorLCAO<std::complex<double>, double>::init(const int ik_in)
{
}

template <>
void hamilt::OperatorLCAO<double, double>::init(const int ik_in)
{
}

template <>
void hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>::contributeHk(const int ik_in)
{
}

template <>
void hamilt::OperatorLCAO<std::complex<double>, double>::contributeHk(const int ik_in)
{
}

template <>
void hamilt::OperatorLCAO<double, double>::contributeHk(const int ik_in)
{
}
// mocks

/************************************************
 *  unit test of class OperatorScLambda
 ***********************************************/

/**
 * - Tested Functions:
 *   -
 */

class ScLambdaLCAOTest : public ::testing::Test
{
  protected:
    std::vector<ModuleBase::Vector3<double>> kvec_d;
    std::vector<int> isk;
    void SetUp()
    {
        kvec_d.push_back(ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
        kvec_d.push_back(ModuleBase::Vector3<double>(0.5, 0.5, 0.5));
        isk.push_back(0);
        isk.push_back(1);
    }
};

TEST_F(ScLambdaLCAOTest, TemplateHelpers)
{
    auto sc_lambda_op
        = hamilt::OperatorScLambda<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>(nullptr,
                                                                                                     this->kvec_d,
                                                                                                     nullptr,
                                                                                                     isk);
    auto sc_lambda_op1 = hamilt::OperatorScLambda<hamilt::OperatorLCAO<std::complex<double>, double>>(nullptr,
                                                                                                      this->kvec_d,
                                                                                                      nullptr,
                                                                                                      isk);
    auto sc_lambda_op2
        = hamilt::OperatorScLambda<hamilt::OperatorLCAO<double, double>>(nullptr, this->kvec_d, nullptr, isk);
    EXPECT_NO_THROW(sc_lambda_op.contributeHR());
    EXPECT_NO_THROW(sc_lambda_op1.contributeHR());
    EXPECT_NO_THROW(sc_lambda_op2.contributeHR());
    EXPECT_NO_THROW(sc_lambda_op2.contributeHk(0));
}
