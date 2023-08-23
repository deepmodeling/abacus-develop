#include "../lambda_lcao.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

/************************************************
 *  unit test of class timer
 ***********************************************/

/**
 * - Tested Functions:
 *   -
 */

template<>
void hamilt::OperatorLCAO<double>::init(const int ik_in)
{
}

template<>
void hamilt::OperatorLCAO<std::complex<double>>::init(const int ik_in)
{
}

template <typename T>
class LambdaLCAOTest : public ::testing::Test
{
  protected:
    std::vector<ModuleBase::Vector3<double>> kvec_d;
    void SetUp()
    {
        kvec_d.push_back(ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
        kvec_d.push_back(ModuleBase::Vector3<double>(0.5, 0.5, 0.5));
    }
};

using MyTypes = ::testing::Types<double, std::complex<double>>;
TYPED_TEST_SUITE(LambdaLCAOTest, MyTypes);

TYPED_TEST(LambdaLCAOTest, SimpleSetters)
{
    hamilt::OperatorLambda<hamilt::OperatorLCAO<TypeParam>> lambda
        = hamilt::OperatorLambda<hamilt::OperatorLCAO<TypeParam>>(nullptr, this->kvec_d, nullptr, nullptr);
    lambda.set_nat(10);
    lambda.set_nloc(200);
    lambda.set_npol(2);
    EXPECT_EQ(lambda.get_nat(), 10);
    EXPECT_EQ(lambda.get_nloc(), 200);
    EXPECT_EQ(lambda.get_npol(), 2);
}