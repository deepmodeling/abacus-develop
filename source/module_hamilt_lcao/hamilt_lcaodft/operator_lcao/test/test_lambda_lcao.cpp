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

TYPED_TEST(LambdaLCAOTest, SetGetters)
{
    hamilt::OperatorLambda<hamilt::OperatorLCAO<TypeParam>> lambda
        = hamilt::OperatorLambda<hamilt::OperatorLCAO<TypeParam>>(nullptr, this->kvec_d, nullptr, nullptr);
    lambda.set_nat(2);
    lambda.set_nloc(200);
    lambda.set_npol(2);
    EXPECT_EQ(lambda.get_nat(), 2);
    EXPECT_EQ(lambda.get_nloc(), 200);
    EXPECT_EQ(lambda.get_npol(), 2);
    std::vector<ModuleBase::Vector3<double>> lambda_in;
    std::vector<ModuleBase::Vector3<double>> lambda_out;
    lambda_in.push_back(ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    lambda_in.push_back(ModuleBase::Vector3<double>(0.5, 0.5, 0.5));
    lambda.set_lambda(lambda_in);
    lambda_out = lambda.get_lambda();
    EXPECT_EQ(lambda_in.size(), lambda_out.size());
    for (int i = 0; i < lambda_in.size(); i++)
    {
        EXPECT_EQ(lambda_out[i].x, lambda_in[i].x);
        EXPECT_EQ(lambda_out[i].y, lambda_in[i].y);
        EXPECT_EQ(lambda_out[i].z, lambda_in[i].z);
    }
}