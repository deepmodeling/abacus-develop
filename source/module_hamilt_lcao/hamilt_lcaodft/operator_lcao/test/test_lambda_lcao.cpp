#include <algorithm>
#include <random>

#include "../lambda_lcao.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

template <typename T>
T get_real_value(const T& val)
{
    return val;
}

template <typename T>
T get_real_value(const std::complex<T>& val)
{
    return val.real();
}

template <typename T>
T get_imag_value(const T& val)
{
    return val;
}

template <typename T>
T get_imag_value(const std::complex<T>& val)
{
    return val.imag();
}

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
    hamilt::OperatorLambda<hamilt::OperatorLCAO<TypeParam>> lambda_op
        = hamilt::OperatorLambda<hamilt::OperatorLCAO<TypeParam>>(nullptr, this->kvec_d, nullptr, nullptr);
    lambda_op.set_nat(2);
    lambda_op.set_nloc(200);
    lambda_op.set_npol(2);
    EXPECT_EQ(lambda_op.get_nat(), 2);
    EXPECT_EQ(lambda_op.get_nloc(), 200);
    EXPECT_EQ(lambda_op.get_npol(), 2);
    // set and get lambda values
    std::vector<ModuleBase::Vector3<double>> lambda_in;
    std::vector<ModuleBase::Vector3<double>> lambda_out;
    lambda_in.push_back(ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    lambda_in.push_back(ModuleBase::Vector3<double>(0.5, 0.5, 0.5));
    lambda_op.set_lambda(lambda_in);
    lambda_out = lambda_op.get_lambda();
    EXPECT_EQ(lambda_in.size(), lambda_out.size());
    for (int i = 0; i < lambda_in.size(); i++)
    {
        EXPECT_EQ(lambda_out[i].x, lambda_in[i].x);
        EXPECT_EQ(lambda_out[i].y, lambda_in[i].y);
        EXPECT_EQ(lambda_out[i].z, lambda_in[i].z);
    }
    // set and get iwt2iat values
    int* iwt2iat_in = new int[400];
    for (int i = 0; i < 200; i++)
    {
        iwt2iat_in[i] = 0;
        iwt2iat_in[i + 200] = 1;
    }
    lambda_op.set_iwt2iat(iwt2iat_in);
    std::vector<int> iwt2iat_out = lambda_op.get_iwt2iat();
    for (int i = 0; i < 200; i++)
    {
        EXPECT_EQ(iwt2iat_out[i], 0);
        EXPECT_EQ(iwt2iat_out[i + 200], 1);
    }
    delete[] iwt2iat_in;
}

TYPED_TEST(LambdaLCAOTest, CalWeightFunc)
{
    hamilt::OperatorLambda<hamilt::OperatorLCAO<TypeParam>> lambda_op
        = hamilt::OperatorLambda<hamilt::OperatorLCAO<TypeParam>>(nullptr, this->kvec_d, nullptr, nullptr);
    lambda_op.set_nat(2);
    lambda_op.set_nloc(200);
    lambda_op.set_npol(2);
    int nat = lambda_op.get_nat();
    int nloc = lambda_op.get_nloc();
    int npol = lambda_op.get_npol();
    // set and get iwt2iat values
    int* iwt2iat_in = new int[400];
    for (int i = 0; i < 200; i++)
    {
        iwt2iat_in[i] = 0;
        iwt2iat_in[i + 200] = 1;
    }
    lambda_op.set_iwt2iat(iwt2iat_in);
    std::vector<int> iwt2iat_out = lambda_op.get_iwt2iat();
    // Generate sample Sloc2 data
    std::vector<TypeParam> Sloc2(nloc * npol * nloc * npol);
    std::generate(Sloc2.begin(), Sloc2.end(), [] { return static_cast<TypeParam>(std::rand() / RAND_MAX); });

    // Run the cal_weight_func method
    lambda_op.cal_weight_func(Sloc2);

    // Verify the results
    for (int i = 0; i < nloc * npol; ++i)
    {
        int iat = iwt2iat_out[i];
        for (int j = i; j < nloc * npol; ++j)
        {
            int jat = iwt2iat_out[j];
            TypeParam expected = (iat == jat) ? Sloc2[i * nloc * npol + j] : Sloc2[i * nloc * npol + j] * 0.5;
            TypeParam actual = lambda_op.W_i_[i * nloc * npol + j];
            EXPECT_DOUBLE_EQ(get_real_value(expected), get_real_value(actual));
            EXPECT_DOUBLE_EQ(get_imag_value(expected), get_imag_value(actual));
        }
    }
}