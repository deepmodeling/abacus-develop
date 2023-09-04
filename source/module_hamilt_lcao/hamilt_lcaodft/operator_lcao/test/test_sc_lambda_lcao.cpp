#include <algorithm>
#include <random>

#include "../sc_lambda_lcao.h"
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
void hamilt::OperatorLCAO<std::complex<double>>::init(const int ik_in)
{
}

class ScLambdaLCAOTest : public ::testing::Test
{
  protected:
    std::vector<ModuleBase::Vector3<double>> kvec_d;
    void SetUp()
    {
        kvec_d.push_back(ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
        kvec_d.push_back(ModuleBase::Vector3<double>(0.5, 0.5, 0.5));
    }
};

TEST_F(ScLambdaLCAOTest, SetGetters)
{
    hamilt::OperatorScLambda<hamilt::OperatorLCAO<std::complex<double>>> sc_lambda_op
        = hamilt::OperatorScLambda<hamilt::OperatorLCAO<std::complex<double>>>(nullptr, this->kvec_d, nullptr, nullptr);
    sc_lambda_op.set_nat(2);
    sc_lambda_op.set_nloc(200*2);
    EXPECT_EQ(sc_lambda_op.get_nat(), 2);
    EXPECT_EQ(sc_lambda_op.get_nloc(), 200*2);
    // set and get lambda values
    std::vector<ModuleBase::Vector3<double>> lambda_in;
    std::vector<ModuleBase::Vector3<double>> lambda_out;
    lambda_in.push_back(ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    lambda_in.push_back(ModuleBase::Vector3<double>(0.5, 0.5, 0.5));
    sc_lambda_op.set_lambda(lambda_in);
    lambda_out = sc_lambda_op.get_lambda();
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
    sc_lambda_op.set_iwt2iat(iwt2iat_in);
    std::vector<int> iwt2iat_out = sc_lambda_op.get_iwt2iat();
    for (int i = 0; i < 200; i++)
    {
        EXPECT_EQ(iwt2iat_out[i], 0);
        EXPECT_EQ(iwt2iat_out[i + 200], 1);
    }
    delete[] iwt2iat_in;
}

TEST_F(ScLambdaLCAOTest, CalWeightFunc)
{
    hamilt::OperatorScLambda<hamilt::OperatorLCAO<std::complex<double>>> sc_lambda_op
        = hamilt::OperatorScLambda<hamilt::OperatorLCAO<std::complex<double>>>(nullptr, this->kvec_d, nullptr, nullptr);
    sc_lambda_op.set_nat(2);
    sc_lambda_op.set_nloc(200*2);
    int nat = sc_lambda_op.get_nat();
    int nloc = sc_lambda_op.get_nloc();
    // set and get iwt2iat values
    int* iwt2iat_in = new int[400];
    for (int i = 0; i < 200; i++)
    {
        iwt2iat_in[i] = 0;
        iwt2iat_in[i + 200] = 1;
    }
    sc_lambda_op.set_iwt2iat(iwt2iat_in);
    std::vector<int> iwt2iat_out = sc_lambda_op.get_iwt2iat();
    // Generate sample Sloc2 data
    std::vector<std::complex<double>> Sloc2(nloc * nloc);
    std::generate(Sloc2.begin(), Sloc2.end(), [] { return static_cast<std::complex<double>>(std::rand() / RAND_MAX); });

    // Run the cal_weight_func method
    sc_lambda_op.cal_weight_func(Sloc2);

    // Verify the results
    for (int i = 0; i < nloc; ++i)
    {
        int iat = iwt2iat_out[i];
        for (int j = i; j < nloc; ++j)
        {
            int jat = iwt2iat_out[j];
            std::complex<double> expected = (iat == jat) ? Sloc2[i * nloc + j] : Sloc2[i * nloc + j] * 0.5;
            std::complex<double> actual = sc_lambda_op.W_i_[i * nloc + j];
            EXPECT_DOUBLE_EQ(expected.real(), actual.real());
            EXPECT_DOUBLE_EQ(expected.imag(), actual.imag());
        }
    }
    delete[] iwt2iat_in;
}