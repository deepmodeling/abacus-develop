#include <algorithm>
#include <random>

#include "../sc_lambda_lcao.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "module_cell/klist.h"

/************************************************
 *  unit test of class timer
 ***********************************************/

/**
 * - Tested Functions:
 *   -
 */
K_Vectors::K_Vectors(){}
K_Vectors::~K_Vectors(){}


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
    EXPECT_EQ(1,1);
}

/*
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

    std::vector<std::complex<double>> Wi_out = sc_lambda_op.get_Wi();

    // Verify the results
    for (int i = 0; i < nloc; ++i)
    {
        int iat = iwt2iat_out[i];
        for (int j = i; j < nloc; ++j)
        {
            int jat = iwt2iat_out[j];
            std::complex<double> expected = (iat == jat) ? Sloc2[i * nloc + j] : Sloc2[i * nloc + j] * 0.5;
            std::complex<double> actual = Wi_out[i * nloc + j];
            EXPECT_DOUBLE_EQ(expected.real(), actual.real());
            EXPECT_DOUBLE_EQ(expected.imag(), actual.imag());
        }
    }
    delete[] iwt2iat_in;
}

TEST_F(ScLambdaLCAOTest, CalHLambda)
{
    hamilt::OperatorScLambda<hamilt::OperatorLCAO<std::complex<double>>> sc_lambda_op
        = hamilt::OperatorScLambda<hamilt::OperatorLCAO<std::complex<double>>>(nullptr, this->kvec_d, nullptr, nullptr);
    const int nat = 3;
    const int nloc = 3*2;
    sc_lambda_op.set_nat(nat);
    sc_lambda_op.set_nloc(nloc);
    std::vector<ModuleBase::Vector3<double>> lambda;
    lambda.push_back(ModuleBase::Vector3<double>(0.1, 0.2, 0.3));
    lambda.push_back(ModuleBase::Vector3<double>(0.4, 0.5, 0.6));
    lambda.push_back(ModuleBase::Vector3<double>(0.7, 0.8, 0.9));
    sc_lambda_op.set_lambda(lambda);
    // Pauli matrices
    // sigma_x = {{0, 1}, {1, 0}}
    // sigma_y = {{0, -i}, {i, 0}}
    // sigma_z = {{1, 0}, {0, -1}}
    // lambda_x * sigma_x + lambda_y * sigma_y + lambda_z * sigma_z
    // = {{lambda_z, lambda_x - i * lambda_y}, {lambda_x + i * lambda_y, -lambda_z}}
    // let sm = lambda_i * sigma
    // then
    // sm_1 = {{0.3, 0.1-0.2i}, {0.1+0.2i, -0.3}}
    // sm_2 = {{0.6, 0.4-0.5i}, {0.4+0.5i, -0.6}}
    // sm_3 = {{0.9, 0.7-0.8i}, {0.7+0.8i, -0.9}}
    const std::vector<int> iwt2iat = {0, 0, 1, 1, 2, 2};
    sc_lambda_op.set_iwt2iat(iwt2iat.data());
    const std::vector<std::complex<double>> Wi = {
        {1.0,  0.0},
        {2.0,  0.0},
        {0.0,  0.0},
        {0.0,  0.0},
        {0.0,  0.0},
        {0.0,  0.0}, // 1st row
        {3.0,  0.0},
        {4.0,  0.0},
        {0.0,  0.0},
        {0.0,  0.0},
        {0.0,  0.0},
        {0.0,  0.0}, // 2nd row
        {0.0,  0.0},
        {0.0,  0.0},
        {5.0,  0.0},
        {6.0,  0.0},
        {0.0,  0.0},
        {0.0,  0.0}, // 3rd row
        {0.0,  0.0},
        {0.0,  0.0},
        {7.0,  0.0},
        {8.0,  0.0},
        {0.0,  0.0},
        {0.0,  0.0}, // 4th row
        {0.0,  0.0},
        {0.0,  0.0},
        {0.0,  0.0},
        {0.0,  0.0},
        {9.0,  0.0},
        {10.0, 0.0}, // 5th row
        {0.0,  0.0},
        {0.0,  0.0},
        {0.0,  0.0},
        {0.0,  0.0},
        {11.0, 0.0},
        {12.0, 0.0}  // 6th row
    };

    sc_lambda_op.set_Wi(Wi);

    // Expected values calculated manually
    std::vector<std::complex<double>> expected_h_lambda = {
        {0.3,   0.0 },
        {0.2,   -0.4},
        {0.0,   0.0 },
        {0.0,   0.0 },
        {0.0,   0.0 },
        {0.0,   0.0 }, // 1st row
        {0.3,   0.6 },
        {-1.2,  0.0 },
        {0.0,   0.0 },
        {0.0,   0.0 },
        {0.0,   0.0 },
        {0.0,   0.0 }, // 2nd row
        {0.0,   0.0 },
        {0.0,   0.0 },
        {3.0,   0.0 },
        {2.4,   -3.0},
        {0.0,   0.0 },
        {0.0,   0.0 }, // 3rd row
        {0.0,   0.0 },
        {0.0,   0.0 },
        {2.8,   3.5 },
        {-4.8,  0.0 },
        {0.0,   0.0 },
        {0.0,   0.0 }, // 4th row
        {0.0,   0.0 },
        {0.0,   0.0 },
        {0.0,   0.0 },
        {0.0,   0.0 },
        {8.1,   0.0 },
        {7.0,   -8.0}, // 5th row
        {0.0,   0.0 },
        {0.0,   0.0 },
        {0.0,   0.0 },
        {0.0,   0.0 },
        {7.7,   8.8 },
        {-10.8, 0.0 }  // 6th row
    };
    const int ik = 0; // Not used in the function
    std::vector<std::complex<double>> h_lambda(nloc * nloc);
    sc_lambda_op.cal_h_lambda(ik, h_lambda.data());
    for (size_t i = 0; i < h_lambda.size(); ++i)
    {
        EXPECT_NEAR(h_lambda[i].real(), expected_h_lambda[i].real(), 1e-9);
        EXPECT_NEAR(h_lambda[i].imag(), expected_h_lambda[i].imag(), 1e-9);
    }
}
*/