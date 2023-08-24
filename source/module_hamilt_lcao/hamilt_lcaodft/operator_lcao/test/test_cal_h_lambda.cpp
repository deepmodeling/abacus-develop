#include <algorithm>
#include <random>

#include "../lambda_lcao.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

template <>
void hamilt::OperatorLCAO<double>::init(const int ik_in)
{
}

template <>
void hamilt::OperatorLCAO<std::complex<double>>::init(const int ik_in)
{
}

class OperatorLambdaTest : public ::testing::Test
{
  protected:
    OperatorLambdaTest() : ::testing::Test()
    {
    }
    std::vector<ModuleBase::Vector3<double>> kvec_d = {(0.0, 0.0, 0.0), (0.5, 0.5, 0.5)};
    // We'll use a mock class to test the function
    class MockOperatorLambda : public hamilt::OperatorLambda<hamilt::OperatorLCAO<std::complex<double>>>
    {
      public:
        MockOperatorLambda()
            : hamilt::OperatorLambda<hamilt::OperatorLCAO<std::complex<double>>>(nullptr, kvec_d, nullptr, nullptr)
        {
        }
        ~MockOperatorLambda()
        {
        }
        void set_parameters(const int nat,
                            const int nloc,
                            const int npol,
                            const std::vector<ModuleBase::Vector3<double>>& lambda,
                            const std::vector<int>& iwt2iat,
                            const std::vector<std::complex<double>>& W_i)
        {
            nat_ = nat;
            nloc_ = nloc;
            npol_ = npol;
            lambda_ = lambda;
            iwt2iat_ = iwt2iat;
            W_i_ = W_i;
        }
    };
    // We'll use a mock class to test the function
    MockOperatorLambda operator_lambda_;
};

TEST_F(OperatorLambdaTest, CalHLambda)
{
    const int nat = 3;
    const int nloc = 3;
    const int npol = 2;
    std::vector<ModuleBase::Vector3<double>> lambda;
    lambda.push_back(ModuleBase::Vector3<double>(0.1, 0.2, 0.3));
    lambda.push_back(ModuleBase::Vector3<double>(0.4, 0.5, 0.6));
    lambda.push_back(ModuleBase::Vector3<double>(0.7, 0.8, 0.9));
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
    const std::vector<std::complex<double>> W_i = {
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

    operator_lambda_.set_parameters(nat, nloc, npol, lambda, iwt2iat, W_i);

    const int ik = 0; // Not used in the function
    std::vector<std::complex<double>> h_lambda(nloc * npol * nloc * npol);
    operator_lambda_.cal_h_lambda(ik, h_lambda.data());

    for (size_t i = 0; i < h_lambda.size(); ++i)
    {
        EXPECT_NEAR(h_lambda[i].real(), expected_h_lambda[i].real(), 1e-9);
        EXPECT_NEAR(h_lambda[i].imag(), expected_h_lambda[i].imag(), 1e-9);
    }
}