// upsi_test.cpp

#include "../upsi.h"

#include <gtest/gtest.h>

TEST(UpsiTest, ScalarProduct)
{
    // Set up test inputs
    int nband = 2;
    int nlocal = 3;
    std::complex<double> U_operator[9] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
    std::complex<double> psi_k_laststep[6] = {7.0, 8.0, 9.0, 10.0, 11.0, 12.0};
    std::complex<double> expected_result[4] = {58.0, 64.0, 139.0, 154.0};
    std::complex<double> actual_result[4];

    Parallel_Orbitals pv;
    pv.desc[0] = 1; // Set dummy values for the parallelization descriptors
    pv.desc_wfc[0] = 1;
    pv.ncol_bands = nband;
    pv.ncol = nlocal;

    // Call the function being tested
    upsi(&pv, nband, nlocal, U_operator, psi_k_laststep, actual_result, false);

    // Check the result
    for (int i = 0; i < 6; i++)
    {
        EXPECT_EQ(actual_result[i], expected_result[i]);
    }
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
