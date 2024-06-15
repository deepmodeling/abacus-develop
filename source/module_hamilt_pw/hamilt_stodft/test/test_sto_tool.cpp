#include "../sto_tool.h"
#include "mpi.h"

#include <gtest/gtest.h>

/************************************************
 *  unit test of sto_tool.cpp
 ***********************************************/

void Stochastic_hchi::hchi_norm(std::complex<double>* chig, std::complex<double>* hchig, const int m)
{
    return;
}

/**
 * - Tested Functions:
 *   - struct parallel_distribution
 *   - void convert_psi(psi_in, psi_out)
 *   - psi::Psi<std::complex<float>>* gatherchi(chi, chi_all, npwx, nrecv_sto, displs_sto, perbands_sto)
 */
class TestStoTool : public ::testing::Test
{
};

TEST_F(TestStoTool, parallel_distribution)
{
    int num_all = 10;
    int np = 4;
    int myrank = 0;
    parallel_distribution pd(num_all, np, myrank);
    EXPECT_EQ(pd.start, 0);
    EXPECT_EQ(pd.num_per, 3);

    myrank = 1;
    parallel_distribution pd1(num_all, np, myrank);
    EXPECT_EQ(pd1.start, 3);
    EXPECT_EQ(pd1.num_per, 3);

    myrank = 2;
    parallel_distribution pd2(num_all, np, myrank);
    EXPECT_EQ(pd2.start, 6);
    EXPECT_EQ(pd2.num_per, 2);

    myrank = 3;
    parallel_distribution pd3(num_all, np, myrank);
    EXPECT_EQ(pd3.start, 8);
    EXPECT_EQ(pd3.num_per, 2);
}

TEST_F(TestStoTool, convert_psi)
{
    psi::Psi<std::complex<double>> psi_in(1, 1, 10);
    psi::Psi<std::complex<float>> psi_out(1, 1, 10);
    for (int i = 0; i < 10; ++i)
    {
        psi_in.get_pointer()[i] = std::complex<double>(i, i);
    }
    convert_psi(psi_in, psi_out);
    for (int i = 0; i < 10; ++i)
    {
        EXPECT_EQ(psi_out.get_pointer()[i], std::complex<float>(i, i));
    }
}

TEST_F(TestStoTool, gatherchi)
{
    psi::Psi<std::complex<float>> chi(1, 1, 10);
    psi::Psi<std::complex<float>> chi_all(1, 1, 10);
    int npwx = 10;
    int nrecv_sto[4] = {1, 2, 3, 4};
    int displs_sto[4] = {0, 1, 3, 6};
    int perbands_sto = 1;
    psi::Psi<std::complex<float>>* p_chi = gatherchi(chi, chi_all, npwx, nrecv_sto, displs_sto, perbands_sto);
    EXPECT_EQ(p_chi, &chi);
}
