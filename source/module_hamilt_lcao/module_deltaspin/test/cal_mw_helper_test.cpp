#include "../spin_constrain.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

K_Vectors::K_Vectors()
{
}
K_Vectors::~K_Vectors()
{
}

/************************************************
 *  unit test of the functions in cal_mw_helper.cpp
 ***********************************************/

/**
 * Tested function:
 *  - SpinConstrain::convert
 *    - convert the data structure for calculation of mw
 *  - SpinConstrain::calculate_MW
 *    - calculate mw from AorbMulP matrix
 */

class SpinConstrainTest : public testing::Test
{
  protected:
    SpinConstrain<std::complex<double>, psi::DEVICE_CPU>& sc
        = SpinConstrain<std::complex<double>, psi::DEVICE_CPU>::getScInstance();
};

TEST_F(SpinConstrainTest, CalculateMW)
{
    std::map<int, int> atomCounts = {
        {0, 1}
    };
    std::map<int, int> orbitalCounts = {
        {0, 1}
    };
    std::map<int, std::map<int, int>> lnchiCounts = {
        {0, {{0, 1}}}
    };
    sc.set_atomCounts(atomCounts);
    sc.set_orbitalCounts(orbitalCounts);
    sc.set_lnchiCounts(lnchiCounts);
    sc.set_nspin(4);
    sc.set_npol(2);
    ModuleBase::matrix orbMulP;
    int nlocal = sc.get_nw() / 2;
    orbMulP.create(sc.get_nspin(), nlocal, true);
    EXPECT_EQ(orbMulP.nc, 1);
    EXPECT_EQ(orbMulP.nr, 4);
    orbMulP(0, 0) = 1.0;
    orbMulP(1, 0) = 2.0;
    orbMulP(2, 0) = 3.0;
    orbMulP(3, 0) = 4.0;
    std::vector<std::vector<std::vector<double>>> AorbMulP = sc.convert(orbMulP);
    EXPECT_EQ(AorbMulP.size(), 4);
    EXPECT_EQ(AorbMulP[0].size(), 1);
    EXPECT_EQ(AorbMulP[0][0].size(), 1);
    EXPECT_EQ(AorbMulP[0][0][0], 1.0);
    EXPECT_EQ(AorbMulP[1][0][0], 2.0);
    EXPECT_EQ(AorbMulP[2][0][0], 3.0);
    EXPECT_EQ(AorbMulP[3][0][0], 4.0);
    // calculate_MW
    sc.calculate_MW(AorbMulP);
    testing::internal::CaptureStdout();
    sc.print_Mi(true);
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("Total Magnetism on atom: 0  (2, 3, 4)"));
}

TEST_F(SpinConstrainTest, CollectMW)
{
    // set paraV
    Parallel_Orbitals paraV;
    int nrow = 2;
    int ncol = 2;
    std::ofstream ofs("test.log");
    paraV.set_global2local(nrow, ncol, false, ofs);
    sc.set_ParaV(&paraV);
    // Prepare the input data
    int nw = sc.get_nw();
    EXPECT_EQ(nw, 2);
    int nlocal = sc.get_nw() / 2;
    ModuleBase::matrix MecMulP(sc.get_nspin(), nlocal, true);
    ModuleBase::ComplexMatrix mud(sc.ParaV->ncol, sc.ParaV->nrow, true);
    mud(0, 0) = std::complex<double>(2.0, 0.0);
    mud(0, 1) = std::complex<double>(3.0, -2.0);
    mud(1, 0) = std::complex<double>(4.0, -3.0);
    mud(1, 1) = std::complex<double>(5.0, 0.0);
    // call the function
    sc.collect_MW(MecMulP, mud, nw, 0);
    // Check the expected results
    ModuleBase::matrix expected_MecMulP(4, 1);
    expected_MecMulP(0, 0) = 7.0;
    expected_MecMulP(1, 0) = 7.0;
    expected_MecMulP(2, 0) = 1.0;
    expected_MecMulP(3, 0) = -3.0;
    // Compare the matrices
    for (size_t i = 0; i < expected_MecMulP.nr; ++i)
    {
        for (size_t j = 0; j < expected_MecMulP.nc; ++j)
        {
            EXPECT_DOUBLE_EQ(MecMulP(i, j), expected_MecMulP(i, j));
        }
    }
    remove("test.log");
}

TEST_F(SpinConstrainTest, CalculateMWS2)
{
    std::map<int, int> atomCounts = {
        {0, 1}
    };
    std::map<int, int> orbitalCounts = {
        {0, 1}
    };
    std::map<int, std::map<int, int>> lnchiCounts = {
        {0, {{0, 1}}}
    };
    sc.set_atomCounts(atomCounts);
    sc.set_orbitalCounts(orbitalCounts);
    sc.set_lnchiCounts(lnchiCounts);
    sc.set_nspin(2);
    sc.set_npol(1);
    ModuleBase::matrix orbMulP;
    int nlocal = sc.get_nw();
    orbMulP.create(sc.get_nspin(), nlocal, true);
    EXPECT_EQ(orbMulP.nc, 1);
    EXPECT_EQ(orbMulP.nr, 2);
    orbMulP(0, 0) = 1.0;
    orbMulP(1, 0) = 2.0;
    std::vector<std::vector<std::vector<double>>> AorbMulP = sc.convert(orbMulP);
    EXPECT_EQ(AorbMulP.size(), 2);
    EXPECT_EQ(AorbMulP[0].size(), 1);
    EXPECT_EQ(AorbMulP[0][0].size(), 1);
    EXPECT_EQ(AorbMulP[0][0][0], 1.0);
    EXPECT_EQ(AorbMulP[1][0][0], 2.0);
    // calculate_MW
    sc.calculate_MW(AorbMulP);
    testing::internal::CaptureStdout();
    sc.print_Mi(true);
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("Total Magnetism on atom: 0  (-1)"));
}

TEST_F(SpinConstrainTest, CollectMWS2)
{
    // set paraV
    Parallel_Orbitals paraV;
    int nrow = 1;
    int ncol = 1;
    std::ofstream ofs("test.log");
    paraV.set_global2local(nrow, ncol, false, ofs);
    sc.set_ParaV(&paraV);
    // Prepare the input data
    int nw = sc.get_nw();
    EXPECT_EQ(nw, 1);
    int nlocal = sc.get_nw();
    ModuleBase::matrix MecMulP(sc.get_nspin(), nlocal, true);
    ModuleBase::ComplexMatrix mud(sc.ParaV->ncol, sc.ParaV->nrow, true);
    mud(0, 0) = std::complex<double>(2.0, 0.0);
    // call the function
    sc.collect_MW(MecMulP, mud, nw, 0);
    sc.collect_MW(MecMulP, mud, nw, 1);
    // Check the expected results
    ModuleBase::matrix expected_MecMulP(2, 1);
    expected_MecMulP(0, 0) = 2.0;
    expected_MecMulP(1, 0) = 2.0;
    // Compare the matrices
    for (size_t i = 0; i < expected_MecMulP.nr; ++i)
    {
        for (size_t j = 0; j < expected_MecMulP.nc; ++j)
        {
            EXPECT_DOUBLE_EQ(MecMulP(i, j), expected_MecMulP(i, j));
        }
    }
    remove("test.log");
}