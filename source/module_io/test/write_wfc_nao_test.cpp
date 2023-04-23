#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "module_io/write_wfc_nao.h"
#include "module_base/global_variable.h"

/************************************************
 *  unit test of write_wfc_nao.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - ModuleIO::write_wfc_nao()
 */

TEST(ModuleIOTest, WriteWfcNaoTest)
{
    // Set up GlobalV
    GlobalV::DRANK = 0;
    GlobalV::NBANDS = 2;
    GlobalV::NLOCAL = 2;
    GlobalV::CURRENT_SPIN = 0;
    // Set up test data
    std::string filename = "test_wfc_nao.txt";
    double** ctot;
    ctot = new double*[2];
    for(int i=0; i<2; i++)
    {
	    ctot[i] = new double[2];
    }
    ctot[0][0] = 0.1;
    ctot[0][1] = 0.2;
    ctot[1][0] = 0.3;
    ctot[1][1] = 0.4;
    //
    ModuleBase::matrix ekb(2, 2);
    ekb(0, 0) = 0.5;
    ekb(1, 0) = 0.6;
    ekb(0, 1) = 0.7;
    ekb(1, 1) = 0.8;
    ModuleBase::matrix wg(2, 2);
    wg(0, 0) = 0.9;
    wg(1, 0) = 1.0;
    wg(0, 1) = 1.1;
    wg(1, 1) = 1.2;

    // Call the function to be tested
    ModuleIO::write_wfc_nao(filename, ctot, ekb, wg);

    for(int i=0; i<2; i++)
    {
	    delete[] ctot[i];
    }
    delete[] ctot;

    // Check the output file
    std::ifstream ifs(filename);
    std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("2 (number of bands)"));
    EXPECT_THAT(str, testing::HasSubstr("2 (number of orbitals)"));
    EXPECT_THAT(str, testing::HasSubstr("1 (band)"));
    EXPECT_THAT(str, testing::HasSubstr("5.00000000e-01 (Ry)"));
    EXPECT_THAT(str, testing::HasSubstr("9.00000000e-01 (Occupations)"));
    EXPECT_THAT(str, testing::HasSubstr("1.00000000e-01 2.00000000e-01"));
    EXPECT_THAT(str, testing::HasSubstr("2 (band)"));
    EXPECT_THAT(str, testing::HasSubstr("7.00000000e-01 (Ry)"));
    EXPECT_THAT(str, testing::HasSubstr("1.10000000e+00 (Occupations)"));
    EXPECT_THAT(str, testing::HasSubstr("3.00000000e-01 4.00000000e-01"));
    std::remove(filename.c_str()); // remove the test file
}
