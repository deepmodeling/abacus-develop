#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "module_io/write_wfc_nao.h"
#include "module_base/global_variable.h"
#include "../binstream.h"

/************************************************
 *  unit test of write_wfc_nao.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - ModuleIO::write_wfc_nao()
 *   - ModuleIO::write_wfc_nao_complex()
 */

TEST(ModuleIOTest, WriteWfcNao)
{
    // Set up GlobalV
    GlobalV::DRANK = 0;
    GlobalV::NBANDS = 2;
    GlobalV::NLOCAL = 2;
    GlobalV::CURRENT_SPIN = 0;
    GlobalV::out_app_flag = 1;

    // Set up test data
    std::string filename = "test_wfc_nao.txt";
    double** ctot = new double*[2];
    ctot[0] = new double[2];
    ctot[1] = new double[2];
    ctot[0][0] = 0.1;
    ctot[0][1] = 0.2;
    ctot[1][0] = 0.3;
    ctot[1][1] = 0.4;
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
    ifs.close();

    //clean up
    delete[] ctot[0];
    delete[] ctot[1];
    delete[] ctot;
    std::remove(filename.c_str()); // remove the test file
}

TEST(ModuleIOTest, WriteWfcNaoBinary)
{
    // Set up GlobalV
    GlobalV::DRANK = 0;
    GlobalV::NBANDS = 2;
    GlobalV::NLOCAL = 2;
    GlobalV::CURRENT_SPIN = 0;
    GlobalV::out_app_flag = 1;

    // Set up test data
    std::string filename = "test_wfc_nao.dat";
    double** ctot = new double*[2];
    ctot[0] = new double[2];
    ctot[1] = new double[2];
    ctot[0][0] = 0.1;
    ctot[0][1] = 0.2;
    ctot[1][0] = 0.3;
    ctot[1][1] = 0.4;
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
    ModuleIO::write_wfc_nao(filename, ctot, ekb, wg, true);

    // Check the output file
    Binstream wfc(filename,"r");
    int nbands, nlocal;
    wfc >> nbands;
    wfc >> nlocal;
    EXPECT_EQ(nbands, 2);
    EXPECT_EQ(nlocal, 2);
    for (int i = 0; i < nbands; i++)
    {
        int band_index;
        double ekb, wg;
        wfc >> band_index;
        wfc >> ekb;
        wfc >> wg;
        EXPECT_EQ(band_index, i+1);
        EXPECT_DOUBLE_EQ(ekb, 0.5 + i*0.2);
        EXPECT_DOUBLE_EQ(wg, 0.9 + i*0.2);
        for (int j = 0; j < nlocal; j++)
        {
            double ctot;
            wfc >> ctot;
            EXPECT_DOUBLE_EQ(ctot, 0.1 + i*0.2 + j*0.1);
        }
    }
    wfc.close();

    //clean up
    delete[] ctot[0];
    delete[] ctot[1];
    delete[] ctot;
    std::remove(filename.c_str()); // remove the test file
}

TEST(ModuleIOTest, WriteWfcNaoComplex)
{
    // Set up GlobalV
    GlobalV::DRANK = 0;
    GlobalV::NBANDS = 2;
    GlobalV::NLOCAL = 3;
    GlobalV::CURRENT_SPIN = 0;
    GlobalV::out_app_flag = 1;
    //set up test data
    std::string name = "test_wfc_nao_complex.txt";
    int ik = 0;
    ModuleBase::Vector3<double> kvec_c {0.0, 0.0, 0.0};
    ModuleBase::matrix ekb(1, 2);
    ekb(0, 0) = 0.9;
    ekb(0, 1) = 1.1;
    ModuleBase::matrix wg(1, 2);
    wg(0, 0) = 0.11;
    wg(0, 1) = 0.22;
    std::complex<double> **ctot = new std::complex<double>*[2];
    ctot[0] = new std::complex<double>[3];
    ctot[1] = new std::complex<double>[3];
    ctot[0][0] = {1.0, 0.0}; ctot[0][1] = {2.0, 0.0}; ctot[0][2] = {3.0, 0.0};
    ctot[1][0] = {0.0, 1.0}; ctot[1][1] = {0.0, 2.0}; ctot[1][2] = {0.0, 3.0};

    // Call the function
    ModuleIO::write_wfc_nao_complex(name, ctot, ik, kvec_c, ekb, wg);
    // Check the output file
    std::ifstream ifs(name);
    std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("1 (index of k points)"));
    EXPECT_THAT(str, testing::HasSubstr("0 0 0"));
    EXPECT_THAT(str, testing::HasSubstr("2 (number of bands)"));
    EXPECT_THAT(str, testing::HasSubstr("3 (number of orbitals)"));
    EXPECT_THAT(str, testing::HasSubstr("1 (band)"));
    EXPECT_THAT(str, testing::HasSubstr("(Ry)"));
    EXPECT_THAT(str, testing::HasSubstr("(Occupations)"));
    EXPECT_THAT(str, testing::HasSubstr("2 (band)"));
    //ifs.close();
    // Clean up
    delete[] ctot[0];
    delete[] ctot[1];
    delete[] ctot;
    std::remove(name.c_str());
}

TEST(ModuleIOTest, WriteWfcNaoComplexBinary)
{
    // Set up GlobalV
    GlobalV::DRANK = 0;
    GlobalV::NBANDS = 2;
    GlobalV::NLOCAL = 3;
    GlobalV::CURRENT_SPIN = 0;
    GlobalV::out_app_flag = 1;
    //set up test data
    std::string name = "test_wfc_nao_complex.dat";
    int ik = 0;
    ModuleBase::Vector3<double> kvec_c {0.0, 0.0, 0.0};
    ModuleBase::matrix ekb(1, 2);
    ekb(0, 0) = 0.9;
    ekb(0, 1) = 1.1;
    ModuleBase::matrix wg(1, 2);
    wg(0, 0) = 0.11;
    wg(0, 1) = 0.22;
    std::complex<double> **ctot = new std::complex<double>*[2];
    ctot[0] = new std::complex<double>[3];
    ctot[1] = new std::complex<double>[3];
    ctot[0][0] = {1.0, 0.0}; ctot[0][1] = {2.0, 2.0}; ctot[0][2] = {3.0, 4.0};
    ctot[1][0] = {4.0, 4.0}; ctot[1][1] = {5.0, 6.0}; ctot[1][2] = {6.0, 8.0};

    // Call the function
    ModuleIO::write_wfc_nao_complex(name, ctot, ik, kvec_c, ekb, wg, true);
    // Check the output file

    Binstream wfc(name,"r");
    int ik1, nbands, nlocal;
    double kx, ky, kz;
    wfc >> ik1;
    wfc >> kx;
    wfc >> ky;
    wfc >> kz;
    wfc >> nbands;
    wfc >> nlocal;
    EXPECT_EQ(ik1, 1);
    EXPECT_DOUBLE_EQ(kx, 0.0);
    EXPECT_DOUBLE_EQ(ky, 0.0);
    EXPECT_DOUBLE_EQ(kz, 0.0);
    EXPECT_EQ(nbands, 2);
    EXPECT_EQ(nlocal, 3);
    for (int i = 0; i < nbands; i++)
    {
        int band_index;
        double ekb, wg;
        wfc >> band_index;
        wfc >> ekb;
        wfc >> wg;
        EXPECT_EQ(band_index, i+1);
        EXPECT_DOUBLE_EQ(ekb, 0.9 + i*0.2);
        EXPECT_DOUBLE_EQ(wg, 0.11 + i*0.11);
        for (int j = 0; j < nlocal; j++)
        {
            std::complex<double> ctot;
            wfc >> ctot;
            EXPECT_DOUBLE_EQ(ctot.real(), 1.0 + i*3.0 + j*1.0);
            EXPECT_DOUBLE_EQ(ctot.imag(), 0.0 + i*4.0 + j*2.0);
        }
    }
    wfc.close();
    // Clean up
    delete[] ctot[0];
    delete[] ctot[1];
    delete[] ctot;
    std::remove(name.c_str());
}
