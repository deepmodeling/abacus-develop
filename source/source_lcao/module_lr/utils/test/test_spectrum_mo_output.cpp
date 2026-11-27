#include <gtest/gtest.h>

#include "../spectrum_mo.h"
#include "source_base/constants.h"

#include <complex>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

TEST(SpectrumMoOutput, WritesLibRpaVelocityMatrixContract)
{
    const std::string filename = testing::TempDir() + "spectrum_mo_librpa_test.txt";
    std::remove(filename.c_str());

    const int nk = 1;
    const int nspin = 1;
    const int nbands = 2;
    const int nbasis = 3;
    const int ks_num = 2;
    std::vector<std::complex<double>> values(nspin * 3 * nk * ks_num * ks_num);
    for (std::size_t index = 0; index < values.size(); ++index)
    {
        values[index] = std::complex<double>(static_cast<double>(index + 1), -static_cast<double>(index + 1) / 2.0);
    }

    LR_Util::output_spectrum_mo_librpa(values, filename, nk, nspin, nbands, nbasis, ks_num);

    std::ifstream input(filename.c_str());
    ASSERT_TRUE(input.good());
    std::string line;
    int header = 0;
    ASSERT_TRUE(static_cast<bool>(std::getline(input, line)));
    std::istringstream(line) >> header;
    EXPECT_EQ(header, nk);
    ASSERT_TRUE(static_cast<bool>(std::getline(input, line)));
    std::istringstream(line) >> header;
    EXPECT_EQ(header, nspin);
    ASSERT_TRUE(static_cast<bool>(std::getline(input, line)));
    std::istringstream(line) >> header;
    EXPECT_EQ(header, nbands);
    ASSERT_TRUE(static_cast<bool>(std::getline(input, line)));
    std::istringstream(line) >> header;
    EXPECT_EQ(header, nbasis);

    const double conversion = ModuleBase::Hartree_to_eV * ModuleBase::BOHR_TO_A;
    for (int direction = 0; direction < 3; ++direction)
    {
        ASSERT_TRUE(static_cast<bool>(std::getline(input, line)));
        std::istringstream block_header(line);
        int actual_direction = 0;
        int actual_kpoint = 0;
        int actual_spin = 0;
        block_header >> actual_direction >> actual_kpoint >> actual_spin;
        EXPECT_EQ(actual_direction, direction + 1);
        EXPECT_EQ(actual_kpoint, 1);
        EXPECT_EQ(actual_spin, 1);

        const std::size_t offset = static_cast<std::size_t>(direction * ks_num * ks_num);
        for (int row = 0; row < ks_num; ++row)
        {
            for (int column = 0; column < ks_num; ++column)
            {
                ASSERT_TRUE(static_cast<bool>(std::getline(input, line)));
                std::istringstream value_line(line);
                double real = 0.0;
                double imag = 0.0;
                value_line >> real >> imag;
                const std::complex<double>& expected = values[offset + row * ks_num + column];
                EXPECT_DOUBLE_EQ(real, expected.real() * conversion);
                EXPECT_DOUBLE_EQ(imag, expected.imag() * conversion);
            }
        }
    }
    EXPECT_FALSE(static_cast<bool>(std::getline(input, line)));
    input.close();
    std::remove(filename.c_str());
}
