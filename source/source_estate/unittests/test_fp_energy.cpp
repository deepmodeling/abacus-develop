#include "../fp_energy.h"

#include <cstring>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
/************************************************
 *  unit test of fp_energy.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - fenergy::calculate_etot()
 *   - fenergy::calculate_harris()
 *   - fenergy::clear_all()
 *   - fenergy::print_all()
 *   - Efermi::set_efval()
 *   - Efermi::get_efval()
 */
class fenergy : public ::testing::Test
{
  protected:
    elecstate::fenergy f_en;
    elecstate::Efermi eferm;
};

TEST_F(fenergy, calculate_etot)
{
    f_en.eband = 1.0;
    f_en.deband = 2.0;
    f_en.calculate_etot();
    EXPECT_EQ(f_en.etot, 3.0);
}

TEST_F(fenergy, calculate_harris)
{
    f_en.eband = 1.0;
    f_en.deband_harris = 2.0;
    f_en.calculate_harris();
    EXPECT_EQ(f_en.etot_harris, 3.0);
}

TEST_F(fenergy, clear_all)
{
    f_en.eband = 1.0;
    f_en.etot = 2.0;
    f_en.clear_all();
    EXPECT_EQ(f_en.eband, 0.0);
    EXPECT_EQ(f_en.etot, 0.0);
}

TEST_F(fenergy, print_all)
{
    testing::internal::CaptureStdout();
    f_en.print_all();
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("total="));
    EXPECT_THAT(output, testing::HasSubstr("entropy(-TS)="));
}

TEST_F(fenergy, eferm_get_ef)
{
    eferm.two_efermi = false;
    eferm.set_efval(0, 0.7);
    EXPECT_EQ(eferm.ef, 0.7);
    eferm.set_efval(2, 0.77);
    EXPECT_EQ(eferm.ef, 0.77);
    eferm.two_efermi = true;
    eferm.set_efval(0, 0.6);
    EXPECT_EQ(eferm.ef_up, 0.6);
    eferm.set_efval(1, -1.0);
    EXPECT_EQ(eferm.ef_dw, -1.0);

    testing::internal::CaptureStdout();
    EXPECT_EXIT(eferm.set_efval(3, 1.0);, ::testing::ExitedWithCode(1), "");
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("Please check NSPIN when TWO_EFERMI is true"));
}

TEST_F(fenergy, eferm_get_efval)
{
    eferm.ef = 0.0;
    eferm.ef_up = 1.0;
    eferm.ef_dw = -1.0;
    eferm.two_efermi = false;
    EXPECT_EQ(eferm.get_efval(0), 0.0);
    eferm.two_efermi = true;
    EXPECT_EQ(eferm.get_efval(0), 1.0);
    EXPECT_EQ(eferm.get_efval(1), -1.0);

    testing::internal::CaptureStdout();
    EXPECT_EXIT(double tmpp = eferm.get_efval(2);, ::testing::ExitedWithCode(1), "");
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("Please check NSPIN when TWO_EFERMI is true"));
}
