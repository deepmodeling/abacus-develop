#include "../clebsch_gordan_coeff.h"

#include <iostream>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

/************************************************
 *  unit test of class Clebsch_Gordan
 ***********************************************/

/**
 * - Tested Functions:
 *   - clebsch_gordan
 *     - computes Clebsch-Gordan coefficient
 *     - functions: gen_rndm_r and compute_ap
 */

namespace ModuleBase
{
void WARNING_QUIT(const std::string& file, const std::string& description)
{
    return;
}
void WARNING(const std::string& file, const std::string& description)
{
    return;
}
} // namespace ModuleBase

TEST(ClebschGordanTest, ClebschGordanExit)
{
    int lmaxkb = -2;
    ModuleBase::realArray ap;
    ModuleBase::IntArray lpx;
    ModuleBase::IntArray lpl;

    std::string output;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(ModuleBase::Clebsch_Gordan::clebsch_gordan(lmaxkb + 1, ap, lpx, lpl), ::testing::ExitedWithCode(1), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("Clebsch_Gordan: lmaxkb + 1 < 0"));
}

TEST(ClebschGordanTest, ClebschGordan)
{
    int lmaxkb = 1;
    ModuleBase::realArray ap;
    ModuleBase::IntArray lpx;
    ModuleBase::IntArray lpl;

    ModuleBase::Clebsch_Gordan::clebsch_gordan(lmaxkb + 1, ap, lpx, lpl);
    EXPECT_DOUBLE_EQ(ap(0, 0, 0), 0.28209479177387564);
    EXPECT_DOUBLE_EQ(ap(8, 0, 0), 0.0);
    EXPECT_EQ(lpx(0, 0), 1);
    EXPECT_EQ(lpx(3, 3), 3);
    EXPECT_EQ(lpl(0, 0, 5), 0);
    EXPECT_EQ(lpl(3, 3, 8), 0);
}