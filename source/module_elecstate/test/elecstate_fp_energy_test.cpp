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
 *   - efermi::get_ef()
 *   - efermi::get_efval()
 */
class fenergy : public ::testing::Test
{
  protected:
    elecstate::fenergy f_en;
    elecstate::efermi eferm;
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
    double& tmp_ef = eferm.get_ef(0);
    tmp_ef = 0.7;
    EXPECT_EQ(eferm.ef, 0.7);
    eferm.two_efermi = true;
    double& tmp_efup = eferm.get_ef(0);
    tmp_efup = 1.0;
    EXPECT_EQ(eferm.ef_up, 1.0);
    double& tmp_efdw = eferm.get_ef(1);
    tmp_efdw = -1.0;
    EXPECT_EQ(eferm.ef_dw, -1.0);

    testing::internal::CaptureStdout();
    EXPECT_EXIT(double& tmpp = eferm.get_ef(2);, ::testing::ExitedWithCode(0), "");
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
    EXPECT_EXIT(double tmpp = eferm.get_efval(2);, ::testing::ExitedWithCode(0), "");
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("Please check NSPIN when TWO_EFERMI is true"));
}

#include "mpi.h"
int main(int argc, char **argv)
{
        MPI_Init(&argc, &argv);
        testing::InitGoogleTest(&argc, argv);

        int nproc;
        int my_rank;

        MPI_Comm_size(MPI_COMM_WORLD,&nproc);
        MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

        int result = RUN_ALL_TESTS();

        MPI_Finalize();

        return result;
}