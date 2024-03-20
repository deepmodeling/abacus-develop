#include <gtest/gtest.h>
#include <iostream>
#include <vector>
#define protected public

#include "module_hsolver/diago_iter_assis.h"
#include "diago_mock.h"

class TestDiagoIterAssist : public ::testing::Test
{
	public:
	using dia_f = hsolver::DiagoIterAssistSolver<float, psi::DEVICE_CPU>;
	using dia_d = hsolver::DiagoIterAssist<std::complex<double>, psi::DEVICE_CPU>;

	hamilt::Hamilt<std::complex<double>> hamilt_test_d;
	hamilt::Hamilt<std::complex<float>> hamilt_test_f;

    DIAGOTEST::hamilt.create(4, 4);

	psi::Psi<std::complex<double>> psi_test_cd;
	psi::Psi<std::complex<float>> psi_test_cf;
	
	elecstate::ElecState elecstate_test;

	std::string method_test = "none";

	std::ofstream temp_ofs;
};

TEST_F(TestDiagoIterAssist, diagH_subspace)
{
    dia_f::diagH_subspace();
    dia_d::diagH_subspace();
    EXPECT_EQ(true);
}

TEST_F(TestDiagoIterAssist, diagH_LAPACK)
{
    EXPECT_EQ(true);
}

TEST_F(TestDiagoIterAssist, test_exit_cond)
{
    EXPECT_EQ(true);
}