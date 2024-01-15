#include <string>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "module_elecstate/elecstate_getters.h"

/***************************************************************
 *  unit test of class Occupy
 ****************************************************************/

/**
 * - Tested functions:
 *   - Occupy::Occupy()
 *   - Occupy::decision()
 */
#define private public
#include "module_elecstate/occupy.h"

class OccupyTest : public ::testing::Test
{
protected:
  Occupy occupy;
  std::string output;
};

TEST_F(OccupyTest, Occupy)
{
  EXPECT_EQ(occupy.use_gaussian_broadening, false);
  EXPECT_FALSE(occupy.gauss());
  EXPECT_EQ(occupy.fixed_occupations, false);
  EXPECT_FALSE(occupy.fix());
  EXPECT_EQ(occupy.gaussian_type, 0);
  EXPECT_DOUBLE_EQ(occupy.gaussian_parameter, 0.01);
}

TEST_F(OccupyTest, DecisionFixed)
{
  EXPECT_NO_THROW(occupy.decision("fixed", "fixed", 0.1));
  EXPECT_EQ(occupy.use_gaussian_broadening, false);
  EXPECT_EQ(occupy.fixed_occupations, false);
  EXPECT_EQ(occupy.gaussian_type, 0);
  EXPECT_DOUBLE_EQ(occupy.gaussian_parameter, 0.0);
}

TEST_F(OccupyTest, DecisionSmearingFixed)
{
  EXPECT_NO_THROW(occupy.decision("smearing", "fixed", 0.1));
  EXPECT_EQ(occupy.use_gaussian_broadening, false);
  EXPECT_EQ(occupy.fixed_occupations, false);
  EXPECT_EQ(occupy.gaussian_type, 0);
  EXPECT_DOUBLE_EQ(occupy.gaussian_parameter, 0.0);
}

TEST_F(OccupyTest, DecisionSmearingGaussian)
{
  EXPECT_NO_THROW(occupy.decision("smearing", "gaussian", 0.1));
  EXPECT_EQ(occupy.use_gaussian_broadening, true);
  EXPECT_EQ(occupy.fixed_occupations, false);
  EXPECT_EQ(occupy.gaussian_type, 0);
  EXPECT_DOUBLE_EQ(occupy.gaussian_parameter, 0.1);
}

TEST_F(OccupyTest, DecisionSmearingGauss)
{
  EXPECT_NO_THROW(occupy.decision("smearing", "gauss", 0.1));
  EXPECT_EQ(occupy.use_gaussian_broadening, true);
  EXPECT_EQ(occupy.fixed_occupations, false);
  EXPECT_EQ(occupy.gaussian_type, 0);
  EXPECT_DOUBLE_EQ(occupy.gaussian_parameter, 0.1);
}

TEST_F(OccupyTest, DecisionSmearingGaussianWarning)
{
  testing::internal::CaptureStdout();
  EXPECT_EXIT(occupy.decision("smearing", "gaussian", 0.0), ::testing::ExitedWithCode(0), "");
  output = testing::internal::GetCapturedStdout();
  // test output on screen
  EXPECT_THAT(
      output,
      testing::HasSubstr("Smearing requires gaussian broadening,but gaussian_parameter = 0(default value = 0.01)"));
}

TEST_F(OccupyTest, DecisionSmearingMethfesselPaxton)
{
  EXPECT_NO_THROW(occupy.decision("smearing", "methfessel-paxton", 0.1));
  EXPECT_EQ(occupy.use_gaussian_broadening, true);
  EXPECT_EQ(occupy.fixed_occupations, false);
  EXPECT_EQ(occupy.gaussian_type, 1);
  EXPECT_DOUBLE_EQ(occupy.gaussian_parameter, 0.1);
}

TEST_F(OccupyTest, DecisionSmearingMP)
{
  EXPECT_NO_THROW(occupy.decision("smearing", "mp", 0.1));
  EXPECT_EQ(occupy.use_gaussian_broadening, true);
  EXPECT_EQ(occupy.fixed_occupations, false);
  EXPECT_EQ(occupy.gaussian_type, 1);
  EXPECT_DOUBLE_EQ(occupy.gaussian_parameter, 0.1);
}

TEST_F(OccupyTest, DecissionMP2)
{
  EXPECT_NO_THROW(occupy.decision("smearing", "mp2", 0.1));
  EXPECT_EQ(occupy.use_gaussian_broadening, true);
  EXPECT_EQ(occupy.fixed_occupations, false);
  EXPECT_EQ(occupy.gaussian_type, 2);
  EXPECT_DOUBLE_EQ(occupy.gaussian_parameter, 0.1);
}

TEST_F(OccupyTest, DecisionMarzariVanderbilt)
{
  EXPECT_NO_THROW(occupy.decision("smearing", "marzari-vanderbilt", 0.1));
  EXPECT_EQ(occupy.use_gaussian_broadening, true);
  EXPECT_EQ(occupy.fixed_occupations, false);
  EXPECT_EQ(occupy.gaussian_type, -1);
  EXPECT_DOUBLE_EQ(occupy.gaussian_parameter, 0.1);
}

TEST_F(OccupyTest, DecisionCold)
{
  EXPECT_NO_THROW(occupy.decision("smearing", "cold", 0.1));
  EXPECT_EQ(occupy.use_gaussian_broadening, true);
  EXPECT_EQ(occupy.fixed_occupations, false);
  EXPECT_EQ(occupy.gaussian_type, -1);
  EXPECT_DOUBLE_EQ(occupy.gaussian_parameter, 0.1);
}

TEST_F(OccupyTest, DecisionMV)
{
  EXPECT_NO_THROW(occupy.decision("smearing", "mv", 0.1));
  EXPECT_EQ(occupy.use_gaussian_broadening, true);
  EXPECT_EQ(occupy.fixed_occupations, false);
  EXPECT_EQ(occupy.gaussian_type, -1);
  EXPECT_DOUBLE_EQ(occupy.gaussian_parameter, 0.1);
}

TEST_F(OccupyTest, DecisionFermiDirac)
{
  EXPECT_NO_THROW(occupy.decision("smearing", "fermi-dirac", 0.1));
  EXPECT_EQ(occupy.use_gaussian_broadening, true);
  EXPECT_EQ(occupy.fixed_occupations, false);
  EXPECT_EQ(occupy.gaussian_type, -99);
  EXPECT_DOUBLE_EQ(occupy.gaussian_parameter, 0.1);
}

TEST_F(OccupyTest, DecisionFD)
{
  EXPECT_NO_THROW(occupy.decision("smearing", "fd", 0.1));
  EXPECT_EQ(occupy.use_gaussian_broadening, true);
  EXPECT_EQ(occupy.fixed_occupations, false);
  EXPECT_EQ(occupy.gaussian_type, -99);
  EXPECT_DOUBLE_EQ(occupy.gaussian_parameter, 0.1);
}

TEST_F(OccupyTest, DecisionTetrahedraWarning)
{
  testing::internal::CaptureStdout();
  EXPECT_EXIT(occupy.decision("tetrahedra", "gaussian", 0.0), ::testing::ExitedWithCode(0), "");
  output = testing::internal::GetCapturedStdout();
  // test output on screen
  EXPECT_THAT(output, testing::HasSubstr("not implemented yet!"));
}

TEST_F(OccupyTest, DecisionFromInput)
{
  EXPECT_NO_THROW(occupy.decision("from_input", "fd", 0.1));
  EXPECT_EQ(occupy.use_gaussian_broadening, false);
  EXPECT_EQ(occupy.fixed_occupations, true);
  EXPECT_EQ(occupy.gaussian_type, 0);
  EXPECT_DOUBLE_EQ(occupy.gaussian_parameter, 0.1);
}

TEST_F(OccupyTest, DecisionArbitrary)
{
  testing::internal::CaptureStdout();
  EXPECT_EXIT(occupy.decision("arbitrary", "gaussian", 0.0), ::testing::ExitedWithCode(0), "");
  output = testing::internal::GetCapturedStdout();
  // test output on screen
  EXPECT_THAT(output, testing::HasSubstr("occupations, not implemented"));
}

TEST_F(OccupyTest, IweightsNOSPIN)
{
  GlobalV::NSPIN = 1;
  double ef = 0.0;
  ModuleBase::matrix wg(1, 1);
  std::vector<double> wk(1, 2.0);
  ModuleBase::matrix ekb(1, 1);
  std::vector<int> isk(1);
  ekb(0, 0) = 0.1;
  occupy.iweights(1, wk, 1, 2.0, ekb, ef, wg, 0, isk);
  EXPECT_DOUBLE_EQ(wg(0, 0), 2.0);
  EXPECT_DOUBLE_EQ(ef, 0.1);
}

TEST_F(OccupyTest, IweightsSPIN)
{
  GlobalV::NSPIN = 2;
  double ef_up = 0.0;
  double ef_dw = 0.0;
  ModuleBase::matrix wg(2, 1);
  std::vector<double> wk(2, 1.0);
  ModuleBase::matrix ekb(2, 1);
  std::vector<int> isk(2);
  isk[0] = 0;
  isk[1] = 1;
  ekb(0, 0) = 0.1;
  ekb(1, 0) = 0.2;
  occupy.iweights(2, wk, 1, 1.0, ekb, ef_up, wg, 0, isk);
  occupy.iweights(2, wk, 1, 1.0, ekb, ef_dw, wg, 1, isk);
  EXPECT_DOUBLE_EQ(wg(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(wg(1, 0), 1.0);
  EXPECT_DOUBLE_EQ(ef_up, 0.1);
  EXPECT_DOUBLE_EQ(ef_dw, 0.2);
}

TEST_F(OccupyTest, IweightsWarning)
{
  GlobalV::NSPIN = 1;
  double ef = 0.0;
  ModuleBase::matrix wg(1, 1);
  std::vector<double> wk(1, 2.0);
  ModuleBase::matrix ekb(1, 1);
  std::vector<int> isk(1);
  ekb(0, 0) = 0.1;

  testing::internal::CaptureStdout();
  EXPECT_EXIT(occupy.iweights(1, wk, 1, 1.0, ekb, ef, wg, -1, isk);, ::testing::ExitedWithCode(0), "");
  output = testing::internal::GetCapturedStdout();
  EXPECT_THAT(output, testing::HasSubstr("It is not a semiconductor or insulator. Please do not set 'smearing_method=fixed', and try other options."));
}

TEST_F(OccupyTest, Wgauss)
{
  EXPECT_DOUBLE_EQ(occupy.wgauss(0.0, 0), 0.5);
  EXPECT_DOUBLE_EQ(occupy.wgauss(0.0, -1), 0.4006259784506005);
  EXPECT_DOUBLE_EQ(occupy.wgauss(0.0, -99), 0.5);
  EXPECT_DOUBLE_EQ(occupy.wgauss(0.0, 1), 0.5);
  EXPECT_DOUBLE_EQ(occupy.wgauss(0.0, 2), 0.5);
  EXPECT_DOUBLE_EQ(occupy.wgauss(10, 0), 1.0);
}

TEST_F(OccupyTest, W1gauss)
{
  EXPECT_DOUBLE_EQ(occupy.w1gauss(0.0, 0), -0.28209479177387814);
  EXPECT_DOUBLE_EQ(occupy.w1gauss(0.0, -1), -0.1710991401561083);
  EXPECT_DOUBLE_EQ(occupy.w1gauss(0.0, -99), -0.69314718055994529);
  EXPECT_DOUBLE_EQ(occupy.w1gauss(0.0, 1), -0.14104739588693907);
  EXPECT_DOUBLE_EQ(occupy.w1gauss(0.0, 2), -0.10578554691520431);
}

TEST_F(OccupyTest, Sumkg)
{
  ModuleBase::matrix ekb(1, 1);
  ekb(0, 0) = -1.0;
  std::vector<double> wk = {1, 1.0};
  double smearing_sigma = 0.1;
  int ngauss = 0;
  double e = 0.0;
  int is = 0;
  std::vector<int> isk = {0, 0};
  EXPECT_DOUBLE_EQ(occupy.sumkg(ekb, 1, 1, wk, smearing_sigma, ngauss, e, is, isk), 1.0);
}

TEST_F(OccupyTest, Efermig)
{
  ModuleBase::matrix ekb(1, 1);
  ekb(0, 0) = -1.0;
  std::vector<double> wk = {1, 1.0};
  double smearing_sigma = 0.1;
  int ngauss = 0;
  double e = 0.0;
  int is = 0;
  std::vector<int> isk = {0, 0};
  double ef = 0.0;
  occupy.efermig(ekb, 1, 1, 1.0, wk, smearing_sigma, ngauss, ef, is, isk);
  EXPECT_NEAR(ef, -0.5, 1e-13);
}

TEST_F(OccupyTest, Gweights)
{
  ModuleBase::matrix ekb(1, 1);
  ekb(0, 0) = -1.0;
  std::vector<double> wk = {1, 1.0};
  double smearing_sigma = 0.1;
  int ngauss = 0;
  double e = 0.0;
  int is = 0;
  std::vector<int> isk = {0, 0};
  double ef = 0.0;
  ModuleBase::matrix wg(1, 1);
  wg(0, 0) = 1.0;
  double demet = 0.0;
  occupy.gweights(1, wk, 1, 1.0, smearing_sigma, ngauss, ekb, ef, demet, wg, is, isk);
  EXPECT_NEAR(ef, -0.5, 1e-13);
  EXPECT_NEAR(demet, 0.0, 1e-13);
  EXPECT_NEAR(wg(0, 0), 1.0, 1e-13);
}

#undef private
