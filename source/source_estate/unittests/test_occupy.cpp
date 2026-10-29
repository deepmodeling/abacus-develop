#include <string>
#include <cmath>
#include <limits>
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "source_io/module_parameter/parameter.h"

/***************************************************************
 *  unit test of class Occupy
 ****************************************************************/

/**
 * - Tested functions:
 *   - Occupy::Occupy()
 *   - Occupy::decision()
 */
#include "source_estate/occupy.h"
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
  EXPECT_EXIT(occupy.decision("smearing", "gaussian", 0.0), ::testing::ExitedWithCode(1), "");
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
  EXPECT_EXIT(occupy.decision("tetrahedra", "gaussian", 0.0), ::testing::ExitedWithCode(1), "");
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
  EXPECT_EXIT(occupy.decision("arbitrary", "gaussian", 0.0), ::testing::ExitedWithCode(1), "");
  output = testing::internal::GetCapturedStdout();
  // test output on screen
  EXPECT_THAT(output, testing::HasSubstr("occupations, not implemented"));
}

TEST_F(OccupyTest, IweightsNOSPIN)
{
  double ef = 0.0;
  ModuleBase::matrix wg(1, 1);
  std::vector<double> wk(1, 2.0);
  ModuleBase::matrix ekb(1, 1);
  std::vector<int> isk(1);
  ekb(0, 0) = 0.1;
  occupy.iweights(1, wk, 1, 0, 2.0, ekb, ef, wg, 0, isk, /*nspin=*/1);
  EXPECT_DOUBLE_EQ(wg(0, 0), 2.0);
  EXPECT_DOUBLE_EQ(ef, 0.1);
}

TEST_F(OccupyTest, IweightsSPIN)
{
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
  occupy.iweights(2, wk, 1, 0, 1.0, ekb, ef_up, wg, 0, isk, /*nspin=*/2);
  occupy.iweights(2, wk, 1, 0, 1.0, ekb, ef_dw, wg, 1, isk, /*nspin=*/2);
  EXPECT_DOUBLE_EQ(wg(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(wg(1, 0), 1.0);
  EXPECT_DOUBLE_EQ(ef_up, 0.1);
  EXPECT_DOUBLE_EQ(ef_dw, 0.2);
}

TEST_F(OccupyTest, IweightsWarning)
{
  double ef = 0.0;
  ModuleBase::matrix wg(1, 1);
  std::vector<double> wk(1, 2.0);
  ModuleBase::matrix ekb(1, 1);
  std::vector<int> isk(1);
  ekb(0, 0) = 0.1;

  testing::internal::CaptureStdout();
  EXPECT_EXIT(occupy.iweights(1, wk, 1, 0, 1.0, ekb, ef, wg, -1, isk, /*nspin=*/1);, ::testing::ExitedWithCode(1), "");
  output = testing::internal::GetCapturedStdout();
  EXPECT_THAT(output, testing::HasSubstr("It is not a semiconductor or insulator. Please do not set 'smearing_method=fixed', and try other options."));
}

TEST_F(OccupyTest, Wgauss)
{
  EXPECT_DOUBLE_EQ(occupy_smearing::wgauss(0.0, 0), 0.5);
  EXPECT_DOUBLE_EQ(occupy_smearing::wgauss(0.0, -1), 0.4006259784506005);
  EXPECT_DOUBLE_EQ(occupy_smearing::wgauss(0.0, -99), 0.5);
  EXPECT_DOUBLE_EQ(occupy_smearing::wgauss(0.0, 1), 0.5);
  EXPECT_DOUBLE_EQ(occupy_smearing::wgauss(0.0, 2), 0.5);
  EXPECT_DOUBLE_EQ(occupy_smearing::wgauss(10, 0), 1.0);
}

TEST_F(OccupyTest, W1gauss)
{
  EXPECT_DOUBLE_EQ(occupy_smearing::w1gauss(0.0, 0), -0.28209479177387814);
  EXPECT_DOUBLE_EQ(occupy_smearing::w1gauss(0.0, -1), -0.1710991401561083);
  EXPECT_DOUBLE_EQ(occupy_smearing::w1gauss(0.0, -99), -0.69314718055994529);
  EXPECT_DOUBLE_EQ(occupy_smearing::w1gauss(0.0, 1), -0.14104739588693907);
  EXPECT_DOUBLE_EQ(occupy_smearing::w1gauss(0.0, 2), -0.10578554691520431);
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
  EXPECT_DOUBLE_EQ(occupy_smearing::sumkg(ekb, 1, 1, wk, smearing_sigma, ngauss, e, is, isk), 1.0);
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
  occupy_smearing::efermig(ekb, 1, 1, 1.0, wk, smearing_sigma, ngauss, ef, is, isk);
  // Full occupation fixes particle number, not a unique finite chemical potential.
  EXPECT_TRUE(std::isfinite(ef));
  EXPECT_GT(ef, ekb(0, 0));
  const double count = occupy_smearing::sumkg(ekb, 1, 1, wk, smearing_sigma, ngauss, ef, is, isk);
  EXPECT_NEAR(count, 1.0, 2e-15);
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
  wg(0, 0) = -123.0;
  double demet = 0.0;
  occupy.gweights(1, wk, 1, 1.0, smearing_sigma, ngauss, ekb, ef, demet, wg, is, isk);
  // Full occupation fixes particle number, not a unique finite chemical potential.
  EXPECT_TRUE(std::isfinite(ef));
  EXPECT_GT(ef, ekb(0, 0));
  const double count = occupy_smearing::sumkg(ekb, 1, 1, wk, smearing_sigma, ngauss, ef, is, isk);
  EXPECT_NEAR(count, 1.0, 2e-15);
  EXPECT_NEAR(demet, 0.0, 1e-13);
  EXPECT_NEAR(wg(0, 0), 1.0, 1e-13);
}

TEST_F(OccupyTest, FractionalFermiDiracOccupation)
{
    ModuleBase::matrix energies(1, 1);
    energies(0, 0) = -0.7;
    const std::vector<double> weights{1.0};
    const std::vector<int> spins{0};
    const double sigma = 0.13;
    const double target = 0.3141592653589793;
    double ef = 0.0;
    occupy_smearing::efermig(energies, 1, 1, target, weights, sigma, -99, ef, -1, spins);
    const double expected = energies(0, 0) + sigma * std::log(target / (1.0 - target));
    EXPECT_NEAR(ef, expected, 2e-15);
    ModuleBase::matrix occupations(1, 1);
    occupations(0, 0) = -123.0;
    double entropy = 0.0;
    occupy.gweights(1, weights, 1, target, sigma, -99, energies, ef, entropy, occupations, -1, spins);
    EXPECT_NEAR(occupations(0, 0), target, 2e-15);
    const double expected_entropy = sigma * (target * std::log(target) + (1.0 - target) * std::log(1.0 - target));
    EXPECT_NEAR(entropy, expected_entropy, 2e-15);
}

TEST_F(OccupyTest, EmptyFullAndGappedSpectra)
{
    ModuleBase::matrix energies(1, 2);
    energies(0, 0) = -1.0;
    energies(0, 1) = 1.0;
    const std::vector<double> weights{1.0};
    const std::vector<int> spins{0};
    for (const double target : {0.0, 1.0, 2.0})
    {
        double ef = 0.0;
        double entropy = 0.0;
        ModuleBase::matrix occupations(1, 2);
        occupy.gweights(1, weights, 2, target, 0.01, 0, energies, ef, entropy, occupations, -1, spins);
        const double count = occupations(0, 0) + occupations(0, 1);
        EXPECT_NEAR(count, target, 4e-15);
        EXPECT_NEAR(entropy, 0.0, 1e-14);
        if (target == 1.0)
        {
            EXPECT_GT(ef, energies(0, 0));
            EXPECT_LT(ef, energies(0, 1));
            EXPECT_NEAR(occupations(0, 0), 1.0, 2e-15);
            EXPECT_NEAR(occupations(0, 1), 0.0, 2e-15);
        }
    }
}

TEST_F(OccupyTest, SelectedSpinAndUnsortedBands)
{
    ModuleBase::matrix energies(2, 3);
    energies(0, 0) = 100.0;
    energies(0, 1) = 200.0;
    energies(0, 2) = 300.0;
    energies(1, 0) = 0.4;
    energies(1, 1) = -0.3;
    energies(1, 2) = 0.1;
    const std::vector<double> weights{1.0, 1.0};
    const std::vector<int> spins{0, 1};
    ModuleBase::matrix occupations(2, 3);
    for (int ib = 0; ib < 3; ++ib)
    {
        occupations(0, ib) = -123.0;
    }
    double ef = 0.0;
    double entropy = 0.0;
    occupy.gweights(2, weights, 3, 1.314159, 0.08, 0, energies, ef, entropy, occupations, 1, spins);
    double count = 0.0;
    for (int ib = 0; ib < 3; ++ib)
    {
        EXPECT_EQ(occupations(0, ib), -123.0);
        count += occupations(1, ib);
    }
    EXPECT_NEAR(count, 1.314159, 3e-15);
}

TEST_F(OccupyTest, CompensatedWeightedOccupationSum)
{
    const int points = 4096;
    ModuleBase::matrix energies(points, 2);
    std::vector<double> weights(points, 1.0 / points);
    const std::vector<int> spins(points, 0);
    long double expected = 0.0L;
    for (int ik = 0; ik < points; ++ik)
    {
        // Repeated degenerate levels expose drift in ordinary accumulation.
        energies(ik, 0) = -0.1;
        energies(ik, 1) = 0.2;
        for (int ib = 0; ib < 2; ++ib)
        {
            const double argument = (0.123 - energies(ik, ib)) / 0.1;
            expected += weights[ik] * 0.5 * std::erfc(-argument);
        }
    }
    const double count = occupy_smearing::sumkg(energies, 2, points, weights, 0.1, 0, 0.123, -1, spins);
    EXPECT_NEAR(count, expected, 5e-16);
}

TEST_F(OccupyTest, FermiSearchAtRepresentableEnergyLimit)
{
    ModuleBase::matrix energies(1, 1);
    energies(0, 0) = 0.0;
    const std::vector<double> weights{1.0};
    const std::vector<int> spins{0};
    double ef = 0.0;
    // A subnormal width leaves too few representable offsets to resolve this
    // occupation, even after moving the energy origin close to the root.
    const double sigma = 1024.0 * std::numeric_limits<double>::denorm_min();
    testing::internal::CaptureStdout();
    EXPECT_EXIT(occupy_smearing::efermig(energies, 1, 1, 0.3, weights, sigma, 0, ef, -1, spins),
                ::testing::ExitedWithCode(1), "");
    const std::string diagnostic = testing::internal::GetCapturedStdout();
    EXPECT_THAT(diagnostic, testing::HasSubstr("Fermi bracket reached floating-point resolution"));
    EXPECT_THAT(diagnostic, testing::HasSubstr("electron residual"));
}

TEST_F(OccupyTest, RejectInvalidWidthAndUnbracketedCharge)
{
    ModuleBase::matrix energies(1, 1);
    energies(0, 0) = 0.0;
    const std::vector<double> weights{1.0};
    const std::vector<int> spins{0};
    double ef = 0.0;
    EXPECT_EXIT(occupy_smearing::efermig(energies, 1, 1, 0.5, weights, 0.0, 0, ef, -1, spins),
                ::testing::ExitedWithCode(1), "");
    EXPECT_EXIT(occupy_smearing::efermig(energies, 1, 1, 1.1, weights, 0.1, 0, ef, -1, spins),
                ::testing::ExitedWithCode(1), "");
}

TEST_F(OccupyTest, FermiParticleNumberAcrossScales)
{
    ModuleBase::matrix energies(3, 64);
    const std::vector<double> weights{0.17, 0.29, 0.54};
    const std::vector<int> spins{0, 0, 0};
    for (int ik = 0; ik < 3; ++ik)
    {
        for (int ib = 0; ib < 64; ++ib)
        {
            energies(ik, ib) = -1.0 + 0.05 * ib + 0.013 * ik;
        }
    }
    // Check the returned occupations across electron-number scales and smearings.
    for (const double target : {0.3141592653589793, 33.0, 63.125})
    {
        for (const int smearing : {0, 1, 2, -1, -99})
        {
            double ef = 0.0;
            double entropy = 0.0;
            ModuleBase::matrix occupations(3, 64);
            occupy.gweights(3, weights, 64, target, 0.1, smearing, energies, ef, entropy, occupations, -1, spins);
            long double count = 0.0L;
            for (int ik = 0; ik < 3; ++ik)
            {
                for (int ib = 0; ib < 64; ++ib)
                {
                    count += occupations(ik, ib);
                }
            }
            const double scale = std::max(1.0, target);
            // Add one epsilon for the independent sum of the returned doubles.
            const double tolerance = 17.0 * std::numeric_limits<double>::epsilon() * scale;
            EXPECT_NEAR(count, target, tolerance) << "N = " << target << ", smearing = " << smearing;
        }
    }
}

TEST_F(OccupyTest, GweightsHonorsExplicitBandRange)
{
    ModuleBase::matrix energies(1, 3);
    energies(0, 0) = -0.7;
    energies(0, 1) = 0.1;
    energies(0, 2) = 0.8;
    ModuleBase::matrix occupations(1, 3);
    for (int ib = 0; ib < 3; ++ib)
    {
        occupations(0, ib) = -123.0;
    }
    const std::vector<double> weights{1.0};
    const std::vector<int> spins{0};
    const double target = 0.3141592653589793;
    double ef = 0.0;
    double entropy = 0.0;
    // The explicit one-band range must govern both root finding and writing.
    occupy.gweights(1, weights, 1, target, 0.13, -99, energies, ef, entropy, occupations, -1, spins);
    EXPECT_NEAR(occupations(0, 0), target, 2e-15);
    EXPECT_DOUBLE_EQ(occupations(0, 1), -123.0);
    EXPECT_DOUBLE_EQ(occupations(0, 2), -123.0);
}

TEST_F(OccupyTest, OccupationAndEntropyEnergyShift)
{
    ModuleBase::matrix energies(1, 3);
    ModuleBase::matrix occupations(1, 3);
    ModuleBase::matrix reference(1, 3);
    const std::vector<double> weights{2.0};
    const std::vector<int> spins{0};
    const double target = 2.3;
    for (const int smearing : {0, 1, 2, -1, -99})
    {
        double reference_entropy = 0.0;
        double reference_fermi = 0.0;
        for (const double shift : {0.0, 1e8})
        {
            for (int ib = 0; ib < 3; ++ib)
            {
                // Binary fractions ensure the shifted input spectrum is exact.
                energies(0, ib) = shift - 0.5 + 0.5 * ib;
            }
            double ef = 0.0;
            double entropy = 0.0;
            occupy.gweights(1, weights, 3, target, 0.125, smearing, energies, ef, entropy, occupations, -1, spins);
            if (shift == 0.0)
            {
                reference_entropy = entropy;
                reference_fermi = ef;
                for (int ib = 0; ib < 3; ++ib)
                {
                    reference(0, ib) = occupations(0, ib);
                }
            }
            double count = 0.0;
            for (int ib = 0; ib < 3; ++ib)
            {
                EXPECT_NEAR(occupations(0, ib), reference(0, ib), 2e-15);
                count += occupations(0, ib);
            }
            EXPECT_NEAR(count, target, 8e-15);
            EXPECT_NEAR(entropy, reference_entropy, 2e-15);
            const double expected_fermi = reference_fermi + shift;
            EXPECT_DOUBLE_EQ(ef, expected_fermi);
        }
    }
}

TEST_F(OccupyTest, NarrowSmearingAtLargeEnergyOffset)
{
    ModuleBase::matrix energies(1, 1);
    energies(0, 0) = 1e8;
    ModuleBase::matrix occupations(1, 1);
    const std::vector<double> weights{1.0};
    const std::vector<int> spins{0};
    double ef = 0.0;
    double entropy = 0.0;
    occupy.gweights(1, weights, 1, 0.3, 1e-8, 0, energies, ef, entropy, occupations, -1, spins);
    const double tolerance = 16.0 * std::numeric_limits<double>::epsilon();
    EXPECT_NEAR(occupations(0, 0), 0.3, tolerance);
    EXPECT_TRUE(std::isfinite(ef));
    EXPECT_TRUE(std::isfinite(entropy));
}

TEST_F(OccupyTest, FermiSurfaceAwayFromSpectralMidpoint)
{
    ModuleBase::matrix energies(1, 7);
    ModuleBase::matrix occupations(1, 7);
    const std::vector<double> weights{2.0};
    const std::vector<int> spins{0};
    const double target = 10.0;
    const double tolerance = 16.0 * std::numeric_limits<double>::epsilon() * target;
    for (int ib = 0; ib < 6; ++ib)
    {
        energies(0, ib) = 0.0;
    }
    energies(0, 6) = 100.0;
    double ef = 0.0;
    double entropy = 0.0;
    occupy.gweights(1, weights, 7, target, 0.01, 0, energies, ef, entropy, occupations, -1, spins);
    long double count = 0.0L;
    for (int ib = 0; ib < 7; ++ib)
    {
        count += occupations(0, ib);
    }
    EXPECT_NEAR(count, target, tolerance);
    EXPECT_DOUBLE_EQ(occupations(0, 6), 0.0);
}

TEST_F(OccupyTest, LowFillingDegenerateSpectrum)
{
    // Gaussian lower-tail cancellation and MP1 cancellation near zero filling.
    struct Case
    {
        int bands;
        int smearing;
        double target;
    };
    const Case cases[] = {{74, 0, 2.0}, {26, 1, 0.3}};
    const std::vector<double> weights{2.0};
    const std::vector<int> spins{0};
    for (const Case& test : cases)
    {
        ModuleBase::matrix energies(1, test.bands);
        ModuleBase::matrix occupations(1, test.bands);
        double ef = 0.0;
        double entropy = 0.0;
        occupy.gweights(1, weights, test.bands, test.target, 0.01, test.smearing,
                        energies, ef, entropy, occupations, -1, spins);
        long double count = 0.0L;
        for (int ib = 0; ib < test.bands; ++ib)
        {
            count += occupations(0, ib);
        }
        const double scale = std::max(1.0, test.target);
        // Include one epsilon for the independent accumulation.
        const double tolerance = 17.0 * std::numeric_limits<double>::epsilon() * scale;
        EXPECT_NEAR(count, test.target, tolerance) << test.smearing;
        EXPECT_TRUE(std::isfinite(entropy));
    }
}

TEST_F(OccupyTest, SmearingHighPrecisionReference)
{
    // Independent 100-decimal-digit erfc/exp evaluations of the closed forms.
    // Include signed MP tails and cold-smearing overshoot; do not clamp to [0,1].
    struct Reference
    {
        int order;
        double x;
        double occupation;
    };
    const Reference references[] = {
        {0, -8.0, 5.61214858649146354e-30},
        {-1, -8.0, 5.1183964952735562742e-34},
        {1, -2.0, -0.007994625186522394009569},
        {2, -0.5, 0.06124734974721450735834},
        {-1, 2.0, 1.041238947228640924942},
        {-99, -0.5, 0.3775406687981454353611}};
    for (const Reference& reference : references)
    {
        const double value = occupy_smearing::wgauss(reference.x, reference.order);
        const double tolerance = 32.0 * std::numeric_limits<double>::epsilon() * std::abs(reference.occupation);
        EXPECT_NEAR(value, reference.occupation, tolerance) << reference.order << ", x = " << reference.x;
    }
}
