#include "../nacv_cs.h"

#include <gtest/gtest.h>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace
{
LR::NacvResult example_result(double gap)
{
    LR::NacvResult result;
    result.gap_rydberg = gap;
    result.ao_connection.create(1, 3);
    result.ao_connection(0, 0) = 0.25;
    result.ao_connection(0, 1) = -0.5;
    return result;
}

ModuleBase::matrix example_force()
{
    ModuleBase::matrix force(1, 3);
    force(0, 0) = -2.0;
    force(0, 1) = 4.0;
    return force;
}
} // namespace

TEST(NacvResult, RequiresCompleteOrbitalWindowWithoutIntegerOverflow)
{
    EXPECT_NO_THROW(LR::validate_nac_window(8, 20, 28));
    EXPECT_THROW(LR::validate_nac_window(7, 20, 28), std::invalid_argument);
    EXPECT_THROW(LR::validate_nac_window(8, 19, 28), std::invalid_argument);
    EXPECT_THROW(LR::validate_nac_window(0, 28, 28), std::invalid_argument);
    EXPECT_THROW(LR::validate_nac_window(28, 0, 28), std::invalid_argument);
    EXPECT_THROW(LR::validate_nac_window(-1, 29, 28), std::invalid_argument);
    const int largest = std::numeric_limits<int>::max();
    EXPECT_THROW(LR::validate_nac_window(largest, largest, 28), std::invalid_argument);
}

TEST(NacvResult, DefersOnlyUnresolvedDefaultsBeforeScf)
{
    EXPECT_NO_THROW(LR::validate_nac_xc_request("default", "PBE"));
    EXPECT_NO_THROW(LR::validate_nac_xc_request("PZ", "default"));
    EXPECT_NO_THROW(LR::validate_nac_xc_request("DEFAULT", "default"));
    EXPECT_NO_THROW(LR::validate_nac_xc_request("NONE", "rpa"));
    EXPECT_NO_THROW(LR::validate_nac_xc_request("HF", "hf"));
    EXPECT_THROW(LR::validate_nac_functionals("default", "pbe"), std::invalid_argument);
    EXPECT_THROW(LR::validate_nac_functionals("PBE", "default"), std::invalid_argument);
    EXPECT_THROW(LR::validate_nac_functionals("default", "default"), std::invalid_argument);
}

TEST(NacvResult, RejectsKnownUnavailableBuiltinsBeforeScf)
{
    const std::vector<std::string> unsupported = {"PBESOL",
        "REVPBE",
        "WC",
        "BLYP",
        "BP",
        "PW91",
        "HCTH",
        "OLYP",
        "SCAN",
        "SCAN0",
        "SCANL",
        "OPT_ORB",
        "MULLER",
        "POWER",
        "WP22",
        "CWP22",
        "BLYP_LR"};
    for (const std::string &name : unsupported)
    {
        EXPECT_THROW(LR::validate_nac_xc_request(name, "pbe"), std::invalid_argument);
        EXPECT_THROW(LR::validate_nac_xc_request("PBE", name), std::invalid_argument);
        EXPECT_THROW(LR::validate_nac_xc_request("default", name), std::invalid_argument);
        EXPECT_THROW(LR::validate_nac_xc_request(name, "default"), std::invalid_argument);
    }
}

TEST(NacvResult, AcceptsIntegerGammaWeightsForBothClosedShellConventions)
{
    for (int spins : {1, 2})
    {
        ModuleBase::matrix weights(spins, 4);
        for (int spin = 0; spin < spins; ++spin)
        {
            weights(spin, 0) = 2.0 / spins;
            weights(spin, 1) = 2.0 / spins;
        }
        EXPECT_NO_THROW(
            LR::validate_nac_occupations(
                weights,
                2,
                4,
                spins));
        weights(0, 2) = std::numeric_limits<double>::epsilon();
        EXPECT_NO_THROW(
            LR::validate_nac_occupations(
                weights,
                2,
                4,
                spins));
    }
}

TEST(NacvResult, RejectsFractionalNonfiniteAndIncompleteOccupations)
{
    ModuleBase::matrix weights(1, 4);
    weights(0, 0) = 2.0;
    weights(0, 1) = 1.9;
    weights(0, 2) = 0.1;
    EXPECT_THROW(
        LR::validate_nac_occupations(
            weights,
            2,
            4,
            1),
        std::invalid_argument);
    weights(0, 1) = 2.0;
    weights(0, 2) = 1e-8;
    EXPECT_THROW(
        LR::validate_nac_occupations(
            weights,
            2,
            4,
            1),
        std::invalid_argument);
    weights(0, 2) = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(
        LR::validate_nac_occupations(
            weights,
            2,
            4,
            1),
        std::invalid_argument);
    EXPECT_THROW(
        LR::validate_nac_occupations(
            weights,
            2,
            5,
            1),
        std::invalid_argument);
    EXPECT_THROW(
        LR::validate_nac_occupations(
            weights,
            0,
            4,
            1),
        std::invalid_argument);
    EXPECT_THROW(
        LR::validate_nac_occupations(
            weights,
            2,
            4,
            4),
        std::invalid_argument);
}

TEST(NacvResult, RequiresImplementedGroundAndExcitedResponseKernels)
{
    EXPECT_NO_THROW(LR::validate_nac_functionals("PBE", "pbe"));
    EXPECT_NO_THROW(LR::validate_nac_functionals("PZ", "lda"));
    EXPECT_NO_THROW(LR::validate_nac_functionals("SLAPWPBXPBC", "pbe"));
    EXPECT_NO_THROW(LR::validate_nac_functionals("PBE0", "pbe"));
    EXPECT_NO_THROW(LR::validate_nac_functionals("HF", "hf"));
    EXPECT_NO_THROW(LR::validate_nac_functionals("NOX+NOC", "rpa"));
    EXPECT_THROW(LR::validate_nac_functionals("SCAN", "pbe"), std::invalid_argument);
    EXPECT_THROW(LR::validate_nac_functionals("pbe", "scan"), std::invalid_argument);
    EXPECT_THROW(LR::validate_nac_functionals("pbe", "bse"), std::invalid_argument);
    EXPECT_THROW(LR::validate_nac_functionals("default", "pbe"), std::invalid_argument);
    EXPECT_THROW(LR::validate_nac_functionals("MULLER", "hf"), std::invalid_argument);
}

TEST(NacvResult, GapResolutionUsesOperatorResidualsAndRoundoff)
{
    const double roundoff = 16.0 * std::numeric_limits<double>::epsilon() * 2.0;
    EXPECT_DOUBLE_EQ(
        LR::nacv_gap_resolution(
            1.0,
            2.0,
            1e-10,
            2e-10),
        3e-10 + roundoff);
    EXPECT_DOUBLE_EQ(
        LR::nacv_gap_resolution(
            0.0,
            0.0,
            0.0,
            0.0),
        0.0);
    EXPECT_DOUBLE_EQ(
        LR::nacv_gap_resolution(
            -2.0,
            1.0,
            0.0,
            0.0),
        roundoff);
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double infinity = std::numeric_limits<double>::infinity();
    EXPECT_THROW(
        LR::nacv_gap_resolution(
            nan,
            1.0,
            0.0,
            0.0),
        std::invalid_argument);
    EXPECT_THROW(
        LR::nacv_gap_resolution(
            0.0,
            infinity,
            0.0,
            0.0),
        std::invalid_argument);
    EXPECT_THROW(
        LR::nacv_gap_resolution(
            0.0,
            1.0,
            -1.0,
            0.0),
        std::invalid_argument);
    EXPECT_THROW(
        LR::nacv_gap_resolution(
            0.0,
            1.0,
            0.0,
            nan),
        std::invalid_argument);
    const double largest = std::numeric_limits<double>::max();
    EXPECT_THROW(
        LR::nacv_gap_resolution(
            0.0,
            1.0,
            largest,
            largest),
        std::overflow_error);
}

TEST(NacvResult, ResolvesSignedGapAndPreservesUnits)
{
    const ModuleBase::matrix force = example_force();
    for (double gap : {2.0, -2.0})
    {
        LR::NacvResult result = example_result(gap);
        LR::assemble_nacv_result(result, force, 1e-8);
        ASSERT_TRUE(result.gap_resolved);
        EXPECT_DOUBLE_EQ(result.numerator_etf(0, 0), 2.0);
        EXPECT_DOUBLE_EQ(result.numerator_raw(0, 0), 2.0 + gap * 0.25);
        EXPECT_DOUBLE_EQ(result.etf(0, 0), 2.0 / gap);
        EXPECT_DOUBLE_EQ(result.raw(0, 0), 2.0 / gap + 0.25);
        EXPECT_DOUBLE_EQ(result.raw(0, 1), -4.0 / gap - 0.5);
        EXPECT_DOUBLE_EQ(result.raw(0, 2), 0.0);
    }
}

TEST(NacvResult, UnresolvedGapRetainsNumeratorsAndClearsDividedVectors)
{
    LR::NacvResult result = example_result(2.0);
    const ModuleBase::matrix force = example_force();
    LR::assemble_nacv_result(result, force, 1e-8);
    for (double gap : {0.0, 1e-8, -1e-8})
    {
        result.gap_rydberg = gap;
        LR::assemble_nacv_result(result, force, 1e-8);
        EXPECT_FALSE(result.gap_resolved);
        EXPECT_DOUBLE_EQ(result.numerator_etf(0, 0), 2.0);
        EXPECT_DOUBLE_EQ(result.numerator_raw(0, 0), 2.0 + gap * 0.25);
        EXPECT_DOUBLE_EQ(result.raw(0, 0), 0.0);
        EXPECT_DOUBLE_EQ(result.etf(0, 1), 0.0);
        EXPECT_DOUBLE_EQ(result.ao_connection(0, 0), 0.25);
    }
}

TEST(NacvResult, ReversePreservesGapWeightedNumerators)
{
    LR::NacvResult result = example_result(2.0);
    const ModuleBase::matrix force = example_force();
    LR::assemble_nacv_result(result, force, 1e-8);
    const LR::NacvResult original = result;
    LR::reverse_nac_pair(result);
    EXPECT_TRUE(result.gap_resolved);
    EXPECT_DOUBLE_EQ(result.gap_rydberg, -original.gap_rydberg);
    EXPECT_DOUBLE_EQ(result.raw(0, 0), -original.raw(0, 0));
    EXPECT_DOUBLE_EQ(result.etf(0, 1), -original.etf(0, 1));
    EXPECT_DOUBLE_EQ(result.ao_connection(0, 0), -original.ao_connection(0, 0));
    EXPECT_DOUBLE_EQ(result.numerator_raw(0, 0), original.numerator_raw(0, 0));
    EXPECT_DOUBLE_EQ(result.numerator_etf(0, 1), original.numerator_etf(0, 1));
    EXPECT_DOUBLE_EQ(result.raw(0, 0) * result.gap_rydberg, result.numerator_raw(0, 0));
    LR::reverse_nac_pair(result);
    EXPECT_DOUBLE_EQ(result.raw(0, 0), original.raw(0, 0));
}

TEST(NacvResult, RejectsInvalidShapeGapAndResolutionInRelease)
{
    LR::NacvResult result = example_result(2.0);
    const ModuleBase::matrix force = example_force();
    ModuleBase::matrix invalid_force(1, 2);
    EXPECT_THROW(LR::assemble_nacv_result(result, invalid_force, 1e-8), std::invalid_argument);
    EXPECT_THROW(LR::assemble_nacv_result(result, force, -1.0), std::invalid_argument);
    const double nan = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(LR::assemble_nacv_result(result, force, nan), std::invalid_argument);
    result.gap_rydberg = nan;
    EXPECT_THROW(LR::assemble_nacv_result(result, force, 1e-8), std::invalid_argument);
}

TEST(NacvResult, RejectsNonfiniteConnectionAndForce)
{
    LR::NacvResult result = example_result(2.0);
    ModuleBase::matrix force = example_force();
    result.ao_connection(0, 0) = std::numeric_limits<double>::infinity();
    EXPECT_THROW(LR::assemble_nacv_result(result, force, 1e-8), std::runtime_error);
    result.ao_connection(0, 0) = 0.0;
    force(0, 0) = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(LR::assemble_nacv_result(result, force, 1e-8), std::runtime_error);
}

TEST(NacvResult, OutputConvertsNumeratorsToHartreeAndRestoresStream)
{
    LR::NacvResult result = example_result(2.0);
    const ModuleBase::matrix force = example_force();
    LR::assemble_nacv_result(result, force, 1e-8);
    std::ostringstream output;
    output << std::fixed << std::setprecision(3);
    const std::ios::fmtflags flags = output.flags();
    LR::print_nac_result(
        output,
        result,
        1,
        0,
        1,
        "singlet");
    EXPECT_EQ(output.flags(), flags);
    EXPECT_EQ(output.precision(), 3);
    std::istringstream input(output.str());
    std::string tag;
    std::string spin;
    int step = 0;
    int bra = 0;
    int ket = 0;
    bool resolved = false;
    double gap = 0.0;
    input >> tag >> step >> bra >> ket >> spin >> resolved >> gap;
    EXPECT_EQ(tag, "NACV_BEGIN");
    EXPECT_EQ(step, 1);
    EXPECT_EQ(ket, 1);
    EXPECT_TRUE(resolved);
    EXPECT_DOUBLE_EQ(gap, 1.0);
    std::string ignored;
    std::getline(input, ignored);
    std::getline(input, ignored);
    int atom = -1;
    int axis = -1;
    double raw = 0.0;
    double etf = 0.0;
    double connection = 0.0;
    double numerator_raw = 0.0;
    double numerator_etf = 0.0;
    input >> tag >> atom >> axis >> raw >> etf >> connection >> numerator_raw >> numerator_etf;
    ASSERT_TRUE(input.good());
    EXPECT_EQ(tag, "NACV_COMPONENT");
    EXPECT_EQ(atom, 0);
    EXPECT_EQ(axis, 0);
    EXPECT_DOUBLE_EQ(raw, 1.25);
    EXPECT_DOUBLE_EQ(etf, 1.0);
    EXPECT_DOUBLE_EQ(connection, 0.25);
    EXPECT_DOUBLE_EQ(numerator_raw, 1.25);
    EXPECT_DOUBLE_EQ(numerator_etf, 1.0);
}

TEST(NacvResult, OutputRejectsInvalidIndicesAndReportsStreamFailure)
{
    LR::NacvResult result = example_result(2.0);
    const ModuleBase::matrix force = example_force();
    LR::assemble_nacv_result(result, force, 1e-8);
    std::ostringstream output;
    EXPECT_THROW(
        LR::print_nac_result(
            output,
            result,
            0,
            0,
            1,
            "singlet"),
        std::invalid_argument);
    EXPECT_THROW(
        LR::print_nac_result(
            output,
            result,
            1,
            1,
            1,
            "singlet"),
        std::invalid_argument);
    output.setstate(std::ios::badbit);
    EXPECT_THROW(
        LR::print_nac_result(
            output,
            result,
            1,
            0,
            1,
            "singlet"),
        std::runtime_error);
}

TEST(NacvResult, EigenpairDiagnosticsKeepPhysicalUnitsAndStreamFormat)
{
    const double energies[] = {2.0, 4.0};
    const std::vector<double> residuals{2e-10, 4e-10};
    std::ostringstream output;
    output << std::fixed << std::setprecision(4);
    const std::ios::fmtflags flags = output.flags();
    LR::print_nac_eigenpairs(
        output,
        3,
        "singlet",
        energies,
        residuals);
    EXPECT_EQ(output.flags(), flags);
    EXPECT_EQ(output.precision(), 4);
    std::istringstream input(output.str());
    for (int state = 1; state <= 2; ++state)
    {
        std::string tag;
        std::string spin;
        int step = 0;
        int index = 0;
        double energy = 0.0;
        double residual = 0.0;
        input >> tag >> step >> spin >> index >> energy >> residual;
        ASSERT_FALSE(input.fail());
        EXPECT_EQ(tag, "LR_EIGENPAIR");
        EXPECT_EQ(step, 3);
        EXPECT_EQ(index, state);
        EXPECT_DOUBLE_EQ(energy, state);
        EXPECT_DOUBLE_EQ(residual, state * 1e-10);
    }
    const std::vector<double> invalid_residuals{-1.0};
    EXPECT_THROW(
        LR::print_nac_eigenpairs(
            output,
            3,
            "singlet",
            energies,
            invalid_residuals),
        std::invalid_argument);
}

TEST(NacvResult, AcceptsUnequalIntegerUnrestrictedSpinSpaces)
{
    ModuleBase::matrix weights(2, 8);
    const std::vector<int> occupied{4, 2};
    for (int spin = 0; spin < 2; ++spin)
    {
        for (int band = 0; band < occupied[spin]; ++band)
        {
            weights(spin, band) = 1.0;
        }
    }
    EXPECT_NO_THROW(LR::validate_nac_occupations(weights, occupied, 8));
    EXPECT_THROW(
        LR::validate_nac_occupations(
            weights,
            4,
            8,
            2),
        std::invalid_argument);
    weights(1, 1) = 0.5;
    EXPECT_THROW(LR::validate_nac_occupations(weights, occupied, 8), std::invalid_argument);
    weights(1, 1) = 1.0;
    weights(0, 7) = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(LR::validate_nac_occupations(weights, occupied, 8), std::invalid_argument);
}

TEST(NacvResult, RejectsMissingUnrestrictedSpinWindows)
{
    ModuleBase::matrix weights(2, 8);
    EXPECT_THROW(LR::validate_nac_occupations(weights, std::vector<int>{4}, 8), std::invalid_argument);
    EXPECT_THROW(LR::validate_nac_occupations(weights, std::vector<int>{4, 0}, 8), std::invalid_argument);
    EXPECT_THROW(LR::validate_nac_occupations(weights, std::vector<int>{4, 8}, 8), std::invalid_argument);
    EXPECT_THROW(LR::validate_nac_occupations(weights, std::vector<int>{4, 2}, 7), std::invalid_argument);
}
