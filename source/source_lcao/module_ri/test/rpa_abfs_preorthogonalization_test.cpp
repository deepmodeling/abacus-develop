#include "source_lcao/module_ri/rpa_abfs_preorthogonalization.h"

#include "source_io/module_parameter/input_parameter.h"

#include <cmath>
#include <gtest/gtest.h>
#include <limits>
#include <string>
#include <vector>

namespace
{

constexpr int kNr = 101;
constexpr double kDr = 0.05;

Numerical_Orbital_Lm make_orbital(const int l,
                                  const int chi,
                                  const std::vector<double>& radial,
                                  const int atom_type = 0)
{
    std::vector<double> r(kNr, 0.0);
    std::vector<double> rab(kNr, kDr);
    for (int ir = 0; ir != kNr; ++ir)
    {
        r[ir] = ir * kDr;
    }

    Numerical_Orbital_Lm orbital;
    orbital.set_orbital_info("N",
                             atom_type,
                             l,
                             chi,
                             kNr,
                             rab.data(),
                             r.data(),
                             Numerical_Orbital_Lm::Psi_Type::Psi,
                             radial.data(),
                             101,
                             0.05,
                             kDr,
                             false,
                             true,
                             false);
    return orbital;
}

std::vector<double> radial_exp(const double power)
{
    std::vector<double> values(kNr, 0.0);
    for (int ir = 0; ir != kNr; ++ir)
    {
        const double r = ir * kDr;
        values[ir] = std::pow(r, power) * std::exp(-r);
    }
    return values;
}

ModuleRI::RpaAbfsOrbitalSet make_set(const int l, const std::vector<std::vector<double>>& radials)
{
    ModuleRI::RpaAbfsOrbitalSet abfs(1);
    abfs[0].resize(l + 1);
    for (std::size_t i = 0; i != radials.size(); ++i)
    {
        abfs[0][l].push_back(make_orbital(l, static_cast<int>(i), radials[i]));
    }
    return abfs;
}

void expect_same_orbital(const Numerical_Orbital_Lm& lhs, const Numerical_Orbital_Lm& rhs)
{
    EXPECT_EQ(lhs.getLabel(), rhs.getLabel());
    EXPECT_EQ(lhs.getType(), rhs.getType());
    EXPECT_EQ(lhs.getL(), rhs.getL());
    EXPECT_EQ(lhs.getChi(), rhs.getChi());
    EXPECT_EQ(lhs.getNr(), rhs.getNr());
    EXPECT_EQ(lhs.getNk(), rhs.getNk());
    for (int ir = 0; ir != lhs.getNr(); ++ir)
    {
        EXPECT_DOUBLE_EQ(lhs.getRadial(ir), rhs.getRadial(ir));
        EXPECT_DOUBLE_EQ(lhs.getRab(ir), rhs.getRab(ir));
        EXPECT_DOUBLE_EQ(lhs.getPsi(ir), rhs.getPsi(ir));
    }
}

} // namespace

TEST(RpaAbfsPreorthogonalization, NoneLeavesOrbitalsExactlyUnchanged)
{
    auto abfs = make_set(0, {radial_exp(0.0), radial_exp(1.0)});
    const auto original = abfs;

    const auto report = ModuleRI::apply_rpa_abfs_preorthogonalization(abfs, "none", 1.0e-2, false);

    ASSERT_EQ(abfs.size(), original.size());
    ASSERT_EQ(abfs[0][0].size(), original[0][0].size());
    for (std::size_t i = 0; i != abfs[0][0].size(); ++i)
    {
        expect_same_orbital(abfs[0][0][i], original[0][0][i]);
    }
    EXPECT_EQ(report.mode, "none");
    EXPECT_EQ(report.input_expanded_size, 2);
    EXPECT_EQ(report.output_expanded_size, 2);
}

TEST(RpaAbfsPreorthogonalization, RejectsDuplicateRadial)
{
    const auto radial = radial_exp(0.0);
    auto abfs = make_set(0, {radial, radial});

    const auto report = ModuleRI::apply_rpa_abfs_preorthogonalization(abfs, "onsite_coulomb", 1.0e-2, false);

    ASSERT_EQ(abfs[0][0].size(), 1U);
    ASSERT_EQ(report.channels.size(), 1U);
    EXPECT_EQ(report.channels[0].input_count, 2);
    EXPECT_EQ(report.channels[0].output_count, 1);
    EXPECT_EQ(report.channels[0].rejected_indices, std::vector<int>({1}));
}

TEST(RpaAbfsPreorthogonalization, RetainedRadialsAreCoulombOrthonormal)
{
    auto abfs = make_set(0, {radial_exp(0.0), radial_exp(2.0)});

    const auto report = ModuleRI::apply_rpa_abfs_preorthogonalization(abfs, "onsite_coulomb", 1.0e-2, false);
    const auto metric = ModuleRI::compute_onsite_coulomb_metric(abfs[0][0]);

    ASSERT_EQ(metric.size(), 2U);
    for (std::size_t i = 0; i != metric.size(); ++i)
    {
        for (std::size_t j = 0; j != metric.size(); ++j)
        {
            EXPECT_NEAR(metric[i][j], i == j ? 1.0 : 0.0, 1.0e-10);
        }
    }
    ASSERT_EQ(report.channels.size(), 1U);
    EXPECT_LE(report.channels[0].maximum_identity_error, 1.0e-10);
}

TEST(RpaAbfsPreorthogonalization, RejectsWholeAngularMultiplet)
{
    const auto radial = radial_exp(2.0);
    auto abfs = make_set(2, {radial, radial});

    const auto report = ModuleRI::apply_rpa_abfs_preorthogonalization(abfs, "onsite_coulomb", 1.0e-2, false);

    EXPECT_EQ(report.input_expanded_size, 10);
    EXPECT_EQ(report.output_expanded_size, 5);
    ASSERT_EQ(abfs[0][2].size(), 1U);
}

TEST(RpaAbfsPreorthogonalization, RejectsInvalidModeThresholdAndNorm)
{
    auto abfs = make_set(0, {radial_exp(0.0)});
    EXPECT_THROW(ModuleRI::apply_rpa_abfs_preorthogonalization(abfs, "unknown", 1.0e-2, false), std::invalid_argument);
    EXPECT_THROW(ModuleRI::apply_rpa_abfs_preorthogonalization(abfs, "onsite_coulomb", 0.0, false),
                 std::invalid_argument);
    EXPECT_THROW(ModuleRI::apply_rpa_abfs_preorthogonalization(abfs, "onsite_coulomb", -1.0, false),
                 std::invalid_argument);
    EXPECT_THROW(ModuleRI::apply_rpa_abfs_preorthogonalization(abfs,
                                                               "onsite_coulomb",
                                                               std::numeric_limits<double>::infinity(),
                                                               false),
                 std::invalid_argument);

    auto zero = make_set(0, {std::vector<double>(kNr, 0.0)});
    EXPECT_THROW(ModuleRI::apply_rpa_abfs_preorthogonalization(zero, "onsite_coulomb", 1.0e-2, false),
                 std::invalid_argument);
}

TEST(RpaAbfsPreorthogonalization, FormatsCompleteProvenance)
{
    const auto radial = radial_exp(0.0);
    auto abfs = make_set(0, {radial, radial});
    const auto report = ModuleRI::apply_rpa_abfs_preorthogonalization(abfs, "onsite_coulomb", 1.0e-2, false);

    const std::string text = ModuleRI::format_rpa_abfs_preorth_report(report);
    EXPECT_NE(text.find("mode=onsite_coulomb"), std::string::npos);
    EXPECT_NE(text.find("input_expanded_size=2"), std::string::npos);
    EXPECT_NE(text.find("output_expanded_size=1"), std::string::npos);
    EXPECT_NE(text.find("rejected_indices=1"), std::string::npos);
    EXPECT_NE(text.find("minimum_residual_norm="), std::string::npos);
    EXPECT_NE(text.find("maximum_identity_error="), std::string::npos);
}

TEST(RpaAbfsPreorthogonalization, InputAdapterUsesOnlyPublicInputs)
{
    const auto radial = radial_exp(0.0);
    Input_para input;
    input.rpa_abfs_preorth = "none";
    input.rpa_abfs_preorth_threshold = 1.0e-2;
    auto legacy = make_set(0, {radial, radial});
    ModuleRI::finalize_rpa_abfs_from_input(legacy, input, false);
    EXPECT_EQ(legacy[0][0].size(), 2U);

    input.rpa_abfs_preorth = "onsite_coulomb";
    auto stable = make_set(0, {radial, radial});
    ModuleRI::finalize_rpa_abfs_from_input(stable, input, false);
    EXPECT_EQ(stable[0][0].size(), 1U);
}
