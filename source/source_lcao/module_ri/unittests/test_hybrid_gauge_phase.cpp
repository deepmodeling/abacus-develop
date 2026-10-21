#include "source_lcao/module_ri/hybrid_gauge_phase.h"
#include "source_lcao/module_ri/mix_dmk_2d.h"
#include "source_base/module_mixing/pulay_mixing.h"
#include "source_base/module_mixing/broyden_mixing.h"

#include "gtest/gtest.h"

#include <cmath>
#include <limits>

namespace
{
void check_borrowed_charge_mixer(Base_Mixing::Mixing& charge_mixer)
{
    const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
    const auto key = std::make_pair(1, std::array<int, 3>{{1, 0, 0}});
    const std::complex<double> physical(0.8, -0.1);
    RI_2D_Comm::HybridGaugeComplexDensity density(1);
    density[0][0][key] = RI::Tensor<std::complex<double>>({1, 1});
    Mix_DMk_2D<std::complex<double>> mixer;
    mixer.set_mixing(&charge_mixer);
    RI_2D_Comm::HybridDensityLayout layout;
    Base_Mixing::Mixing_Data charge_history;
    charge_mixer.init_mixing_data(charge_history, 2, sizeof(double));
    const auto inner_product = [](double* first, double* second)
    {
        return first[0] * second[0] + first[1] * second[1];
    };
    for (int iteration = 1; iteration <= 6; ++iteration)
    {
        // Also exercise a reset inside the same time step. Source time remains
        // the endpoint even though both charge and EXX histories are restarted.
        const bool restart = iteration == 1 || iteration == 4;
        if (restart)
        {
            charge_mixer.reset();
            charge_mixer.init_mixing_data(charge_history, 2, sizeof(double));
        }
        const auto time = RI_2D_Comm::hybrid_density_time(2, iteration);
        const double left_value = 0.1;
        const double right_value = 0.2 * iteration;
        const double source_value = time == RI_2D_Comm::HybridDensityTime::Left ? left_value : right_value;
        const ModuleBase::Vector3<double> source_A(source_value, 0.0, 0.0);
        const double raw_source_value = iteration == 1 ? left_value : right_value;
        density[0][0][key](0, 0) = physical * std::exp(std::complex<double>(0.0, raw_source_value));
        RI_2D_Comm::mix_hybrid_density(lattice, 1.0, source_A, restart, mixer, layout, density);
        EXPECT_NEAR(std::abs(density[0][0][key](0, 0) - physical), 0.0, 1e-12);
        // Charge mixing supplies the coefficients consumed by EXX next iteration.
        const double input[2] = {0.2, 0.4};
        const double output[2] = {0.3 + 0.1 * iteration, 0.8 - 0.01 * iteration * iteration};
        charge_mixer.push_data(charge_history, input, output, nullptr, true);
        charge_mixer.cal_coef(charge_history, inner_product);
    }
}
} // namespace

TEST(HybridGaugePhase, BorrowsPulayCoefficientsAcrossRestart)
{
    Base_Mixing::Pulay_Mixing charge_mixer(2, 0.4);
    check_borrowed_charge_mixer(charge_mixer);
}

TEST(HybridGaugePhase, BorrowsBroydenCoefficientsAcrossRestart)
{
    Base_Mixing::Broyden_Mixing charge_mixer(2, 0.4);
    check_borrowed_charge_mixer(charge_mixer);
}

TEST(HybridGaugePhase, NegativeRIConventionAndRoundTrip)
{
    const ModuleBase::Vector3<double> source_A(0.3, -0.1, 0.2);
    const ModuleBase::Vector3<double> translation(2.0, 3.0, -1.0);
    const std::complex<double> phase = RI_2D_Comm::hybrid_input_phase(source_A, translation);
    EXPECT_NEAR(phase.real(), std::cos(0.1), 1e-15);
    EXPECT_NEAR(phase.imag(), -std::sin(0.1), 1e-15);
    const ModuleBase::Vector3<double> reverse_A(-0.3, 0.1, -0.2);
    const std::complex<double> reverse = RI_2D_Comm::hybrid_input_phase(reverse_A, translation);
    EXPECT_NEAR(std::abs(phase * reverse - std::complex<double>(1.0, 0.0)), 0.0, 1e-15);
}

TEST(HybridGaugePhase, SelectsExplicitDensityTimeLayer)
{
    using RI_2D_Comm::HybridDensityTime;
    EXPECT_EQ(RI_2D_Comm::hybrid_density_time(0, 1), HybridDensityTime::Right);
    EXPECT_EQ(RI_2D_Comm::hybrid_density_time(1, 1), HybridDensityTime::Left);
    EXPECT_EQ(RI_2D_Comm::hybrid_density_time(1, 2), HybridDensityTime::Right);
    EXPECT_EQ(RI_2D_Comm::hybrid_density_time(2, 5), HybridDensityTime::Right);
}

TEST(HybridGaugePhase, ExplicitSourceFrameAndUnchangedCache)
{
    const std::array<int, 3> cell = {{1, -1, 0}};
    const auto key = std::make_pair(2, cell);
    RI_2D_Comm::HybridGaugeComplexDensity original(2);
    original[1][0][key] = RI::Tensor<std::complex<double>>({1, 2});
    original[1][0][key](0, 0) = std::complex<double>(0.7, 0.4);
    original[1][0][key](0, 1) = std::complex<double>(-0.2, 0.1);
    auto transformed = original;
    const ModuleBase::Matrix3 lattice(2.0, 0.5, 0.0,
                                     0.0, 3.0, 0.2,
                                     0.0, 0.0, 4.0);
    const ModuleBase::Vector3<double> source_A(0.2, 0.1, -0.3);
    RI_2D_Comm::apply_hybrid_input_phase(lattice, 2.0, source_A, transformed);
    // R=(a1-a2)*lat0=(4,-5,-0.4), A.R=0.42.
    const std::complex<double> expected_phase(std::cos(0.42), -std::sin(0.42));
    const std::complex<double> source(0.7, 0.4);
    EXPECT_NEAR(std::abs(transformed[1][0][key](0, 0) - source * expected_phase), 0.0, 1e-15);
    EXPECT_EQ(original[1][0][key](0, 0), source);
    EXPECT_NE(original[1][0][key].ptr(), transformed[1][0][key].ptr());
    EXPECT_TRUE(transformed[0].empty());
}

TEST(HybridGaugePhase, EqualPhysicalDensityAcrossTwoSourceFrames)
{
    const ModuleBase::Vector3<double> translation(1.1, -0.5, 0.3);
    const ModuleBase::Vector3<double> left_A(0.1, 0.2, -0.1);
    const ModuleBase::Vector3<double> right_A(0.4, -0.1, 0.2);
    const std::complex<double> physical(0.7, -0.2);
    const auto left_phase = RI_2D_Comm::hybrid_input_phase(left_A, translation);
    const auto right_phase = RI_2D_Comm::hybrid_input_phase(right_A, translation);
    const auto left_raw = physical / left_phase;
    const auto right_raw = physical / right_phase;
    const auto common_mixed = 0.4 * left_raw * left_phase + 0.6 * right_raw * right_phase;
    EXPECT_NEAR(std::abs(common_mixed - physical), 0.0, 1e-15);
    const auto wrong_stage_mixed = (0.4 * left_raw + 0.6 * right_raw) * right_phase;
    EXPECT_GT(std::abs(wrong_stage_mixed - physical), 0.01);
}

TEST(HybridGaugePhase, ImageRelabelCommutesWithOrdinaryFrameConversion)
{
    // Algebraic oracle for the inspected conj-first DMK / RI transpose
    // conventions. It does not enable transport of production snapshots.
    const ModuleBase::Matrix3 lattice(2.0, 0.5, 0.0,
                                     0.0, 3.0, 0.2,
                                     0.1, 0.0, 4.0);
    const double lat0 = 1.7;
    const ModuleBase::Vector3<double> a(0.2, -0.1, 0.07);
    const ModuleBase::Vector3<double> k_phase(0.31, 0.19, -0.23);
    const ModuleBase::Vector3<int> bi(-1, 0, 1);
    const ModuleBase::Vector3<int> bj(0, 1, -1);
    const ModuleBase::Vector3<int> r(1, -2, 0);
    const auto delta = bi - bj;
    const auto rn = r + delta;
    const auto cartesian = [&](const ModuleBase::Vector3<int>& image)
    {
        return ModuleBase::Vector3<double>(image.x, image.y, image.z) * lattice * lat0;
    };
    const auto coefficient_phase = [&](const ModuleBase::Vector3<int>& image)
    {
        const double argument = k_phase.x * image.x + k_phase.y * image.y
            + k_phase.z * image.z + a * cartesian(image);
        return std::exp(std::complex<double>(0.0, argument));
    };
    const std::complex<double> ci(0.7, -0.2);
    const std::complex<double> cj(-0.1, 0.6);
    const std::complex<double> ui = coefficient_phase(bi);
    const std::complex<double> uj = coefficient_phase(bj);
    // dm_from_psi stores conj(C_i)*C_j; RI's row-major tensor view
    // transposes its column-major buffer before Fourier conversion.
    const auto stored_old = std::conj(ci) * cj;
    const auto stored_new = std::conj(ui * ci) * (uj * cj);
    EXPECT_NEAR(std::abs(stored_new - std::conj(ui) * stored_old * uj), 0.0, 1e-14);
    EXPECT_GT(std::abs(stored_new - ui * stored_old * std::conj(uj)), 1e-3);

    const double k_delta = k_phase.x * delta.x + k_phase.y * delta.y + k_phase.z * delta.z;
    const auto transposed_new_ratio = ui * std::conj(uj);
    const auto fourier_relabel = std::exp(std::complex<double>(0.0, -k_delta));
    const auto raw_ri_ratio = transposed_new_ratio * fourier_relabel;
    EXPECT_NEAR(std::abs(raw_ri_ratio - std::exp(std::complex<double>(0.0, a * cartesian(delta)))),
                0.0, 1e-14);

    const auto old_key = std::make_pair(1, std::array<int, 3>{{r.x, r.y, r.z}});
    const auto new_key = std::make_pair(1, std::array<int, 3>{{rn.x, rn.y, rn.z}});
    RI_2D_Comm::HybridGaugeComplexDensity old_density(1);
    RI_2D_Comm::HybridGaugeComplexDensity new_density(1);
    old_density[0][0][old_key] = RI::Tensor<std::complex<double>>({1, 1});
    new_density[0][0][new_key] = RI::Tensor<std::complex<double>>({1, 1});
    const std::complex<double> raw(0.8, -0.3);
    old_density[0][0][old_key](0, 0) = raw;
    new_density[0][0][new_key](0, 0) = raw * raw_ri_ratio;
    const auto old_source = old_density;
    const auto new_source = new_density;
    RI_2D_Comm::apply_hybrid_input_phase(lattice, lat0, a, old_density);
    RI_2D_Comm::apply_hybrid_input_phase(lattice, lat0, a, new_density);
    EXPECT_NEAR(std::abs(old_density[0][0][old_key](0, 0) - new_density[0][0][new_key](0, 0)),
                0.0, 1e-14);
    EXPECT_EQ(old_source[0].at(0).at(old_key)(0, 0), raw);
    EXPECT_EQ(new_source[0].at(0).at(new_key)(0, 0), raw * raw_ri_ratio);
}

TEST(HybridGaugePhase, ProductionPlainMixerAndTransactionalTopologyGuard)
{
    const std::array<int, 3> cell = {{1, 0, 0}};
    const auto key = std::make_pair(1, cell);
    RI_2D_Comm::HybridGaugeComplexDensity density(1);
    density[0][0][key] = RI::Tensor<std::complex<double>>({1, 2});
    density[0][0][key](0, 0) = std::complex<double>(0.5, 0.2);
    density[0][0][key](0, 1) = std::complex<double>(0.3, -0.1);
    const auto accepted = density;
    const auto layout = RI_2D_Comm::hybrid_density_layout(density);
    auto first = RI_2D_Comm::flatten_hybrid_density(density);
    auto second = first;
    second[0][0] += std::complex<double>(0.4, -0.2);
    Mix_DMk_2D<std::complex<double>> mixer;
    mixer.set_nks(1);
    mixer.set_mixing_plain(0.25);
    mixer.mix(first, true);
    mixer.mix(second, false);
    const auto output = mixer.get_DMk_out();
    RI_2D_Comm::assign_flat_hybrid_density(output, layout, density);
    EXPECT_NEAR(std::abs(density[0][0][key](0, 0) - std::complex<double>(0.6, 0.15)), 0.0, 1e-15);
    EXPECT_EQ(accepted[0].at(0).at(key)(0, 0), first[0][0]);

    auto changed = density;
    const std::array<int, 3> other_cell = {{2, 0, 0}};
    const auto other_key = std::make_pair(1, other_cell);
    changed[0][0][other_key] = changed[0][0][key];
    changed[0][0].erase(key);
    EXPECT_THROW(RI_2D_Comm::assign_flat_hybrid_density(output, layout, changed), std::invalid_argument);
    EXPECT_EQ(changed[0][0][other_key](0, 0), density[0][0][key](0, 0));
    changed = density;
    changed[0][0][key] = RI::Tensor<std::complex<double>>({2, 1});
    EXPECT_THROW(RI_2D_Comm::assign_flat_hybrid_density(output, layout, changed), std::invalid_argument);
    const std::vector<const std::vector<std::complex<double>>*> invalid(1, nullptr);
    EXPECT_THROW(RI_2D_Comm::assign_flat_hybrid_density(invalid, layout, density), std::invalid_argument);
}

TEST(HybridGaugePhase, ProductionMixerReceivesCommonFrameNotRawSourceFrames)
{
    const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0,
                                     0.0, 1.0, 0.0,
                                     0.0, 0.0, 1.0);
    const ModuleBase::Vector3<double> left_A(0.2, 0.0, 0.0);
    const ModuleBase::Vector3<double> right_A(0.5, 0.0, 0.0);
    const std::array<int, 3> cell = {{1, 0, 0}};
    const auto key = std::make_pair(1, cell);
    const std::complex<double> physical(0.8, -0.1);
    RI_2D_Comm::HybridGaugeComplexDensity left(1);
    RI_2D_Comm::HybridGaugeComplexDensity right(1);
    left[0][0][key] = RI::Tensor<std::complex<double>>({1, 1});
    right[0][0][key] = RI::Tensor<std::complex<double>>({1, 1});
    left[0][0][key](0, 0) = physical * std::exp(std::complex<double>(0.0, 0.2));
    right[0][0][key](0, 0) = physical * std::exp(std::complex<double>(0.0, 0.5));
    Mix_DMk_2D<std::complex<double>> mixer;
    mixer.set_mixing_plain(0.6);
    RI_2D_Comm::HybridDensityLayout layout;
    RI_2D_Comm::mix_hybrid_density(lattice, 1.0, left_A, true, mixer, layout, left);
    RI_2D_Comm::mix_hybrid_density(lattice, 1.0, right_A, false, mixer, layout, right);
    EXPECT_NEAR(std::abs(right[0][0][key](0, 0) - physical), 0.0, 1e-15);

    const auto before = *mixer.get_DMk_out()[0];
    auto changed = right;
    const auto other_key = std::make_pair(1, std::array<int, 3>{{2, 0, 0}});
    changed[0][0][other_key] = changed[0][0][key];
    changed[0][0].erase(key);
    EXPECT_THROW(RI_2D_Comm::mix_hybrid_density(lattice, 1.0, right_A, false,
                                              mixer, layout, changed), std::invalid_argument);
    EXPECT_EQ(*mixer.get_DMk_out()[0], before);
    EXPECT_EQ(changed[0][0][other_key](0, 0), right[0][0][key](0, 0));
    const ModuleBase::Vector3<double> invalid_A(std::numeric_limits<double>::infinity(), 0.0, 0.0);
    EXPECT_THROW(RI_2D_Comm::mix_hybrid_density(lattice, 1.0, invalid_A, true,
                                              mixer, layout, changed), std::invalid_argument);
    EXPECT_EQ(*mixer.get_DMk_out()[0], before);
    RI_2D_Comm::mix_hybrid_density(lattice, 1.0, right_A, true, mixer, layout, changed);
    EXPECT_EQ(layout, RI_2D_Comm::hybrid_density_layout(changed));
}
