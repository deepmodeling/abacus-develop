#include "source_base/constants.h"
#include "source_estate/module_pot/td_field_manager.h"
#include "source_estate/module_pot/td_field_profiles.h"
#include "source_io/module_parameter/input_parameter.h"

#include "gtest/gtest.h"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>

namespace
{

TEST(TDFieldProfileTest, RepresentativeSamples)
{
    elecstate::TDGaussianProfile gaussian(0.5, 0.25, 2.0, 0.0, 2.0, 1.0);
    elecstate::TDTrapezoidProfile trapezoid(0.4, 0.4, 1.0, 2.0, 4.0, 3.0);
    elecstate::TDTrigonometricProfile trigonometric(0.7, 0.3, 0.2, 0.4, 4.0);
    elecstate::TDHeavisideProfile heaviside(3.0, 5.0);
    elecstate::TDSupersineProfile supersine(-0.4, 0.3, 0.75, 0, 8, 6.0, 1.0);

    struct ProfileCase
    {
        const char* name;
        double value;
        double expected;
    };
    const std::vector<ProfileCase> cases
        = {{"Gaussian", gaussian.electric_field(elecstate::TDFieldSample(0, 3, 4, 0.75)), 1.5118014898049421},
           {"Trapezoid", trapezoid.electric_field(elecstate::TDFieldSample(2, 1, 2, 4.0)), -0.93633038223107046},
           {"Trigonometric", trigonometric.electric_field(elecstate::TDFieldSample(1, 1, 2, 1.25)), 0.93167894632639381},
           {"HeavisideBefore", heaviside.electric_field(elecstate::TDFieldSample(2, 1, 2, 2.5)), 5.0},
           {"HeavisideSwitch", heaviside.electric_field(elecstate::TDFieldSample(3, 0, 2, 3.0)), 0.0},
           {"Supersine", supersine.electric_field(elecstate::TDFieldSample(2, 0, 2, 2.0)), -3.6184740935466699}};

    for (const ProfileCase& profile_case: cases)
    {
        EXPECT_NEAR(profile_case.value, profile_case.expected, 1.0e-12) << profile_case.name;
    }
}

TEST(TDFieldManagerTest, MixedFieldsSumAndAccumulate)
{
    Input_para input;
    input.td_vext = true;
    input.td_stype = 2;
    input.td_tstart = 0;
    input.td_tend = 2;
    input.td_dt = ModuleBase::AU_to_FS;
    input.td_ttype = {2, 3};
    input.td_vext_dire = {1, 1};
    input.td_trigo_freq1 = {0.0};
    input.td_trigo_freq2 = {0.0};
    input.td_trigo_phase1 = {0.0};
    input.td_trigo_phase2 = {ModuleBase::PI / 2.0};
    input.td_trigo_amp = {2.0};
    input.td_heavi_t0 = {10.0};
    input.td_heavi_amp = {3.0};

    const double amplitude_conversion = ModuleBase::BOHR_TO_A / ModuleBase::Hartree_to_eV;
    const double total_field = 5.0 * amplitude_conversion;
    std::shared_ptr<elecstate::TDFieldManager> manager = elecstate::create_td_field_manager(input);

    ASSERT_EQ(manager->fields().size(), 2U);
    EXPECT_EQ(manager->fields()[0].direction(), 0);
    EXPECT_EQ(manager->fields()[1].direction(), 0);
    EXPECT_EQ(manager->fields()[0].subdivisions(), 2);
    EXPECT_EQ(manager->fields()[1].subdivisions(), 2);

    manager->prepare_sample(0);
    manager->prepare_interval(0, 1);
    ASSERT_EQ(manager->field_vals_ha().size(), 2U);
    EXPECT_NEAR(manager->field_vals_ha()[0], 2.0 * amplitude_conversion, 1.0e-14);
    EXPECT_NEAR(manager->field_vals_ha()[1], 3.0 * amplitude_conversion, 1.0e-14);
    EXPECT_NEAR(manager->efield_ha()[0], total_field, 1.0e-14);
    EXPECT_NEAR(manager->efield_ha()[0], total_field, 1.0e-14);
    EXPECT_NEAR(manager->A_right_ha()[0], -total_field, 1.0e-14);
    EXPECT_NEAR(manager->A_prop_ha()[0], -0.5 * total_field, 1.0e-14);

    manager->prepare_interval(1, 2);
    EXPECT_NEAR(manager->A_prop_ha()[0], -1.5 * total_field, 1.0e-14);

    Input_para oscillatory_input;
    oscillatory_input.td_vext = true;
    oscillatory_input.td_stype = 1;
    oscillatory_input.td_dt = 0.05;
    oscillatory_input.td_ttype = {0, 0};
    oscillatory_input.td_vext_dire = {1, 2};
    oscillatory_input.td_gauss_freq = {1.0, -1.0};
    oscillatory_input.td_gauss_phase = {0.0, 0.0};
    oscillatory_input.td_gauss_sigma = {1.0, 1.0};
    oscillatory_input.td_gauss_t0 = {0.0, 0.0};
    oscillatory_input.td_gauss_amp = {1.0, 1.0};

    std::shared_ptr<elecstate::TDFieldManager> oscillatory_manager = elecstate::create_td_field_manager(oscillatory_input);
    ASSERT_EQ(oscillatory_manager->fields().size(), 2U);
    EXPECT_EQ(oscillatory_manager->fields()[0].subdivisions(), 12);
    EXPECT_EQ(oscillatory_manager->fields()[1].subdivisions(), 12);
}

TEST(TDFieldManagerTest, SubintervalsPreserveFieldAndOriginalGrid)
{
    Input_para input;
    input.td_vext = true;
    input.td_stype = 1;
    input.td_dt = ModuleBase::AU_to_FS;
    input.td_tstart = 0;
    input.td_tend = 10;
    input.td_ttype = std::vector<int>{2};
    input.td_vext_dire = std::vector<int>{1};
    const double omega = 0.02;
    input.td_trigo_freq1 = {omega / (2 * ModuleBase::PI * ModuleBase::AU_to_FS)};
    input.td_trigo_freq2 = {0.0};
    input.td_trigo_phase1 = {0.0};
    input.td_trigo_phase2 = {ModuleBase::PI / 2};
    input.td_trigo_amp = {1.0};
    const double amp = ModuleBase::BOHR_TO_A / ModuleBase::Hartree_to_eV;

    for (bool from_file: {false, true})
    {
        SCOPED_TRACE(from_file);
        std::shared_ptr<elecstate::TDFieldManager> manager = elecstate::create_td_field_manager(input);
        if (from_file)
        {
            const std::vector<ModuleBase::Vector3<double>> samples = {{0, 0, 0}, {2, 0, 0}, {5, 0, 0}};
            manager->set_A_samples(samples);
        }
        manager->prepare_sample(0);
        manager->prepare_interval(0, 1);
        const elecstate::TDFieldState full = manager->field_state();
        const elecstate::TDFieldState left = manager->subinterval(full, 0, 2);
        manager->bind_field(left);
        const elecstate::TDFieldState right = manager->subinterval(full, 1, 2);
        EXPECT_DOUBLE_EQ(left.dt, 0.5);
        EXPECT_DOUBLE_EQ(right.dt, 0.5);
        EXPECT_DOUBLE_EQ(left.left.x, full.left.x);
        EXPECT_DOUBLE_EQ(left.right.x, right.left.x);
        EXPECT_DOUBLE_EQ(right.right.x, full.right.x);
        const double mid_a = from_file ? 1.0 : -amp * std::sin(omega * 0.5) / omega;
        const double mid_e = from_file ? -2.0 : amp * std::cos(omega * 0.5);
        EXPECT_NEAR(left.right.x, mid_a, 1e-11);
        EXPECT_NEAR(left.electric.x, mid_e, 1e-14);
        EXPECT_DOUBLE_EQ(left.propagation.x, (left.left.x + left.right.x) * 0.5);
        EXPECT_DOUBLE_EQ(right.propagation.x, (right.left.x + right.right.x) * 0.5);
        manager->bind_field(right);
        EXPECT_EQ(manager->current_step(), 1);
        EXPECT_DOUBLE_EQ(manager->interval_dt_ha(), 0.5);
        manager->bind_field(full);
        manager->prepare_interval(0, 1);
        EXPECT_EQ(manager->current_step(), 1);
        EXPECT_DOUBLE_EQ(manager->interval_dt_ha(), 1.0);
        EXPECT_DOUBLE_EQ(manager->A_right_ha().x, full.right.x);
        EXPECT_DOUBLE_EQ(manager->efield_ha().x, full.electric.x);
        manager->prepare_interval(1, 2);
        const double end_a = from_file ? 5.0 : -amp * std::sin(2 * omega) / omega;
        EXPECT_EQ(manager->current_step(), 2);
        EXPECT_DOUBLE_EQ(manager->A_left_ha().x, full.right.x);
        EXPECT_NEAR(manager->A_right_ha().x, end_a, 1e-10);
    }
}

TEST(TDFieldManagerTest, RestartRequiresCompleteState)
{
    Input_para input;
    input.td_stype = 1;
    input.td_dt = ModuleBase::AU_to_FS;
    std::shared_ptr<elecstate::TDFieldManager> manager = elecstate::create_td_field_manager(input);
    const std::string restart_prefix = "td_field_manager_test_";
    const std::string restart_path = restart_prefix + "Restart_td.txt";

    {
        std::ofstream output(restart_path.c_str());
        output << "# Hartree units\n\n7 6 1 1 0\n1 2 3 # A_left\n4 5 6 # A_right\n";
    }
    manager->read_restart(restart_prefix);
    EXPECT_EQ(manager->current_step(), 7);
    EXPECT_DOUBLE_EQ(manager->A_left_ha()[0], 1.0);
    EXPECT_DOUBLE_EQ(manager->A_left_ha()[1], 2.0);
    EXPECT_DOUBLE_EQ(manager->A_left_ha()[2], 3.0);
    EXPECT_DOUBLE_EQ(manager->A_right_ha()[0], 4.0);
    EXPECT_DOUBLE_EQ(manager->A_right_ha()[1], 5.0);
    EXPECT_DOUBLE_EQ(manager->A_right_ha()[2], 6.0);

    {
        std::ofstream output(restart_path.c_str(), std::ofstream::out);
        output << "7\n1 2 3\n4 5\n";
    }
    EXPECT_EXIT(manager->read_restart(restart_prefix), testing::ExitedWithCode(1), "");
    std::remove(restart_path.c_str());
}

} // namespace
