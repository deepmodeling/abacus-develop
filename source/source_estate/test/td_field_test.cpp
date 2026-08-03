#include "source_estate/module_pot/td_field.h"

#include "source_base/constants.h"
#include "source_io/module_parameter/input_parameter.h"

#include "gtest/gtest.h"

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

    const double amplitude_conversion = ModuleBase::BOHR_TO_A / ModuleBase::Ry_to_eV;
    const double total_field = 5.0 * amplitude_conversion;
    std::shared_ptr<elecstate::TDFieldManager> manager = elecstate::create_td_field_manager(input);

    ASSERT_EQ(manager->fields().size(), 2U);
    EXPECT_EQ(manager->fields()[0].direction(), 0);
    EXPECT_EQ(manager->fields()[1].direction(), 0);
    EXPECT_EQ(manager->fields()[0].subdivisions(), 2);
    EXPECT_EQ(manager->fields()[1].subdivisions(), 2);

    manager->advance_vector_gauge();
    ASSERT_EQ(manager->field_values().size(), 2U);
    EXPECT_NEAR(manager->field_values()[0], 2.0 * amplitude_conversion, 1.0e-14);
    EXPECT_NEAR(manager->field_values()[1], 3.0 * amplitude_conversion, 1.0e-14);
    EXPECT_NEAR(manager->electric_field()[0], total_field, 1.0e-14);
    EXPECT_NEAR(manager->total_electric_field()[0], total_field, 1.0e-14);
    EXPECT_NEAR(manager->vector_potential_laststep()[0], -total_field, 1.0e-14);
    EXPECT_NEAR(manager->vector_potential()[0], -0.5 * total_field, 1.0e-14);

    manager->advance_vector_gauge();
    EXPECT_NEAR(manager->vector_potential()[0], -1.5 * total_field, 1.0e-14);

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

} // namespace
