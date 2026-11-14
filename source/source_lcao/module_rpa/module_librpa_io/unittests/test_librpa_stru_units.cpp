#include "../librpa_stru_units.h"

#include "source_base/constants.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace
{
void require(const bool condition, const char* message)
{
    if (!condition)
    {
        std::cerr << "FAILED: " << message << '\n';
        std::exit(1);
    }
}

void require_near(const double actual, const double expected, const double tolerance, const char* message)
{
    require(std::abs(actual - expected) <= tolerance, message);
}

void test_general_lat0()
{
    const double lat0_bohr = 5.3;
    const auto scales = RpaLriDetail::librpa_stru_unit_scales(lat0_bohr);

    require(scales.real_space_bohr == lat0_bohr,
            "real-space stru_out scale must use the actual lattice constant");
    require_near(scales.reciprocal_space_bohr_inv,
                 ModuleBase::TWO_PI / lat0_bohr,
                 0.0,
                 "reciprocal-space stru_out scale must use 2*pi/lat0");
    require_near(scales.real_space_bohr * scales.reciprocal_space_bohr_inv,
                 ModuleBase::TWO_PI,
                 4.0 * std::numeric_limits<double>::epsilon(),
                 "real and reciprocal stru_out scales must be mutually consistent");
}

void test_exact_one_angstrom_compatibility()
{
    const double lat0_bohr = 1.0 / ModuleBase::BOHR_TO_A;
    const auto scales = RpaLriDetail::librpa_stru_unit_scales(lat0_bohr);

    require(scales.real_space_bohr == 1.0 / ModuleBase::BOHR_TO_A,
            "exact one-Angstrom real-space scale must match the legacy conversion");
    require(scales.reciprocal_space_bohr_inv == ModuleBase::TWO_PI * ModuleBase::BOHR_TO_A,
            "exact one-Angstrom reciprocal scale must match the legacy conversion");
}

void test_common_truncated_one_angstrom_constant()
{
    const double lat0_bohr = 1.8897162;
    const auto scales = RpaLriDetail::librpa_stru_unit_scales(lat0_bohr);
    const double legacy_real_scale = 1.0 / ModuleBase::BOHR_TO_A;
    const double legacy_reciprocal_scale = ModuleBase::TWO_PI * ModuleBase::BOHR_TO_A;
    const double real_relative_change = std::abs(scales.real_space_bohr / legacy_real_scale - 1.0);
    const double reciprocal_relative_change
        = std::abs(scales.reciprocal_space_bohr_inv / legacy_reciprocal_scale - 1.0);
    const double expected_real_relative_change
        = std::abs(lat0_bohr * ModuleBase::BOHR_TO_A - 1.0);
    const double expected_reciprocal_relative_change
        = std::abs(1.0 / (lat0_bohr * ModuleBase::BOHR_TO_A) - 1.0);

    require_near(real_relative_change,
                 expected_real_relative_change,
                 std::numeric_limits<double>::epsilon(),
                 "the common-case real-space change must equal the input conversion error");
    require_near(reciprocal_relative_change,
                 expected_reciprocal_relative_change,
                 std::numeric_limits<double>::epsilon(),
                 "the common-case reciprocal change must equal the input conversion error");
    require(real_relative_change < 5.7e-6,
            "the common truncated one-Angstrom constant must stay within its conversion error");
    require(reciprocal_relative_change < 5.7e-6,
            "the reciprocal common-case change must stay within its conversion error");
    require_near(scales.real_space_bohr * scales.reciprocal_space_bohr_inv,
                 ModuleBase::TWO_PI,
                 4.0 * std::numeric_limits<double>::epsilon(),
                 "the common truncated constant must still produce reciprocal geometry");

    const double tau_component = 0.573;
    const double legacy_direct_position = tau_component * lat0_bohr;
    const double new_direct_position = tau_component * scales.real_space_bohr;
    require(new_direct_position == legacy_direct_position,
            "Direct-coordinate atom output must remain exactly unchanged");
}

void test_invalid_lat0()
{
    for (const double value : {0.0, -1.0, std::numeric_limits<double>::infinity(),
                               std::numeric_limits<double>::quiet_NaN()})
    {
        bool rejected = false;
        try
        {
            (void)RpaLriDetail::librpa_stru_unit_scales(value);
        }
        catch (const std::invalid_argument&)
        {
            rejected = true;
        }
        require(rejected, "invalid lattice constants must be rejected");
    }
}
} // namespace

int main()
{
    test_general_lat0();
    test_exact_one_angstrom_compatibility();
    test_common_truncated_one_angstrom_constant();
    test_invalid_lat0();
    std::cout << "LibRPA stru_out unit tests passed\n";
    return 0;
}
