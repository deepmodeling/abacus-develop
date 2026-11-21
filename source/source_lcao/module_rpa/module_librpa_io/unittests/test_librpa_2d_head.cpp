#include "../librpa_2d_coulomb_head.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

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

void require_near(const double actual, const double expected, const double tolerance,
                  const char* message)
{
    require(std::abs(actual - expected) <= tolerance, message);
}

void test_mos2_rotated_shrink_basis_coefficient()
{
    constexpr double area_parallel_bohr2 = 31.347506240476704;
    const std::vector<std::vector<double>> moments{{11.09877004933396}, {6.073507331002011}};
    const std::vector<int> atoms_per_type{1, 2};
    const auto normalization = RpaLriDetail::strict_2d_coulomb_head_normalization(
        area_parallel_bohr2, moments, atoms_per_type);

    require_near(normalization.multipole_norm_squared, 196.9576792074629, 2.0e-13,
                 "MoS2 shrink-basis monopole norm is wrong");
    require_near(normalization.raw_head_coefficient, 496.0890637033998, 5.0e-13,
                 "MoS2 raw 2D Coulomb head coefficient is wrong");
    require_near(normalization.sheet_to_raw_scale, 8.885664111490273, 2.0e-14,
                 "MoS2 sheet-to-raw scale is wrong");
}

void test_invalid_inputs_are_rejected()
{
    for (const double area : {0.0, -1.0, std::numeric_limits<double>::infinity()})
    {
        bool rejected = false;
        try
        {
            (void)RpaLriDetail::strict_2d_coulomb_head_normalization(
                area, {{1.0}}, {1});
        }
        catch (const std::invalid_argument&)
        {
            rejected = true;
        }
        require(rejected, "invalid in-plane area was accepted");
    }
}

void test_sidecar_uses_librpa_key_value_syntax()
{
    const auto normalization = RpaLriDetail::strict_2d_coulomb_head_normalization(
        31.347506240476704, {{11.09877004933396}, {6.073507331002011}}, {1, 2});
    const std::string sidecar =
        RpaLriDetail::format_strict_2d_coulomb_head_sidecar(normalization);

    require(sidecar.find("version = 1\n") != std::string::npos,
            "sidecar version is not in LibRPA key-value syntax");
    require(sidecar.find("strict_2d_coulomb_head_coefficient = ") != std::string::npos,
            "sidecar coefficient is not in LibRPA key-value syntax");
    require(sidecar.find("strict_2d_sheet_to_raw_scale = ") != std::string::npos,
            "sidecar scale is not in LibRPA key-value syntax");
}
} // namespace

int main()
{
    test_mos2_rotated_shrink_basis_coefficient();
    test_invalid_inputs_are_rejected();
    test_sidecar_uses_librpa_key_value_syntax();
    std::cout << "LibRPA strict 2D Coulomb head normalization tests passed\n";
    return 0;
}
