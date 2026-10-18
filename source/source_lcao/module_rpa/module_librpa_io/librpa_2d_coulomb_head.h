#pragma once

#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Strict-2D Coulomb head normalization for the reader-v1 output.
//
// In a strict 2D (slab) cell the q_parallel -> 0 limit of the full Coulomb
// does not vanish: the auxiliary monopoles form a charged sheet whose
// interaction is the "head" of the separable Coulomb representation. With
// s_multipoles the per-type s-channel (l=0) moments, the head coefficient is
//     8 pi^2 * ||s||^2 / A_parallel,   ||s||^2 = sum_I nat_I * sum_n s_{I,n}^2,
// and the sheet-to-raw scale sqrt(coefficient / 2 pi) converts between the
// sheet representation LibRPA consumes and the raw coefficient. Both are
// written to the librpa_2d_coulomb_head.txt sidecar next to the V files.
namespace RpaLriDetail
{
struct Strict2dCoulombHeadNormalization
{
    double area_parallel_bohr2 = 0.0;
    double multipole_norm_squared = 0.0;
    double raw_head_coefficient = 0.0;
    double sheet_to_raw_scale = 0.0;
};

inline Strict2dCoulombHeadNormalization strict_2d_coulomb_head_normalization(
    const double area_parallel_bohr2,
    const std::vector<std::vector<double>>& s_multipoles_by_type,
    const std::vector<int>& atoms_per_type)
{
    if (!(area_parallel_bohr2 > 0.0) || !std::isfinite(area_parallel_bohr2)
        || s_multipoles_by_type.size() != atoms_per_type.size()
        || s_multipoles_by_type.empty())
    {
        throw std::invalid_argument("invalid strict 2D Coulomb head normalization inputs");
    }

    double multipole_norm_squared = 0.0;
    for (std::size_t it = 0; it != s_multipoles_by_type.size(); ++it)
    {
        if (atoms_per_type[it] < 0)
        {
            throw std::invalid_argument("negative atom count in strict 2D Coulomb head normalization");
        }
        double type_norm_squared = 0.0;
        for (const double moment : s_multipoles_by_type[it])
        {
            if (!std::isfinite(moment))
            {
                throw std::invalid_argument("non-finite strict 2D auxiliary multipole");
            }
            type_norm_squared += moment * moment;
        }
        multipole_norm_squared += atoms_per_type[it] * type_norm_squared;
    }
    if (!(multipole_norm_squared > 0.0) || !std::isfinite(multipole_norm_squared))
    {
        throw std::invalid_argument("strict 2D auxiliary basis has no finite monopole channel");
    }

    constexpr double pi = 3.141592653589793238462643383279502884;
    const double raw_head_coefficient =
        8.0 * pi * pi * multipole_norm_squared / area_parallel_bohr2;
    return {area_parallel_bohr2,
            multipole_norm_squared,
            raw_head_coefficient,
            std::sqrt(raw_head_coefficient / (2.0 * pi))};
}

inline std::string format_strict_2d_coulomb_head_sidecar(
    const Strict2dCoulombHeadNormalization& normalization)
{
    std::ostringstream output;
    output << "# ABACUS reader-v1 strict 2D Coulomb head normalization\n"
           << "version = 1\n" << std::setprecision(17)
           << "area_parallel_bohr2 = " << normalization.area_parallel_bohr2 << '\n'
           << "multipole_norm_squared = " << normalization.multipole_norm_squared << '\n'
           << "strict_2d_coulomb_head_coefficient = "
           << normalization.raw_head_coefficient << '\n'
           << "strict_2d_sheet_to_raw_scale = " << normalization.sheet_to_raw_scale << '\n';
    return output.str();
}
} // namespace RpaLriDetail
