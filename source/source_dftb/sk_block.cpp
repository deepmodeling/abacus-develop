#include "source_dftb/sk_block.h"

#include <cmath>
#include <stdexcept>

namespace ModuleDFTB
{
namespace
{
constexpr std::size_t index(SkfChannel channel)
{
    return static_cast<std::size_t>(channel);
}

void fill_block(const SkfIntegralEvaluation& forward,
                const SkfIntegralEvaluation& reverse,
                const std::array<double, 3>& displacement,
                double distance,
                bool overlap,
                std::array<double, 16>& block,
                std::array<std::array<double, 16>, 3>& derivative)
{
    const auto& value_ab = overlap ? forward.overlap : forward.hamiltonian;
    const auto& deriv_ab = overlap ? forward.d_overlap_dr : forward.d_hamiltonian_dr;
    const auto& value_ba = overlap ? reverse.overlap : reverse.hamiltonian;
    const auto& deriv_ba = overlap ? reverse.d_overlap_dr : reverse.d_hamiltonian_dr;
    const std::size_t ss = index(SkfChannel::ss_sigma);
    const std::size_t sp = index(SkfChannel::sp_sigma);
    const std::size_t pp_sigma = index(SkfChannel::pp_sigma);
    const std::size_t pp_pi = index(SkfChannel::pp_pi);

    const double direction[3] = {displacement[0] / distance,
                                 displacement[1] / distance,
                                 displacement[2] / distance};
    // DFTB+ stores p orbitals in the order p_y, p_z, p_x.
    constexpr std::size_t cartesian_axis[3] = {1, 2, 0};

    block[0] = value_ab[ss];
    for (std::size_t cart = 0; cart < 3; ++cart)
    {
        derivative[cart][0] = deriv_ab[ss] * direction[cart];
    }

    for (std::size_t p = 0; p < 3; ++p)
    {
        const std::size_t axis = cartesian_axis[p];
        const std::size_t row = p + 1;
        const std::size_t col = p + 1;
        const double cosine = direction[axis];
        // <p_B|H|s_A> uses the A-B table. <s_B|H|p_A> uses the
        // reverse B-A table and the odd-parity sign from SK rotation.
        block[row * 4] = cosine * value_ab[sp];
        block[col] = -cosine * value_ba[sp];
        for (std::size_t cart = 0; cart < 3; ++cart)
        {
            const double dcosine = ((axis == cart ? 1.0 : 0.0)
                                    - cosine * direction[cart]) / distance;
            derivative[cart][row * 4] = direction[cart] * cosine * deriv_ab[sp]
                                          + dcosine * value_ab[sp];
            derivative[cart][col] = -(direction[cart] * cosine * deriv_ba[sp]
                                      + dcosine * value_ba[sp]);
        }
    }

    for (std::size_t p_row = 0; p_row < 3; ++p_row)
    {
        const std::size_t axis_row = cartesian_axis[p_row];
        const double l_row = direction[axis_row];
        for (std::size_t p_col = 0; p_col < 3; ++p_col)
        {
            const std::size_t axis_col = cartesian_axis[p_col];
            const double l_col = direction[axis_col];
            const double delta = p_row == p_col ? 1.0 : 0.0;
            const double sigma_minus_pi = value_ab[pp_sigma] - value_ab[pp_pi];
            const std::size_t element = (p_row + 1) * 4 + p_col + 1;
            block[element] = delta * value_ab[pp_pi] + l_row * l_col * sigma_minus_pi;
            for (std::size_t cart = 0; cart < 3; ++cart)
            {
                const double dl_row = ((axis_row == cart ? 1.0 : 0.0)
                                       - l_row * direction[cart]) / distance;
                const double dl_col = ((axis_col == cart ? 1.0 : 0.0)
                                       - l_col * direction[cart]) / distance;
                derivative[cart][element] = direction[cart]
                                                * (delta * deriv_ab[pp_pi]
                                                   + l_row * l_col * (deriv_ab[pp_sigma] - deriv_ab[pp_pi]))
                                            + (dl_row * l_col + l_row * dl_col) * sigma_minus_pi;
            }
        }
    }
}
} // namespace

SkfSpBlockEvaluation evaluate_sp_pair_block(const SkfData& ab,
                                             const SkfData& ba,
                                             const std::array<double, 3>& displacement_bohr)
{
    if (ab.homonuclear != ba.homonuclear
        || ab.declared_grid_points != ba.declared_grid_points
        || ab.hamiltonian.size() != ba.hamiltonian.size()
        || std::abs(ab.grid_spacing_bohr - ba.grid_spacing_bohr) > 1.0e-8)
    {
        throw std::runtime_error("Incompatible directed SKF integral grids");
    }
    if (ab.homonuclear && ab.filename != ba.filename)
    {
        throw std::runtime_error("A homonuclear SK block must use the same SKF file in both directions");
    }
    double distance_squared = 0.0;
    for (const double component : displacement_bohr)
    {
        if (!std::isfinite(component))
        {
            throw std::invalid_argument("SK pair displacement must be finite");
        }
        distance_squared += component * component;
    }
    const double distance = std::sqrt(distance_squared);
    if (!(distance > 0.0))
    {
        throw std::invalid_argument("SK pair displacement must have positive length");
    }

    const auto ab_values = ab.evaluate_integrals(distance);
    const auto ba_values = ba.evaluate_integrals(distance);
    SkfSpBlockEvaluation result;
    fill_block(ab_values, ba_values, displacement_bohr, distance, false,
               result.hamiltonian, result.d_hamiltonian_dR);
    fill_block(ab_values, ba_values, displacement_bohr, distance, true,
               result.overlap, result.d_overlap_dR);
    return result;
}

} // namespace ModuleDFTB
