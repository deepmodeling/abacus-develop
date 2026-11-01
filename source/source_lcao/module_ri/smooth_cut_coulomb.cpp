#include "smooth_cut_coulomb.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace Conv_Coulomb_Pot_K
{
double log_erfc_cut_factor(const double wave_number,
                           const double cutoff_radius,
                           const double width_factor)
{
    if (!std::isfinite(wave_number) || wave_number < 0.0
        || !std::isfinite(cutoff_radius) || cutoff_radius <= 0.0
        || !std::isfinite(width_factor) || width_factor <= 1.0 || width_factor > 1.25)
    {
        throw std::invalid_argument("log-erfc cutoff requires k>=0, Rc>0 and 1<width<=1.25, all finite");
    }
    const double cut_log_width = std::log(width_factor);
    const double cut_argument = wave_number * cutoff_radius;
    if (cut_argument == 0.0)
    {
        return 0.0;
    }
    // f(r) is a positive Gaussian mixture of sharp spheres with radius
    // Rc*exp(log(width)*t). Thus k^2*v(k)/(4*pi) is the Gaussian average
    // of 1-cos(k*Rc*exp(log(width)*t)). This also preserves positivity.
    // Integrate t in [-8,8]; omitted Gaussian mass is < 1.2e-29.
    // Bound the phase advance by pi/2 per panel, even at the upper endpoint.
    // A fixed small Gauss-Hermite rule would alias at large k*Rc.
    const double cut_pi = std::acos(-1.0);
    const double cut_panel_estimate = std::ceil(32.0 * cut_argument * cut_log_width
                                               * std::exp(8.0 * cut_log_width) / cut_pi);
    if (!std::isfinite(cut_panel_estimate) || cut_panel_estimate > 1000000.0)
    {
        throw std::runtime_error("log-erfc cutoff quadrature exceeds one million panels");
    }
    const int cut_panel_count = std::max(128, static_cast<int>(cut_panel_estimate));
    const double cut_step = 16.0 / cut_panel_count;
    double cut_integral = 0.0;
    for (int cut_panel = 0; cut_panel < cut_panel_count; ++cut_panel)
    {
        const double cut_coordinate = -8.0 + (cut_panel + 0.5) * cut_step;
        const double cut_phase_half = 0.5 * cut_argument * std::exp(cut_log_width * cut_coordinate);
        const double cut_sine = std::sin(cut_phase_half);
        cut_integral += std::exp(-cut_coordinate * cut_coordinate) * 2.0 * cut_sine * cut_sine;
    }
    return cut_integral * cut_step / std::sqrt(cut_pi);
}

double bvk_inscribed_radius(const std::array<std::array<double, 3>, 3>& bvk_lattice)
{
    std::array<std::array<double, 3>, 3> bvk_cofactor{};
    double bvk_shortest_squared = std::numeric_limits<double>::max();
    for (int bvk_axis = 0; bvk_axis < 3; ++bvk_axis)
    {
        double bvk_length_squared = 0.0;
        for (int cart_axis = 0; cart_axis < 3; ++cart_axis)
        {
            if (!std::isfinite(bvk_lattice[bvk_axis][cart_axis]))
            {
                throw std::invalid_argument("non-finite BvK lattice");
            }
            bvk_length_squared += bvk_lattice[bvk_axis][cart_axis] * bvk_lattice[bvk_axis][cart_axis];
            const int next_row = (bvk_axis + 1) % 3;
            const int last_row = (bvk_axis + 2) % 3;
            const int next_col = (cart_axis + 1) % 3;
            const int last_col = (cart_axis + 2) % 3;
            bvk_cofactor[bvk_axis][cart_axis]
                = bvk_lattice[next_row][next_col] * bvk_lattice[last_row][last_col]
                - bvk_lattice[next_row][last_col] * bvk_lattice[last_row][next_col];
        }
        bvk_shortest_squared = std::min(bvk_shortest_squared, bvk_length_squared);
    }
    double bvk_determinant = 0.0;
    for (int cart_axis = 0; cart_axis < 3; ++cart_axis)
    {
        bvk_determinant += bvk_lattice[0][cart_axis] * bvk_cofactor[0][cart_axis];
    }
    if (bvk_determinant == 0.0 || !std::isfinite(bvk_determinant) || bvk_shortest_squared <= 0.0)
    {
        throw std::invalid_argument("singular BvK lattice");
    }
    // For d=n*B, |n_i| <= |d|*||(B^-1)[:,i]||. These bounds include
    // every vector shorter than the shortest supplied basis vector.
    std::array<int, 3> bvk_bounds{};
    double bvk_candidate_count = 1.0;
    for (int bvk_axis = 0; bvk_axis < 3; ++bvk_axis)
    {
        double bvk_inverse_norm_squared = 0.0;
        for (int cart_axis = 0; cart_axis < 3; ++cart_axis)
        {
            const double bvk_inverse_entry = bvk_cofactor[bvk_axis][cart_axis] / bvk_determinant;
            bvk_inverse_norm_squared += bvk_inverse_entry * bvk_inverse_entry;
        }
        const double bvk_bound = std::floor(std::sqrt(bvk_shortest_squared * bvk_inverse_norm_squared) + 1.e-10);
        if (!std::isfinite(bvk_bound) || bvk_bound > 1000000.0)
        {
            throw std::runtime_error("BvK shortest-vector search needs a reduced lattice");
        }
        bvk_bounds[bvk_axis] = static_cast<int>(bvk_bound);
        bvk_candidate_count *= 2.0 * bvk_bounds[bvk_axis] + 1.0;
    }
    if (bvk_candidate_count > 10000000.0)
    {
        throw std::runtime_error("BvK shortest-vector search exceeds ten million candidates");
    }
    for (int first_cell = -bvk_bounds[0]; first_cell <= bvk_bounds[0]; ++first_cell)
    {
        for (int second_cell = -bvk_bounds[1]; second_cell <= bvk_bounds[1]; ++second_cell)
        {
            for (int third_cell = -bvk_bounds[2]; third_cell <= bvk_bounds[2]; ++third_cell)
            {
                if (first_cell == 0 && second_cell == 0 && third_cell == 0) { continue; }
                double bvk_distance_squared = 0.0;
                for (int cart_axis = 0; cart_axis < 3; ++cart_axis)
                {
                    const double bvk_displacement = first_cell * bvk_lattice[0][cart_axis]
                        + second_cell * bvk_lattice[1][cart_axis] + third_cell * bvk_lattice[2][cart_axis];
                    bvk_distance_squared += bvk_displacement * bvk_displacement;
                }
                bvk_shortest_squared = std::min(bvk_shortest_squared, bvk_distance_squared);
            }
        }
    }
    return 0.5 * std::sqrt(bvk_shortest_squared);
}
}
