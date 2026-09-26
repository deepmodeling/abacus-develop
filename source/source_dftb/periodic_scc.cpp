#include "source_dftb/periodic_scc.h"
#include "source_dftb/charge_mixing.h"

#include "source_base/parallel_common.h"
#include "source_base/parallel_reduce.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>

namespace ModuleDFTB
{
namespace
{
using Vec3 = std::array<double, 3>;
using Mat3 = std::array<Vec3, 3>;
constexpr double pi = 3.141592653589793238462643383279502884;
constexpr double gamma_tolerance = 1.0e-10;
constexpr double ewald_tolerance = 1.0e-9;
constexpr double same_distance_tolerance = 1.0e-5;

double dot(const Vec3& a, const Vec3& b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Vec3 add(const Vec3& a, const Vec3& b)
{
    return {{a[0] + b[0], a[1] + b[1], a[2] + b[2]}};
}

Vec3 subtract(const Vec3& a, const Vec3& b)
{
    return {{a[0] - b[0], a[1] - b[1], a[2] - b[2]}};
}

Vec3 scale(const Vec3& a, double factor)
{
    return {{factor * a[0], factor * a[1], factor * a[2]}};
}

double norm(const Vec3& a)
{
    return std::sqrt(dot(a, a));
}

Vec3 cross(const Vec3& a, const Vec3& b)
{
    return {{a[1] * b[2] - a[2] * b[1],
             a[2] * b[0] - a[0] * b[2],
             a[0] * b[1] - a[1] * b[0]}};
}

double determinant(const Mat3& a)
{
    return dot(a[0], cross(a[1], a[2]));
}

Mat3 reciprocal_lattice(const Mat3& lattice, double volume)
{
    Mat3 reciprocal{};
    reciprocal[0] = scale(cross(lattice[1], lattice[2]), 2.0 * pi / volume);
    reciprocal[1] = scale(cross(lattice[2], lattice[0]), 2.0 * pi / volume);
    reciprocal[2] = scale(cross(lattice[0], lattice[1]), 2.0 * pi / volume);
    return reciprocal;
}

Vec3 lattice_translation(const Mat3& lattice, int i, int j, int k)
{
    Vec3 result{};
    result = add(result, scale(lattice[0], static_cast<double>(i)));
    result = add(result, scale(lattice[1], static_cast<double>(j)));
    result = add(result, scale(lattice[2], static_cast<double>(k)));
    return result;
}

double short_gamma(double distance, double ua, double ub)
{
    const double taua = 3.2 * ua;
    const double taub = 3.2 * ub;
    if (!(ua > 0.0) || !(ub > 0.0) || distance < 0.0)
    {
        throw std::invalid_argument("DFTB Hubbard U and pair distance must be non-negative");
    }
    if (distance < same_distance_tolerance)
    {
        if (std::abs(ua - ub) < 3.125e-6)
        {
            return -0.5 * (ua + ub);
        }
        const double sum = taua + taub;
        return -0.5 * ((taua * taub) / sum + std::pow(taua * taub, 2) / std::pow(sum, 3));
    }
    if (std::abs(ua - ub) < 3.125e-6)
    {
        const double tau = 0.5 * (taua + taub);
        const double r = distance;
        return std::exp(-tau * r)
               * (1.0 / r + 0.6875 * tau + 0.1875 * r * tau * tau
                  + 0.0208333333333333333333333333333333 * r * r * tau * tau * tau);
    }
    const auto gamma_sub = [distance](double tau1, double tau2) {
        const double denominator = tau1 * tau1 - tau2 * tau2;
        return std::exp(-tau1 * distance)
               * (0.5 * std::pow(tau2, 4) * tau1 / (denominator * denominator)
                  - (std::pow(tau2, 6) - 3.0 * std::pow(tau2, 4) * tau1 * tau1)
                        / (distance * std::pow(denominator, 3)));
    };
    if (std::abs(taua - taub) < 1.0e-5)
    {
        // Match DFTB+'s near-degenerate limiting branch.
        const double tau = 0.5 * (taua + taub);
        const double r = distance;
        return std::exp(-tau * r)
               * (1.0 / r + 0.6875 * tau + 0.1875 * r * tau * tau
                  + 0.0208333333333333333333333333333333 * r * r * tau * tau * tau);
    }
    return gamma_sub(taua, taub) + gamma_sub(taub, taua);
}

double damped_gamma_correction(double distance, double ua, double ub)
{
    // DFTB+ HCorrection=Damping only damps pairs involving hydrogen-like species.
    // The C2N target contains no hydrogen, so this is exactly the undamped 3ob gamma.
    return -short_gamma(distance, ua, ub);
}

double gamma_cutoff(double ua, double ub)
{
    double high = 1.0;
    while (std::abs(short_gamma(high, ua, ub)) > gamma_tolerance && high < 256.0)
    {
        high *= 2.0;
    }
    if (high >= 256.0)
    {
        throw std::runtime_error("Could not determine the DFTB short-gamma cutoff");
    }
    double low = 0.0;
    for (int i = 0; i < 100; ++i)
    {
        const double middle = 0.5 * (low + high);
        if (std::abs(short_gamma(middle, ua, ub)) > gamma_tolerance) low = middle;
        else high = middle;
    }
    return high;
}

double gamma_derivative_wrt_ua(double distance, double ua, double ub)
{
    // Port of DFTB+'s gamma2pU() for the undamped C/N 3ob case. Using the
    // same equal-U limit is important: differentiating a numerically selected
    // near-degenerate branch gives a different derivative at Ua == Ub.
    if (!(ua > 0.0) || !(ub > 0.0) || distance < 0.0)
    {
        throw std::invalid_argument("DFTB Hubbard U and pair distance must be non-negative");
    }
    const double taua = 3.2 * ua;
    const double taub = 3.2 * ub;
    constexpr double min_hubbard_difference = 3.125e-6;
    if (distance < same_distance_tolerance)
    {
        if (std::abs(ua - ub) < min_hubbard_difference)
        {
            return 0.5;
        }
        const double inverse_sum = 1.0 / (taua + taub);
        return 1.6 * inverse_sum
               * (taub + inverse_sum
                                   * (-taua * taub
                                      + inverse_sum
                                            * (2.0 * taua * taub * taub
                                               + inverse_sum * (-3.0 * taua * taua * taub * taub))));
    }
    if (std::abs(ua - ub) < min_hubbard_difference)
    {
        const double tau = 0.5 * (taua + taub);
        const double tau_r = tau * distance;
        const double g = (48.0 + 33.0 * tau_r + 9.0 * tau_r * tau_r
                          + tau_r * tau_r * tau_r)
                         / (48.0 * distance);
        const double gp = 33.0 / 48.0 + 18.0 * tau_r / 48.0 + 3.0 * tau_r * tau_r / 48.0;
        return -3.2 * std::exp(-tau_r) * (gp - distance * g);
    }

    const double denominator = taua * taua - taub * taub;
    const double denominator2 = denominator * denominator;
    const double denominator3 = denominator2 * denominator;
    const double denominator4 = denominator3 * denominator;
    const double f = 0.5 * taua * std::pow(taub, 4) / denominator2
                     - (std::pow(taub, 6) - 3.0 * taua * taua * std::pow(taub, 4))
                           / (denominator3 * distance);
    const double fpt1 = -0.5 * (std::pow(taub, 6) + 3.0 * taua * taua * std::pow(taub, 4))
                            / denominator3
                        - 12.0 * std::pow(taua, 3) * std::pow(taub, 4)
                              / (denominator4 * distance);
    const double reversed_denominator = -denominator;
    const double fpt2 = 2.0 * std::pow(taub, 3) * std::pow(taua, 3)
                            / std::pow(reversed_denominator, 3)
                        + 12.0 * std::pow(taub, 4) * std::pow(taua, 3)
                              / (std::pow(reversed_denominator, 4) * distance);
    const double shortp_t = std::exp(-taua * distance) * (fpt1 - distance * f)
                            + std::exp(-taub * distance) * fpt2;
    return -3.2 * shortp_t;
}

const DftbPairParameters& pair_parameters(std::size_t species_a,
                                           std::size_t species_b,
                                           const std::vector<DftbPairParameters>& parameters)
{
    for (const auto& pair : parameters)
    {
        if (pair.species_a == species_a && pair.species_b == species_b
            && pair.ab != nullptr && pair.ba != nullptr)
        {
            return pair;
        }
    }
    throw std::runtime_error("No directed SKF pair is available for DFTB species indices "
                             + std::to_string(species_a) + "-" + std::to_string(species_b));
}

double pair_interaction_cutoff(const DftbPairParameters& pair)
{
    const double electronic = static_cast<double>(pair.ab->declared_grid_points)
                                  * pair.ab->grid_spacing_bohr
                              + 1.0;
    const double electronic_ba = static_cast<double>(pair.ba->declared_grid_points)
                                     * pair.ba->grid_spacing_bohr
                                 + 1.0;
    double cutoff = std::max(electronic, electronic_ba);
    if (pair.ab->has_repulsive_spline) cutoff = std::max(cutoff, pair.ab->repulsive.cutoff_bohr);
    return cutoff;
}

std::vector<DftbPairImage> make_pair_images(const DftbPeriodicInput& input,
                                            const Mat3& lattice,
                                            const Mat3& reciprocal)
{
    double cutoff = 0.0;
    for (std::size_t a = 0; a < input.atoms.size(); ++a)
    {
        for (std::size_t b = a; b < input.atoms.size(); ++b)
        {
            const auto& pair = pair_parameters(input.atoms[a].species, input.atoms[b].species,
                                               input.pair_parameters);
            cutoff = std::max(cutoff, pair_interaction_cutoff(pair));
            const double ua = input.atoms[a].homonuclear_data->hubbard_u_hartree[0];
            const double ub = input.atoms[b].homonuclear_data->hubbard_u_hartree[0];
            cutoff = std::max(cutoff, gamma_cutoff(ua, ub));
        }
    }

    std::array<int, 3> bounds{};
    for (int axis = 0; axis < 3; ++axis)
    {
        bounds[axis] = static_cast<int>(std::ceil(cutoff * norm(reciprocal[axis]) / (2.0 * pi))) + 2;
    }

    std::vector<DftbPairImage> result;
    for (std::size_t a = 0; a < input.atoms.size(); ++a)
    {
        for (std::size_t b = a; b < input.atoms.size(); ++b)
        {
            for (int i = -bounds[0]; i <= bounds[0]; ++i)
            {
                for (int j = -bounds[1]; j <= bounds[1]; ++j)
                {
                    for (int k = -bounds[2]; k <= bounds[2]; ++k)
                    {
                        if (a == b)
                        {
                            if (i == 0 && j == 0 && k == 0) continue;
                        }
                        const Vec3 translation = lattice_translation(lattice, i, j, k);
                        if (a == b)
                        {
                            bool canonical_translation = false;
                            for (const double component : translation)
                            {
                                if (component > 1.0e-12)
                                {
                                    canonical_translation = true;
                                    break;
                                }
                                if (component < -1.0e-12) break;
                            }
                            if (!canonical_translation) continue;
                        }
                        const Vec3 displacement = subtract(add(input.atoms[b].position_bohr, translation),
                                                            input.atoms[a].position_bohr);
                        const double distance = norm(displacement);
                        if (distance <= cutoff && distance > same_distance_tolerance)
                        {
                            DftbPairImage pair;
                            pair.atom_a = a;
                            pair.atom_b = b;
                            pair.translation_bohr = translation;
                            result.push_back(pair);
                        }
                    }
                }
            }
        }
    }
    return result;
}

double ewald_real_cutoff(double alpha)
{
    double high = 1.0;
    while (std::erfc(alpha * high) / high > ewald_tolerance) high *= 2.0;
    double low = 0.0;
    for (int i = 0; i < 100; ++i)
    {
        const double middle = 0.5 * (low + high);
        if (std::erfc(alpha * middle) / middle > ewald_tolerance) low = middle;
        else high = middle;
    }
    return high;
}

double ewald_reciprocal_cutoff(double alpha, double volume)
{
    const auto term = [alpha, volume](double g) {
        return 4.0 * pi * std::exp(-0.25 * g * g / (alpha * alpha)) / (volume * g * g);
    };
    double high = 1.0e-8;
    while (term(high) > ewald_tolerance) high *= 2.0;
    double low = 0.0;
    for (int i = 0; i < 100; ++i)
    {
        const double middle = 0.5 * (low + high);
        if (middle > 0.0 && term(middle) > ewald_tolerance) low = middle;
        else high = middle;
    }
    return high;
}

double ewald_alpha(const Mat3& lattice, const Mat3& reciprocal, double volume)
{
    const double min_r = std::min({norm(lattice[0]), norm(lattice[1]), norm(lattice[2])});
    const double min_g = std::min({norm(reciprocal[0]), norm(reciprocal[1]), norm(reciprocal[2])});
    const auto difference = [min_r, min_g, volume](double alpha) {
        const auto g_term = [alpha, volume](double g) {
            return 4.0 * pi * std::exp(-0.25 * g * g / (alpha * alpha)) / (volume * g * g);
        };
        const auto r_term = [alpha](double r) { return std::erfc(alpha * r) / r; };
        return (g_term(4.0 * min_g) - g_term(5.0 * min_g))
               - (r_term(2.0 * min_r) - r_term(3.0 * min_r));
    };
    double low = 1.0e-8;
    double high = low;
    while (difference(high) <= 0.0 && high < 64.0) high *= 2.0;
    if (high >= 64.0) throw std::runtime_error("Could not determine Ewald alpha");
    low = 0.5 * high;
    for (int i = 0; i < 120; ++i)
    {
        const double middle = 0.5 * (low + high);
        if (difference(middle) < 0.0) low = middle;
        else high = middle;
    }
    return 0.5 * (low + high);
}

std::vector<Vec3> make_real_translations(const Mat3& lattice,
                                          const Mat3& reciprocal,
                                          double cutoff)
{
    std::array<int, 3> bounds{};
    for (int axis = 0; axis < 3; ++axis)
    {
        bounds[axis] = static_cast<int>(std::ceil(cutoff * norm(reciprocal[axis]) / (2.0 * pi))) + 2;
    }
    std::vector<Vec3> result;
    for (int i = -bounds[0]; i <= bounds[0]; ++i)
    {
        for (int j = -bounds[1]; j <= bounds[1]; ++j)
        {
            for (int k = -bounds[2]; k <= bounds[2]; ++k)
            {
                const Vec3 translation = lattice_translation(lattice, i, j, k);
                result.push_back(translation);
            }
        }
    }
    return result;
}

std::vector<Vec3> make_reciprocal_vectors(const Mat3& lattice,
                                           const Mat3& reciprocal,
                                           double cutoff)
{
    std::array<int, 3> bounds{};
    for (int axis = 0; axis < 3; ++axis)
    {
        bounds[axis] = static_cast<int>(std::ceil(cutoff * norm(lattice[axis]) / (2.0 * pi))) + 2;
    }
    std::vector<Vec3> result;
    for (int i = -bounds[0]; i <= bounds[0]; ++i)
    {
        for (int j = -bounds[1]; j <= bounds[1]; ++j)
        {
            for (int k = -bounds[2]; k <= bounds[2]; ++k)
            {
                if (i == 0 && j == 0 && k == 0) continue;
                const int first_nonzero = i != 0 ? i : (j != 0 ? j : k);
                if (first_nonzero < 0) continue;
                Vec3 g{};
                g = add(g, scale(reciprocal[0], static_cast<double>(i)));
                g = add(g, scale(reciprocal[1], static_cast<double>(j)));
                g = add(g, scale(reciprocal[2], static_cast<double>(k)));
                if (norm(g) <= cutoff + 1.0e-10) result.push_back(g);
            }
        }
    }
    return result;
}

double ewald_kernel(const Vec3& displacement,
                    const std::vector<Vec3>& translations,
                    const std::vector<Vec3>& reciprocal_vectors,
                    double alpha,
                    double volume,
                    double real_cutoff)
{
    double result = 0.0;
    for (const Vec3& translation : translations)
    {
        const Vec3 dr = add(displacement, translation);
        const double distance = norm(dr);
        if (distance > same_distance_tolerance && distance <= real_cutoff + 1.0e-10)
        {
            result += std::erfc(alpha * distance) / distance;
        }
    }
    for (const Vec3& g : reciprocal_vectors)
    {
        const double g2 = dot(g, g);
        result += 8.0 * pi / volume * std::exp(-g2 / (4.0 * alpha * alpha))
                  * std::cos(dot(g, displacement)) / g2;
    }
    result -= pi / (volume * alpha * alpha);
    if (norm(displacement) <= same_distance_tolerance)
    {
        result -= 2.0 * alpha / std::sqrt(pi);
    }
    return result;
}

std::vector<double> build_gamma_matrix(const DftbPeriodicInput& input,
                                        const std::vector<DftbPairImage>& pairs,
                                        const Mat3& lattice,
                                        const Mat3& reciprocal,
                                        double alpha,
                                        double volume)
{
    const double real_cutoff = ewald_real_cutoff(alpha);
    const double reciprocal_cutoff = ewald_reciprocal_cutoff(alpha, volume);
    const auto translations = make_real_translations(lattice, reciprocal, real_cutoff);
    const auto g_vectors = make_reciprocal_vectors(lattice, reciprocal, reciprocal_cutoff);
    const std::size_t n = input.atoms.size();
    std::vector<double> gamma(n * n, 0.0);
    const auto at = [n, &gamma](std::size_t i, std::size_t j) -> double& { return gamma[i * n + j]; };
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j <= i; ++j)
        {
            const Vec3 displacement = subtract(input.atoms[i].position_bohr, input.atoms[j].position_bohr);
            const double value = ewald_kernel(displacement, translations, g_vectors, alpha, volume,
                                              real_cutoff);
            at(i, j) = value;
            at(j, i) = value;
        }
    }

    // The short-range part is -expGamma. Add the origin term (which is +U)
    // and periodic image contributions. Self-image records are half listed.
    for (std::size_t atom = 0; atom < n; ++atom)
    {
        const double u = input.atoms[atom].homonuclear_data->hubbard_u_hartree[0];
        at(atom, atom) += damped_gamma_correction(0.0, u, u);
    }
    for (const auto& pair : pairs)
    {
        const Vec3 displacement = subtract(add(input.atoms[pair.atom_b].position_bohr,
                                                pair.translation_bohr),
                                           input.atoms[pair.atom_a].position_bohr);
        const double ua = input.atoms[pair.atom_a].homonuclear_data->hubbard_u_hartree[0];
        const double ub = input.atoms[pair.atom_b].homonuclear_data->hubbard_u_hartree[0];
        const double correction = damped_gamma_correction(norm(displacement), ua, ub);
        if (pair.atom_a == pair.atom_b)
        {
            at(pair.atom_a, pair.atom_a) += 2.0 * correction;
        }
        else
        {
            at(pair.atom_a, pair.atom_b) += correction;
            at(pair.atom_b, pair.atom_a) += correction;
        }
    }
    return gamma;
}

std::vector<double> third_order_potential(const DftbPeriodicInput& input,
                                           const std::vector<DftbPairImage>& pairs,
                                           const std::vector<double>& charges)
{
    std::vector<double> potential(input.atoms.size(), 0.0);
    if (!input.third_order) return potential;
    for (const auto& atom : input.atoms)
    {
        if (atom.species >= input.hubbard_derivative.size())
        {
            throw std::invalid_argument("Missing DFTB3 Hubbard derivative for a species");
        }
    }

    const auto add_self_event = [&](std::size_t atom, double distance, double multiplicity) {
        const double ua = input.atoms[atom].homonuclear_data->hubbard_u_hartree[0];
        const double ub = ua;
        const double dgamma_dqa = gamma_derivative_wrt_ua(distance, ua, ub)
                                  * input.hubbard_derivative[input.atoms[atom].species];
        const double q = charges[atom];
        // DFTB+'s three full-order terms are each divided by three before
        // being combined; for one self event their sum equals dGamma/dq*q^2.
        potential[atom] += multiplicity * dgamma_dqa * q * q;
    };

    for (std::size_t atom = 0; atom < input.atoms.size(); ++atom)
    {
        add_self_event(atom, 0.0, 1.0);
    }
    for (const auto& pair : pairs)
    {
        const Vec3 displacement = subtract(add(input.atoms[pair.atom_b].position_bohr,
                                                pair.translation_bohr),
                                           input.atoms[pair.atom_a].position_bohr);
        const double distance = norm(displacement);
        const double ua = input.atoms[pair.atom_a].homonuclear_data->hubbard_u_hartree[0];
        const double ub = input.atoms[pair.atom_b].homonuclear_data->hubbard_u_hartree[0];
        const double gp_a = gamma_derivative_wrt_ua(distance, ua, ub)
                            * input.hubbard_derivative[input.atoms[pair.atom_a].species];
        const double gp_b = gamma_derivative_wrt_ua(distance, ub, ua)
                            * input.hubbard_derivative[input.atoms[pair.atom_b].species];
        const double qa = charges[pair.atom_a];
        const double qb = charges[pair.atom_b];
        if (pair.atom_a == pair.atom_b)
        {
            // The neighbor list contains both +R and -R self images, whereas
            // this canonical pair list contains one representative.
            potential[pair.atom_a] += 2.0 * gp_a * qa * qa;
        }
        else
        {
            potential[pair.atom_a] += (2.0 * gp_a * qa * qb + gp_b * qb * qb) / 3.0;
            potential[pair.atom_b] += (2.0 * gp_b * qa * qb + gp_a * qa * qa) / 3.0;
        }
    }
    return potential;
}

struct KPointResult
{
    std::vector<DftbKPointSpectrum> spectra;
    DftbFermiFilling filling;
    std::vector<double> populations;
};

KPointResult solve_for_charges(const DftbPeriodicInput& input,
                               const std::vector<DftbPairImage>& pairs,
                               const Mat3& reciprocal,
                               const std::vector<double>& gamma,
                               const std::vector<double>& charges)
{
    const std::size_t n_atoms = input.atoms.size();
    const std::size_t n_orbitals = 4 * n_atoms;
    std::vector<double> potential(n_atoms, 0.0);
    std::vector<double> third = third_order_potential(input, pairs, charges);
    for (std::size_t i = 0; i < n_atoms; ++i)
    {
        for (std::size_t j = 0; j < n_atoms; ++j)
        {
            potential[i] += gamma[i * n_atoms + j] * charges[j];
        }
        potential[i] += third[i];
    }

    KPointResult result;
    result.spectra.resize(input.kpoints.size());
    const int mpi_rank = Parallel_Common::get_rank();
    const int mpi_size = Parallel_Common::get_size();
    for (std::size_t ik = 0; ik < input.kpoints.size(); ++ik)
    {
        const auto& kp = input.kpoints[ik];
        result.spectra[ik].weight = kp.weight;
        result.spectra[ik].solution.dimension = n_orbitals;
        result.spectra[ik].solution.eigenvalues_hartree.resize(n_orbitals);
        const int owner = static_cast<int>(ik % static_cast<std::size_t>(mpi_size));
        int solve_status = 0;
        char error_message[1024] = {};
        if (mpi_rank == owner)
        {
            try
            {
                Vec3 k{};
                for (int axis = 0; axis < 3; ++axis)
                {
                    k = add(k, scale(reciprocal[axis], kp.fractional[axis]));
                }
                auto matrices = assemble_sp_bloch_matrices(input.atoms, pairs, input.pair_parameters, k);
                for (std::size_t row = 0; row < n_orbitals; ++row)
                {
                    const std::size_t atom_row = row / 4;
                    for (std::size_t col = 0; col < n_orbitals; ++col)
                    {
                        const std::size_t atom_col = col / 4;
                        matrices.h(row, col) += matrices.s(row, col)
                                                * 0.5 * (potential[atom_row] + potential[atom_col]);
                    }
                }
                result.spectra[ik].solution = solve_generalized_hermitian(matrices);
            }
            catch (const std::exception& error)
            {
                solve_status = 1;
                std::strncpy(error_message, error.what(), sizeof(error_message) - 1);
            }
        }
        int all_ranks_succeeded = solve_status == 0 ? 1 : 0;
        Parallel_Reduce::reduce_min(all_ranks_succeeded);
        solve_status = all_ranks_succeeded == 0 ? 1 : 0;
        std::vector<int> error_codes(sizeof(error_message), 0);
        if (mpi_rank == owner)
        {
            for (std::size_t i = 0; i < sizeof(error_message); ++i)
                error_codes[i] = static_cast<unsigned char>(error_message[i]);
        }
        Parallel_Reduce::reduce_all(error_codes.data(), static_cast<int>(error_codes.size()));
        for (std::size_t i = 0; i < sizeof(error_message); ++i)
            error_message[i] = static_cast<char>(error_codes[i]);
        if (solve_status != 0)
        {
            throw std::runtime_error(std::string("Native DFTB k-point eigensolver failed: ") + error_message);
        }
        Parallel_Reduce::reduce_all(result.spectra[ik].solution.eigenvalues_hartree.data(),
                                    static_cast<int>(n_orbitals));
    }

    result.filling = fill_fermi_occupations(result.spectra, input.total_electrons,
                                             input.thermal_energy_hartree);

    // Reassemble only S(k) for the Mulliken analysis; this avoids carrying a
    // second full H/S matrix set through the SCC loop.
    result.populations.assign(n_atoms, 0.0);
    for (std::size_t ik = 0; ik < input.kpoints.size(); ++ik)
    {
        if (static_cast<int>(ik % static_cast<std::size_t>(mpi_size)) != mpi_rank) continue;
        const auto& kp = input.kpoints[ik];
        Vec3 k{};
        for (int axis = 0; axis < 3; ++axis)
        {
            k = add(k, scale(reciprocal[axis], kp.fractional[axis]));
        }
        const auto matrices = assemble_sp_bloch_matrices(input.atoms, pairs, input.pair_parameters, k);
        const auto& eig = result.spectra[ik].solution;
        const auto& occupations = result.filling.occupations[ik];
        for (std::size_t mu = 0; mu < n_orbitals; ++mu)
        {
            const std::size_t atom = mu / 4;
            double mulliken = 0.0;
            for (std::size_t nu = 0; nu < n_orbitals; ++nu)
            {
                std::complex<double> density{};
                for (std::size_t band = 0; band < n_orbitals; ++band)
                {
                    const double occupation = occupations[band];
                    if (occupation == 0.0) continue;
                    density += occupation * eig.coefficient(mu, band)
                               * std::conj(eig.coefficient(nu, band));
                }
                mulliken += std::real(density * matrices.s(nu, mu));
            }
            result.populations[atom] += kp.weight * mulliken;
        }
    }
    Parallel_Reduce::reduce_all(result.populations.data(), static_cast<int>(n_atoms));
    return result;
}

std::vector<DftbBandPoint> solve_band_path(const DftbPeriodicInput& input,
                                              const std::vector<DftbPairImage>& pairs,
                                              const Mat3& reciprocal,
                                              const std::vector<double>& gamma,
                                              const std::vector<double>& charges)
{
    const std::size_t n_atoms = input.atoms.size();
    const std::size_t n_orbitals = 4 * n_atoms;
    std::vector<double> potential(n_atoms, 0.0);
    const std::vector<double> third = third_order_potential(input, pairs, charges);
    for (std::size_t i = 0; i < n_atoms; ++i)
    {
        for (std::size_t j = 0; j < n_atoms; ++j) potential[i] += gamma[i * n_atoms + j] * charges[j];
        potential[i] += third[i];
    }

    const int mpi_rank = Parallel_Common::get_rank();
    const int mpi_size = Parallel_Common::get_size();
    std::vector<DftbBandPoint> result(input.band_kpoints.size());
    for (std::size_t ik = 0; ik < input.band_kpoints.size(); ++ik)
    {
        result[ik].fractional = input.band_kpoints[ik].fractional;
        result[ik].label = input.band_kpoints[ik].label;
        result[ik].eigenvalues_hartree.resize(n_orbitals);
        const int owner = static_cast<int>(ik % static_cast<std::size_t>(mpi_size));
        int solve_status = 0;
        char error_message[1024] = {};
        if (mpi_rank == owner)
        {
            try
            {
                Vec3 k{};
                for (int axis = 0; axis < 3; ++axis)
                    k = add(k, scale(reciprocal[axis], input.band_kpoints[ik].fractional[axis]));
                auto matrices = assemble_sp_bloch_matrices(input.atoms, pairs, input.pair_parameters, k);
                for (std::size_t row = 0; row < n_orbitals; ++row)
                {
                    const std::size_t atom_row = row / 4;
                    for (std::size_t col = 0; col < n_orbitals; ++col)
                    {
                        const std::size_t atom_col = col / 4;
                        matrices.h(row, col) += matrices.s(row, col)
                                                * 0.5 * (potential[atom_row] + potential[atom_col]);
                    }
                }
                result[ik].eigenvalues_hartree = solve_generalized_hermitian(matrices).eigenvalues_hartree;
            }
            catch (const std::exception& error)
            {
                solve_status = 1;
                std::strncpy(error_message, error.what(), sizeof(error_message) - 1);
            }
        }
        int all_ranks_succeeded = solve_status == 0 ? 1 : 0;
        Parallel_Reduce::reduce_min(all_ranks_succeeded);
        solve_status = all_ranks_succeeded == 0 ? 1 : 0;
        std::vector<int> error_codes(sizeof(error_message), 0);
        if (mpi_rank == owner)
        {
            for (std::size_t i = 0; i < sizeof(error_message); ++i)
                error_codes[i] = static_cast<unsigned char>(error_message[i]);
        }
        Parallel_Reduce::reduce_all(error_codes.data(), static_cast<int>(error_codes.size()));
        for (std::size_t i = 0; i < sizeof(error_message); ++i)
            error_message[i] = static_cast<char>(error_codes[i]);
        if (solve_status != 0)
            throw std::runtime_error(std::string("Native DFTB band eigensolver failed: ") + error_message);
        Parallel_Reduce::reduce_all(result[ik].eigenvalues_hartree.data(), static_cast<int>(n_orbitals));
        if (ik > 0)
        {
            Vec3 delta{};
            for (int axis = 0; axis < 3; ++axis)
                delta = add(delta, scale(reciprocal[axis], input.band_kpoints[ik].fractional[axis]
                                                           - input.band_kpoints[ik - 1].fractional[axis]));
            result[ik].distance_inverse_bohr = result[ik - 1].distance_inverse_bohr + std::sqrt(dot(delta, delta));
        }
    }
    return result;
}

double quadratic_form(const std::vector<double>& matrix,
                      const std::vector<double>& x,
                      std::vector<double>* product = nullptr)
{
    const std::size_t n = x.size();
    std::vector<double> local(n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < n; ++j) local[i] += matrix[i * n + j] * x[j];
    }
    double value = 0.0;
    for (std::size_t i = 0; i < n; ++i) value += x[i] * local[i];
    if (product != nullptr) *product = std::move(local);
    return value;
}


} // namespace

DftbPeriodicResult solve_periodic_dftb(const DftbPeriodicInput& input)
{
    return solve_periodic_dftb(input, std::function<void(const DftbSccIteration&)>());
}

DftbPeriodicResult solve_periodic_dftb(
    const DftbPeriodicInput& input,
    const std::function<void(const DftbSccIteration&)>& on_iteration)
{
    if (input.atoms.empty() || input.kpoints.empty() || input.pair_parameters.empty()
        || input.maximum_scc_iterations <= 0 || !(input.scc_tolerance > 0.0)
        || !(input.mixing_parameter > 0.0 && input.mixing_parameter <= 1.0)
        || (input.mixing_method != "linear" && input.mixing_method != "pulay"
            && input.mixing_method != "broyden")
        || input.mixing_history < 2 || input.mixing_history > 20
        || !(input.broyden_inverse_jacobi_weight > 0.0)
        || !(input.broyden_minimal_weight > 0.0)
        || !(input.broyden_maximal_weight >= input.broyden_minimal_weight)
        || !(input.broyden_weight_factor > 0.0)
        || !(input.total_electrons >= 0.0) || input.hubbard_derivative.empty())
    {
        throw std::invalid_argument("Incomplete or invalid native periodic DFTB input");
    }
    for (const auto& atom : input.atoms)
    {
        if (atom.homonuclear_data == nullptr || !atom.homonuclear_data->has_atomic_data
            || atom.species >= input.hubbard_derivative.size())
        {
            throw std::invalid_argument("Each DFTB atom needs homonuclear SKF data and Hubbard parameters");
        }
    }
    double weight_sum = 0.0;
    for (const auto& kpoint : input.kpoints)
    {
        if (!(kpoint.weight >= 0.0) || !std::isfinite(kpoint.weight))
            throw std::invalid_argument("DFTB k-point weights must be finite and non-negative");
        weight_sum += kpoint.weight;
    }
    if (!std::isfinite(weight_sum) || std::abs(weight_sum - 1.0) > 1.0e-10)
    {
        throw std::invalid_argument("DFTB k-point weights must sum to one");
    }

    Mat3 lattice = input.lattice_bohr;
    const double volume = std::abs(determinant(lattice));
    if (!(volume > 1.0e-10)) throw std::invalid_argument("DFTB lattice must have non-zero volume");
    const Mat3 reciprocal = reciprocal_lattice(lattice, determinant(lattice));
    const auto pairs = make_pair_images(input, lattice, reciprocal);
    const double alpha = ewald_alpha(lattice, reciprocal, volume);
    const auto gamma = build_gamma_matrix(input, pairs, lattice, reciprocal, alpha, volume);

    double neutral_electrons = 0.0;
    for (const auto& atom : input.atoms)
    {
        neutral_electrons += atom.homonuclear_data->valence_electron_count();
    }
    const double target_charge = input.total_electrons - neutral_electrons;
    std::vector<double> charges(input.atoms.size(), target_charge / static_cast<double>(input.atoms.size()));
    DftbChargeMixerParameters mixer_parameters;
    mixer_parameters.method = input.mixing_method;
    mixer_parameters.mixing_parameter = input.mixing_parameter;
    mixer_parameters.history = input.mixing_history;
    mixer_parameters.inverse_jacobi_weight = input.broyden_inverse_jacobi_weight;
    mixer_parameters.minimal_weight = input.broyden_minimal_weight;
    mixer_parameters.maximal_weight = input.broyden_maximal_weight;
    mixer_parameters.weight_factor = input.broyden_weight_factor;
    DftbChargeMixer charge_mixer(mixer_parameters);
    const double repulsive_energy_hartree = calculate_repulsive_energy_hartree(input.atoms, pairs,
                                                                                 input.pair_parameters);

    DftbPeriodicResult result;
    KPointResult state;
    double previous_band_free_energy = 0.0;
    double previous_electronic_energy = 0.0;
    bool has_previous_energy = false;
    for (int iteration = 1; iteration <= input.maximum_scc_iterations; ++iteration)
    {
        state = solve_for_charges(input, pairs, reciprocal, gamma, charges);
        std::vector<double> output_charges(input.atoms.size(), 0.0);
        double charge_sum = 0.0;
        for (std::size_t atom = 0; atom < input.atoms.size(); ++atom)
        {
            output_charges[atom] = state.populations[atom]
                                   - input.atoms[atom].homonuclear_data->valence_electron_count();
            charge_sum += output_charges[atom];
        }
        const double charge_correction = (target_charge - charge_sum)
                                         / static_cast<double>(output_charges.size());
        for (double& charge : output_charges) charge += charge_correction;
        charge_sum = 0.0;
        for (const double charge : output_charges) charge_sum += charge;

        std::vector<double> residual(charges.size(), 0.0);
        result.maximum_charge_residual = 0.0;
        for (std::size_t atom = 0; atom < charges.size(); ++atom)
        {
            residual[atom] = output_charges[atom] - charges[atom];
            result.maximum_charge_residual = std::max(result.maximum_charge_residual, std::abs(residual[atom]));
        }
        result.scc_iterations = iteration;
        result.scc_residual_history.push_back(result.maximum_charge_residual);

        const bool converged = result.maximum_charge_residual <= input.scc_tolerance;
        std::vector<double> iteration_gamma_charge;
        const double iteration_scc_quadratic = quadratic_form(gamma, output_charges, &iteration_gamma_charge);
        const std::vector<double> iteration_third = third_order_potential(input, pairs, output_charges);
        double iteration_third_energy = 0.0;
        double iteration_potential_expectation = 0.0;
        for (std::size_t atom = 0; atom < output_charges.size(); ++atom)
        {
            iteration_third_energy += output_charges[atom] * iteration_third[atom] / 3.0;
            iteration_potential_expectation += state.populations[atom]
                                                * (iteration_gamma_charge[atom] + iteration_third[atom]);
        }
        const double iteration_h0_energy = state.filling.band_free_energy_hartree
                                           - iteration_potential_expectation;
        const double iteration_electronic_energy = iteration_h0_energy + 0.5 * iteration_scc_quadratic
                                                   + iteration_third_energy;
        std::vector<double> next_charges;
        DftbSccIteration iteration_output;
        iteration_output.iteration = iteration;
        iteration_output.maximum_charge_residual = result.maximum_charge_residual;
        iteration_output.net_electron_excess = charge_sum;
        iteration_output.band_energy_hartree = state.filling.band_energy_hartree;
        iteration_output.band_free_energy_hartree = state.filling.band_free_energy_hartree;
        iteration_output.fermi_energy_hartree = state.filling.fermi_energy_hartree;
        iteration_output.has_previous_energy = has_previous_energy;
        iteration_output.electronic_energy_hartree = iteration_electronic_energy;
        iteration_output.total_free_energy_hartree = iteration_electronic_energy + repulsive_energy_hartree;
        if (has_previous_energy)
        {
            iteration_output.band_free_energy_change_hartree = state.filling.band_free_energy_hartree
                                                                - previous_band_free_energy;
            iteration_output.electronic_energy_change_hartree = iteration_electronic_energy
                                                                - previous_electronic_energy;
        }

        if (converged)
        {
            next_charges = output_charges;
            iteration_output.mixing_step = "converged";
            result.converged = true;
        }
        else
        {
            next_charges = charge_mixer.mix(charges, residual, target_charge);
            iteration_output.mixing_step = charge_mixer.last_step();
        }
        if (on_iteration) on_iteration(iteration_output);
        previous_band_free_energy = state.filling.band_free_energy_hartree;
        previous_electronic_energy = iteration_electronic_energy;
        has_previous_energy = true;
        charges = std::move(next_charges);
        if (result.converged) break;
    }

    if (!result.converged)
    {
        throw std::runtime_error("Native periodic DFTB SCC did not converge in "
                                 + std::to_string(input.maximum_scc_iterations) + " iterations; max charge residual="
                                 + std::to_string(result.maximum_charge_residual));
    }

    // Re-evaluate at the converged density so reported eigenvalues, populations,
    // and potential energies are mutually consistent.
    state = solve_for_charges(input, pairs, reciprocal, gamma, charges);
    std::vector<double> final_charges(input.atoms.size(), 0.0);
    for (std::size_t atom = 0; atom < input.atoms.size(); ++atom)
    {
        final_charges[atom] = state.populations[atom]
                              - input.atoms[atom].homonuclear_data->valence_electron_count();
    }
    std::vector<double> gamma_charges;
    const double scc_quadratic = quadratic_form(gamma, final_charges, &gamma_charges);
    const auto third = third_order_potential(input, pairs, final_charges);
    double third_energy = 0.0;
    for (std::size_t atom = 0; atom < final_charges.size(); ++atom)
    {
        third_energy += final_charges[atom] * third[atom] / 3.0;
    }
    std::vector<double> total_potential = gamma_charges;
    double potential_expectation = 0.0;
    for (std::size_t atom = 0; atom < final_charges.size(); ++atom)
    {
        total_potential[atom] += third[atom];
        potential_expectation += state.populations[atom] * total_potential[atom];
    }

    result.electron_excess_charges = std::move(final_charges);
    result.atomic_electron_populations = state.populations;
    result.fermi_energy_hartree = state.filling.fermi_energy_hartree;
    result.band_energy_hartree = state.filling.band_energy_hartree;
    result.band_free_energy_hartree = state.filling.band_free_energy_hartree;
    result.kpoint_eigenvalues.resize(state.spectra.size());
    for (std::size_t ik = 0; ik < state.spectra.size(); ++ik)
    {
        result.kpoint_eigenvalues[ik].fractional = input.kpoints[ik].fractional;
        result.kpoint_eigenvalues[ik].weight = input.kpoints[ik].weight;
        result.kpoint_eigenvalues[ik].eigenvalues_hartree = state.spectra[ik].solution.eigenvalues_hartree;
        result.kpoint_eigenvalues[ik].occupations = state.filling.occupations[ik];
    }
    result.h0_energy_hartree = state.filling.band_free_energy_hartree - potential_expectation;
    result.scc_energy_hartree = 0.5 * scc_quadratic;
    result.third_order_energy_hartree = third_energy;
    result.repulsive_energy_hartree = repulsive_energy_hartree;
    result.total_free_energy_hartree = result.h0_energy_hartree + result.scc_energy_hartree
                                       + result.third_order_energy_hartree
                                       + result.repulsive_energy_hartree;
    if (!input.band_kpoints.empty())
        result.band_structure = solve_band_path(input, pairs, reciprocal, gamma, result.electron_excess_charges);
    return result;
}

} // namespace ModuleDFTB
