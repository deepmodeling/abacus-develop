#include "vdwd3_evaluator.h"

#include "vdwd3_data.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

namespace vdw
{
namespace d3
{
namespace
{

constexpr double kCoordinationSteepness = 16.0;
constexpr double kReferenceWeight = 4.0;
constexpr double kAtmRadiusScale = 4.0 / 3.0;

double dot(const Vec3& lhs, const Vec3& rhs)
{
    return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

Vec3 cross(const Vec3& lhs, const Vec3& rhs)
{
    return Vec3(lhs.y * rhs.z - rhs.y * lhs.z,
                lhs.z * rhs.x - rhs.z * lhs.x,
                lhs.x * rhs.y - rhs.x * lhs.y);
}

double norm(const Vec3& value)
{
    return std::sqrt(dot(value, value));
}

void add_outer(Matrix3& matrix, const Vec3& lhs, const Vec3& rhs, double scale = 1.0)
{
    const double left[3] = {lhs.x, lhs.y, lhs.z};
    const double right[3] = {rhs.x, rhs.y, rhs.z};
    for (int row = 0; row < 3; ++row)
    {
        for (int column = 0; column < 3; ++column)
        {
            matrix.value[row][column] += scale * left[row] * right[column];
        }
    }
}

bool validate(const Structure& structure, const Cutoffs& cutoffs, std::string& error)
{
    if (structure.atomic_numbers.size() != structure.positions.size())
    {
        error = "atomic-number and position arrays have different sizes";
        return false;
    }
    for (int atomic_number : structure.atomic_numbers)
    {
        if (atomic_number <= 0 || atomic_number > data::max_element)
        {
            error = "DFT-D3 supports atomic numbers 1 through 103";
            return false;
        }
    }
    if (!std::isfinite(cutoffs.disp2) || !std::isfinite(cutoffs.disp3)
        || !std::isfinite(cutoffs.cn) || !(cutoffs.disp2 > 0.0)
        || !(cutoffs.disp3 > 0.0) || !(cutoffs.cn > 0.0))
    {
        error = "DFT-D3 real-space cutoffs must be positive";
        return false;
    }
    if (!std::isfinite(cutoffs.width2) || cutoffs.width2 < 0.0
        || cutoffs.width2 > cutoffs.disp2)
    {
        error = "two-body smooth width must satisfy 0 <= width <= cutoff";
        return false;
    }
    if (!std::isfinite(cutoffs.width3) || cutoffs.width3 < 0.0
        || cutoffs.width3 > cutoffs.disp3)
    {
        error = "three-body smooth width must satisfy 0 <= width <= cutoff";
        return false;
    }
    return true;
}

std::vector<Vec3> lattice_points(const Structure& structure, double cutoff)
{
    if (!structure.periodic[0] && !structure.periodic[1] && !structure.periodic[2])
    {
        return std::vector<Vec3>(1, Vec3());
    }

    std::array<Vec3, 3> normal = {{cross(structure.lattice[1], structure.lattice[2]),
                                   cross(structure.lattice[2], structure.lattice[0]),
                                   cross(structure.lattice[0], structure.lattice[1])}};
    std::array<int, 3> repeat;
    for (int direction = 0; direction < 3; ++direction)
    {
        normal[direction] = normal[direction] / norm(normal[direction]);
        repeat[direction] = static_cast<int>(
            std::ceil(std::abs(cutoff / dot(normal[direction], structure.lattice[direction]))));
    }

    std::vector<Vec3> translations;
    translations.reserve(static_cast<std::size_t>((2 * repeat[0] + 1)
                                                   * (2 * repeat[1] + 1)
                                                   * (2 * repeat[2] + 1)));
    for (int ix = -repeat[0]; ix <= repeat[0]; ++ix)
    {
        for (int iy = -repeat[1]; iy <= repeat[1]; ++iy)
        {
            for (int iz = -repeat[2]; iz <= repeat[2]; ++iz)
            {
                translations.push_back(structure.lattice[0] * ix
                                       + structure.lattice[1] * iy
                                       + structure.lattice[2] * iz);
            }
        }
    }
    return translations;
}

void smooth_cutoff(double distance, double cutoff, double width, double& value, double& derivative)
{
    // Match s-dftd3 exactly: width == cutoff selects the legacy sharp path.
    if (width <= 0.0 || width >= cutoff)
    {
        value = 1.0;
        derivative = 0.0;
        return;
    }
    const double inner = cutoff - width;
    if (distance <= inner)
    {
        value = 1.0;
        derivative = 0.0;
    }
    else if (distance >= cutoff)
    {
        value = 0.0;
        derivative = 0.0;
    }
    else
    {
        const double x = (cutoff - distance) / width;
        value = x * x * x * (10.0 + x * (-15.0 + 6.0 * x));
        derivative = -30.0 * x * x * (1.0 - x) * (1.0 - x) / width;
    }
}

double coordination_count(int atomic_number_i, int atomic_number_j, double distance)
{
    const double radius = data::covalent_radius(atomic_number_i)
                          + data::covalent_radius(atomic_number_j);
    return 1.0 / (1.0 + std::exp(-kCoordinationSteepness * (radius / distance - 1.0)));
}

double coordination_derivative(int atomic_number_i, int atomic_number_j, double distance)
{
    const double radius = data::covalent_radius(atomic_number_i)
                          + data::covalent_radius(atomic_number_j);
    const double exponential = std::exp(-kCoordinationSteepness * (radius / distance - 1.0));
    const double denominator = exponential + 1.0;
    return -kCoordinationSteepness * radius * exponential
           / (distance * distance * denominator * denominator);
}

std::vector<double> coordination_numbers(const Structure& structure,
                                         const std::vector<Vec3>& translations,
                                         double cutoff)
{
    const std::size_t atoms = structure.positions.size();
    const double cutoff2 = cutoff * cutoff;
    std::vector<double> coordination(atoms, 0.0);
    for (std::size_t i = 0; i < atoms; ++i)
    {
        for (std::size_t j = 0; j <= i; ++j)
        {
            for (const Vec3& translation : translations)
            {
                const Vec3 vector = structure.positions[i] - (structure.positions[j] + translation);
                const double distance2 = dot(vector, vector);
                if (distance2 > cutoff2 || distance2 < 1.0e-12)
                {
                    continue;
                }
                const double count = coordination_count(structure.atomic_numbers[i],
                                                        structure.atomic_numbers[j],
                                                        std::sqrt(distance2));
                coordination[i] += count;
                if (i != j)
                {
                    coordination[j] += count;
                }
            }
        }
    }
    return coordination;
}

void reference_weights(const Structure& structure,
                       const std::vector<double>& coordination,
                       bool derivatives,
                       std::vector<double>& weights,
                       std::vector<double>& weight_derivatives)
{
    const std::size_t atoms = structure.positions.size();
    weights.assign(atoms * data::max_reference, 0.0);
    if (derivatives)
    {
        weight_derivatives.assign(atoms * data::max_reference, 0.0);
    }

    for (std::size_t atom = 0; atom < atoms; ++atom)
    {
        const int atomic_number = structure.atomic_numbers[atom];
        const int references = data::reference_count(atomic_number);
        double normalization = 0.0;
        double normalization_derivative = 0.0;
        double raw[data::max_reference] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        double maximum_cn = data::reference_cn(atomic_number, 0);
        for (int reference = 0; reference < references; ++reference)
        {
            const double reference_cn = data::reference_cn(atomic_number, reference);
            const double delta = reference_cn - coordination[atom];
            raw[reference] = std::exp(-kReferenceWeight * delta * delta);
            normalization += raw[reference];
            normalization_derivative += 2.0 * kReferenceWeight * delta * raw[reference];
            maximum_cn = std::max(maximum_cn, reference_cn);
        }
        normalization = 1.0 / normalization;

        for (int reference = 0; reference < references; ++reference)
        {
            const std::size_t index = atom * data::max_reference + reference;
            const double reference_cn = data::reference_cn(atomic_number, reference);
            double weight = raw[reference] * normalization;
            if (!std::isfinite(weight))
            {
                weight = reference_cn == maximum_cn ? 1.0 : 0.0;
            }
            weights[index] = weight;

            if (derivatives)
            {
                const double raw_derivative = 2.0 * kReferenceWeight
                                              * (reference_cn - coordination[atom]) * raw[reference];
                double derivative = raw_derivative * normalization
                                    - raw[reference] * normalization_derivative
                                          * normalization * normalization;
                if (!std::isfinite(derivative))
                {
                    derivative = 0.0;
                }
                weight_derivatives[index] = derivative;
            }
        }
    }
}

void atomic_c6(const Structure& structure,
               const std::vector<double>& weights,
               const std::vector<double>& weight_derivatives,
               bool derivatives,
               std::vector<double>& c6,
               std::vector<double>& dc6dcn)
{
    const std::size_t atoms = structure.positions.size();
    c6.assign(atoms * atoms, 0.0);
    if (derivatives)
    {
        dc6dcn.assign(atoms * atoms, 0.0);
    }

    for (std::size_t i = 0; i < atoms; ++i)
    {
        const int references_i = data::reference_count(structure.atomic_numbers[i]);
        for (std::size_t j = 0; j <= i; ++j)
        {
            const int references_j = data::reference_count(structure.atomic_numbers[j]);
            double coefficient = 0.0;
            double derivative_i = 0.0;
            double derivative_j = 0.0;
            for (int reference_i = 0; reference_i < references_i; ++reference_i)
            {
                for (int reference_j = 0; reference_j < references_j; ++reference_j)
                {
                    const double reference_c6 = data::reference_c6(structure.atomic_numbers[i],
                                                                  reference_i,
                                                                  structure.atomic_numbers[j],
                                                                  reference_j);
                    const double weight_i = weights[i * data::max_reference + reference_i];
                    const double weight_j = weights[j * data::max_reference + reference_j];
                    coefficient += weight_i * weight_j * reference_c6;
                    if (derivatives)
                    {
                        derivative_i += weight_derivatives[i * data::max_reference + reference_i]
                                        * weight_j * reference_c6;
                        derivative_j += weight_i
                                        * weight_derivatives[j * data::max_reference + reference_j]
                                        * reference_c6;
                    }
                }
            }
            c6[i * atoms + j] = coefficient;
            c6[j * atoms + i] = coefficient;
            if (derivatives)
            {
                dc6dcn[i * atoms + j] = derivative_i;
                dc6dcn[j * atoms + i] = derivative_j;
            }
        }
    }
}

void pairwise_dispersion(const Structure& structure,
                         const Parameters& parameters,
                         const Cutoffs& cutoffs,
                         const std::vector<Vec3>& translations,
                         const std::vector<double>& c6,
                         const std::vector<double>& dc6dcn,
                         bool derivatives,
                         std::vector<double>& atomic_energy,
                         std::vector<double>& dEdcn,
                         Result& result)
{
    const std::size_t atoms = structure.positions.size();
    const double cutoff2 = cutoffs.disp2 * cutoffs.disp2;
    const double epsilon = std::numeric_limits<double>::epsilon();

    for (std::size_t i = 0; i < atoms; ++i)
    {
        const int atomic_number_i = structure.atomic_numbers[i];
        for (std::size_t j = 0; j <= i; ++j)
        {
            const int atomic_number_j = structure.atomic_numbers[j];
            const double rrij = 3.0 * data::r4r2(atomic_number_i) * data::r4r2(atomic_number_j);
            const double coefficient = c6[j * atoms + i];

            double r0 = 0.0;
            double r0_sixth = 0.0;
            double r0_eighth = 0.0;
            if (parameters.damping == Damping::Rational)
            {
                r0 = parameters.a1 * std::sqrt(rrij) + parameters.a2;
                const double r0_squared = r0 * r0;
                r0_sixth = r0_squared * r0_squared * r0_squared;
                r0_eighth = r0_sixth * r0_squared;
            }
            else
            {
                r0 = data::vdw_radius(atomic_number_i, atomic_number_j);
            }

            for (const Vec3& translation : translations)
            {
                const Vec3 vector = structure.positions[i] - (structure.positions[j] + translation);
                const double distance2 = dot(vector, vector);
                if (distance2 > cutoff2 || distance2 < epsilon)
                {
                    continue;
                }
                const double distance = std::sqrt(distance2);
                double switching = 1.0;
                double switching_derivative = 0.0;
                smooth_cutoff(distance,
                              cutoffs.disp2,
                              cutoffs.width2,
                              switching,
                              switching_derivative);

                double bare_energy = 0.0;
                double bare_radial = 0.0;
                if (parameters.damping == Damping::Rational)
                {
                    const double fourth = distance2 * distance2;
                    const double inverse_sixth = 1.0 / (fourth * distance2 + r0_sixth);
                    const double inverse_eighth = 1.0 / (fourth * fourth + r0_eighth);
                    const double derivative_sixth = -6.0 * fourth * inverse_sixth * inverse_sixth;
                    const double derivative_eighth = -8.0 * fourth * distance2
                                                      * inverse_eighth * inverse_eighth;
                    bare_energy = parameters.s6 * inverse_sixth
                                  + parameters.s8 * rrij * inverse_eighth;
                    bare_radial = parameters.s6 * derivative_sixth
                                  + parameters.s8 * rrij * derivative_eighth;
                }
                else
                {
                    const double sixth = distance2 * distance2 * distance2;
                    const double eighth = sixth * distance2;
                    const double exponent6 = parameters.alp;
                    const double exponent8 = parameters.alp + 2.0;
                    const double t6 = std::pow(parameters.rs6 * r0 / distance, exponent6);
                    const double t8 = std::pow(parameters.rs8 * r0 / distance, exponent8);
                    const double f6 = 1.0 / (1.0 + 6.0 * t6);
                    const double f8 = 1.0 / (1.0 + 6.0 * t8);
                    const double derivative6 = -6.0 * f6 / distance2
                                               + 6.0 * exponent6 * t6 * f6 * f6 / distance2;
                    const double derivative8 = -8.0 * f8 / distance2
                                               + 6.0 * exponent8 * t8 * f8 * f8 / distance2;
                    bare_energy = parameters.s6 * f6 / sixth
                                  + parameters.s8 * rrij * f8 / eighth;
                    bare_radial = parameters.s6 * derivative6 / sixth
                                  + parameters.s8 * rrij * derivative8 / eighth;
                }

                const double dispersion = switching * bare_energy;
                const double energy = -coefficient * dispersion * 0.5;
                atomic_energy[i] += energy;
                if (i != j)
                {
                    atomic_energy[j] += energy;
                }

                if (derivatives)
                {
                    const double radial = switching * bare_radial
                                          + switching_derivative * bare_energy / distance;
                    const Vec3 gradient = vector * (-coefficient * radial);
                    dEdcn[i] -= dc6dcn[i * atoms + j] * dispersion;
                    add_outer(result.virial, gradient, vector, 0.5);
                    if (i != j)
                    {
                        dEdcn[j] -= dc6dcn[j * atoms + i] * dispersion;
                        result.gradient[i] += gradient;
                        result.gradient[j] -= gradient;
                        add_outer(result.virial, gradient, vector, 0.5);
                    }
                }
            }
        }
    }
}

double triple_scale(std::size_t i, std::size_t j, std::size_t k)
{
    if (i == j)
    {
        return i == k ? 1.0 / 6.0 : 0.5;
    }
    return i != k && j != k ? 1.0 : 0.5;
}

void atm_dispersion(const Structure& structure,
                    const Parameters& parameters,
                    const Cutoffs& cutoffs,
                    const std::vector<Vec3>& translations,
                    const std::vector<double>& c6,
                    const std::vector<double>& dc6dcn,
                    bool derivatives,
                    std::vector<double>& atomic_energy,
                    std::vector<double>& dEdcn,
                    Result& result)
{
    if (std::abs(parameters.s9) < std::numeric_limits<double>::epsilon())
    {
        return;
    }

    const std::size_t atoms = structure.positions.size();
    const double cutoff2 = cutoffs.disp3 * cutoffs.disp3;
    const double epsilon = std::numeric_limits<double>::epsilon();
    const double alpha = parameters.alp + 2.0;
    const double alpha_third = alpha / 3.0;

    for (std::size_t i = 0; i < atoms; ++i)
    {
        const int atomic_number_i = structure.atomic_numbers[i];
        for (std::size_t j = 0; j <= i; ++j)
        {
            const int atomic_number_j = structure.atomic_numbers[j];
            const double c6ij = c6[j * atoms + i];
            const double r0ij = kAtmRadiusScale * data::vdw_radius(atomic_number_j, atomic_number_i);
            for (const Vec3& translation_j : translations)
            {
                const Vec3 vij = structure.positions[j] + translation_j - structure.positions[i];
                const double r2ij = dot(vij, vij);
                if (r2ij > cutoff2 || r2ij < epsilon)
                {
                    continue;
                }
                const double rij = std::sqrt(r2ij);
                double swij = 1.0;
                double dswij = 0.0;
                smooth_cutoff(rij, cutoffs.disp3, cutoffs.width3, swij, dswij);

                for (std::size_t k = 0; k <= j; ++k)
                {
                    const int atomic_number_k = structure.atomic_numbers[k];
                    const double c6ik = c6[k * atoms + i];
                    const double c6jk = c6[k * atoms + j];
                    const double c9 = -parameters.s9 * std::sqrt(std::abs(c6ij * c6ik * c6jk));
                    const double r0ik = kAtmRadiusScale
                                        * data::vdw_radius(atomic_number_k, atomic_number_i);
                    const double r0jk = kAtmRadiusScale
                                        * data::vdw_radius(atomic_number_k, atomic_number_j);
                    const double r0 = r0ij * r0ik * r0jk;
                    const double scale = triple_scale(i, j, k);

                    for (const Vec3& translation_k : translations)
                    {
                        const Vec3 vik = structure.positions[k] + translation_k - structure.positions[i];
                        const double r2ik = dot(vik, vik);
                        if (r2ik > cutoff2 || r2ik < epsilon)
                        {
                            continue;
                        }
                        const double rik = std::sqrt(r2ik);
                        double swik = 1.0;
                        double dswik = 0.0;
                        smooth_cutoff(rik, cutoffs.disp3, cutoffs.width3, swik, dswik);

                        const Vec3 vjk = vik - vij;
                        const double r2jk = dot(vjk, vjk);
                        if (r2jk > cutoff2 || r2jk < epsilon)
                        {
                            continue;
                        }
                        const double rjk = std::sqrt(r2jk);
                        double swjk = 1.0;
                        double dswjk = 0.0;
                        smooth_cutoff(rjk, cutoffs.disp3, cutoffs.width3, swjk, dswjk);

                        const double switching = swij * swik * swjk;
                        const double r2 = r2ij * r2ik * r2jk;
                        const double r1 = std::sqrt(r2);
                        const double r3 = r2 * r1;
                        const double r5 = r3 * r2;
                        const double damping_power = std::pow(r0 / r1, alpha_third);
                        const double damping = 1.0 / (1.0 + 6.0 * damping_power);
                        const double angular = 0.375 * (r2ij + r2jk - r2ik)
                                               * (r2ij - r2jk + r2ik)
                                               * (-r2ij + r2jk + r2ik) / r5
                                               + 1.0 / r3;
                        const double rr = angular * damping;
                        const double base_energy = rr * c9;
                        const double energy = base_energy * scale * switching;
                        atomic_energy[i] -= energy / 3.0;
                        atomic_energy[j] -= energy / 3.0;
                        atomic_energy[k] -= energy / 3.0;

                        if (!derivatives)
                        {
                            continue;
                        }

                        const double damping_derivative = -2.0 * alpha * damping_power
                                                          * damping * damping;
                        double angular_derivative = -0.375
                            * (r2ij * r2ij * r2ij
                               + r2ij * r2ij * (r2jk + r2ik)
                               + r2ij * (3.0 * r2jk * r2jk + 2.0 * r2jk * r2ik
                                          + 3.0 * r2ik * r2ik)
                               - 5.0 * (r2jk - r2ik) * (r2jk - r2ik) * (r2jk + r2ik))
                            / r5;
                        const Vec3 gradient_ij
                            = vij * (switching * c9
                                     * (-angular_derivative * damping
                                        + angular * damping_derivative) / r2ij)
                              - vij * (base_energy * dswij / rij * swik * swjk);

                        angular_derivative = -0.375
                            * (r2ik * r2ik * r2ik
                               + r2ik * r2ik * (r2jk + r2ij)
                               + r2ik * (3.0 * r2jk * r2jk + 2.0 * r2jk * r2ij
                                          + 3.0 * r2ij * r2ij)
                               - 5.0 * (r2jk - r2ij) * (r2jk - r2ij) * (r2jk + r2ij))
                            / r5;
                        const Vec3 gradient_ik
                            = vik * (switching * c9
                                     * (-angular_derivative * damping
                                        + angular * damping_derivative) / r2ik)
                              - vik * (base_energy * dswik / rik * swij * swjk);

                        angular_derivative = -0.375
                            * (r2jk * r2jk * r2jk
                               + r2jk * r2jk * (r2ik + r2ij)
                               + r2jk * (3.0 * r2ik * r2ik + 2.0 * r2ik * r2ij
                                          + 3.0 * r2ij * r2ij)
                               - 5.0 * (r2ik - r2ij) * (r2ik - r2ij) * (r2ik + r2ij))
                            / r5;
                        const Vec3 gradient_jk
                            = vjk * (switching * c9
                                     * (-angular_derivative * damping
                                        + angular * damping_derivative) / r2jk)
                              - vjk * (base_energy * dswjk / rjk * swij * swik);

                        result.gradient[i] -= (gradient_ij + gradient_ik) * scale;
                        result.gradient[j] += (gradient_ij - gradient_jk) * scale;
                        result.gradient[k] += (gradient_ik + gradient_jk) * scale;
                        add_outer(result.virial, gradient_ij, vij, scale);
                        add_outer(result.virial, gradient_ik, vik, scale);
                        add_outer(result.virial, gradient_jk, vjk, scale);

                        dEdcn[i] -= energy * 0.5
                                    * (dc6dcn[i * atoms + j] / c6ij
                                       + dc6dcn[i * atoms + k] / c6ik);
                        dEdcn[j] -= energy * 0.5
                                    * (dc6dcn[j * atoms + i] / c6ij
                                       + dc6dcn[j * atoms + k] / c6jk);
                        dEdcn[k] -= energy * 0.5
                                    * (dc6dcn[k * atoms + i] / c6ik
                                       + dc6dcn[k * atoms + j] / c6jk);
                    }
                }
            }
        }
    }
}

void add_coordination_derivatives(const Structure& structure,
                                  const std::vector<Vec3>& translations,
                                  double cutoff,
                                  const std::vector<double>& dEdcn,
                                  Result& result)
{
    const std::size_t atoms = structure.positions.size();
    const double cutoff2 = cutoff * cutoff;
    for (std::size_t i = 0; i < atoms; ++i)
    {
        for (std::size_t j = 0; j <= i; ++j)
        {
            for (const Vec3& translation : translations)
            {
                const Vec3 vector = structure.positions[i] - (structure.positions[j] + translation);
                const double distance2 = dot(vector, vector);
                if (distance2 > cutoff2 || distance2 < 1.0e-12)
                {
                    continue;
                }
                const double distance = std::sqrt(distance2);
                const Vec3 count_derivative
                    = vector * (coordination_derivative(structure.atomic_numbers[i],
                                                        structure.atomic_numbers[j],
                                                        distance)
                                / distance);
                const double pair_derivative = dEdcn[i] + dEdcn[j];
                result.gradient[i] += count_derivative * pair_derivative;
                result.gradient[j] -= count_derivative * pair_derivative;
                const double strain_derivative = dEdcn[i] + (i == j ? 0.0 : dEdcn[j]);
                add_outer(result.virial, count_derivative, vector, strain_derivative);
            }
        }
    }
}

} // namespace

bool evaluate(const Structure& structure,
              const Parameters& parameters,
              const Cutoffs& cutoffs,
              bool derivatives,
              Result& result,
              std::string& error)
{
    error.clear();
    result = Result();
    if (!validate(structure, cutoffs, error))
    {
        return false;
    }

    const std::size_t atoms = structure.positions.size();
    if (derivatives)
    {
        result.gradient.assign(atoms, Vec3());
    }
    if (atoms == 0)
    {
        return true;
    }

    const std::vector<Vec3> cn_translations = lattice_points(structure, cutoffs.cn);
    const std::vector<double> coordination
        = coordination_numbers(structure, cn_translations, cutoffs.cn);

    std::vector<double> weights;
    std::vector<double> weight_derivatives;
    reference_weights(structure, coordination, derivatives, weights, weight_derivatives);

    std::vector<double> c6;
    std::vector<double> dc6dcn;
    atomic_c6(structure, weights, weight_derivatives, derivatives, c6, dc6dcn);

    std::vector<double> atomic_energy(atoms, 0.0);
    std::vector<double> dEdcn(derivatives ? atoms : 0, 0.0);
    const std::vector<Vec3> pair_translations = lattice_points(structure, cutoffs.disp2);
    pairwise_dispersion(structure,
                        parameters,
                        cutoffs,
                        pair_translations,
                        c6,
                        dc6dcn,
                        derivatives,
                        atomic_energy,
                        dEdcn,
                        result);

    const std::vector<Vec3> atm_translations = lattice_points(structure, cutoffs.disp3);
    atm_dispersion(structure,
                   parameters,
                   cutoffs,
                   atm_translations,
                   c6,
                   dc6dcn,
                   derivatives,
                   atomic_energy,
                   dEdcn,
                   result);

    if (derivatives)
    {
        add_coordination_derivatives(structure, cn_translations, cutoffs.cn, dEdcn, result);
    }
    result.energy = std::accumulate(atomic_energy.begin(), atomic_energy.end(), 0.0);
    return true;
}

} // namespace d3
} // namespace vdw
