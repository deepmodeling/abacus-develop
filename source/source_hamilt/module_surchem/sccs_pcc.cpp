#include "sccs_pcc.h"

#include "source_base/constants.h"
#include "source_base/matrix3.h"

#include <cmath>
#include <stdexcept>

namespace ModuleSccs
{
namespace
{

double dot(const ModuleBase::Vector3<double>& left, const ModuleBase::Vector3<double>& right)
{
    return left.x * right.x + left.y * right.y + left.z * right.z;
}

double norm_squared(const ModuleBase::Vector3<double>& value)
{
    return dot(value, value);
}

double norm(const ModuleBase::Vector3<double>& value)
{
    return std::sqrt(norm_squared(value));
}

ModuleBase::Vector3<double> lattice_row(const ModuleBase::Matrix3& lattice,
                                        const int row,
                                        const double scale)
{
    if (row == 0)
    {
        return ModuleBase::Vector3<double>(scale * lattice.e11,
                                           scale * lattice.e12,
                                           scale * lattice.e13);
    }
    if (row == 1)
    {
        return ModuleBase::Vector3<double>(scale * lattice.e21,
                                           scale * lattice.e22,
                                           scale * lattice.e23);
    }
    return ModuleBase::Vector3<double>(scale * lattice.e31,
                                       scale * lattice.e32,
                                       scale * lattice.e33);
}

ModuleBase::Vector3<double> scaled(const ModuleBase::Vector3<double>& value,
                                   const double factor)
{
    return ModuleBase::Vector3<double>(factor * value.x,
                                       factor * value.y,
                                       factor * value.z);
}

ModuleBase::Vector3<double> add(const ModuleBase::Vector3<double>& left,
                                const ModuleBase::Vector3<double>& right)
{
    return ModuleBase::Vector3<double>(left.x + right.x,
                                       left.y + right.y,
                                       left.z + right.z);
}

bool finite_vector(const ModuleBase::Vector3<double>& value)
{
    return std::isfinite(value.x) && std::isfinite(value.y)
           && std::isfinite(value.z);
}

double minimum_image(const double displacement, const double length)
{
    return displacement - length * std::floor(displacement / length + 0.5);
}

MultipoleMoments moments_from_relative_positions(
    const std::vector<PointCharge>& charges,
    const std::vector<ModuleBase::Vector3<double>>& relative_positions)
{
    if (charges.size() != relative_positions.size())
    {
        throw std::invalid_argument(
            "PCC charges and relative positions must have the same size");
    }
    MultipoleMoments moments;
    for (std::size_t index = 0; index < charges.size(); ++index)
    {
        const PointCharge& point = charges[index];
        const ModuleBase::Vector3<double>& relative = relative_positions[index];
        if (!std::isfinite(point.charge) || !finite_vector(point.position)
            || !finite_vector(relative))
        {
            throw std::domain_error("PCC point charges and positions must be finite");
        }
        moments.charge += point.charge;
        moments.dipole.x += point.charge * relative.x;
        moments.dipole.y += point.charge * relative.y;
        moments.dipole.z += point.charge * relative.z;
        moments.quadrupole_trace += point.charge * norm_squared(relative);
    }
    return moments;
}

} // namespace

PccGeometry pcc_geometry(const ModuleBase::Matrix3& lattice_vectors,
                         const double lattice_scale,
                         const double relative_tolerance)
{
    if (!std::isfinite(lattice_scale) || lattice_scale <= 0.0
        || !std::isfinite(relative_tolerance) || relative_tolerance <= 0.0)
    {
        throw std::invalid_argument(
            "zero-dimensional PCC geometry requires positive finite scale and tolerance");
    }
    const ModuleBase::Vector3<double> a
        = lattice_row(lattice_vectors, 0, lattice_scale);
    const ModuleBase::Vector3<double> b
        = lattice_row(lattice_vectors, 1, lattice_scale);
    const ModuleBase::Vector3<double> c
        = lattice_row(lattice_vectors, 2, lattice_scale);
    const double length_a = norm(a);
    const double length_b = norm(b);
    const double length_c = norm(c);
    if (!finite_vector(a) || !finite_vector(b) || !finite_vector(c)
        || !std::isfinite(length_a) || !std::isfinite(length_b)
        || !std::isfinite(length_c) || length_a <= 0.0 || length_b <= 0.0
        || length_c <= 0.0)
    {
        throw std::invalid_argument(
            "zero-dimensional PCC lattice vectors must have positive finite lengths");
    }
    const double length = (length_a + length_b + length_c) / 3.0;
    if (std::abs(length_a - length) > relative_tolerance * length
        || std::abs(length_b - length) > relative_tolerance * length
        || std::abs(length_c - length) > relative_tolerance * length
        || std::abs(dot(a, b)) > relative_tolerance * length * length
        || std::abs(dot(a, c)) > relative_tolerance * length * length
        || std::abs(dot(b, c)) > relative_tolerance * length * length)
    {
        throw std::invalid_argument(
            "zero-dimensional SCCS PCC requires an orthogonal equal-edge cubic cell");
    }

    PccGeometry geometry;
    geometry.parameters.cube_length = length;
    geometry.origin = scaled(add(add(a, b), c), 0.5);
    geometry.axis_a = scaled(a, 1.0 / length_a);
    geometry.axis_b = scaled(b, 1.0 / length_b);
    geometry.axis_c = scaled(c, 1.0 / length_c);
    validate_pcc_geometry(geometry);
    return geometry;
}

void validate_pcc_parameters(const PccParameters& parameters)
{
    if (!std::isfinite(parameters.cube_length) || parameters.cube_length <= 0.0)
    {
        throw std::invalid_argument("zero-dimensional PCC requires a positive finite cube length");
    }
    if (!std::isfinite(parameters.madelung) || parameters.madelung <= 0.0)
    {
        throw std::invalid_argument("zero-dimensional PCC requires a positive finite Madelung constant");
    }
}

void validate_pcc_geometry(const PccGeometry& geometry)
{
    validate_pcc_parameters(geometry.parameters);
    if (!finite_vector(geometry.origin) || !finite_vector(geometry.axis_a)
        || !finite_vector(geometry.axis_b) || !finite_vector(geometry.axis_c))
    {
        throw std::invalid_argument("zero-dimensional PCC geometry must be finite");
    }
    const double tolerance = 1.0e-10;
    if (std::abs(norm_squared(geometry.axis_a) - 1.0) > tolerance
        || std::abs(norm_squared(geometry.axis_b) - 1.0) > tolerance
        || std::abs(norm_squared(geometry.axis_c) - 1.0) > tolerance
        || std::abs(dot(geometry.axis_a, geometry.axis_b)) > tolerance
        || std::abs(dot(geometry.axis_a, geometry.axis_c)) > tolerance
        || std::abs(dot(geometry.axis_b, geometry.axis_c)) > tolerance)
    {
        throw std::invalid_argument(
            "zero-dimensional PCC geometry requires orthonormal lattice axes");
    }
}

ModuleBase::Vector3<double> pcc_relative_position(
    const ModuleBase::Vector3<double>& position,
    const PccGeometry& geometry)
{
    validate_pcc_geometry(geometry);
    if (!finite_vector(position))
    {
        throw std::domain_error("zero-dimensional PCC positions must be finite");
    }
    const ModuleBase::Vector3<double> displacement(position.x - geometry.origin.x,
                                                    position.y - geometry.origin.y,
                                                    position.z - geometry.origin.z);
    const double length = geometry.parameters.cube_length;
    const double a = minimum_image(dot(displacement, geometry.axis_a), length);
    const double b = minimum_image(dot(displacement, geometry.axis_b), length);
    const double c = minimum_image(dot(displacement, geometry.axis_c), length);
    return add(add(scaled(geometry.axis_a, a), scaled(geometry.axis_b, b)),
               scaled(geometry.axis_c, c));
}

ModuleBase::Vector3<double> pcc_system_center(
    const std::vector<ModuleBase::Vector3<double>>& positions,
    const std::vector<double>& weights,
    const PccGeometry& geometry)
{
    validate_pcc_geometry(geometry);
    if (positions.empty() || positions.size() != weights.size())
    {
        throw std::invalid_argument(
            "zero-dimensional PCC system center requires matching non-empty positions and weights");
    }
    if (!finite_vector(positions[0]))
    {
        throw std::domain_error(
            "zero-dimensional PCC system-center positions must be finite");
    }
    PccGeometry reference_geometry = geometry;
    reference_geometry.origin = positions[0];
    double total_weight = 0.0;
    ModuleBase::Vector3<double> weighted_displacement;
    for (std::size_t index = 0; index < positions.size(); ++index)
    {
        if (!finite_vector(positions[index]) || !std::isfinite(weights[index])
            || weights[index] <= 0.0)
        {
            throw std::domain_error(
                "zero-dimensional PCC system-center positions and weights must be finite and positive");
        }
        const ModuleBase::Vector3<double> relative
            = pcc_relative_position(positions[index], reference_geometry);
        total_weight += weights[index];
        weighted_displacement.x += weights[index] * relative.x;
        weighted_displacement.y += weights[index] * relative.y;
        weighted_displacement.z += weights[index] * relative.z;
    }
    if (!std::isfinite(total_weight) || !finite_vector(weighted_displacement))
    {
        throw std::domain_error("zero-dimensional PCC system center must be finite");
    }
    const ModuleBase::Vector3<double> unwrapped_center(
        positions[0].x + weighted_displacement.x / total_weight,
        positions[0].y + weighted_displacement.y / total_weight,
        positions[0].z + weighted_displacement.z / total_weight);
    const ModuleBase::Vector3<double> centered
        = pcc_relative_position(unwrapped_center, geometry);
    return ModuleBase::Vector3<double>(geometry.origin.x + centered.x,
                                       geometry.origin.y + centered.y,
                                       geometry.origin.z + centered.z);
}

MultipoleMoments point_charge_moments(const std::vector<PointCharge>& charges,
                                      const ModuleBase::Vector3<double>& origin)
{
    std::vector<ModuleBase::Vector3<double>> relative_positions(charges.size());
    for (std::size_t index = 0; index < charges.size(); ++index)
    {
        const PointCharge& point = charges[index];
        if (!std::isfinite(point.charge) || !std::isfinite(point.position.x)
            || !std::isfinite(point.position.y) || !std::isfinite(point.position.z))
        {
            throw std::domain_error("PCC point charges and positions must be finite");
        }
        relative_positions[index]
            = ModuleBase::Vector3<double>(point.position.x - origin.x,
                                          point.position.y - origin.y,
                                          point.position.z - origin.z);
    }
    return moments_from_relative_positions(charges, relative_positions);
}

MultipoleMoments point_charge_moments(const std::vector<PointCharge>& charges,
                                      const PccGeometry& geometry)
{
    std::vector<ModuleBase::Vector3<double>> relative_positions(charges.size());
    for (std::size_t index = 0; index < charges.size(); ++index)
    {
        relative_positions[index] = pcc_relative_position(charges[index].position,
                                                           geometry);
    }
    return moments_from_relative_positions(charges, relative_positions);
}

MultipoleMoments density_moments(const std::vector<double>& density,
                                 const std::vector<ModuleBase::Vector3<double>>& positions,
                                 const double volume_element,
                                 const ModuleBase::Vector3<double>& origin)
{
    if (density.size() != positions.size())
    {
        throw std::invalid_argument("PCC density and position arrays must have the same size");
    }
    if (!std::isfinite(volume_element) || volume_element <= 0.0)
    {
        throw std::invalid_argument("PCC density integration requires a positive finite volume element");
    }

    std::vector<PointCharge> charges(density.size());
    for (std::size_t index = 0; index < density.size(); ++index)
    {
        charges[index].charge = density[index] * volume_element;
        charges[index].position = positions[index];
    }
    return point_charge_moments(charges, origin);
}

MultipoleMoments density_moments(const std::vector<double>& density,
                                 const std::vector<ModuleBase::Vector3<double>>& positions,
                                 const double volume_element,
                                 const PccGeometry& geometry)
{
    if (density.size() != positions.size())
    {
        throw std::invalid_argument("PCC density and position arrays must have the same size");
    }
    std::vector<ModuleBase::Vector3<double>> relative_positions(positions.size());
    for (std::size_t index = 0; index < positions.size(); ++index)
    {
        relative_positions[index] = pcc_relative_position(positions[index], geometry);
    }
    return density_moments_from_relative_positions(density,
                                                    relative_positions,
                                                    volume_element);
}

MultipoleMoments density_moments_from_relative_positions(
    const std::vector<double>& density,
    const std::vector<ModuleBase::Vector3<double>>& relative_positions,
    const double volume_element)
{
    if (density.size() != relative_positions.size())
    {
        throw std::invalid_argument(
            "PCC density and relative-position arrays must have the same size");
    }
    if (!std::isfinite(volume_element) || volume_element <= 0.0)
    {
        throw std::invalid_argument(
            "PCC density integration requires a positive finite volume element");
    }
    std::vector<PointCharge> charges(density.size());
    for (std::size_t index = 0; index < density.size(); ++index)
    {
        charges[index].charge = density[index] * volume_element;
        charges[index].position = relative_positions[index];
    }
    return moments_from_relative_positions(charges, relative_positions);
}

double pcc_potential(const MultipoleMoments& moments,
                     const ModuleBase::Vector3<double>& position,
                     const PccParameters& parameters)
{
    validate_pcc_parameters(parameters);
    const double volume = parameters.cube_length * parameters.cube_length * parameters.cube_length;
    const double parabolic = moments.charge * norm_squared(position)
                             - 2.0 * dot(moments.dipole, position)
                             + moments.quadrupole_trace;
    return parameters.madelung * moments.charge / parameters.cube_length
           - 2.0 * ModuleBase::PI * parabolic / (3.0 * volume);
}

ModuleBase::Vector3<double> pcc_potential_gradient(const MultipoleMoments& moments,
                                                   const ModuleBase::Vector3<double>& position,
                                                   const PccParameters& parameters)
{
    validate_pcc_parameters(parameters);
    const double volume = parameters.cube_length * parameters.cube_length * parameters.cube_length;
    const double factor = -4.0 * ModuleBase::PI / (3.0 * volume);
    return ModuleBase::Vector3<double>(factor * (moments.charge * position.x - moments.dipole.x),
                                       factor * (moments.charge * position.y - moments.dipole.y),
                                       factor * (moments.charge * position.z - moments.dipole.z));
}

ModuleBase::Vector3<double> pcc_point_charge_force(
    const MultipoleMoments& total_moments,
    const PointCharge& point,
    const ModuleBase::Vector3<double>& origin,
    const PccParameters& parameters)
{
    if (!std::isfinite(point.charge) || !std::isfinite(point.position.x)
        || !std::isfinite(point.position.y) || !std::isfinite(point.position.z)
        || !std::isfinite(origin.x) || !std::isfinite(origin.y) || !std::isfinite(origin.z))
    {
        throw std::domain_error("PCC point charge force requires finite charge and coordinates");
    }
    const ModuleBase::Vector3<double> relative(point.position.x - origin.x,
                                                point.position.y - origin.y,
                                                point.position.z - origin.z);
    const ModuleBase::Vector3<double> gradient
        = pcc_potential_gradient(total_moments, relative, parameters);
    return ModuleBase::Vector3<double>(-point.charge * gradient.x,
                                       -point.charge * gradient.y,
                                       -point.charge * gradient.z);
}

ModuleBase::Vector3<double> pcc_point_charge_force(
    const MultipoleMoments& total_moments,
    const PointCharge& point,
    const PccGeometry& geometry)
{
    if (!std::isfinite(point.charge) || !finite_vector(point.position))
    {
        throw std::domain_error("PCC point charge force requires finite charge and coordinates");
    }
    const ModuleBase::Vector3<double> gradient
        = pcc_potential_gradient(total_moments,
                                 pcc_relative_position(point.position, geometry),
                                 geometry.parameters);
    return ModuleBase::Vector3<double>(-point.charge * gradient.x,
                                       -point.charge * gradient.y,
                                       -point.charge * gradient.z);
}

double pcc_bilinear_energy(const MultipoleMoments& left,
                           const MultipoleMoments& right,
                           const PccParameters& parameters)
{
    validate_pcc_parameters(parameters);
    const double volume = parameters.cube_length * parameters.cube_length * parameters.cube_length;
    const double monopole = parameters.madelung * left.charge * right.charge / parameters.cube_length;
    const double multipole = left.quadrupole_trace * right.charge
                             + left.charge * right.quadrupole_trace
                             - 2.0 * dot(left.dipole, right.dipole);
    return monopole - 2.0 * ModuleBase::PI * multipole / (3.0 * volume);
}

double pcc_self_energy(const MultipoleMoments& moments, const PccParameters& parameters)
{
    return 0.5 * pcc_bilinear_energy(moments, moments, parameters);
}

} // namespace ModuleSccs
