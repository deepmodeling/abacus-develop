#include "pcc_2d.h"
#include "../common/lattice_row.h"

#include "source_base/constants.h"
#include "source_base/matrix3.h"

#include <cmath>
#include <stdexcept>

namespace ModulePcc
{
namespace
{

bool finite_position(const ModuleBase::Vector3<double>& position)
{
    return std::isfinite(position.x) && std::isfinite(position.y)
           && std::isfinite(position.z);
}

double norm(const ModuleBase::Vector3<double>& vector)
{
    const double squared_length = vector.x * vector.x + vector.y * vector.y + vector.z * vector.z;
    return std::sqrt(squared_length);
}

ModuleBase::Vector3<double> cross(const ModuleBase::Vector3<double>& left,
                                  const ModuleBase::Vector3<double>& right)
{
    const double component_x = left.y * right.z - left.z * right.y;
    const double component_y = left.z * right.x - left.x * right.z;
    const double component_z = left.x * right.y - left.y * right.x;
    return ModuleBase::Vector3<double>(component_x, component_y, component_z);
}

double dot(const ModuleBase::Vector3<double>& left, const ModuleBase::Vector3<double>& right)
{
    return left.x * right.x + left.y * right.y + left.z * right.z;
}

double potential_constant(const Pcc2dParameters& parameters)
{
    // Open planar kernel -2*pi*|u|/A, zero on the plane of the source, minus the
    // zero-mean periodic kernel, at u = 0. With this gauge the energy of a
    // charged slab converges with the vacuum size. ENVIRON core_1da
    // (calc_1da_vperiodic, dim=2) writes -pi*q/(3*L), which agrees only for
    // A = L^2 and adds a term linear in L to charged-slab energies.
    return -ModuleBase::PI * parameters.cell_length / (3.0 * parameters.periodic_area);
}

double inverse_volume_factor(const Pcc2dParameters& parameters)
{
    return 2.0 * ModuleBase::PI
           / (parameters.periodic_area * parameters.cell_length);
}

double minimum_image(const double displacement, const double cell_length)
{
    const double fractional_image = displacement / cell_length + 0.5;
    const double image = std::floor(fractional_image);
    return displacement - cell_length * image;
}

double wrap(const double coordinate, const double cell_length)
{
    const double fractional_coordinate = coordinate / cell_length;
    const double image = std::floor(fractional_coordinate);
    return coordinate - cell_length * image;
}

} // namespace

Pcc2dGeometry pcc_2d_geometry(const ModuleBase::Matrix3& lattice_vectors,
                              const double lattice_scale,
                              const int axis,
                              const double relative_tolerance)
{
    if (!std::isfinite(lattice_scale) || lattice_scale <= 0.0
        || !std::isfinite(relative_tolerance) || relative_tolerance <= 0.0)
    {
        throw std::invalid_argument(
            "two-dimensional PCC geometry requires positive finite scale and tolerance");
    }
    if (axis < 0 || axis > 2)
    {
        throw std::invalid_argument("two-dimensional PCC open axis must be 0, 1 or 2");
    }

    // The periodic vectors follow the open one cyclically, so their cross
    // product points along the open vector in a right-handed cell.
    const ModuleBase::Vector3<double> open
        = ModuleSurchem::lattice_row(lattice_vectors, axis, lattice_scale);
    const int first_axis = (axis + 1) % 3;
    const int second_axis = (axis + 2) % 3;
    const ModuleBase::Vector3<double> first
        = ModuleSurchem::lattice_row(lattice_vectors, first_axis, lattice_scale);
    const ModuleBase::Vector3<double> second
        = ModuleSurchem::lattice_row(lattice_vectors, second_axis, lattice_scale);
    if (!finite_position(open) || !finite_position(first) || !finite_position(second))
    {
        throw std::invalid_argument("two-dimensional PCC lattice vectors must be finite");
    }

    const double length_open = norm(open);
    const double length_first = norm(first);
    const double length_second = norm(second);
    if (!std::isfinite(length_open) || !std::isfinite(length_first) || !std::isfinite(length_second)
        || length_open <= 0.0 || length_first <= 0.0 || length_second <= 0.0)
    {
        throw std::invalid_argument(
            "two-dimensional PCC lattice vectors must have positive lengths");
    }

    const ModuleBase::Vector3<double> plane_normal = cross(first, second);
    const double periodic_area = norm(plane_normal);
    if (!std::isfinite(periodic_area)
        || periodic_area <= relative_tolerance * length_first * length_second)
    {
        throw std::invalid_argument(
            "two-dimensional PCC requires independent periodic lattice vectors");
    }
    const double orientation = dot(plane_normal, open) < 0.0 ? -1.0 : 1.0;
    const double normal_x = orientation * plane_normal.x / periodic_area;
    const double normal_y = orientation * plane_normal.y / periodic_area;
    const double normal_z = orientation * plane_normal.z / periodic_area;
    const ModuleBase::Vector3<double> normal(normal_x, normal_y, normal_z);
    const ModuleBase::Vector3<double> alignment = cross(open, normal);
    if (norm(alignment) > relative_tolerance * length_open)
    {
        throw std::invalid_argument(
            "two-dimensional PCC requires the open lattice vector perpendicular to the periodic ones");
    }

    Pcc2dGeometry geometry;
    geometry.parameters.periodic_area = periodic_area;
    geometry.parameters.cell_length = length_open;
    geometry.axis = axis;
    geometry.normal = normal;
    const double cell_center_x = 0.5 * (open.x + first.x + second.x);
    const double cell_center_y = 0.5 * (open.y + first.y + second.y);
    const double cell_center_z = 0.5 * (open.z + first.z + second.z);
    const ModuleBase::Vector3<double> cell_center(cell_center_x, cell_center_y, cell_center_z);
    geometry.origin = dot(cell_center, normal);
    validate_pcc_2d_parameters(geometry.parameters);
    return geometry;
}

void validate_pcc_2d_parameters(const Pcc2dParameters& parameters)
{
    if (!std::isfinite(parameters.periodic_area) || parameters.periodic_area <= 0.0)
    {
        throw std::invalid_argument("two-dimensional PCC requires a positive finite periodic area");
    }
    if (!std::isfinite(parameters.cell_length) || parameters.cell_length <= 0.0)
    {
        throw std::invalid_argument("two-dimensional PCC requires a positive finite open cell length");
    }
}

void validate_pcc_2d_geometry(const Pcc2dGeometry& geometry)
{
    validate_pcc_2d_parameters(geometry.parameters);
    if (geometry.axis < 0 || geometry.axis > 2)
    {
        throw std::invalid_argument("two-dimensional PCC open axis must be 0, 1 or 2");
    }
    const double normal_length_error = norm(geometry.normal) - 1.0;
    if (!finite_position(geometry.normal) || std::abs(normal_length_error) > 1.0e-10)
    {
        throw std::invalid_argument("two-dimensional PCC requires a finite unit normal");
    }
    if (!std::isfinite(geometry.origin))
    {
        throw std::invalid_argument("two-dimensional PCC requires a finite origin");
    }
}

double pcc_2d_coordinate(const ModuleBase::Vector3<double>& position,
                         const Pcc2dGeometry& geometry)
{
    if (!finite_position(position))
    {
        throw std::domain_error("two-dimensional PCC positions must be finite");
    }
    return dot(position, geometry.normal);
}

double pcc_2d_relative_coordinate(const ModuleBase::Vector3<double>& position,
                                  const Pcc2dGeometry& geometry)
{
    const double displacement = pcc_2d_coordinate(position, geometry) - geometry.origin;
    return minimum_image(displacement, geometry.parameters.cell_length);
}

double pcc_2d_system_center(const std::vector<double>& coordinates,
                            const std::vector<double>& weights,
                            const double cell_length)
{
    if (!std::isfinite(cell_length) || cell_length <= 0.0)
    {
        throw std::invalid_argument(
            "two-dimensional PCC system center requires a positive finite cell length");
    }
    if (coordinates.empty() || coordinates.size() != weights.size())
    {
        throw std::invalid_argument(
            "two-dimensional PCC system center requires matching non-empty positions and weights");
    }
    const double reference = coordinates[0];
    if (!std::isfinite(reference))
    {
        throw std::domain_error("two-dimensional PCC system-center positions must be finite");
    }
    double total_weight = 0.0;
    double weighted_displacement = 0.0;
    for (std::size_t index = 0; index < coordinates.size(); ++index)
    {
        if (!std::isfinite(coordinates[index]) || !std::isfinite(weights[index])
            || weights[index] <= 0.0)
        {
            throw std::domain_error(
                "two-dimensional PCC system-center positions and weights must be finite and positive");
        }
        total_weight += weights[index];
        weighted_displacement
            += weights[index] * minimum_image(coordinates[index] - reference, cell_length);
    }
    if (!std::isfinite(total_weight) || !std::isfinite(weighted_displacement))
    {
        throw std::domain_error("two-dimensional PCC system center must be finite");
    }
    return wrap(reference + weighted_displacement / total_weight, cell_length);
}

Pcc2dMoments pcc_2d_point_charge_moments(const std::vector<PointCharge>& charges,
                                          const Pcc2dGeometry& geometry)
{
    validate_pcc_2d_geometry(geometry);
    Pcc2dMoments moments;
    for (std::size_t index = 0; index < charges.size(); ++index)
    {
        const PointCharge& point = charges[index];
        if (!std::isfinite(point.charge))
        {
            throw std::domain_error("two-dimensional PCC point charges must be finite");
        }
        const double relative = pcc_2d_relative_coordinate(point.position, geometry);
        moments.charge += point.charge;
        moments.dipole += point.charge * relative;
        moments.quadrupole += point.charge * relative * relative;
    }
    return moments;
}

Pcc2dMoments pcc_2d_density_moments_from_relative_coordinates(
    const std::vector<double>& density,
    const std::vector<double>& relative_coordinates,
    const double volume_element)
{
    if (density.size() != relative_coordinates.size())
    {
        throw std::invalid_argument(
            "two-dimensional PCC density and relative coordinates must have the same size");
    }
    if (!std::isfinite(volume_element) || volume_element <= 0.0)
    {
        throw std::invalid_argument("two-dimensional PCC requires a positive finite volume element");
    }
    Pcc2dMoments moments;
    for (std::size_t index = 0; index < density.size(); ++index)
    {
        if (!std::isfinite(density[index]) || !std::isfinite(relative_coordinates[index]))
        {
            throw std::domain_error(
                "two-dimensional PCC density values and relative coordinates must be finite");
        }
        const double charge = density[index] * volume_element;
        moments.charge += charge;
        moments.dipole += charge * relative_coordinates[index];
        moments.quadrupole += charge * relative_coordinates[index] * relative_coordinates[index];
    }
    return moments;
}

Pcc2dMoments pcc_2d_sum_moments(const Pcc2dMoments& left, const Pcc2dMoments& right)
{
    Pcc2dMoments result;
    result.charge = left.charge + right.charge;
    result.dipole = left.dipole + right.dipole;
    result.quadrupole = left.quadrupole + right.quadrupole;
    return result;
}

Pcc2dMoments pcc_2d_density_moments(
    const std::vector<double>& density,
    const std::vector<ModuleBase::Vector3<double>>& positions,
    const double volume_element,
    const Pcc2dGeometry& geometry)
{
    if (density.size() != positions.size())
    {
        throw std::invalid_argument("two-dimensional PCC density and positions must have the same size");
    }
    validate_pcc_2d_geometry(geometry);
    std::vector<double> relative_coordinates(positions.size());
    for (std::size_t index = 0; index < positions.size(); ++index)
    {
        relative_coordinates[index] = pcc_2d_relative_coordinate(positions[index], geometry);
    }
    return pcc_2d_density_moments_from_relative_coordinates(density,
                                                            relative_coordinates,
                                                            volume_element);
}

double pcc_2d_potential(const Pcc2dMoments& moments,
                        const double relative_coordinate,
                        const Pcc2dParameters& parameters)
{
    const double parabolic = moments.charge * relative_coordinate * relative_coordinate
                             - 2.0 * moments.dipole * relative_coordinate
                             + moments.quadrupole;
    return potential_constant(parameters) * moments.charge
           - inverse_volume_factor(parameters) * parabolic;
}

ModuleBase::Vector3<double> pcc_2d_potential_gradient(
    const Pcc2dMoments& moments,
    const double relative_coordinate,
    const Pcc2dGeometry& geometry)
{
    const double derivative = -2.0 * inverse_volume_factor(geometry.parameters)
                              * (moments.charge * relative_coordinate - moments.dipole);
    const double component_x = derivative * geometry.normal.x;
    const double component_y = derivative * geometry.normal.y;
    const double component_z = derivative * geometry.normal.z;
    return ModuleBase::Vector3<double>(component_x, component_y, component_z);
}

ModuleBase::Vector3<double> pcc_2d_point_charge_force(
    const Pcc2dMoments& total_moments,
    const PointCharge& point,
    const Pcc2dGeometry& geometry)
{
    if (!std::isfinite(point.charge) || !finite_position(point.position))
    {
        throw std::domain_error("two-dimensional PCC point charge force inputs must be finite");
    }
    validate_pcc_2d_geometry(geometry);
    if (!std::isfinite(total_moments.charge) || !std::isfinite(total_moments.dipole)
        || !std::isfinite(total_moments.quadrupole))
    {
        throw std::domain_error("two-dimensional PCC force moments must be finite");
    }
    const double relative = pcc_2d_relative_coordinate(point.position, geometry);
    const ModuleBase::Vector3<double> gradient
        = pcc_2d_potential_gradient(total_moments, relative, geometry);
    const double component_x = -point.charge * gradient.x;
    const double component_y = -point.charge * gradient.y;
    const double component_z = -point.charge * gradient.z;
    return ModuleBase::Vector3<double>(component_x, component_y, component_z);
}

double pcc_2d_bilinear_energy(const Pcc2dMoments& left,
                              const Pcc2dMoments& right,
                              const Pcc2dParameters& parameters)
{
    validate_pcc_2d_parameters(parameters);
    if (!std::isfinite(left.charge) || !std::isfinite(left.dipole)
        || !std::isfinite(left.quadrupole) || !std::isfinite(right.charge)
        || !std::isfinite(right.dipole) || !std::isfinite(right.quadrupole))
    {
        throw std::domain_error("two-dimensional PCC energy moments must be finite");
    }

    const double multipole = left.charge * right.quadrupole
                             + left.quadrupole * right.charge
                             - 2.0 * left.dipole * right.dipole;
    return potential_constant(parameters) * left.charge * right.charge
           - inverse_volume_factor(parameters) * multipole;
}

double pcc_2d_self_energy(const Pcc2dMoments& moments,
                          const Pcc2dParameters& parameters)
{
    return 0.5 * pcc_2d_bilinear_energy(moments, moments, parameters);
}

} // namespace ModulePcc
