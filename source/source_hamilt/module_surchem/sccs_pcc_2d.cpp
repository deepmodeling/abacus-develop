#include "sccs_pcc_2d.h"

#include "source_base/constants.h"
#include "source_base/matrix3.h"

#include <cmath>
#include <stdexcept>

namespace ModuleSccs
{
namespace
{

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

bool finite_position(const ModuleBase::Vector3<double>& position)
{
    return std::isfinite(position.x) && std::isfinite(position.y)
           && std::isfinite(position.z);
}

double norm(const ModuleBase::Vector3<double>& vector)
{
    return std::sqrt(vector.x * vector.x + vector.y * vector.y + vector.z * vector.z);
}

ModuleBase::Vector3<double> cross(const ModuleBase::Vector3<double>& left,
                                  const ModuleBase::Vector3<double>& right)
{
    return ModuleBase::Vector3<double>(left.y * right.z - left.z * right.y,
                                       left.z * right.x - left.x * right.z,
                                       left.x * right.y - left.y * right.x);
}

double potential_constant(const Pcc2dParameters& parameters)
{
    // Andreussi-Marzari, Phys. Rev. B 90, 245101 (2014), Eq. (88).
    return ModuleBase::PI / (3.0 * parameters.cell_length_y);
}

double inverse_volume_factor(const Pcc2dParameters& parameters)
{
    return 2.0 * ModuleBase::PI
           / (parameters.periodic_area * parameters.cell_length_y);
}

double minimum_image_y(const double displacement_y, const double cell_length_y)
{
    return displacement_y
           - cell_length_y * std::floor(displacement_y / cell_length_y + 0.5);
}

double wrap_y(const double position_y, const double cell_length_y)
{
    return position_y - cell_length_y * std::floor(position_y / cell_length_y);
}

void validate_pcc_2d_geometry(const Pcc2dGeometry& geometry)
{
    validate_pcc_2d_parameters(geometry.parameters);
    if (!std::isfinite(geometry.origin_y))
    {
        throw std::invalid_argument("two-dimensional PCC requires a finite y origin");
    }
}

} // namespace

Pcc2dGeometry pcc_2d_geometry(const ModuleBase::Matrix3& lattice_vectors,
                              const double lattice_scale,
                              const double relative_tolerance)
{
    if (!std::isfinite(lattice_scale) || lattice_scale <= 0.0
        || !std::isfinite(relative_tolerance) || relative_tolerance <= 0.0)
    {
        throw std::invalid_argument(
            "two-dimensional PCC geometry requires positive finite scale and tolerance");
    }

    const ModuleBase::Vector3<double> a
        = lattice_row(lattice_vectors, 0, lattice_scale);
    const ModuleBase::Vector3<double> b
        = lattice_row(lattice_vectors, 1, lattice_scale);
    const ModuleBase::Vector3<double> c
        = lattice_row(lattice_vectors, 2, lattice_scale);
    if (!finite_position(a) || !finite_position(b) || !finite_position(c))
    {
        throw std::invalid_argument("two-dimensional PCC lattice vectors must be finite");
    }

    const double length_a = norm(a);
    const double length_b = norm(b);
    const double length_c = norm(c);
    if (!std::isfinite(length_a) || !std::isfinite(length_b) || !std::isfinite(length_c)
        || length_a <= 0.0 || length_b <= 0.0 || length_c <= 0.0)
    {
        throw std::invalid_argument(
            "two-dimensional PCC lattice vectors must have positive lengths");
    }
    if (b.y <= 0.0 || std::abs(a.y) > relative_tolerance * length_a
        || std::abs(c.y) > relative_tolerance * length_c
        || std::abs(b.x) > relative_tolerance * length_b
        || std::abs(b.z) > relative_tolerance * length_b)
    {
        throw std::invalid_argument(
            "two-dimensional PCC requires a and c in the x-z plane and b along +y");
    }

    const double periodic_area = norm(cross(a, c));
    if (!std::isfinite(periodic_area)
        || periodic_area <= relative_tolerance * length_a * length_c)
    {
        throw std::invalid_argument(
            "two-dimensional PCC requires independent in-plane lattice vectors");
    }

    Pcc2dGeometry geometry;
    geometry.parameters.periodic_area = periodic_area;
    geometry.parameters.cell_length_y = length_b;
    geometry.origin_y = 0.5 * (a.y + b.y + c.y);
    validate_pcc_2d_parameters(geometry.parameters);
    return geometry;
}

void validate_pcc_2d_parameters(const Pcc2dParameters& parameters)
{
    if (!std::isfinite(parameters.periodic_area) || parameters.periodic_area <= 0.0)
    {
        throw std::invalid_argument("two-dimensional PCC requires a positive finite periodic area");
    }
    if (!std::isfinite(parameters.cell_length_y) || parameters.cell_length_y <= 0.0)
    {
        throw std::invalid_argument("two-dimensional PCC requires a positive finite y cell length");
    }
}

double pcc_2d_relative_y(const double position_y, const Pcc2dGeometry& geometry)
{
    validate_pcc_2d_geometry(geometry);
    if (!std::isfinite(position_y))
    {
        throw std::domain_error("two-dimensional PCC y positions must be finite");
    }
    return minimum_image_y(position_y - geometry.origin_y,
                           geometry.parameters.cell_length_y);
}

double pcc_2d_system_center_y(const std::vector<double>& positions_y,
                              const std::vector<double>& weights,
                              const double cell_length_y)
{
    if (!std::isfinite(cell_length_y) || cell_length_y <= 0.0)
    {
        throw std::invalid_argument(
            "two-dimensional PCC system center requires a positive finite cell length");
    }
    if (positions_y.empty() || positions_y.size() != weights.size())
    {
        throw std::invalid_argument(
            "two-dimensional PCC system center requires matching non-empty positions and weights");
    }
    const double reference_y = positions_y[0];
    if (!std::isfinite(reference_y))
    {
        throw std::domain_error("two-dimensional PCC system-center positions must be finite");
    }
    double total_weight = 0.0;
    double weighted_displacement = 0.0;
    for (std::size_t index = 0; index < positions_y.size(); ++index)
    {
        if (!std::isfinite(positions_y[index]) || !std::isfinite(weights[index])
            || weights[index] <= 0.0)
        {
            throw std::domain_error(
                "two-dimensional PCC system-center positions and weights must be finite and positive");
        }
        total_weight += weights[index];
        weighted_displacement
            += weights[index]
               * minimum_image_y(positions_y[index] - reference_y, cell_length_y);
    }
    if (!std::isfinite(total_weight) || !std::isfinite(weighted_displacement))
    {
        throw std::domain_error("two-dimensional PCC system center must be finite");
    }
    return wrap_y(reference_y + weighted_displacement / total_weight,
                  cell_length_y);
}

Pcc2dMoments pcc_2d_point_charge_moments(const std::vector<PointCharge>& charges,
                                          const Pcc2dGeometry& geometry)
{
    validate_pcc_2d_geometry(geometry);
    Pcc2dMoments moments;
    for (std::size_t index = 0; index < charges.size(); ++index)
    {
        const PointCharge& point = charges[index];
        if (!std::isfinite(point.charge) || !finite_position(point.position))
        {
            throw std::domain_error("two-dimensional PCC point charges must be finite");
        }
        const double relative_y
            = minimum_image_y(point.position.y - geometry.origin_y,
                              geometry.parameters.cell_length_y);
        moments.charge += point.charge;
        moments.dipole_y += point.charge * relative_y;
        moments.quadrupole_yy += point.charge * relative_y * relative_y;
    }
    return moments;
}

Pcc2dMoments pcc_2d_density_moments_from_relative_y(
    const std::vector<double>& density,
    const std::vector<double>& relative_y,
    const double volume_element)
{
    if (density.size() != relative_y.size())
    {
        throw std::invalid_argument(
            "two-dimensional PCC density and relative-y coordinates must have the same size");
    }
    if (!std::isfinite(volume_element) || volume_element <= 0.0)
    {
        throw std::invalid_argument("two-dimensional PCC requires a positive finite volume element");
    }
    Pcc2dMoments moments;
    for (std::size_t index = 0; index < density.size(); ++index)
    {
        if (!std::isfinite(density[index]) || !std::isfinite(relative_y[index]))
        {
            throw std::domain_error(
                "two-dimensional PCC density values and relative-y coordinates must be finite");
        }
        const double charge = density[index] * volume_element;
        moments.charge += charge;
        moments.dipole_y += charge * relative_y[index];
        moments.quadrupole_yy += charge * relative_y[index] * relative_y[index];
    }
    return moments;
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
    std::vector<double> relative_y(positions.size());
    for (std::size_t index = 0; index < positions.size(); ++index)
    {
        if (!finite_position(positions[index]))
        {
            throw std::domain_error("two-dimensional PCC positions must be finite");
        }
        relative_y[index]
            = minimum_image_y(positions[index].y - geometry.origin_y,
                              geometry.parameters.cell_length_y);
    }
    return pcc_2d_density_moments_from_relative_y(density,
                                                   relative_y,
                                                   volume_element);
}

double pcc_2d_potential(const Pcc2dMoments& moments,
                        const double relative_y,
                        const Pcc2dParameters& parameters)
{
    validate_pcc_2d_parameters(parameters);
    if (!std::isfinite(moments.charge) || !std::isfinite(moments.dipole_y)
        || !std::isfinite(moments.quadrupole_yy) || !std::isfinite(relative_y))
    {
        throw std::domain_error("two-dimensional PCC potential inputs must be finite");
    }

    const double parabolic = moments.charge * relative_y * relative_y
                             - 2.0 * moments.dipole_y * relative_y
                             + moments.quadrupole_yy;
    return potential_constant(parameters) * moments.charge
           - inverse_volume_factor(parameters) * parabolic;
}

ModuleBase::Vector3<double> pcc_2d_potential_gradient(
    const Pcc2dMoments& moments,
    const double relative_y,
    const Pcc2dParameters& parameters)
{
    validate_pcc_2d_parameters(parameters);
    if (!std::isfinite(moments.charge) || !std::isfinite(moments.dipole_y)
        || !std::isfinite(moments.quadrupole_yy) || !std::isfinite(relative_y))
    {
        throw std::domain_error("two-dimensional PCC gradient inputs must be finite");
    }

    const double gradient_y = -2.0 * inverse_volume_factor(parameters)
                              * (moments.charge * relative_y - moments.dipole_y);
    return ModuleBase::Vector3<double>(0.0, gradient_y, 0.0);
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
    const ModuleBase::Vector3<double> gradient
        = pcc_2d_potential_gradient(total_moments,
                                    pcc_2d_relative_y(point.position.y, geometry),
                                    geometry.parameters);
    return ModuleBase::Vector3<double>(-point.charge * gradient.x,
                                       -point.charge * gradient.y,
                                       -point.charge * gradient.z);
}

double pcc_2d_bilinear_energy(const Pcc2dMoments& left,
                              const Pcc2dMoments& right,
                              const Pcc2dParameters& parameters)
{
    validate_pcc_2d_parameters(parameters);
    if (!std::isfinite(left.charge) || !std::isfinite(left.dipole_y)
        || !std::isfinite(left.quadrupole_yy) || !std::isfinite(right.charge)
        || !std::isfinite(right.dipole_y) || !std::isfinite(right.quadrupole_yy))
    {
        throw std::domain_error("two-dimensional PCC energy moments must be finite");
    }

    const double multipole = left.charge * right.quadrupole_yy
                             + left.quadrupole_yy * right.charge
                             - 2.0 * left.dipole_y * right.dipole_y;
    return potential_constant(parameters) * left.charge * right.charge
           - inverse_volume_factor(parameters) * multipole;
}

double pcc_2d_self_energy(const Pcc2dMoments& moments,
                          const Pcc2dParameters& parameters)
{
    return 0.5 * pcc_2d_bilinear_energy(moments, moments, parameters);
}

double pcc_2d_ionic_shape_energy(const double polarization_charge,
                                 const Pcc2dMoments& smooth_ionic_moments,
                                 const Pcc2dMoments& point_ionic_moments,
                                 const Pcc2dParameters& parameters)
{
    validate_pcc_2d_parameters(parameters);
    if (!std::isfinite(polarization_charge)
        || !std::isfinite(smooth_ionic_moments.charge)
        || !std::isfinite(smooth_ionic_moments.dipole_y)
        || !std::isfinite(smooth_ionic_moments.quadrupole_yy)
        || !std::isfinite(point_ionic_moments.charge)
        || !std::isfinite(point_ionic_moments.dipole_y)
        || !std::isfinite(point_ionic_moments.quadrupole_yy))
    {
        throw std::domain_error("two-dimensional PCC ionic-shape energy inputs must be finite");
    }
    if (point_ionic_moments.charge == 0.0)
    {
        throw std::domain_error(
            "two-dimensional PCC ionic-shape energy requires non-zero ionic charge");
    }
    const double ionic_center_y
        = point_ionic_moments.dipole_y / point_ionic_moments.charge;
    const double charge_difference
        = smooth_ionic_moments.charge - point_ionic_moments.charge;
    const double dipole_difference
        = smooth_ionic_moments.dipole_y - point_ionic_moments.dipole_y;
    const double quadrupole_difference
        = smooth_ionic_moments.quadrupole_yy
          - point_ionic_moments.quadrupole_yy;
    const double centered_quadrupole_difference
        = quadrupole_difference - 2.0 * ionic_center_y * dipole_difference
          + ionic_center_y * ionic_center_y * charge_difference;
    // Andreussi-Marzari, Phys. Rev. B 90, 245101 (2014), Eq. (A2).
    // Its coordinate origin is the center of ionic charge.
    return ModuleBase::PI * polarization_charge
           * centered_quadrupole_difference
           / (parameters.periodic_area * parameters.cell_length_y);
}

} // namespace ModuleSccs
