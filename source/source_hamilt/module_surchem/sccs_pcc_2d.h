#ifndef SCCS_PCC_2D_H
#define SCCS_PCC_2D_H

#include "sccs_pcc.h"

#include <vector>

namespace ModuleBase
{
class Matrix3;
}

namespace ModuleSccs
{

struct Pcc2dMoments
{
    double charge = 0.0;
    double dipole_y = 0.0;
    double quadrupole_yy = 0.0;
};

struct Pcc2dParameters
{
    double periodic_area = 0.0;
    double cell_length_y = 0.0;
};

struct Pcc2dGeometry
{
    Pcc2dParameters parameters;
    double origin_y = 0.0;
};

Pcc2dGeometry pcc_2d_geometry(const ModuleBase::Matrix3& lattice_vectors,
                              double lattice_scale,
                              double relative_tolerance);

void validate_pcc_2d_parameters(const Pcc2dParameters& parameters);

double pcc_2d_relative_y(double position_y, const Pcc2dGeometry& geometry);

double pcc_2d_system_center_y(const std::vector<double>& positions_y,
                              const std::vector<double>& weights,
                              double cell_length_y);

Pcc2dMoments pcc_2d_point_charge_moments(const std::vector<PointCharge>& charges,
                                          const Pcc2dGeometry& geometry);

Pcc2dMoments pcc_2d_density_moments(const std::vector<double>& density,
                                    const std::vector<ModuleBase::Vector3<double>>& positions,
                                    double volume_element,
                                    const Pcc2dGeometry& geometry);

Pcc2dMoments pcc_2d_density_moments_from_relative_y(
    const std::vector<double>& density,
    const std::vector<double>& relative_y,
    double volume_element);

double pcc_2d_potential(const Pcc2dMoments& moments,
                        double relative_y,
                        const Pcc2dParameters& parameters);

ModuleBase::Vector3<double> pcc_2d_potential_gradient(
    const Pcc2dMoments& moments,
    double relative_y,
    const Pcc2dParameters& parameters);

ModuleBase::Vector3<double> pcc_2d_point_charge_force(
    const Pcc2dMoments& total_moments,
    const PointCharge& point,
    const Pcc2dGeometry& geometry);

double pcc_2d_bilinear_energy(const Pcc2dMoments& left,
                              const Pcc2dMoments& right,
                              const Pcc2dParameters& parameters);

double pcc_2d_self_energy(const Pcc2dMoments& moments,
                          const Pcc2dParameters& parameters);

double pcc_2d_ionic_shape_energy(double polarization_charge,
                                 const Pcc2dMoments& smooth_ionic_moments,
                                 const Pcc2dMoments& point_ionic_moments,
                                 const Pcc2dParameters& parameters);

} // namespace ModuleSccs

#endif
