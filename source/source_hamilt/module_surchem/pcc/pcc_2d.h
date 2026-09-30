#ifndef PCC_2D_H
#define PCC_2D_H

#include "pcc_0d.h"

#include <vector>

namespace ModuleBase
{
class Matrix3;
}

namespace ModulePcc
{

// Moments along the slab normal: charge, dipole and quadrupole of the
// coordinate measured from the geometry origin.
struct Pcc2dMoments
{
    double charge = 0.0;
    double dipole = 0.0;
    double quadrupole = 0.0;
};

struct Pcc2dParameters
{
    double periodic_area = 0.0;
    // Cell length along the slab normal.
    double cell_length = 0.0;
};

// The slab is open along lattice vector `axis` (0, 1 or 2), which must be
// perpendicular to the two periodic lattice vectors; `normal` is its unit
// vector and `origin` the reference coordinate along it. The multipoles are
// taken about the origin, and coordinates wrap half a cell length from it, so
// that plane must lie in the vacuum. The defaults match pcc_2d_axis 2.
struct Pcc2dGeometry
{
    Pcc2dParameters parameters;
    int axis = 2;
    ModuleBase::Vector3<double> normal = ModuleBase::Vector3<double>(0.0, 0.0, 1.0);
    double origin = 0.0;
};

Pcc2dGeometry pcc_2d_geometry(const ModuleBase::Matrix3& lattice_vectors,
                              double lattice_scale,
                              int axis,
                              double relative_tolerance);

void validate_pcc_2d_parameters(const Pcc2dParameters& parameters);

void validate_pcc_2d_geometry(const Pcc2dGeometry& geometry);

// Coordinate of a Cartesian position along the slab normal, not wrapped.
double pcc_2d_coordinate(const ModuleBase::Vector3<double>& position,
                         const Pcc2dGeometry& geometry);

// Per-point functions (pcc_2d_relative_coordinate, pcc_2d_potential,
// pcc_2d_potential_gradient) do not validate the geometry, parameters or
// moments; callers validate them once per grid.
double pcc_2d_relative_coordinate(const ModuleBase::Vector3<double>& position,
                                  const Pcc2dGeometry& geometry);

// Weighted center of coordinates along the normal, wrapped into one period.
double pcc_2d_system_center(const std::vector<double>& coordinates,
                            const std::vector<double>& weights,
                            double cell_length);

Pcc2dMoments pcc_2d_point_charge_moments(const std::vector<PointCharge>& charges,
                                          const Pcc2dGeometry& geometry);

Pcc2dMoments pcc_2d_density_moments(const std::vector<double>& density,
                                    const std::vector<ModuleBase::Vector3<double>>& positions,
                                    double volume_element,
                                    const Pcc2dGeometry& geometry);

Pcc2dMoments pcc_2d_density_moments_from_relative_coordinates(
    const std::vector<double>& density,
    const std::vector<double>& relative_coordinates,
    double volume_element);

// Moments of two charge distributions taken about the same plane, added.
Pcc2dMoments pcc_2d_sum_moments(const Pcc2dMoments& left, const Pcc2dMoments& right);

double pcc_2d_potential(const Pcc2dMoments& moments,
                        double relative_coordinate,
                        const Pcc2dParameters& parameters);

ModuleBase::Vector3<double> pcc_2d_potential_gradient(
    const Pcc2dMoments& moments,
    double relative_coordinate,
    const Pcc2dGeometry& geometry);

ModuleBase::Vector3<double> pcc_2d_point_charge_force(
    const Pcc2dMoments& total_moments,
    const PointCharge& point,
    const Pcc2dGeometry& geometry);

double pcc_2d_bilinear_energy(const Pcc2dMoments& left,
                              const Pcc2dMoments& right,
                              const Pcc2dParameters& parameters);

double pcc_2d_self_energy(const Pcc2dMoments& moments,
                          const Pcc2dParameters& parameters);

} // namespace ModulePcc

#endif
