#ifndef CELL_GEOMETRY_H
#define CELL_GEOMETRY_H

#include "source_base/vector3.h"

#include <array>
#include <string>
#include <vector>

namespace ModuleBase
{
class Matrix3;
}

namespace unitcell
{

/// Orthogonal cell in Cartesian coordinates; lengths and origin are in Bohr.
/// Construct with make_orthogonal_cell before using relative_position.
struct OrthogonalCell
{
    std::array<ModuleBase::Vector3<double>, 3> axes;
    std::array<double, 3> lengths;
    ModuleBase::Vector3<double> origin;
};

/// Accept rectangular cells, including rotated cells. On failure, leave cell
/// unchanged and explain the invalid input in error. No model parameters are read.
bool make_orthogonal_cell(const ModuleBase::Matrix3& lattice,
                          double lattice_scale,
                          double relative_tolerance,
                          OrthogonalCell& cell,
                          std::string& error);

/// Minimum-image displacement from origin, with projections in [-L/2, L/2).
/// Requires a cell constructed above and a finite Cartesian position.
ModuleBase::Vector3<double> relative_position(const ModuleBase::Vector3<double>& position,
                                             const OrthogonalCell& cell);

/// Unwrap positions about the first position, take their positive-weight center,
/// then wrap the center about cell.origin. Intended for localized systems whose
/// extent about the first position is less than half a cell in each direction.
/// Requires a cell constructed by make_orthogonal_cell.
/// On failure, leave center unchanged and set error.
bool weighted_center(const std::vector<ModuleBase::Vector3<double>>& positions,
                      const std::vector<double>& weights,
                      const OrthogonalCell& cell,
                      ModuleBase::Vector3<double>& center,
                      std::string& error);

/// A periodic plane and its perpendicular open lattice vector, in Bohr.
/// The two in-plane lattice vectors need not be orthogonal.
struct SlabCell
{
    ModuleBase::Vector3<double> normal;
    double length = 0.0;
    double area = 0.0;
    double origin = 0.0;
};

bool make_slab_cell(const ModuleBase::Matrix3& lattice,
                    double lattice_scale,
                    int open_axis,
                    double relative_tolerance,
                    SlabCell& cell,
                    std::string& error);

/// Normal displacement in [-L/2, L/2); requires a validated cell and position.
double relative_coordinate(const ModuleBase::Vector3<double>& position, const SlabCell& cell);

/// Mass/weight center along the normal, unwrapped about the first position.
/// Requires a localized slab of extent less than half the normal period.
/// On failure, leave center unchanged and set error.
bool weighted_center(const std::vector<ModuleBase::Vector3<double>>& positions,
                      const std::vector<double>& weights,
                      const SlabCell& cell,
                      double& center,
                      std::string& error);

} // namespace unitcell

#endif
