#ifndef PW_GRID_GEOMETRY_H
#define PW_GRID_GEOMETRY_H

#include "source_base/vector3.h"

#include <string>
#include <vector>

namespace ModuleBase
{
class Matrix3;
}

namespace ModulePW
{
class PW_Basis;

/// Cartesian coordinates in Bohr of this rank's integer FFT grid nodes.
/// Supports an empty local slab; on failure, leave positions unchanged.
bool grid_positions(const PW_Basis& basis,
                     const ModuleBase::Matrix3& lattice,
                     double lattice_scale,
                     std::vector<ModuleBase::Vector3<double>>& positions,
                     std::string& error);
} // namespace ModulePW

#endif
