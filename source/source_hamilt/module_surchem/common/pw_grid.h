#ifndef SURCHEM_PW_GRID_H
#define SURCHEM_PW_GRID_H

#include "source_base/vector3.h"

#include <vector>

namespace ModuleBase
{
class Matrix3;
}

namespace ModulePW
{
class PW_Basis;
}

namespace ModuleSurchem
{

std::vector<ModuleBase::Vector3<double>> pw_grid_positions(
    const ModulePW::PW_Basis& basis,
    const ModuleBase::Matrix3& lattice_vectors,
    double lattice_scale);

} // namespace ModuleSurchem

#endif
