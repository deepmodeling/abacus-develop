#ifndef SURCHEM_PW_GRID_H
#define SURCHEM_PW_GRID_H

#include "source_base/matrix3.h"
#include "source_base/vector3.h"

#include <vector>

namespace ModulePW
{
class PW_Basis;
}

namespace ModuleSurchem
{

ModuleBase::Vector3<double> cell_center(const ModuleBase::Matrix3& lattice_vectors,
                                       double lattice_scale);

std::vector<ModuleBase::Vector3<double>> pw_grid_positions(
    const ModulePW::PW_Basis& basis,
    const ModuleBase::Matrix3& lattice_vectors,
    double lattice_scale);

} // namespace ModuleSurchem

#endif
