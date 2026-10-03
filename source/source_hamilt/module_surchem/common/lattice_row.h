#ifndef SURCHEM_LATTICE_ROW_H
#define SURCHEM_LATTICE_ROW_H

#include "source_base/matrix3.h"
#include "source_base/vector3.h"

namespace ModuleSurchem
{

// Cartesian lattice vector row (0, 1 or 2) of lattice_vectors, times scale.
inline ModuleBase::Vector3<double> lattice_row(const ModuleBase::Matrix3& lattice_vectors,
                                               const int row,
                                               const double scale)
{
    if (row == 0)
    {
        const double component_x = scale * lattice_vectors.e11;
        const double component_y = scale * lattice_vectors.e12;
        const double component_z = scale * lattice_vectors.e13;
        return ModuleBase::Vector3<double>(component_x, component_y, component_z);
    }
    if (row == 1)
    {
        const double component_x = scale * lattice_vectors.e21;
        const double component_y = scale * lattice_vectors.e22;
        const double component_z = scale * lattice_vectors.e23;
        return ModuleBase::Vector3<double>(component_x, component_y, component_z);
    }
    const double component_x = scale * lattice_vectors.e31;
    const double component_y = scale * lattice_vectors.e32;
    const double component_z = scale * lattice_vectors.e33;
    return ModuleBase::Vector3<double>(component_x, component_y, component_z);
}

} // namespace ModuleSurchem

#endif
