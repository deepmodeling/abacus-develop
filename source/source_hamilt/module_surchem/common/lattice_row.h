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
        return ModuleBase::Vector3<double>(scale * lattice_vectors.e11,
                                           scale * lattice_vectors.e12,
                                           scale * lattice_vectors.e13);
    }
    if (row == 1)
    {
        return ModuleBase::Vector3<double>(scale * lattice_vectors.e21,
                                           scale * lattice_vectors.e22,
                                           scale * lattice_vectors.e23);
    }
    return ModuleBase::Vector3<double>(scale * lattice_vectors.e31,
                                       scale * lattice_vectors.e32,
                                       scale * lattice_vectors.e33);
}

} // namespace ModuleSurchem

#endif
