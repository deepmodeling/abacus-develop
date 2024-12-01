#ifndef GINT_HELPER_H
#define GINT_HELPER_H

#include "module_base/vector3.h"

namespace Gint
{
    // Define the alias for the Vector3<double> and Matrix3
    using Vector3 = ModuleBase::Vector3<double>;
    using Matrix3 = ModuleBase::Matrix3;

    inline int index3Dto1D(const int id_x, const int id_y, const int id_z,
                           const int dim_x, const int dim_y, const int dim_z)
    {
        return id_z + id_y * dim_z + id_x * dim_y * dim_z;
    };

    inline ModuleBase::Vector3<int> index1Dto3D(const int index_1d,
                                                const int dim_x, const int dim_y, const int dim_z)
    {
        int id_x = index_1d / (dim_y * dim_z);
        int id_y = (index_1d - id_x * dim_y * dim_z) / dim_z;
        int id_z = index_1d % dim_z;
        return ModuleBase::Vector3<int>(id_x, id_y, id_z);
    };
}

#endif