#ifndef MESH_GRID_INFO_H
#define MESH_GRID_INFO_H

#include "module_base/vector3.h"
#include "module_base/matrix3.h"
#include "module_cell/unitcell.h"

namespace Gint
{

class MeshGridInfo
{
    public:
        // constructor
        MeshGridInfo(
            const Vector3& meshgrid_vec1,
            const Vector3& meshgrid_vec2,
            const Vector3& meshgrid_vec3)
            : meshgrid_vec1_(meshgrid_vec1),
              meshgrid_vec2_(meshgrid_vec2),
              meshgrid_vec3_(meshgrid_vec3)
            {       
                // initialize the meshgrid_latvec0_
                meshgrid_latvec0_.e11 = meshgrid_vec1_.x;
                meshgrid_latvec0_.e12 = meshgrid_vec1_.y;
                meshgrid_latvec0_.e13 = meshgrid_vec1_.z;

                meshgrid_latvec0_.e21 = meshgrid_vec2_.x;
                meshgrid_latvec0_.e22 = meshgrid_vec2_.y;
                meshgrid_latvec0_.e23 = meshgrid_vec2_.z;

                meshgrid_latvec0_.e31 = meshgrid_vec3_.x;
                meshgrid_latvec0_.e32 = meshgrid_vec3_.y;
                meshgrid_latvec0_.e33 = meshgrid_vec3_.z;

                // initialize the GT matrix
                meshgrid_GT_ = meshgrid_latvec0_.Inverse();
            };
        
        // getter functions
        const Vector3 &get_vec1() const { return meshgrid_vec1_; };
        const Vector3 &get_vec2() const { return meshgrid_vec2_; };
        const Vector3 &get_vec3() const { return meshgrid_vec3_; };
        const Matrix3 &get_latvec0() const { return meshgrid_latvec0_; };
        const Matrix3 &get_GT() const { return meshgrid_GT_; };

    private:
        // basis vectors of meshgrid
        Vector3 meshgrid_vec1_;
        Vector3 meshgrid_vec2_;
        Vector3 meshgrid_vec3_;

        // used to convert the (i, j, k) index of the meshgrid to the Cartesian coordinate
        // if meshrid_vec1_ is row vector,
        // then meshgrid_latvec0_ = [meshgrid_vec1_; meshgrid_vec2_; meshgrid_vec3_],
        // (i, j, k) * meshgrid_latvec0_ = (x, y, z)
        Matrix3 meshgrid_latvec0_;

        // used to convert the Cartesian coordinate to the (i, j, k) index of the mesh grid
        // meshgrid_GT_ = meshgrid_latvec0_.Inverse()
        // (x, y, z) * meshgrid_GT_ = (i, j, k)
        Matrix3 meshgrid_GT_;
};

}
#endif