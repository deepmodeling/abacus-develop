#pragma once

#include "gint_type.h"

namespace ModuleGint
{

class MeshGridInfo
{
    public:
        // constructor
        MeshGridInfo(
            Vec3d meshgrid_vec1,
            Vec3d meshgrid_vec2,
            Vec3d meshgrid_vec3)
            : meshgrid_vec1_(meshgrid_vec1),
              meshgrid_vec2_(meshgrid_vec2),
              meshgrid_vec3_(meshgrid_vec3)
            {       
                // initialize the mgrid_latvec0_
                mgrid_latvec0_.e11 = meshgrid_vec1_.x;
                mgrid_latvec0_.e12 = meshgrid_vec1_.y;
                mgrid_latvec0_.e13 = meshgrid_vec1_.z;

                mgrid_latvec0_.e21 = meshgrid_vec2_.x;
                mgrid_latvec0_.e22 = meshgrid_vec2_.y;
                mgrid_latvec0_.e23 = meshgrid_vec2_.z;

                mgrid_latvec0_.e31 = meshgrid_vec3_.x;
                mgrid_latvec0_.e32 = meshgrid_vec3_.y;
                mgrid_latvec0_.e33 = meshgrid_vec3_.z;

                // initialize the GT matrix
                meshgrid_GT_ = mgrid_latvec0_.Inverse();

                mgrid_vol_ = std::abs(mgrid_latvec0_.Det());
            }

        double get_volume() const { return mgrid_vol_; }
        Vec3d get_cartesian_coord(const Vec3i& index_3d) const { return index_3d * mgrid_latvec0_; }
        Vec3d get_direct_coord(const Vec3d& cart_coord) const { return cart_coord * meshgrid_GT_; }

    private:
        // basis vectors of meshgrid
        Vec3d meshgrid_vec1_;
        Vec3d meshgrid_vec2_;
        Vec3d meshgrid_vec3_;

        // used to convert the (i, j, k) index of the meshgrid to the Cartesian coordinate
        // if meshrid_vec1_ is row vector,
        // then mgrid_latvec0_ = [meshgrid_vec1_; meshgrid_vec2_; meshgrid_vec3_],
        // (i, j, k) * mgrid_latvec0_ = (x, y, z)
        Matrix3 mgrid_latvec0_;

        // used to convert the Cartesian coordinate to the (i, j, k) index of the mesh grid
        // meshgrid_GT_ = mgrid_latvec0_.Inverse()
        // (x, y, z) * meshgrid_GT_ = (i, j, k)
        Matrix3 meshgrid_GT_;

        double mgrid_vol_;
};

} // namespace ModuleGint