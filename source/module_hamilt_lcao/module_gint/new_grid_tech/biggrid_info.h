#ifndef BIGGIRD_INFO_H
#define BIGGIRD_INFO_H

#include <memory>
#include "module_base/vector3.h"
#include "module_base/matrix3.h"
#include "meshgrid_info.h"

namespace Gint
{

/**
 * @class BigGridInfo
 * @brief This class stores some basic properties common to all big grids.
 */
class BigGridInfo
{
    public:
        // constructor
        BigGridInfo(
            const Vec3d& biggrid_vec1,
            const Vec3d& biggrid_vec2,
            const Vec3d& biggrid_vec3,
            const int nmx, const int nmy, const int nmz);
        
        // getter functions
        const Vec3d &get_vec1() const { return biggrid_vec1_; };
        const Vec3d &get_vec2() const { return biggrid_vec2_; };
        const Vec3d &get_vec3() const { return biggrid_vec3_; };
        const Matrix3 &get_latvec0() const { return biggrid_latvec0_; };
        const Matrix3 &get_GT() const { return biggrid_GT_; };
        int get_nmx() const { return nmx_; };
        int get_nmy() const { return nmy_; };
        int get_nmz() const { return nmz_; };
        int get_nmxyz() const { return nmxyz_; };
        const std::vector<Vec3d>& get_meshgrid_coords() const { return meshgrid_coords_; };
        std::shared_ptr<const MeshGridInfo> get_meshgrid_info() const { return meshgrid_info_; };

        // get the 3D index of a meshgrid in the big grid from the 1D index
        Vec3i meshgrid_idx_1Dto3D(int index_1d) const
        {
            return Gint::index1Dto3D(index_1d, nmx_, nmy_, nmz_);
        };

        // get the 1D index of a meshgrid in the big grid from the 3D index
        int meshgrid_idx_3Dto1D(const Vec3i index_3d) const
        {
            return Gint::index3Dto1D(index_3d.x, index_3d.y, index_3d.z, nmx_, nmy_, nmz_);
        };
        
        // get the cartesian coordinate of a meshgrid in the big grid from the 3D index
        Vec3d get_meshgrid_coord(Vec3i index_3d) const
        {
            return index_3d * meshgrid_info_->get_latvec0();
        };

        // get the cartesian coordinate of a meshgrid in the big grid from the 1D index
        Vec3d get_meshgrid_coord(int index_1d) const
        {
            return get_meshgrid_coord(meshgrid_idx_1Dto3D(index_1d));
        };

    private:
        // basis vectors of the big grid
        Vec3d biggrid_vec1_;
        Vec3d biggrid_vec2_;
        Vec3d biggrid_vec3_;

        // used to convert the (i, j, k) index of the big grid to the Cartesian coordinate
        // if biggrid_vec1_ is row vector,
        // then biggrid_latvec0_ = [biggrid_vec1_; biggrid_vec2_; biggrid_vec3_],
        // (i, j, k) * biggrid_latvec0_ = (x, y, z)
        Matrix3 biggrid_latvec0_;

        // used to convert the Cartesian coordinate to the (i, j, k) index of the big grid
        // biggrid_GT_ = biggrid_latvec0_.Inverse()
        // (x, y, z) * biggrid_GT_ = (i, j, k)
        Matrix3 biggrid_GT_;

        //-------------------------------------------
        // some member variables related to meshgrid 
        //-------------------------------------------

        // basic attributes of meshgrid
        std::shared_ptr<MeshGridInfo> meshgrid_info_;

        // the number of meshgrids of a biggrid along the first basis vector
        // nmx may be a confusing name, because it is not the number of meshgrids along x axis
        // but it's used in the original code, so I keep it, maybe it will be changed later
        int nmx_;

        // the number of meshgrids of a biggrid along the second basis vector
        int nmy_;

        // the number of meshgrids of a biggrid along the third basis vector
        int nmz_;

        // total number of meshgrids in the biggrid
        int nmxyz_;

        // store the relative Cartesian coordinates of all meshgrids in the biggrid
        // the size of vector is nbxyz_
        std::vector<Vec3d> meshgrid_coords_;
};

}
#endif