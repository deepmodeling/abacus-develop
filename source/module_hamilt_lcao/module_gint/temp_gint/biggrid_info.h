#pragma once

#include <memory>
#include "gint_type.h"
#include "gint_helper.h"
#include "meshgrid_info.h"

namespace ModuleGint
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
            Vec3d biggrid_vec1,
            Vec3d biggrid_vec2,
            Vec3d biggrid_vec3,
            int nmx, int nmy, int nmz);
        
        Vec3d get_cartesian_coord(const Vec3d& index_3d) const { return index_3d * biggrid_latvec0_; };
        Vec3d get_cartesian_coord(const Vec3i& index_3d) const { return index_3d * biggrid_latvec0_; };
        const Vec3d get_direct_coord(const Vec3d& cart_coord) const { return cart_coord * biggrid_GT_; };

        // Return the maximum number of big grids that can fit inside a sphere of radius r,
        // along the three lattice vector directions.
        Vec3i max_ext_bgrid_num(double r) const;

        // get number of meshgrids along three lattice directions
        int get_nmx() const { return nmx_; };
        int get_nmy() const { return nmy_; };
        int get_nmz() const { return nmz_; };
        int get_mgrids_num() const { return nmxyz_; };

        const std::vector<Vec3d>& get_mgrids_coord() const { return meshgrid_coords_; };
        const Vec3d& get_mgrid_coord(int index_1d) const { return meshgrid_coords_[index_1d]; };

        std::shared_ptr<const MeshGridInfo> get_mgrid_info() const { return meshgrid_info_; };

        // get the 3D index of a meshgrid in the big grid from the 1D index
        Vec3i mgrid_idx_1Dto3D(int index_1d) const
        {
            return index1Dto3D(index_1d, nmx_, nmy_, nmz_);
        };

        // get the 1D index of a meshgrid in the big grid from the 3D index
        int mgrid_idx_3Dto1D(const Vec3i index_3d) const
        {
            return index3Dto1D(index_3d.x, index_3d.y, index_3d.z, nmx_, nmy_, nmz_);
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

        //======================================================
        // some member variables related to meshgrid 
        //======================================================

        // basic attributes of meshgrid
        std::shared_ptr<const MeshGridInfo> meshgrid_info_;

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

} // namespace ModuleGint