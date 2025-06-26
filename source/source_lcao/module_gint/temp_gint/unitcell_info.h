#pragma once

#include <memory>
#include <cmath>
#include "biggrid_info.h"
#include "gint_helper.h"
#include "gint_type.h"

namespace ModuleGint
{

class UnitCellInfo
{
    public:
        // constructor
        UnitCellInfo(
            const Vec3d& unitcell_vec1,
            const Vec3d& unitcell_vec2,
            const Vec3d& unitcell_vec3,
            int nbx, int nby, int nbz,
            int nmx, int nmy, int nmz);
        
        // getter functions
        int get_nbx() const { return nbx_; };
        int get_nby() const { return nby_; };
        int get_nbz() const { return nbz_; };
        int get_bgrids_num() const { return nbxyz_; };
        int get_nmx() const { return nmx_; };
        int get_nmy() const { return nmy_; };
        int get_nmz() const { return nmz_; };
        int get_mgrids_num() const { return nmxyz_; };
        std::shared_ptr<const BigGridInfo> get_bgrid_info() const { return biggrid_info_; };
        std::shared_ptr<const MeshGridInfo> get_mgrid_info() const { return meshgrid_info_; };

        //====================================================================
        // functions related to the big grid
        //====================================================================

        // transform the 1D index of a big grid in the unit cell to the 3D index
        Vec3i bgrid_idx_1Dto3D(const int index_1d) const
        {
            return index1Dto3D(index_1d, nbx_, nby_, nbz_);
        };

        // transform the 3D index of a biggrid in the unit cell to the 1D index
        int bgrid_idx_3Dto1D(const Vec3i index_3d) const
        {
            return index3Dto1D(index_3d.x, index_3d.y, index_3d.z, nbx_, nby_, nbz_);
        };

        // get the cartesian coordinate of a big grid in the unit cell from the 3D index
        Vec3d get_bgrid_coord(Vec3i index_3d) const
        {
            return biggrid_info_->get_cartesian_coord(index_3d);
        };

        // get the cartesian coordinate of a big grid in the unit cell from the 1D index
        Vec3d get_bgrid_coord(int index_1d) const
        {
            return get_bgrid_coord(bgrid_idx_1Dto3D(index_1d));
        };

        // get the 3D index of a big grid in the unit cell from the cartesian coordinate
        Vec3i get_bgrid_idx_3d(const Vec3d coord) const
        {   
            Vec3d direct_coord = biggrid_info_->get_direct_coord(coord);
            return Vec3i(
                static_cast<int>(floor(direct_coord.x)),
                static_cast<int>(floor(direct_coord.y)),
                static_cast<int>(floor(direct_coord.z)));
        };

        // Get the relative Cartesian coordinates of big grid A relative to big grid B
        // returned vector = coordinates of point A - coordinates of point B
        // this function is more efficient than obtaining two 3D coordinates separately 
        // through two 3D indices and then subtracting them
        Vec3d get_relative_coord(Vec3i index_3d_a, Vec3i index_3d_b) const
        {
            return get_bgrid_coord(index_3d_a - index_3d_b);
        };

        // get the extended unitcell index of a big grid
        Vec3i get_unitcell_idx(const Vec3i index_3d) const
        {
            return Vec3i(floor_div(index_3d.x, nbx_),
                         floor_div(index_3d.y, nby_),
                         floor_div(index_3d.z, nbz_));
        };

        // map the extended big grid index to the big grid index in unitcell
        Vec3i map_ext_idx_to_ucell(const Vec3i index_3d) const
        {
            return Vec3i(index_3d.x - floor_div(index_3d.x, nbx_) * nbx_,
                         index_3d.y - floor_div(index_3d.y, nby_) * nby_,
                         index_3d.z - floor_div(index_3d.z, nbz_) * nbz_);
        };


        //====================================================================
        // functions related to the meshgrid
        //====================================================================

        // transform the 1D index of a meshgrid in the unit cell to the 3D index
        Vec3i mgrid_idx_1Dto3D(const int index_1d) const
        {
            return index1Dto3D(index_1d, nmx_, nmy_, nmz_);
        }

        // transform the 3D index of a meshgrid in the unit cell to the 1D index
        int mgrid_idx_3Dto1D(const Vec3i index_3d) const
        {
            return index3Dto1D(index_3d.x, index_3d.y, index_3d.z, nmx_, nmy_, nmz_);
        }

        // get the cartesian coordinate of a meshgrid in the unit cell from the 3D index
        Vec3d get_mgrid_coord(Vec3i index_3d) const
        {
            return meshgrid_info_->get_cartesian_coord(index_3d);
        };

        // get the cartesian coordinate of a meshgrid in the unit cell from the 1D index
        Vec3d get_mgrid_coord(int index_1d) const
        {
            return get_mgrid_coord(mgrid_idx_1Dto3D(index_1d));
        }
        
    private:
        // basis vectors of the unit cell
        Vec3d unitcell_vec1_;
        Vec3d unitcell_vec2_;
        Vec3d unitcell_vec3_;

        //====================================================================
        // member variables related to the Big Grid
        //====================================================================

        // the number of big cells along the first lattice vector
        int nbx_;

        // the number of big cells along the second lattice vector
        int nby_;

        // the number of big cells along the third lattice vector
        int nbz_;

        // the total number of big cells
        int nbxyz_;

        // basic attributes of the big grid
        std::shared_ptr<const BigGridInfo> biggrid_info_;

        std::shared_ptr<const MeshGridInfo> meshgrid_info_;

        //====================================================================
        // member variables related to meshgrid
        //====================================================================

        // the number of meshgrids along the first lattice vector
        int nmx_;

        // the number of meshgrids along the second lattice vector
        int nmy_;

        // the number of meshgrids along the third lattice vector
        int nmz_;

        // the total number of meshgrids in the unitcell
        int nmxyz_;

};

} // namespace ModuleGint