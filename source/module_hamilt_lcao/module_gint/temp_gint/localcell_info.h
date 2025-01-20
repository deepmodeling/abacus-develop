#pragma once

#include <memory>
#include "gint_type.h"
#include "unitcell_info.h"

namespace ModuleGint
{

class LocalCellInfo
{
    public:
        // constructor
        LocalCellInfo(
            int startidx_x, int startidx_y, int startidx_z,
            int nbx, int nby, int nbz,
            std::shared_ptr<const UnitCellInfo> unitcell_info);

        // getter functions
        const int get_startidx_bx() const { return startidx_bx_; };
        const int get_startidx_by() const { return startidx_by_; };
        const int get_startidx_bz() const { return startidx_bz_; };
        const int get_nbx() const { return nbx_; };
        const int get_nby() const { return nby_; };
        const int get_nbz() const { return nbz_; };
        const int get_bgrids_num() const { return nbxyz_; };
        const int get_mgrids_num() const { return nmxyz_; };
        std::shared_ptr<const UnitCellInfo> get_unitcell_info() const { return unitcell_info_; };
        std::shared_ptr<const BigGridInfo> get_bgrid_info() const { return unitcell_info_->get_bgrid_info(); };

        //====================================================================
        // functions related to the big grid
        //====================================================================

        // transform the 3D index of a big grid in the local cell to the 3D index in the local cell
        int bgrid_idx_3Dto1D(const Vec3i index_3d) const;

        // transform the 1D index of a big grid in the local cell to the 1D index in the local cell
        Vec3i bgrid_idx_1Dto3D(const int index_1d) const;

        // transform the 3D index of a big grid in the local cell to the 3D index in the unit cell
        Vec3i get_bgrid_global_idx_3D(const Vec3i index_3d) const;

        // transform the 1D index of a big grid in the local cell to the 3D index in the unit cell
        Vec3i get_bgrid_global_idx_3D(const int index_1d) const;

        // transform the 1D index of a big grid in the local cell to the 1D index in the unit cell
        int get_bgrid_global_idx_1D(const int index_1d) const;

        // transform the 3D index of a big grid in the unit cell to the 3D index in the local cell
        Vec3i get_bgrid_local_idx_3D(const Vec3i index_3d) const;

        // transform the 1D index of a big grid in the unit cell to the 1D index in the local cell
        int get_bgrid_local_idx_1D(const int index_1d) const;

        // transform the 3D index of a big grid in the unit cell to the 1D index in the local cell
        int get_bgrid_local_idx_1D(const Vec3i index_3d) const;

        // get the cartesian coordinate of a big grid in the unit cell from the 1D index
        Vec3d get_bgrid_global_coord_3D(const int index_1d) const;

        // the input is the 3D index of a big grid in the unitcell
        // return true if the big grid is in the local cell
        bool is_bgrid_in_lcell(const Vec3i index_3d) const;


        //====================================================================
        // functions related to the meshgrid
        //====================================================================

        // transform the 3D index of a meshgrid in the local cell to the 3D index in the local cell
        int mgrid_idx_3Dto1D(const Vec3i index_3d) const;

        // transform the 1D index of a meshgrid in the local cell to the 1D index in the local cell
        Vec3i mgrid_idx_1Dto3D(const int index_1d) const;

        // transform the 3D index of a meshgrid in the local cell to the 3D index in the unit cell
        Vec3i get_mgrid_global_idx_3D(const Vec3i index_3d) const;

        // transform the 1D index of a meshgrid in the local cell to the 1D index in the unit cell
        int get_mgrid_global_idx_1D(const int index_1d) const;

    private:
        //====================================================================
        // information about the big grid
        //====================================================================

        // 3D index of the first big grid in the local cell within the unit cell
        int startidx_bx_;
        int startidx_by_;
        int startidx_bz_;

        // Number of big grids in the local cell along the three basis vectors of the local cell
        int nbx_;
        int nby_;
        int nbz_;

        // Total number of big grids in the local cell
        int nbxyz_;

        //====================================================================
        // information about the meshgrid
        //====================================================================

        // 3D index of the first meshgrid in the local cell within the unit cell
        int startidx_mx_;
        int startidx_my_;
        int startidx_mz_;

        // Number of meshgrids in the local cell along the three basis vectors of the local cell
        int nmx_;
        int nmy_;
        int nmz_;

        // Total number of meshgrids in the local cell
        int nmxyz_;

        // information about the Unitcell
        std::shared_ptr<const UnitCellInfo> unitcell_info_;

        // information about the big grid
        std::shared_ptr<const BigGridInfo> biggrid_info_;
        
};

} // namespace ModuleGint