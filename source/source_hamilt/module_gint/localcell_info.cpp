#include "localcell_info.h"

namespace ModuleGint
{
    LocalCellInfo::LocalCellInfo(
        int startidx_bx, int startidx_by, int startidx_bz,
        int nbx, int nby, int nbz,
        std::shared_ptr<const UnitCellInfo> unitcell_info)
        : startidx_bx_(startidx_bx), startidx_by_(startidx_by), startidx_bz_(startidx_bz),
          nbx_(nbx), nby_(nby), nbz_(nbz), nbxyz_(nbx*nby*nbz),
          unitcell_info_(unitcell_info), biggrid_info_(unitcell_info->get_bgrid_info())
    {
        startidx_mx_ = startidx_bx_ * biggrid_info_->get_nmx();
        startidx_my_ = startidx_by_ * biggrid_info_->get_nmy();
        startidx_mz_ = startidx_bz_ * biggrid_info_->get_nmz();
        nmx_ = nbx_ * biggrid_info_->get_nmx();
        nmy_ = nby_ * biggrid_info_->get_nmy();
        nmz_ = nbz_ * biggrid_info_->get_nmz();
        nmxyz_ = nmx_ * nmy_ * nmz_;
    }

    //====================================================================
    // functions related to the big grid
    //====================================================================

    int LocalCellInfo::bgrid_idx_3Dto1D(const Vec3i index_3d) const
    {
        return index3Dto1D(index_3d.x, index_3d.y, index_3d.z, nbx_, nby_, nbz_);
    }

    Vec3i LocalCellInfo::bgrid_idx_1Dto3D(const int index_1d) const
    {
        return index1Dto3D(index_1d, nbx_, nby_, nbz_);
    }

    Vec3i LocalCellInfo::get_bgrid_global_idx_3D(const Vec3i index_3d) const
    {
        return Vec3i(
            startidx_bx_ + index_3d.x,
            startidx_by_ + index_3d.y,
            startidx_bz_ + index_3d.z);
    }

    Vec3i LocalCellInfo::get_bgrid_global_idx_3D(const int index_1d) const
    {
        return get_bgrid_global_idx_3D(bgrid_idx_1Dto3D(index_1d));
    }

    int LocalCellInfo::get_bgrid_global_idx_1D(const int index_1d) const
    {
        Vec3i ucell_idx_3d = get_bgrid_global_idx_3D(bgrid_idx_1Dto3D(index_1d));
        return unitcell_info_->bgrid_idx_3Dto1D(ucell_idx_3d);
    }


    Vec3i LocalCellInfo::get_bgrid_local_idx_3D(const Vec3i index_3d) const
    {
        int x = index_3d.x - startidx_bx_;
        int y = index_3d.y - startidx_by_;
        int z = index_3d.z - startidx_bz_;
        assert(x >= 0 && x < nbx_);
        assert(y >= 0 && y < nby_);
        assert(z >= 0 && z < nbz_);
        return Vec3i(x, y, z);
    }

    int LocalCellInfo::get_bgrid_local_idx_1D(const Vec3i index_3d) const
    {
        return bgrid_idx_3Dto1D(get_bgrid_local_idx_3D(index_3d));
    }

    int LocalCellInfo::get_bgrid_local_idx_1D(const int index_1d) const
    {
        Vec3i idx_3d = unitcell_info_->bgrid_idx_1Dto3D(index_1d);
        return bgrid_idx_3Dto1D(get_bgrid_local_idx_3D(idx_3d));
    }

    Vec3d LocalCellInfo::get_bgrid_global_coord_3D(const int index_1d) const
    {
        Vec3i ucell_idx_3d = get_bgrid_global_idx_3D(index_1d);
        return unitcell_info_->get_bgrid_coord(ucell_idx_3d);
    }

    bool LocalCellInfo::is_bgrid_in_lcell(const Vec3i index_3d) const
    {
        return (index_3d.x >= startidx_bx_ && index_3d.x < startidx_bx_ + nbx_ &&
                index_3d.y >= startidx_by_ && index_3d.y < startidx_by_ + nby_ &&
                index_3d.z >= startidx_bz_ && index_3d.z < startidx_bz_ + nbz_);
    }

    //====================================================================
    // functions related to the meshgrid
    //====================================================================

    int LocalCellInfo::mgrid_idx_3Dto1D(const Vec3i index_3d) const
    {
        return index3Dto1D(index_3d.x, index_3d.y, index_3d.z, nmx_, nmy_, nmz_);
    }

    Vec3i LocalCellInfo::mgrid_idx_1Dto3D(const int index_1d) const
    {
        return index1Dto3D(index_1d, nmx_, nmy_, nmz_);
    }

    Vec3i LocalCellInfo::get_mgrid_global_idx_3D(const Vec3i index_3d) const
    {
        return Vec3i(
            startidx_mx_ + index_3d.x,
            startidx_my_ + index_3d.y,
            startidx_mz_ + index_3d.z);
    }

    int LocalCellInfo::get_mgrid_global_idx_1D(const int index_1d) const
    {
        Vec3i ucell_idx_3d = get_mgrid_global_idx_3D(mgrid_idx_1Dto3D(index_1d));
        return unitcell_info_->mgrid_idx_3Dto1D(ucell_idx_3d);
    }

}