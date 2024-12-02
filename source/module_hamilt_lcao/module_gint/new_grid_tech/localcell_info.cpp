#include "localcell_info.h"

namespace Gint
{
    LocalCellInfo::LocalCellInfo(
        const int startind_bx, const int startind_by, const int startind_bz,
        const int nbx, const int nby, const int nbz,
        std::shared_ptr<UnitCellInfo> unitcell_info)
        : startind_bx_(startind_bx), startind_by_(startind_by), startind_bz_(startind_bz),
          nbx_(nbx), nby_(nby), nbz_(nbz), nbxyz_(nbx*nby*nbz),
          unitcell_info_(unitcell_info)
    {
        startind_mx_ = startind_bx_ * unitcell_info_->get_biggrid_info()->get_nmx();
        startind_my_ = startind_by_ * unitcell_info_->get_biggrid_info()->get_nmy();
        startind_mz_ = startind_bz_ * unitcell_info_->get_biggrid_info()->get_nmz();
        nmx_ = nbx_ * unitcell_info_->get_biggrid_info()->get_nmx();
        nmy_ = nby_ * unitcell_info_->get_biggrid_info()->get_nmy();
        nmz_ = nbz_ * unitcell_info_->get_biggrid_info()->get_nmz();
        nmxyz_ = nmx_ * nmy_ * nmz_;
    }

    //----------------------------------
    // functions related to the big grid
    //----------------------------------

    int LocalCellInfo::biggrid_idx_3Dto1D(const Vec3i index_3d) const
    {
        return Gint::index3Dto1D(index_3d.x, index_3d.y, index_3d.z, nbx_, nby_, nbz_);
    }

    Vec3i LocalCellInfo::biggrid_idx_1Dto3D(const int index_1d) const
    {
        return Gint::index1Dto3D(index_1d, nbx_, nby_, nbz_);
    }

    Vec3i LocalCellInfo::get_biggrid_ucell_idx(const Vec3i index_3d) const
    {
        return Vec3i(
            startind_bx_ + index_3d.x,
            startind_by_ + index_3d.y,
            startind_bz_ + index_3d.z);
    }

    int LocalCellInfo::get_biggrid_ucell_idx(const int index_1d) const
    {
        Vec3i ucell_idx_3d = get_biggrid_ucell_idx(biggrid_idx_1Dto3D(index_1d));
        return unitcell_info_->biggrid_idx_3Dto1D(ucell_idx_3d);
    }

    //----------------------------------
    // functions related to the meshgrid
    //----------------------------------

    int LocalCellInfo::meshgrid_idx_3Dto1D(const Vec3i index_3d) const
    {
        return Gint::index3Dto1D(index_3d.x, index_3d.y, index_3d.z, nmx_, nmy_, nmz_);
    }

    Vec3i LocalCellInfo::meshgrid_idx_1Dto3D(const int index_1d) const
    {
        return Gint::index1Dto3D(index_1d, nmx_, nmy_, nmz_);
    }

    Vec3i LocalCellInfo::get_meshgrid_ucell_idx(const Vec3i index_3d) const
    {
        return Vec3i(
            startind_mx_ + index_3d.x,
            startind_my_ + index_3d.y,
            startind_mz_ + index_3d.z);
    }

    int LocalCellInfo::get_meshgrid_ucell_idx(const int index_1d) const
    {
        Vec3i ucell_idx_3d = get_meshgrid_ucell_idx(meshgrid_idx_1Dto3D(index_1d));
        return unitcell_info_->meshgrid_idx_3Dto1D(ucell_idx_3d);
    }


}