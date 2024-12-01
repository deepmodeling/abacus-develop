#include "unitcell_info.h"
#include "gint_helper.h"

namespace Gint
{
    UnitCellInfo::UnitCellInfo(
        const Vector3& unitcell_vec1,
        const Vector3& unitcell_vec2,
        const Vector3& unitcell_vec3,
        const int nbx, const int nby, const int nbz,
        const int nmx, const int nmy, const int nmz)
        : unitcell_vec1_(unitcell_vec1),
            unitcell_vec2_(unitcell_vec2),
            unitcell_vec3_(unitcell_vec3),
            nbx_(nbx), nby_(nby), nbz_(nbz), nbxyz_(nbx*nby*nbz),
            nmx_(nmx), nmy_(nmy), nmz_(nmz), nmxyz_(nmx*nmy*nmz)
        {
            // initialize the biggrid_info_
            biggrid_info_ = std::make_shared<BigGridInfo>(
                unitcell_vec1_ / static_cast<double>(nbx),
                unitcell_vec2_ / static_cast<double>(nby),
                unitcell_vec3_ / static_cast<double>(nbz),
                nmx/nbx, nmy/nby, nmz/nbz);
        }
    
    //----------------------------------
    // functions related to the big grid
    //----------------------------------

    ModuleBase::Vector3<int> UnitCellInfo::biggrid_idx_1Dto3D(const int index_1d) const
    {
        return Gint::index1Dto3D(index_1d, nbx_, nby_, nbz_);
    }

    int UnitCellInfo::biggrid_idx_3Dto1D(const ModuleBase::Vector3<int> index_3d) const
    {
        return Gint::index3Dto1D(index_3d.x, index_3d.y, index_3d.z, nbx_, nby_, nbz_);
    }

    Vector3 UnitCellInfo::get_biggrid_coord(ModuleBase::Vector3<int> index_3d) const
    {
        return index_3d * biggrid_info_->get_latvec0();
    }

    Vector3 UnitCellInfo::get_biggrid_coord(int index_1d) const
    {
        return get_biggrid_coord(biggrid_idx_1Dto3D(index_1d));
    }

    Vector3 UnitCellInfo::get_relative_coord(int index_1d_a, int index_1d_b) const
    {
        return get_biggrid_coord(index_1d_a) - get_biggrid_coord(index_1d_b);
    }
    
    //----------------------------------
    // functions related to the meshgrid
    //----------------------------------

    ModuleBase::Vector3<int> UnitCellInfo::meshgrid_idx_1Dto3D(const int index_1d) const
    {
        return Gint::index1Dto3D(index_1d, nmx_, nmy_, nmz_);
    }

    int UnitCellInfo::meshgrid_idx_3Dto1D(const ModuleBase::Vector3<int> index_3d) const
    {
        return Gint::index3Dto1D(index_3d.x, index_3d.y, index_3d.z, nmx_, nmy_, nmz_);
    }

    Vector3 UnitCellInfo::get_meshgrid_coord(ModuleBase::Vector3<int> index_3d) const
    {
        return index_3d * biggrid_info_->get_meshgrid_info()->get_latvec0();
    }

    Vector3 UnitCellInfo::get_meshgrid_coord(int index_1d) const
    {
        return get_meshgrid_coord(meshgrid_idx_1Dto3D(index_1d));
    }

}