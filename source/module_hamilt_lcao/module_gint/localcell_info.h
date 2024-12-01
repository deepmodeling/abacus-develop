#ifndef LOCALCELL_INFO_H
#define LOCALCELL_INFO_H

#include <memory>
#include "module_base/vector3.h"
#include "unitcell_info.h"

namespace Gint
{

class LocalCellInfo
{
    public:
        // constructor
        LocalCellInfo(
            const int startind_x, const int startind_y, const int startind_z,
            const int nbx, const int nby, const int nbz,
            std::shared_ptr<UnitCellInfo> unitcell_info);

        // getter functions
        const int get_startind_bx() const { return startind_bx_; };
        const int get_startind_by() const { return startind_by_; };
        const int get_startind_bz() const { return startind_bz_; };
        const int get_nbx() const { return nbx_; };
        const int get_nby() const { return nby_; };
        const int get_nbz() const { return nbz_; };
        const int get_nbxyz() const { return nbxyz_; };
        std::shared_ptr<const UnitCellInfo> get_unitcell_info() const { return unitcell_info_; };

        //----------------------------------
        // functions related to the big grid
        //----------------------------------

        // transform the 3D index of a big grid in the local cell to the 3D index in the local cell
        int biggrid_idx_3Dto1D(const ModuleBase::Vector3<int> index_3d) const;

        // transform the 1D index of a big grid in the local cell to the 1D index in the local cell
        ModuleBase::Vector3<int> biggrid_idx_1Dto3D(const int index_1d) const;

        // transform the 3D index of a big grid in the local cell to the 3D index in the unit cell
        ModuleBase::Vector3<int> get_biggrid_global_idx(const ModuleBase::Vector3<int> index_3d) const;

        // transform the 1D index of a big grid in the local cell to the 1D index in the unit cell
        int get_biggrid_global_idx(const int index_1d) const;


        //-----------------------------------
        // functions related to the meshgrid
        //-----------------------------------

        // transform the 3D index of a meshgrid in the local cell to the 3D index in the local cell
        int meshgrid_idx_3Dto1D(const ModuleBase::Vector3<int> index_3d) const;

        // transform the 1D index of a meshgrid in the local cell to the 1D index in the local cell
        ModuleBase::Vector3<int> meshgrid_idx_1Dto3D(const int index_1d) const;

        // transform the 3D index of a meshgrid in the local cell to the 3D index in the unit cell
        ModuleBase::Vector3<int> get_meshgrid_global_idx(const ModuleBase::Vector3<int> index_3d) const;

        // transform the 1D index of a meshgrid in the local cell to the 1D index in the unit cell
        int get_meshgrid_global_idx(const int index_1d) const;

    private:
        //-------------------------------
        // information about the big grid
        //-------------------------------

        // 3D index of the first big grid in the local cell within the unit cell
        int startind_bx_;
        int startind_by_;
        int startind_bz_;

        // Number of big grids in the local cell along the three basis vectors of the local cell
        int nbx_;
        int nby_;
        int nbz_;

        // Total number of big grids in the local cell
        int nbxyz_;

        //--------------------------------
        // information about the meshgrid
        //--------------------------------

        // 3D index of the first meshgrid in the local cell within the unit cell
        int startind_mx_;
        int startind_my_;
        int startind_mz_;

        // Number of meshgrids in the local cell along the three basis vectors of the local cell
        int nmx_;
        int nmy_;
        int nmz_;

        // Total number of meshgrids in the local cell
        int nmxyz_;

        //--------------------------------
        // information about the Unitcell
        //--------------------------------
        std::shared_ptr<UnitCellInfo> unitcell_info_;
        
}

}
#endif // LOCALCELL_INFO_H