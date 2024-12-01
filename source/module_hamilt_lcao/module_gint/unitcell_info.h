#ifndef UNITCELL_INFO_H
#define UNITCELL_INFO_H

#include <memory>
#include "module_base/vector3.h"
#include "module_base/matrix3.h"
#include "biggrid_info.h"
#include "gint_helper.h"

namespace Gint
{

class UnitCellInfo
{
    public:
        // constructor
        UnitCellInfo(
            const Vector3& unitcell_vec1,
            const Vector3& unitcell_vec2,
            const Vector3& unitcell_vec3,
            const int nbx, const int nby, const int nbz,
            const int nmx, const int nmy, const int nmz);
        
        // getter functions
        const Vector3 &get_vec1() const { return unitcell_vec1_; };
        const Vector3 &get_vec2() const { return unitcell_vec2_; };
        const Vector3 &get_vec3() const { return unitcell_vec3_; };
        int get_nbx() const { return nbx_; };
        int get_nby() const { return nby_; };
        int get_nbz() const { return nbz_; };
        int get_nbxyz() const { return nbxyz_; };

        // get the number of meshcells along the first lattice vector of the unit cell
        int get_nmx() const { return nbx_ * biggrid_info_->get_nmx(); };

        // get the number of meshcells along the second lattice vector of the unit cell
        int get_nmy() const { return nby_ * biggrid_info_->get_nmy(); };

        // get the number of meshcells along the third lattice vector of the unit cell
        int get_nmz() const { return nbz_ * biggrid_info_->get_nmz(); };

        // get the total number of meshcells in the unit cell
        int get_nmxyz() const { return nbxyz_ * biggrid_info_->get_nmxyz(); };

        std::shared_ptr<const BigGridInfo> get_biggrid_info() const { return biggrid_info_; };
        std::shared_ptr<const MeshGridInfo> get_meshgrid_info() const { return biggrid_info_->get_meshgrid_info(); };

        //----------------------------------
        // functions related to the big grid
        //----------------------------------

        // transform the 1D index of a big grid in the unit cell to the 3D index
        ModuleBase::Vector3<int> biggrid_idx_1Dto3D(const int index_1d) const
        {
            return Gint::index1Dto3D(index_1d, nbx_, nby_, nbz_);
        };

        // transform the 3D index of a biggrid in the unit cell to the 1D index
        int biggrid_idx_3Dto1D(const ModuleBase::Vector3<int> index_3d) const
        {
            return Gint::index3Dto1D(index_3d.x, index_3d.y, index_3d.z, nbx_, nby_, nbz_);
        };

        // get the cartesian coordinate of a big grid in the unit cell from the 3D index
        Vector3 get_biggrid_coord(ModuleBase::Vector3<int> index_3d) const
        {
            return index_3d * biggrid_info_->get_latvec0();
        };

        // get the cartesian coordinate of a big grid in the unit cell from the 1D index
        Vector3 get_biggrid_coord(int index_1d) const
        {
            return get_biggrid_coord(biggrid_idx_1Dto3D(index_1d));
        };

        // Get the relative Cartesian coordinates of big grid A relative to big grid B
        // returned vector = coordinates of point A - coordinates of point B
        Vector3 get_relative_coord(int index_1d_a, int index_1d_b) const
        {
            return get_biggrid_coord(index_1d_a) - get_biggrid_coord(index_1d_b);
        };

        //----------------------------------
        // functions related to the meshgrid
        //----------------------------------

        // transform the 1D index of a meshgrid in the unit cell to the 3D index
        ModuleBase::Vector3<int> meshgrid_idx_1Dto3D(const int index_1d) const
        {
            return Gint::index1Dto3D(index_1d, nmx_, nmy_, nmz_);
        }

        // transform the 3D index of a meshgrid in the unit cell to the 1D index
        int meshgrid_idx_3Dto1D(const ModuleBase::Vector3<int> index_3d) const
        {
            return Gint::index3Dto1D(index_3d.x, index_3d.y, index_3d.z, nmx_, nmy_, nmz_);
        }

        // get the cartesian coordinate of a meshgrid in the unit cell from the 3D index
        Vector3 get_meshgrid_coord(ModuleBase::Vector3<int> index_3d) const
        {
            return index_3d * biggrid_info_->get_meshgrid_info()->get_latvec0();
        };

        // get the cartesian coordinate of a meshgrid in the unit cell from the 1D index
        Vector3 get_meshgrid_coord(int index_1d) const
        {
            return get_meshgrid_coord(meshgrid_idx_1Dto3D(index_1d));
        }
        
    private:
        // basis vectors of the unit cell
        Vector3 unitcell_vec1_;
        Vector3 unitcell_vec2_;
        Vector3 unitcell_vec3_;

        //----------------------------------------------
        // member variables related to the Big Grid
        //----------------------------------------------

        // the number of big cells along the first lattice vector
        int nbx_;

        // the number of big cells along the second lattice vector
        int nby_;

        // the number of big cells along the third lattice vector
        int nbz_;

        // the total number of big cells
        int nbxyz_;

        // basic attributes of the big grid
        std::shared_ptr<BigGridInfo> biggrid_info_;

        //-------------------------------------------
        // member variables related to meshgrid
        //-------------------------------------------

        // the number of meshgrids along the first lattice vector
        int nmx_;

        // the number of meshgrids along the second lattice vector
        int nmy_;

        // the number of meshgrids along the third lattice vector
        int nmz_;

        // the total number of meshgrids in the unitcell
        int nmxyz_;

}

}

#endif