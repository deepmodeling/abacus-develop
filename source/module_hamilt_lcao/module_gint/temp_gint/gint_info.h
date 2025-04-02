#pragma once

#include <memory>
#include <vector>
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_cell/unitcell.h"
#include "module_cell/atom_spec.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "gint_type.h"
#include "big_grid.h"
#include "gint_atom.h"
#include "unitcell_info.h"
#include "localcell_info.h"
#include "divide_info.h"

namespace ModuleGint
{

class GintInfo
{
    public:
    // constructor
    GintInfo(
        int nbx, int nby, int nbz,
        int nmx, int nmy, int nmz,
        int startidx_bx, int startidx_by, int startidx_bz,
        int nbx_local, int nby_local, int nbz_local,
        const Numerical_Orbital* Phi,
        const UnitCell& ucell, Grid_Driver& gd);

    // getter functions
    std::vector<std::shared_ptr<BigGrid>> get_biggrids() const { return biggrids_; };
    double get_local_mgrid_num() const { return localcell_info_->get_mgrids_num(); };
    double get_mgrid_volume() const { return meshgrid_info_->get_volume(); };

    //=========================================
    // functions about hcontainer
    //=========================================
    template <typename T>
    std::shared_ptr<HContainer<T>> get_hr(int npol = 1) const;
    
    private:
    // initialize the atoms
    void init_atoms_(int ntype, const Atom* atoms, const Numerical_Orbital* Phi);

    // initialize the ijr_info
    void init_ijr_info_(const UnitCell& ucell, Grid_Driver& gd);

    const UnitCell* ucell_;

    // the unitcell information
    std::shared_ptr<const UnitCellInfo> unitcell_info_;

    // the biggrid information
    std::shared_ptr<const BigGridInfo> biggrid_info_;

    // the meshgrid information
    std::shared_ptr<const MeshGridInfo> meshgrid_info_;

    // the divide information
    std::shared_ptr<const DivideInfo> divide_info_;

    // the localcell information
    std::shared_ptr<const LocalCellInfo> localcell_info_;

    // the big grids on this processor
    std::vector<std::shared_ptr<BigGrid>> biggrids_;

    // the total atoms in the unitcell(include extended unitcell) on this processor
    // atoms[iat][Vec3i] is the atom with index iat in the unitcell with index Vec3i
    // Note: Since GintAtom does not implement a default constructor,
    // the map should not be accessed using [], but rather using the at function
    std::vector<std::map<Vec3i, GintAtom>> atoms_;

    // if the iat-th(global index) atom is in this processor, return true
    std::vector<bool> is_atom_in_proc_;

    // format for storing atomic pair information in hcontainer, used for initializing hcontainer
    std::vector<int> ijr_info_;
};

} // namespace ModuleGint
