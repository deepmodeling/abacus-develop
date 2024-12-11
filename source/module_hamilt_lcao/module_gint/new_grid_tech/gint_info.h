#pragma once

#include <memory>
#include <vector>
#include "module_cell/atom_spec.h"
#include "gint_type.h"
#include "biggrid.h"
#include "gint_atom.h"
#include "unitcell_info.h"
#include "localcell_info.h"
#include "divide_info.h"

namespace Gint
{

class GintInfo
{
    public:
    // constructor
    GintInfo(
        Vec3d unitcell_vec1,
        Vec3d unitcell_vec2,
        Vec3d unitcell_vec3,
        int nbx, int nby, int nbz,
        int nmx, int nmy, int nmz,
        int startidx_bx, int startidx_by, int startidx_bz,
        int nbx_local, int nby_local, int nbz_local,
        int ntype, Atom* atoms, Numerical_Orbital* Phi);

    private:
    // initialize the atoms
    void init_atoms_(int ntype, Atom* atoms, Numerical_Orbital* Phi);

    // the unitcell information
    std::shared_ptr<const UnitCellInfo> unitcell_info_;

    // the divide information
    std::shared_ptr<const DivideInfo> divide_info_;

    // the localcell information
    std::shared_ptr<const LocalCellInfo> localcell_info_;

    // the big grids on this processor
    std::vector<std::shared_ptr<BigGrid>> biggrids_;

    // the total atoms in the unitcell(include extended unitcell)
    std::vector<std::shared_ptr<GintAtom>> atoms_;
}

} // namespace Gint
