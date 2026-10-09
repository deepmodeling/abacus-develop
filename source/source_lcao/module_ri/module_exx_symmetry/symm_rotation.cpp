#include "symm_rotation.h"
#include "source_base/constants.h"
#include "source_io/module_parameter/parameter.h"
#include <cmath>
#include "source_base/parallel_reduce.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_base/tool_title.h"
#include "source_base/timer.h"
#include "source_base/mathzone.h"
#include "source_lcao/module_ri/ri_util.h"

namespace ModuleSymmetry
{
    void Symmetry_rotation::set_Cs_rotation(const std::vector<std::vector<int>>& abfs_l_nchi)
    {
        this->reduce_Cs_ = true;
        this->abfs_l_nchi_ = abfs_l_nchi;
        for (auto& abfs_T : abfs_l_nchi) { this->abfs_Lmax_ = std::max(this->abfs_Lmax_, static_cast<int>(abfs_T.size()) - 1); }
    }

    std::vector<TC> Symmetry_rotation::get_Rs_from_adjacent_list(const UnitCell& ucell,
                                                                 const Grid_Driver& gd,
                                                                 const Parallel_Orbitals& pv) const
    {
        // find the union set of Rs for all the atom pairs
        std::set<TC> Rs_set;
        for (int iat1 = 0;iat1 < ucell.nat;++iat1)
        {
            auto tau1 = ucell.get_tau(iat1);
            int it1 = ucell.iat2it[iat1], ia1 = ucell.iat2ia[iat1];
            AdjacentAtomInfo adjs;
            gd.Find_atom(ucell, tau1, it1, ia1, &adjs);
            for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
            {
                const int it2 = adjs.ntype[ad];
                const int ia2 = adjs.natom[ad];
                int iat2 = ucell.itia2iat(it2, ia2);
                if (pv.get_nrow_atom(iat1) && pv.get_ncol_atom(iat2))
                {
                    const ModuleBase::Vector3<int>& R_index = adjs.box[ad];
                    if (ucell.cal_dtau(iat1, iat2, R_index).norm() * ucell.lat0
                        < ucell.atoms[it1].Rcut + ucell.atoms[it2].Rcut) {
                        Rs_set.insert({ R_index.x, R_index.y, R_index.z });
}
                }
            }
        }
        // set to vector
        std::vector<TC> Rs(Rs_set.size());
        for (auto& R : Rs_set) { Rs.push_back(R);
}
        return Rs;
    }

    std::vector<TC> Symmetry_rotation::get_Rs_from_BvK(const K_Vectors& kv) const
    {
        const TC& period = RI_Util::get_Born_vonKarmen_period(kv);
        return RI_Util::get_Born_von_Karmen_cells(period);
    }

}
