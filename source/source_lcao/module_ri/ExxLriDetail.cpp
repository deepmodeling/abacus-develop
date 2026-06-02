//=======================
// AUTHOR : Huanjing Gong
// DATE :   2026-03-17
//=======================


#include "ExxLriDetail.h"

#include "RI_Util.h"
#include "source_cell/klist.h"
#include "source_cell/unitcell.h"
#include "source_base/global_function.h"

namespace ExxLriDetail
{

double default_spencer_rcut(const UnitCell& ucell, const K_Vectors& kv)
{
    return std::pow(0.75 * kv.get_nkstot_full() * ucell.omega / (ModuleBase::PI), 1.0 / 3.0);
}

CoulombParam build_center2_cut_coulomb_param(const CoulombParam& coulomb_param,
                                             const UnitCell& ucell,
                                             const K_Vectors& kv,
                                             bool* synthesized_rcut)
{
    CoulombParam center2_param = RI_Util::update_coulomb_param(coulomb_param, ucell, &kv);
    const double fallback_rcut = default_spencer_rcut(ucell, kv);
    bool used_fallback_rcut = false;

    for (auto& param_list: center2_param)
    {
        if (param_list.first != Conv_Coulomb_Pot_K::Coulomb_Type::Fock)
        {
            continue;
        }
        for (auto& param: param_list.second)
        {
            auto rcut_it = param.find("Rcut");
            if (rcut_it == param.end() || rcut_it->second.empty())
            {
                param["Rcut"] = ModuleBase::GlobalFunc::TO_STRING(fallback_rcut);
                used_fallback_rcut = true;
            }
        }
    }

    if (synthesized_rcut != nullptr)
    {
        *synthesized_rcut = used_fallback_rcut;
    }
    return center2_param;
}
}
