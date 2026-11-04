//=======================
// AUTHOR : Huanjing Gong
// DATE :   2026-03-17
//=======================


#include "exx_lri_detail.h"
#include "smooth_cut_coulomb.h"
#include <iostream>
#include <sstream>

#include "ri_util.h"
#include "source_cell/klist.h"
#include "source_cell/unitcell.h"
#include "source_base/global_function.h"

namespace ExxLriDetail
{

double default_spencer_rcut(const UnitCell& ucell, const K_Vectors& kv)
{
    return std::pow(0.75 * kv.get_nkstot_nospin() * ucell.omega / (ModuleBase::PI), 1.0 / 3.0);
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
            const auto cut_radius_policy = param.find("rpa_cut_radius");
            if (cut_radius_policy != param.end() && cut_radius_policy->second == "inscribed")
            {
                const std::array<std::array<double, 3>, 3> cut_bvk_lattice = {{
                    {{ucell.lat0 * kv.nmp[0] * ucell.a1.x, ucell.lat0 * kv.nmp[0] * ucell.a1.y, ucell.lat0 * kv.nmp[0] * ucell.a1.z}},
                    {{ucell.lat0 * kv.nmp[1] * ucell.a2.x, ucell.lat0 * kv.nmp[1] * ucell.a2.y, ucell.lat0 * kv.nmp[1] * ucell.a2.z}},
                    {{ucell.lat0 * kv.nmp[2] * ucell.a3.x, ucell.lat0 * kv.nmp[2] * ucell.a3.y, ucell.lat0 * kv.nmp[2] * ucell.a3.z}}
                }};
                std::ostringstream cut_radius_text;
                cut_radius_text.precision(17);
                cut_radius_text << Conv_Coulomb_Pot_K::bvk_inscribed_radius(cut_bvk_lattice);
                param["Rcut"] = cut_radius_text.str();
            }
            auto rcut_it = param.find("Rcut");
            if (rcut_it == param.end() || rcut_it->second.empty())
            {
                param["Rcut"] = ModuleBase::GlobalFunc::TO_STRING(fallback_rcut);
                used_fallback_rcut = true;
            }
            const auto cut_mode = param.find("rpa_cut_coulomb");
            if (cut_mode != param.end())
            {
                std::cout << "RPA cut Coulomb: mode=" << cut_mode->second
                          << " radius=" << param.at("rpa_cut_radius")
                          << " Rcut_bohr=" << param.at("Rcut")
                          << " width=" << param.at("rpa_cut_width") << std::endl;
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
