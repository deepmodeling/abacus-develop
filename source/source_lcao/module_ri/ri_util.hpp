//=======================
// AUTHOR : Peize Lin
// DATE :   2022-08-17
//=======================

#ifndef RI_UTIL_HPP
#define RI_UTIL_HPP

#include "ri_util.h"
#include "source_base/global_function.h"
#include "source_io/module_parameter/parameter.h"

namespace RI_Util
{
	inline std::map<Conv_Coulomb_Pot_K::Coulomb_Type, std::vector<std::map<std::string,std::string>>>
	update_coulomb_param(
		const std::map<Conv_Coulomb_Pot_K::Coulomb_Type, std::vector<std::map<std::string,std::string>>> &coulomb_param,
		const UnitCell &ucell,
		const K_Vectors *p_kv)
	{
		std::map<Conv_Coulomb_Pot_K::Coulomb_Type, std::vector<std::map<std::string,std::string>>> coulomb_param_updated = coulomb_param;
		for(auto &param_list : coulomb_param_updated)
		{
			for(auto &param : param_list.second)
			{
				if(param.at("singularity_correction") == "spencer")
				{
					// 4/3 * pi * Rcut^3 = V_{supercell} = V_{unitcell} * Nk
					const int nspin0 = (PARAM.inp.nspin==2) ? 2 : 1;
					const double Rcut = std::pow(0.75 * p_kv->get_nkstot_nospin() * ucell.omega / (ModuleBase::PI), 1.0/3.0);
					param["Rcut"] = ModuleBase::GlobalFunc::TO_STRING(Rcut);
				}
                else if(param.at("singularity_correction") == "revised_spencer")
				{
					const double bvk_a1 = ucell.a1.norm() * p_kv->nmp[0];
                    const double bvk_a2 = ucell.a2.norm() * p_kv->nmp[1];
                    const double bvk_a3 = ucell.a3.norm() * p_kv->nmp[2];
                    const double Rcut = 0.5 * std::min({bvk_a1, bvk_a2, bvk_a3});
                    param["Rcut"] = ModuleBase::GlobalFunc::TO_STRING(Rcut);
				}
			}
		}
		return coulomb_param_updated;
	}

	inline std::map<Conv_Coulomb_Pot_K::Coulomb_Method,
            std::pair<bool,
                std::map<Conv_Coulomb_Pot_K::Coulomb_Type,
                    std::vector<std::map<std::string,std::string>>>>>
	update_coulomb_settings(
		const std::map<Conv_Coulomb_Pot_K::Coulomb_Type, std::vector<std::map<std::string,std::string>>> &coulomb_param,
		const UnitCell &ucell,
		const K_Vectors *p_kv)
	{
		const std::map<Conv_Coulomb_Pot_K::Coulomb_Type, std::vector<std::map<std::string,std::string>>>
			coulomb_param_updated = update_coulomb_param(coulomb_param, ucell, p_kv);

		// Separate the parameters into Center2 and Ewald methods
		std::map<Conv_Coulomb_Pot_K::Coulomb_Type, std::vector<std::map<std::string,std::string>>> coulomb_param_center2;
		std::map<Conv_Coulomb_Pot_K::Coulomb_Type, std::vector<std::map<std::string,std::string>>> coulomb_param_ewald;
		for(auto &param_list : coulomb_param_updated)
		{
			switch(param_list.first)
			{
				case Conv_Coulomb_Pot_K::Coulomb_Type::Fock:
				{
					for(auto &param : param_list.second)
					{
						if(param.at("singularity_correction") == "spencer" || param.at("singularity_correction") == "limits"
							|| param.at("singularity_correction") == "revised_spencer")
						{
							coulomb_param_center2[param_list.first].push_back(param);
						}
						else if (param.at("singularity_correction") == "massidda" || param.at("singularity_correction") == "carrier" )
						{
							coulomb_param_ewald[param_list.first].push_back(param);
						}
					}
					break;
				}
				case Conv_Coulomb_Pot_K::Coulomb_Type::Erfc:
				{
					coulomb_param_center2[param_list.first] = param_list.second; // Erfc is always calculated with Center2 method.
					break;
				}
				default:
				{
					throw std::invalid_argument( std::string(__FILE__) + " line " + std::to_string(__LINE__) );
				}
			}
		}

		std::map<Conv_Coulomb_Pot_K::Coulomb_Method,
            std::pair<bool,
                std::map<Conv_Coulomb_Pot_K::Coulomb_Type,
                    std::vector<std::map<std::string,std::string>>>>> coulomb_settings;

		const bool cal_center = !coulomb_param_center2.empty();
		const bool cal_ewald = !coulomb_param_ewald.empty();
		if(cal_center)
		{
			coulomb_settings[Conv_Coulomb_Pot_K::Coulomb_Method::Center2] = std::make_pair(cal_center, coulomb_param_center2);
		}
		if(cal_ewald)
		{
			coulomb_settings[Conv_Coulomb_Pot_K::Coulomb_Method::Ewald] = std::make_pair(cal_ewald, coulomb_param_ewald);
		}
		if(cal_center && cal_ewald)
		{
			coulomb_settings[Conv_Coulomb_Pot_K::Coulomb_Method::Center2].first = false; // If both methods are available, only HF for C is needed.
		}

		return coulomb_settings;
	}
}

#endif
