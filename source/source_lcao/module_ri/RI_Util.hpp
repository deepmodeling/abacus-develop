//=======================
// AUTHOR : Peize Lin
// DATE :   2022-08-17
//=======================

#ifndef RI_UTIL_HPP
#define RI_UTIL_HPP

#include "RI_Util.h"
#include "source_pw/module_pwdft/global.h"
#include "source_base/global_function.h"

namespace RI_Util
{
	inline std::array<int,3>
	get_Born_vonKarmen_period(const K_Vectors &kv)
	{
		return std::array<int,3>{kv.nmp[0], kv.nmp[1], kv.nmp[2]};
	}

	template<typename Tcell>
	std::vector<std::array<Tcell,1>>
	get_Born_von_Karmen_cells( const std::array<Tcell,1> &Born_von_Karman_period )
	{
		using namespace RI::Array_Operator;
		std::vector<std::array<Tcell,1>> Born_von_Karman_cells;
		for( int c=0; c<Born_von_Karman_period[0]; ++c )
			Born_von_Karman_cells.emplace_back( std::array<Tcell,1>{c} % Born_von_Karman_period );
		return Born_von_Karman_cells;
	}

	template<typename Tcell, size_t Ndim>
	std::vector<std::array<Tcell,Ndim>>
	get_Born_von_Karmen_cells( const std::array<Tcell,Ndim> &Born_von_Karman_period )
	{
		using namespace RI::Array_Operator;

		std::array<Tcell,Ndim-1> sub_Born_von_Karman_period;
		for(int i=0; i<Ndim-1; ++i)
			sub_Born_von_Karman_period[i] = Born_von_Karman_period[i];

		std::vector<std::array<Tcell,Ndim>> Born_von_Karman_cells;
		for( const std::array<Tcell,Ndim-1> &sub_cell : get_Born_von_Karmen_cells(sub_Born_von_Karman_period) )
			for( Tcell c=0; c<Born_von_Karman_period.back(); ++c )
			{
				std::array<Tcell,Ndim> cell;
				for(int i=0; i<Ndim-1; ++i)
					cell[i] = sub_cell[i];
				cell.back() = (std::array<Tcell,1>{c} % std::array<Tcell,1>{Born_von_Karman_period.back()})[0];
				Born_von_Karman_cells.emplace_back(std::move(cell));
			}
		return Born_von_Karman_cells;
	}

	/* example for Ndim=3:
	template<typename Tcell, size_t Ndim>
	std::vector<std::array<Tcell,Ndim>>
	get_Born_von_Karmen_cells( const std::array<Tcell,Ndim> &Born_von_Karman_period )
	{
		using namespace Array_Operator;
		std::vector<std::array<Tcell,Ndim>> Born_von_Karman_cells;
		for( int ix=0; ix<Born_von_Karman_period[0]; ++ix )
			for( int iy=0; iy<Born_von_Karman_period[1]; ++iy )
				for( int iz=0; iz<Born_von_Karman_period[2]; ++iz )
					Born_von_Karman_cells.push_back( std::array<Tcell,Ndim>{ix,iy,iz} % Born_von_Karman_period );
		return Born_von_Karman_cells;
	}
	*/

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
					const double Rcut = std::pow(0.75 * p_kv->get_nkstot_full()/nspin0 * ucell.omega / (ModuleBase::PI), 1.0/3.0);
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
		}

		std::map<Conv_Coulomb_Pot_K::Coulomb_Method, 
            std::pair<bool, 
                std::map<Conv_Coulomb_Pot_K::Coulomb_Type, 
                    std::vector<std::map<std::string,std::string>>>>> coulomb_settings;

		if(!coulomb_param_center2.empty())
		{
			coulomb_settings[Conv_Coulomb_Pot_K::Coulomb_Method::Center2] = std::make_pair(true, coulomb_param_center2);
		}
		if (!coulomb_param_ewald.empty())
		{
			coulomb_settings[Conv_Coulomb_Pot_K::Coulomb_Method::Ewald] = std::make_pair(true, coulomb_param_ewald);
		}

		return coulomb_settings;
	}
}

#endif