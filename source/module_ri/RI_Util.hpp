//=======================
// AUTHOR : Peize Lin
// DATE :   2022-08-17
//=======================

#ifndef RI_UTIL_HPP
#define RI_UTIL_HPP

#include "RI_Util.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

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
}

#endif