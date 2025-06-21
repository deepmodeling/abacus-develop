#ifndef CONV_COULOMB_POT_K_HPP
#define CONV_COULOMB_POT_K_HPP

#include "conv_coulomb_pot_k.h"
#include <cmath>
#include <cassert>

namespace Conv_Coulomb_Pot_K
{

	template< typename T >
	std::vector<T> cal_orbs_ccp(
		const std::vector<T> & orbs,
		const std::map<Conv_Coulomb_Pot_K::Coulomb_Type, std::vector<std::map<std::string,std::string>>> &coulomb_param,
		const double rmesh_times)
	{
		std::vector<T> orbs_ccp(orbs.size());
		for( size_t i=0; i!=orbs.size(); ++i )
			orbs_ccp[i] = cal_orbs_ccp(orbs[i], coulomb_param, rmesh_times);
		return orbs_ccp;
	}

	template< typename T >
	double get_rmesh_proportion(
		const std::vector<T> & orbs,
		const double psi_threshold)
	{
		double rmesh_proportion=0;
		for( const auto &orb : orbs )
			rmesh_proportion = std::max(rmesh_proportion, get_rmesh_proportion(orb,psi_threshold));
		return rmesh_proportion;
	}

	// for cal_orbs_ccp()
	template<typename T>
	std::vector<T> operator*(const T &s, const std::vector<T> &v_in)
	{
		std::vector<T> v(v_in.size());
		for(std::size_t i=0; i<v.size(); ++i)
			{ v[i] = s * v_in[i]; }
		return v;
	}
	template<typename T>
	std::vector<T> operator+ (const std::vector<T> &v1, const std::vector<T> &v2)
	{
		assert(v1.size()==v2.size());
		std::vector<T> v(v1.size());
		for(std::size_t i=0; i<v.size(); ++i)
			{ v[i] = v1[i] + v2[i]; }
		return v;
	}
}

#endif