#ifndef CONV_COULOMB_POT_K_HPP
#define CONV_COULOMB_POT_K_HPP

#include "conv_coulomb_pot_k.h"
#include <vector>
#include <cmath>

namespace Conv_Coulomb_Pot_K
{

	template< typename T >
	std::vector<T> cal_orbs_ccp(
		const std::vector<T> & orbs,
		const Ccp_Type &ccp_type,
		const std::map<std::string,double> &parameter,
		const double rmesh_times)
	{
		std::vector<T> orbs_ccp(orbs.size());
		for( size_t i=0; i!=orbs.size(); ++i )
			orbs_ccp[i] = cal_orbs_ccp(orbs[i], ccp_type, parameter, rmesh_times);
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

}

#endif