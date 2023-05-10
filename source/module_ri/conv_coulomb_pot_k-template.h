#ifndef CONV_COULOMB_POT_K_TEMPLATE_H
#define CONV_COULOMB_POT_K_TEMPLATE_H

#include "conv_coulomb_pot_k.h"
#include <vector>
#include <cmath>

#include "../module_ri/test_code/exx_abfs-construct_orbs-test.h"


template< typename T >
T Conv_Coulomb_Pot_K::cal_orbs_ccp(
	const T & orbs,
	const Ccp_Type &ccp_type,
	const std::map<std::string,double> &parameter,
	const double rmesh_times, 
    const int& nks)
{
	T orbs_ccp(orbs.size());
	for( size_t i=0; i!=orbs.size(); ++i )
		orbs_ccp[i] = cal_orbs_ccp(orbs[i], ccp_type, parameter, rmesh_times, nks );
	return orbs_ccp;
}

extern template
Numerical_Orbital_Lm Conv_Coulomb_Pot_K::cal_orbs_ccp<Numerical_Orbital_Lm>(
	const Numerical_Orbital_Lm & orbs,
	const Ccp_Type &ccp_type,
	const std::map<std::string,double> &parameter,
    const double rmesh_times,
    const int& nks);

	
	
template< typename T >
double Conv_Coulomb_Pot_K::get_rmesh_proportion(
	const T & orbs,
	const double psi_threshold)
{
	double rmesh_proportion=0;
	for( const auto &orb : orbs )
		rmesh_proportion = std::max(rmesh_proportion, get_rmesh_proportion(orb,psi_threshold));
	return rmesh_proportion;
}

extern template
double Conv_Coulomb_Pot_K::get_rmesh_proportion(
	const Numerical_Orbital_Lm & orbs,
	const double psi_threshold);
	
#endif