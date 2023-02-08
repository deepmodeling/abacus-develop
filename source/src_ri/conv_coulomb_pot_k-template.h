#ifndef CONV_COULOMB_POT_K_TEMPLATE_H
#define CONV_COULOMB_POT_K_TEMPLATE_H

#include "conv_coulomb_pot_k.h"
#include <vector>
#include <cmath>

#include "../src_ri/test_code/exx_abfs-construct_orbs-test.h"

template<typename T>
T Conv_Coulomb_Pot_K::cal_orbs_ccp_rmesh(
	const T & orbs,
	const Ccp_Type &ccp_type,
	const std::map<std::string,double> &parameter,
	const double psi_threshold,
	double &rmesh_times)
{
	const double rmesh_times_trial = 10;						// Peize Lin test
	const T orbs_ccp_trial = cal_orbs_ccp(orbs, ccp_type, parameter, rmesh_times_trial);
	rmesh_times = std::max(1.0, get_rmesh_proportion(orbs_ccp_trial, psi_threshold) * rmesh_times_trial);
	const T orbs_ccp = cal_orbs_ccp(orbs, ccp_type, parameter, rmesh_times);
	return orbs_ccp;
}



template< typename T >
T Conv_Coulomb_Pot_K::cal_orbs_ccp(
	const T & orbs,
	const Ccp_Type &ccp_type,
	const std::map<std::string,double> &parameter,
	const double rmesh_times)
{
	T orbs_ccp(orbs.size());
	for( size_t i=0; i!=orbs.size(); ++i )
		orbs_ccp[i] = cal_orbs_ccp( orbs[i], ccp_type, parameter, rmesh_times );
	return orbs_ccp;
}

extern template
Numerical_Orbital_Lm Conv_Coulomb_Pot_K::cal_orbs_ccp<Numerical_Orbital_Lm>(
	const Numerical_Orbital_Lm & orbs,
	const Ccp_Type &ccp_type,
	const std::map<std::string,double> &parameter,
	const double rmesh_times);

	
	
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