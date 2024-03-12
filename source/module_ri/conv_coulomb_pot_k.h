#ifndef CONV_COULOMB_POT_K_H
#define CONV_COULOMB_POT_K_H

#include <vector>
#include <map>
#include <string>

namespace Conv_Coulomb_Pot_K
{
	enum class Ccp_Type{		//	parameter:
		Ccp,					//
		Hf,						//		"hf_Rcut"
		Hse};					//		"hse_omega"

	template<typename T> T cal_orbs_ccp(
		const T &orbs,
		const Ccp_Type &ccp_type,
		const std::map<std::string,double> &parameter,
        const double rmesh_times);

  //private:
	template< typename T > double get_rmesh_proportion(
		const T &orbs,
		const double psi_threshold);

  //private:
	std::vector<double> cal_psi_ccp(
		const std::vector<double> & psif);
	std::vector<double> cal_psi_hf(
		const std::vector<double> &psif,
		const std::vector<double> &k_radial,
		const double hf_Rcut);
	std::vector<double> cal_psi_hse(
		const std::vector<double> & psif,
		const std::vector<double> & k_radial,
		const double hse_omega);
}

#include "conv_coulomb_pot_k.hpp"

#endif