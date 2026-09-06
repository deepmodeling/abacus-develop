#ifndef CONV_COULOMB_POT_K_H
#define CONV_COULOMB_POT_K_H

#include "source_hamilt/module_xc/coulomb_config.h"

namespace Conv_Coulomb_Pot_K
{
	template<typename T> extern T cal_orbs_ccp(
		const T &orbs,
		const CoulombParam &coulomb_param,
		const double rmesh_times);

	template<typename T> extern T cal_orbs_ccp_spencer(
		const T &orbs,
		const CoulombParam &coulomb_param,
		const double rmesh_times);

  //private:
	template< typename T > extern double get_rmesh_proportion(
		const T &orbs,
		const double psi_threshold);

  //private:
	extern std::vector<double> cal_psi_fock_limits(
		const std::vector<double> & psif);
	extern std::vector<double> cal_psi_fock_spencer(
		const std::vector<double> &psif,
		const std::vector<double> &k_radial,
		const double rcut);
	extern std::vector<double> cal_psi_erfc_limits(
		const std::vector<double> & psif,
		const std::vector<double> & k_radial,
		const double erfc_omega);
	extern std::vector<double> cal_psi_erfc_spencer(
		const std::vector<double> & psif,
		const std::vector<double> & k_radial,
		const double erfc_omega,
		const double rcut);
}

#include "conv_coulomb_pot_k.hpp"

#endif