#ifndef CONV_COULOMB_POT_K_H
#define CONV_COULOMB_POT_K_H

#include <vector>
#include <map>
#include <string>

namespace Conv_Coulomb_Pot_K
{
	enum class Coulomb_Type{Fock, Erfc};
	enum class Ccp_Type{		//	parameter:
		Ccp,					//
		Hf,						//		"hf_Rcut"
		Erfc,					//		"hse_omega"
		Erf};					//		"hse_omega", "hf_Rcut"
	enum class Coulomb_Method{Center2, Ewald}; // Different methods for constructing the Coulomb matrix.

	template<typename T> extern T cal_orbs_ccp(
		const T &orbs,
		const std::map<Conv_Coulomb_Pot_K::Coulomb_Type, std::vector<std::map<std::string,std::string>>> &coulomb_param,
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