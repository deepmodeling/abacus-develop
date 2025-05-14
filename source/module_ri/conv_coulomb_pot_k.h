#ifndef CONV_COULOMB_POT_K_H
#define CONV_COULOMB_POT_K_H

#include <map>
#include <string>
#include <vector>

#include "module_cell/klist.h"

namespace Conv_Coulomb_Pot_K
{
	enum class Ccp_Type{		//	parameter:
		Ccp,					//
		Hf,						//		"hf_Rcut"
		Erfc,					//		"hse_omega"
		Erf,    //  	"hse_omega", "hf_Rcut"
        Cam,    //  	"hse_omega", "hybrid_alpha", "hybrid_beta", "hf_Rcut"
		Ccp_Cam // "hse_omega", "hybrid_alpha", "hybrid_beta"
    };

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
		const int Rcut_type,
		const double Rc);
	std::vector<double> cal_psi_hse(
		const std::vector<double> & psif,
		const std::vector<double> & k_radial,
		const double hse_omega);
  std::vector<double> cal_psi_cam(
                                           const std::vector<double>& psif,
                                           const std::vector<double>& k_radial,
                                           const double omega,
                                           const double hybrid_alpha,
                                           const double hybrid_beta,
										   const int Rcut_type,
                                           const double Rc);

  std::vector<double> cal_psi_ccp_cam(const std::vector<double>& psif,
                                               const std::vector<double>& k_radial,
                                               const double omega,
                                               const double hybrid_alpha,
                                               const double hybrid_beta);
}

#include "conv_coulomb_pot_k.hpp"

#endif