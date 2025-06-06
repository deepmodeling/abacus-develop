#ifndef CONV_COULOMB_POT_K_H
#define CONV_COULOMB_POT_K_H

#include "module_cell/klist.h"

#include <map>
#include <string>
#include <vector>

namespace Conv_Coulomb_Pot_K
{
enum class Ccp_Type
{   //	parameter:
    // bare Coulomb
    // 1/r
    Ccp, // "hybrid_alpha"
    // Spencer-Alavi truncated bare Coulomb
    Hf, //	"hf_Rcut", "hybrid_alpha"
    // short-range Coulomb
    // erfc(omega * r) / r
    Erfc, //	"hse_omega", "hybrid_beta"
    // long-range Coulomb
    // erf(omega * r) / r
	Erf, //"hse_omega"
    // Coulomb-attenuating method for truncated range-separated Coulomb
    Cam,  //  "hse_omega", "hybrid_alpha", "hybrid_beta", "hf_Rcut"
    // Coulomb-attenuating method for range-separated Coulomb
    // (alpha + beta * erfc(omega * r)) / r
    Ccp_Cam // "hse_omega", "hybrid_alpha", "hybrid_beta"
};

template <typename T>
T cal_orbs_ccp(const T& orbs,
               const Ccp_Type& ccp_type,
               const std::map<std::string, double>& parameter,
               const double rmesh_times);

// private:
template <typename T>
double get_rmesh_proportion(const T& orbs, const double psi_threshold);

// private:
std::vector<double> cal_psi_ccp(const std::vector<double>& psif, const double hybrid_alpha);
std::vector<double> cal_psi_hf(const std::vector<double>& psif,
                               const std::vector<double>& k_radial,
                               const int Rcut_type,
                               const double hybrid_alpha,
                               const double Rc);
std::vector<double> cal_psi_erfc(const std::vector<double>& psif,
                                 const std::vector<double>& k_radial,
                                 const double hse_omega,
                                 const double hybrid_beta);
std::vector<double> cal_psi_cam(const std::vector<double>& psif,
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
} // namespace Conv_Coulomb_Pot_K

#include "conv_coulomb_pot_k.hpp"

#endif