#ifndef EXX_INFO_GLOBAL_H
#define EXX_INFO_GLOBAL_H

#include "source_lcao/module_ri/conv_coulomb_pot_k.h"

#include <vector>
#include <map>
#include <string>
#include <cstddef>

struct Exx_Info_Global
{
    bool cal_exx = false;

    std::map<Conv_Coulomb_Pot_K::Coulomb_Type, std::vector<std::map<std::string, std::string>>> coulomb_param;

    // Fock:
    //      "alpha":        "0"
    //      "singularity_correction":   "limits" / "spencer" / "revised_spencer" / "massidda" / "carrier"
    //      "lambda":       "0.3"
    //      "Rcut"
    // Erfc:
    //      "alpha":        "0"
    //      "omega":        "0.11"
    //      "singularity_correction":   "limits" / "spencer" / "revised_spencer"
    //      "Rcut"

    Conv_Coulomb_Pot_K::Ccp_Type ccp_type;
    double hybrid_alpha = 0.25;
    double hse_omega = 0.11;
    double mixing_beta_for_loop1 = 1.0;

    bool separate_loop = true;
    size_t hybrid_step = 1;
};

#endif