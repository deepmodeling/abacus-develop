#ifndef GENERAL_EXX_INFO_H
#define GENERAL_EXX_INFO_H

#include "coulomb_config.h"

#include <cstddef>

/// General EXX configuration, independent of basis type.
/// Contains the fields needed by both PW and LCAO EXX calculations.
struct General_Exx_Info
{
    bool cal_exx = false;

    CoulombParam coulomb_param;

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

/// Forward declaration for Input_para (full definition in input_parameter.h)
struct Input_para;

/// Initialize General_Exx_Info from input parameters.
/// Returns true if opt_orb mode is requested (generate_opt_orb).
bool init_general_exx_info(General_Exx_Info& info, const Input_para& inp);

#endif // GENERAL_EXX_INFO_H
