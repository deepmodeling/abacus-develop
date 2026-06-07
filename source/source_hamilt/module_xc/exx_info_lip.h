#ifndef EXX_INFO_LIP_H
#define EXX_INFO_LIP_H

#include "exx_info_global.h"

struct Exx_Info_Lip
{
    const Conv_Coulomb_Pot_K::Ccp_Type& ccp_type;  // reference to info_global.ccp_type
    const double& hse_omega;                        // reference to info_global.hse_omega
    double lambda = 0.3;

    Exx_Info_Lip(const Exx_Info_Global& info_global)
        : ccp_type(info_global.ccp_type),
          hse_omega(info_global.hse_omega) {}
};

#endif