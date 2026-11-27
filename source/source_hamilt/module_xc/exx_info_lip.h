#ifndef EXX_INFO_LIP_H
#define EXX_INFO_LIP_H

#include "source_hamilt/module_xc/coulomb_config.h"

struct Exx_Info_Lip
{
    Conv_Coulomb_Pot_K::Ccp_Type ccp_type;
    double hse_omega = 0.11;
    double lambda = 0.3;
};

#endif