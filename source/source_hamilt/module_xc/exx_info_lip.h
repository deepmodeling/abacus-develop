#ifndef EXX_INFO_LIP_H
#define EXX_INFO_LIP_H

#include "source_lcao/module_ri/conv_coulomb_pot_k.h"

struct Exx_Info_Lip
{
    Conv_Coulomb_Pot_K::Ccp_Type ccp_type;
    double hse_omega = 0.11;
    double lambda = 0.3;
};

#endif