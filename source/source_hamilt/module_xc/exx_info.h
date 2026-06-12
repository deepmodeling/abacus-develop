#ifndef EXX_INFO_H
#define EXX_INFO_H

#include "exx_info_global.h"
#include "exx_info_lip.h"
#include "exx_info_ri.h"
#include "exx_info_opt_abfs.h"

struct Exx_Info
{
    Exx_Info_Global info_global;
    Exx_Info_Lip info_lip;
    Exx_Info_RI info_ri;
    Exx_Info_Opt_ABFs info_opt_abfs;

    void sync_from_global()
    {
        info_lip.ccp_type = info_global.ccp_type;
        info_lip.hse_omega = info_global.hse_omega;
        info_ri.coulomb_param = info_global.coulomb_param;
    }
};

namespace GlobalC
{
    extern Exx_Info exx_info;
}

#endif