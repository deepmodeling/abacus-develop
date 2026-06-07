#ifndef EXX_INFO_H
#define EXX_INFO_H

#include "exx_info_global.h"
#include "exx_info_lip.h"
#include "exx_info_ri.h"
#include "exx_info_opt_abfs.h"

struct Exx_Info
{
    // WARNING: Lifetime dependency
    // Exx_Info_Lip and Exx_Info_RI hold const references to members of Exx_Info_Global.
    // Their lifetimes MUST NOT exceed the Exx_Info_Global they were initialized from.
    // Currently this is safe because Exx_Info owns both info_global and info_lip/info_ri,
    // and GlobalC::exx_info is a global variable with program-lifetime scope.
    // However, if Exx_Info is ever copied or moved, these references will become dangling.
    // Do NOT add copy/move constructors to Exx_Info without addressing this.

    Exx_Info_Global info_global;
    Exx_Info_Lip info_lip;
    Exx_Info_RI info_ri;
    Exx_Info_Opt_ABFs info_opt_abfs;

    Exx_Info() : info_lip(this->info_global), info_ri(this->info_global)
    {
    }
};

namespace GlobalC
{
    extern Exx_Info exx_info;
}

#endif