#ifndef EXX_INFO_GLOBAL_H
#define EXX_INFO_GLOBAL_H

#include "general_exx_info.h"

/// Backward-compatible alias: Exx_Info_Global inherits all members
/// from General_Exx_Info without adding anything.
struct Exx_Info_Global : public General_Exx_Info
{
};

#endif