#ifndef WRITE_ORB_INFO
#define WRITE_ORB_INFO
#include "source_cell/unitcell.h"

#include <iosfwd>

namespace ModuleIO
{
    void write_orb_info(const UnitCell* ucell, std::ofstream& ofs_running);
}

#endif
