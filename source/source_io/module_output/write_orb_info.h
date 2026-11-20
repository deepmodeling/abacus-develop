#ifndef WRITE_ORB_INFO
#define WRITE_ORB_INFO

#include <iosfwd>

class UnitCell;

namespace ModuleIO
{
    void write_orb_info(const UnitCell* ucell, std::ofstream& ofs_running);
}

#endif
