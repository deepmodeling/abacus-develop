#ifndef READ_ORB_H
#define READ_ORB_H

#include "module_cell/unitcell.h"

namespace elecstate 
{

    /**
     * @brief read number of numerical orbitals for each angular momentum
     * @param it index of atom type
     * @param orb_file orbital filename
     * @param ofs_running ofstream
     * @param atom Atom instance stored in UnitCell
    */
    void read_orb_file(int it,
                       std::string& orb_file,
                       std::ofstream& ofs_running,
                       Atom* atom);

}

#endif