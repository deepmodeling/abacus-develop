#ifndef READ_PSEUDO_H
#define READ_PSEUDO_H

#include "module_cell/unitcell.h"

namespace elecstate {

    void read_pseudo(std::ofstream& ofs, UnitCell& ucell);

    // read in pseudopotential from files for each type of atom
    void read_cell_pseudopots(const std::string& fn, std::ofstream& log, UnitCell& ucell);

    void print_unitcell_pseudo(const std::string& fn, UnitCell& ucell);

}

#endif