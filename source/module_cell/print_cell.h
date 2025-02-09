#ifndef PRINT_CELL_H
#define PRINT_CELL_H

#include "atom_spec.h"
#include "module_cell/unitcell.h"
namespace unitcell
{
    void print_tau(Atom* atoms,
                   const std::string& Coordinate,
                   const int ntype,
                   const double lat0);
}

#endif