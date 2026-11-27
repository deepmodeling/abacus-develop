/**
 * @file unitcell.cpp
 * @brief Implementation of UnitCell class: constructor and destructor.
 *        Indexing helpers live in unitcell_index.cpp, per-species counting
 *        helpers in unitcell_stats.cpp, and setup routines in unitcell_setup.cpp.
 */
#include "unitcell.h"

UnitCell::UnitCell()
{
    itia2iat.create(1, 1);
}

UnitCell::~UnitCell()
{
    if (set_atom_flag)
    {
        delete[] atoms;
    }
}
