#ifndef MOCK_REMAKE_CELL_H
#define MOCK_REMAKE_CELL_H

#include "source_cell/unitcell_data.h"
#include <string>

namespace unitcell
{
    // Mock state tracking
    extern bool remake_cell_called;
    extern std::string remake_cell_latName;
    extern ModuleBase::Matrix3 remake_cell_latvec;

    // Mock implementation of remake_cell
    void remake_cell(Lattice& lat);

    // Helper functions for testing
    void reset_remake_cell_mock();
    bool was_remake_cell_called();
}

#endif // MOCK_REMAKE_CELL_H
