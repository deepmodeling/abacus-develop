#include "mock_remake_cell.h"
#include <cmath>

namespace unitcell
{
    // Mock state variables
    bool remake_cell_called = false;
    std::string remake_cell_latName = "";
    ModuleBase::Matrix3 remake_cell_latvec;

    void remake_cell(Lattice& lat)
    {
        remake_cell_called = true;
        remake_cell_latName = lat.latName;
        remake_cell_latvec = lat.latvec;

        // Mock implementation: enforce simple cubic structure for "sc"
        if (lat.latName == "sc")
        {
            double celldm = std::sqrt(lat.latvec.e11 * lat.latvec.e11 +
                                     lat.latvec.e12 * lat.latvec.e12 +
                                     lat.latvec.e13 * lat.latvec.e13);
            lat.latvec.Zero();
            lat.latvec.e11 = celldm;
            lat.latvec.e22 = celldm;
            lat.latvec.e33 = celldm;
        }
        // Mock implementation: enforce FCC structure for "fcc"
        else if (lat.latName == "fcc")
        {
            double celldm = std::sqrt(lat.latvec.e11 * lat.latvec.e11 +
                                     lat.latvec.e12 * lat.latvec.e12 +
                                     lat.latvec.e13 * lat.latvec.e13) / std::sqrt(2.0);
            lat.latvec.e11 = -celldm;
            lat.latvec.e12 = 0.0;
            lat.latvec.e13 = celldm;
            lat.latvec.e21 = 0.0;
            lat.latvec.e22 = celldm;
            lat.latvec.e23 = celldm;
            lat.latvec.e31 = -celldm;
            lat.latvec.e32 = celldm;
            lat.latvec.e33 = 0.0;
        }
        // For other lattice types, do nothing (just track the call)
    }

    void reset_remake_cell_mock()
    {
        remake_cell_called = false;
        remake_cell_latName = "";
        remake_cell_latvec.Zero();
    }

    bool was_remake_cell_called()
    {
        return remake_cell_called;
    }
}
