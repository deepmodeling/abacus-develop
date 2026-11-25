#ifndef SCCS_CELL_TEST_H
#define SCCS_CELL_TEST_H

#include "source_cell/unitcell.h"
#include "source_estate/module_charge/charge.h"
#include "source_hamilt/module_sccs/test/sccs_test.h"

#include <vector>

namespace SccsTest
{
// One unit-charge atom at fractional (1/4, 1/4, 1/4) with 0.6 uniform
// electrons in one spin channel, on the cell of the PW basis.
class CellTest : public PwTest
{
protected:
    void set_up_cell()
    {
        cell.lat0 = length;
        cell.tpiba = tpiba;
        cell.omega = basis.omega;
        cell.ntype = 1;
        cell.nat = 1;
        cell.atoms = &atom;
        atom.na = 1;
        atom.mass = 1.0;
        atom.ncpp.zv = 1.0;
        atom.tau = {ModuleBase::Vector3<double>(0.25, 0.25, 0.25)};
        const double density_value = 0.6 / basis.omega;
        density.assign(basis.nrxx, density_value);
        channel = density.data();
        charge.nspin = 1;
        charge.rho = &channel;
    }

    UnitCell cell;
    Atom atom;
    Charge charge;
    std::vector<double> density;
    double* channel = nullptr;
};
} // namespace SccsTest

#endif
