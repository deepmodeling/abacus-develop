#ifndef TEST_ESOLVER_H
#define TEST_ESOLVER_H

#include "module_cell/unitcell.h"

UnitCell::UnitCell()
{
    Coordinate = "Direct";
    latName = "none";
    lat0_angstrom = 10.0;

    latvec.e11 = 1.0;
    latvec.e12 = 0.0;
    latvec.e13 = 0.0;
    latvec.e21 = 0.0;
    latvec.e22 = 1.0;
    latvec.e23 = 0.0;
    latvec.e31 = 0.0;
    latvec.e32 = 0.0;
    latvec.e33 = 1.0;

    ntype = 2;
    nat = 2;

    atom_label = new std::string[ntype];
    atom_label[0] = "Al";
    atom_label[1] = "Cu";

    atoms = new Atom[ntype];
    set_atom_flag = true;

    for (int it = 0; it < ntype; it++)
    {
        Atom* atom = &atoms[it];
        for (int ia = 0; ia < atom->na; ia++)
        {
            for (int ik = 0; ik < 3; ++ik)
            {
                atom->tau[ia][ik] = 3.0 * ia + ik;
            }
        }
    }
}
UnitCell::~UnitCell()
{
}
Magnetism::Magnetism()
{
}
Magnetism::~Magnetism()
{
}
Atom::Atom()
{
    na = 1;
    tau = new ModuleBase::Vector3<double>[na];
    mbl = new ModuleBase::Vector3<int>[na];
}
Atom::~Atom()
{
    delete[] tau;
    delete[] mbl;
}
Atom_pseudo::Atom_pseudo()
{
}
Atom_pseudo::~Atom_pseudo()
{
}
pseudo_nc::pseudo_nc()
{
}
pseudo_nc::~pseudo_nc()
{
}
#endif