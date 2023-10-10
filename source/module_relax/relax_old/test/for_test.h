#ifndef TEST_RELAX_OLD_H
#define TEST_RELAX_OLD_H

#include "module_cell/unitcell.h"

UnitCell::UnitCell()
{
    Coordinate = "Direct";
    latName = "none";
    lat0 = 10.0;

    latvec.e11 = 1.0;
    latvec.e12 = 0.0;
    latvec.e13 = 0.0;
    latvec.e21 = 0.0;
    latvec.e22 = 1.0;
    latvec.e23 = 0.0;
    latvec.e31 = 0.0;
    latvec.e32 = 0.0;
    latvec.e33 = 1.0;

    ntype = 1;
    nat = 2;
    namax = 0;
    nwmax = 0;

    iat2it = nullptr;
    iat2ia = nullptr;
    iwt2iat = nullptr;
    iwt2iw = nullptr;

    itia2iat.create(1, 1);
    lc = new int[3];

    latvec = ModuleBase::Matrix3();
    latvec_supercell = ModuleBase::Matrix3();
    G = ModuleBase::Matrix3();
    GT = ModuleBase::Matrix3();
    GGT = ModuleBase::Matrix3();
    invGGT = ModuleBase::Matrix3();

    tpiba = 0.0;
    tpiba2 = 0.0;
    omega = 0.0;

    atom_label = new std::string[1];
    atom_mass = nullptr;
    pseudo_fn = new std::string[1];
    pseudo_type = new std::string[1];
    orbital_fn = new std::string[1];

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
                atom->mbl[ia][ik] = 1;
            }
        }
    }
}
UnitCell::~UnitCell()
{
}
void UnitCell::update_pos_tau(const double* pos)
{
}
void UnitCell::print_tau(void) const
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
    na = 2;
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
pseudo::pseudo()
{
}
pseudo::~pseudo()
{
}

#endif