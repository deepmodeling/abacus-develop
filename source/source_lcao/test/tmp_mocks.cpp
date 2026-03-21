
#include "source_cell/unitcell.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"

// constructor of Atom
Atom::Atom() {}
Atom::~Atom() {}

Atom_pseudo::Atom_pseudo() {}
Atom_pseudo::~Atom_pseudo() {}

Magnetism::Magnetism() {}
Magnetism::~Magnetism() {}

#ifdef __LCAO
InfoNonlocal::InfoNonlocal() {}
InfoNonlocal::~InfoNonlocal() {}
LCAO_Orbitals::LCAO_Orbitals() {}
LCAO_Orbitals::~LCAO_Orbitals() {}
#endif

pseudo::pseudo() {}
pseudo::~pseudo() {}
SepPot::SepPot(){}
SepPot::~SepPot(){}
Sep_Cell::Sep_Cell() noexcept {}
Sep_Cell::~Sep_Cell() noexcept {}

// constructor of UnitCell
UnitCell::UnitCell() {}
UnitCell::~UnitCell() {}

void UnitCell::set_iat2iwt(const int& npol_in)
{
    this->iat2iwt.resize(this->nat);
    this->npol = npol_in;
    int iat = 0;
    int iwt = 0;
    for (int it = 0; it < this->ntype; it++)
    {
        for (int ia = 0; ia < atoms[it].na; ia++)
        {
            this->iat2iwt[iat] = iwt;
            iwt += atoms[it].nw * this->npol;
            ++iat;
        }
    }
    return;
}

// stub for Grid_Driver::Find_atom (used by density_matrix_io.cpp but not exercised in test)
void Grid_Driver::Find_atom(const UnitCell& ucell,
                            const ModuleBase::Vector3<double>& tau,
                            const int& T,
                            const int& I,
                            AdjacentAtomInfo* adjs) const
{
}
