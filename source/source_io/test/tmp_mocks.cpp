
#include "source_cell/unitcell.h"

// constructor of Atom
Atom::Atom()
{
}
Atom::~Atom()
{
}

Atom_pseudo::Atom_pseudo()
{
}
Atom_pseudo::~Atom_pseudo()
{
}

Magnetism::Magnetism()
{
}
Magnetism::~Magnetism()
{
}

InfoNonlocal::InfoNonlocal()
{
}
InfoNonlocal::~InfoNonlocal()
{
}

pseudo::pseudo()
{
}
pseudo::~pseudo()
{
}

SepPot::SepPot()
{
}
SepPot::~SepPot()
{
}

Sep_Cell::Sep_Cell() noexcept : ntype(0), omega(0.0), tpiba2(0.0)
{
}
Sep_Cell::~Sep_Cell() noexcept = default;

// constructor of UnitCell
UnitCell::UnitCell()
{
}
UnitCell::~UnitCell()
{
}
