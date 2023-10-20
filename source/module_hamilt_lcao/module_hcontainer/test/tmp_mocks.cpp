
#include "module_cell/unitcell.h"

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

// constructor of UnitCell
UnitCell::UnitCell()
{
}
UnitCell::~UnitCell()
{
}

void UnitCell::set_iat2iwt(const int& npol_in)
{
	this->iat2iwt.resize(this->nat);
	this->npol = npol_in;
	int iat=0;
	int iwt=0;
	for(int it = 0;it < this->ntype; it++)
	{
		for(int ia=0; ia<atoms[it].na; ia++)
		{
			this->iat2iwt[iat] = iwt;
			iwt += atoms[it].nw * this->npol;
			++iat;
		}	
	}
	return;
}