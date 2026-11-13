
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
SepPot::SepPot(){}
SepPot::~SepPot(){}
Sep_Cell::Sep_Cell() noexcept {}
Sep_Cell::~Sep_Cell() noexcept {}

void Atom::set_index()
{
    // Minimal mock: every orbital has (iw, l, m) = identity index.  The
    // angular-momentum kernel only inspects iw2l (l of each orbital block),
    // which assemble_angular_momentum computes from l_nchi directly, so the
    // other maps need to exist but can be zero-filled.
    iw2l.assign(nw, 0);
    iw2m.assign(nw, 0);
    iw2n.assign(nw, 0);
    iw2_ylm.assign(nw, 0);
    iw2_new.assign(nw, false);
}

void UnitCell::set_iat2itia()
{
    // Minimal mock: populate the forward (it, ia) -> iat lookup together with
    // the reverse iat -> it / iat -> ia vectors.  Used by itiaiw2iwt inside
    // assemble_angular_momentum and by get_tau / cal_dtau.
    iat2it.resize(nat);
    iat2ia.resize(nat);
    itia2iat.create(ntype, namax);
    int iat = 0;
    for (int it = 0; it < ntype; ++it)
    {
        for (int ia = 0; ia < atoms[it].na; ++ia)
        {
            itia2iat(it, ia) = iat;
            iat2it[iat] = it;
            iat2ia[iat] = ia;
            ++iat;
        }
    }
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
