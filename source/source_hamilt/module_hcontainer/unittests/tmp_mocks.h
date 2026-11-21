#ifndef MODULE_HCONTAINER_TMP_MOCKS_H
#define MODULE_HCONTAINER_TMP_MOCKS_H

#include "source_cell/unitcell.h"

// Minimal definitions of the UnitCell/Atom members that the hcontainer tests
// need. They live in a header and are `inline` so that every test target can
// include them without dragging in source_cell; the real implementations are
// linked by the production binary.

inline Atom::Atom()
{
}
inline Atom::~Atom()
{
}

inline Atom_pseudo::Atom_pseudo()
{
}
inline Atom_pseudo::~Atom_pseudo()
{
}

inline Magnetism::Magnetism()
{
}
inline Magnetism::~Magnetism()
{
}

inline pseudo::pseudo()
{
}
inline pseudo::~pseudo()
{
}

inline UnitCell::UnitCell()
{
}
inline UnitCell::~UnitCell()
{
}

inline SepPot::SepPot() {}
inline SepPot::~SepPot() {}
inline Sep_Cell::Sep_Cell() noexcept {}
inline Sep_Cell::~Sep_Cell() noexcept {}

inline void Atom::set_index()
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

inline void UnitCell::set_iat2itia()
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

inline void UnitCell::set_iat2iwt(const int& npol_in)
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

#endif
