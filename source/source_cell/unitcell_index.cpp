/**
 * @file unitcell_index.cpp
 * @brief Atom/orbital indexing helpers of UnitCell: set_iat2itia, set_iat2iwt.
 */
#include "unitcell.h"

#include <cassert>

void UnitCell::set_iat2itia() {
    assert(nat > 0);
    this->iat2it.resize(nat);
    this->iat2ia.resize(nat);
    int iat = 0;
    for (int it = 0; it < ntype; it++) {
        for (int ia = 0; ia < atoms[it].na; ia++) {
            this->iat2it[iat] = it;
            this->iat2ia[iat] = ia;
            ++iat;
        }
    }
    return;
}

void UnitCell::set_iat2iwt(const int& npol_in)
{
#ifdef __DEBUG
    assert(npol_in == 1 || npol_in == 2);
    assert(this->nat > 0);
    assert(this->ntype > 0);
#endif
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
