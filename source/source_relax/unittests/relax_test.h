#include "source_base/matrix.h"
#include "source_cell/unitcell.h"

UnitCell::UnitCell() {};
UnitCell::~UnitCell() {};

Magnetism::Magnetism() {};
Magnetism::~Magnetism() {};

Atom::Atom() {};
Atom::~Atom() {};
Atom_pseudo::Atom_pseudo() {};
Atom_pseudo::~Atom_pseudo() {};
pseudo::pseudo() {};
pseudo::~pseudo() {};
SepPot::SepPot() {}
SepPot::~SepPot() {}
Sep_Cell::Sep_Cell() noexcept {}
Sep_Cell::~Sep_Cell() noexcept {}
int ModuleSymmetry::Symmetry::symm_flag = 0;
void ModuleSymmetry::Symmetry::symmetrize_mat3(ModuleBase::matrix &sigma, const Lattice &) const
{
    // Average over the identity and a twofold rotation about z.
    sigma(0, 2) = sigma(2, 0) = sigma(1, 2) = sigma(2, 1) = 0.0;
}
void ModuleSymmetry::Symmetry::symmetrize_vec3_nat(double *v) const {};
