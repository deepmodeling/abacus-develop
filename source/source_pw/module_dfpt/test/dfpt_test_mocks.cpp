// ============================================================
// Minimal test-support ctor/dtor stubs shared by all DFPT unit
// test binaries (test/ and test_serial/), mirroring the
// tmp_mocks.cpp convention of the other module test suites.
//
// In production these symbols live in full link closures the DFPT
// tests do not want to pull in:
//   - cell/spepot/stru_fac/charge_mixing closures (via UnitCell),
//   - the pwdft DFT+U closure (Plus_U_Base, see the former
//     dftu_test_support.cpp this file absorbed).
// The DFPT tests only need default-constructible objects, so the
// empty definitions are replicated here once instead of in every
// test translation unit. Keep the signatures in sync with the
// production sources.
// ============================================================

#include "source_cell/atom_pseudo.h"
#include "source_cell/atom_spec.h"
#include "source_cell/magnetism.h"
#include "source_cell/pseudo.h"
#include "source_cell/unitcell.h"
#include "source_estate/module_charge/charge_mixing.h"
#include "source_pw/module_pwdft/dftu_base.h"
#include "source_pw/module_pwdft/stru_fac.h"

pseudo::pseudo()
{
}
pseudo::~pseudo()
{
}
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
SepPot::SepPot()
{
}
SepPot::~SepPot()
{
}
Sep_Cell::Sep_Cell() noexcept
{
}
Sep_Cell::~Sep_Cell() noexcept
{
}
UnitCell::UnitCell()
{
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

Structure_Factor::Structure_Factor()
{
}
Structure_Factor::~Structure_Factor()
{
}

Charge_Mixing::~Charge_Mixing()
{
}

Plus_U_Base::Plus_U_Base()
{
}
Plus_U_Base::~Plus_U_Base()
{
}
