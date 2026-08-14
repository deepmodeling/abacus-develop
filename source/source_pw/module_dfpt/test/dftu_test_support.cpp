// ============================================================
// Minimal test-support definitions for constructing Plus_U in
// the DFPT unit tests (DFT+U interface reservation, U0).
//
// In production these symbols live in module_dftu/dftu.cpp,
// which pulls a large link closure (init -> dftu_io/occup ->
// scalapack ...). The DFPT tests only need to *construct* a
// Plus_U and read the public inline accessors, so the ctor,
// dtor and static data members are replicated here instead.
// This keeps the DFPT (PW) unit tests free of the LCAO-side
// DFT+U dependency. Keep in sync with dftu.cpp.
// ============================================================

#include "source_lcao/module_dftu/dftu.h"

#include <vector>

double Plus_U::energy_u = 0.0;

std::vector<double> Plus_U::U = {};

std::vector<double> Plus_U::U0 = {};

std::vector<int> Plus_U::orbital_corr = {};

double Plus_U::uramping = 0.0;

int Plus_U::omc = 0;

int Plus_U::mixing_dftu = 0;

int Plus_U::nspin = 0;

bool Plus_U::Yukawa = false;

Plus_U::Plus_U()
{}

Plus_U::~Plus_U()
{}