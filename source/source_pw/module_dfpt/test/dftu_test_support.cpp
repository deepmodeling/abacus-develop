// ============================================================
// Minimal test-support definitions for constructing the DFT+U
// provider in the DFPT unit tests (DFT+U interface reservation, U0).
//
// In production these symbols live in
// source_pw/module_pwdft/dftu_base.cpp, which is part of the pwdft
// object library and pulls the PW-side DFT+U link closure. The DFPT
// tests only need a default-constructed Plus_U_Base (occupation
// matrices not initialized -> u_active() false), so the ctor and dtor
// are replicated here instead. This keeps the DFPT unit tests free of
// that link closure. Keep in sync with dftu_base.cpp.
// ============================================================

#include "source_pw/module_pwdft/dftu_base.h"

Plus_U_Base::Plus_U_Base()
{}

Plus_U_Base::~Plus_U_Base()
{}
