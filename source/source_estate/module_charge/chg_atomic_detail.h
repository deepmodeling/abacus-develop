#ifndef CHG_ATOMIC_DETAIL_H
#define CHG_ATOMIC_DETAIL_H

// Internal helpers for atomic_rho (chg_atomic.cpp).
// Not part of the public module_charge API: only chg_atomic.cpp and
// chg_atomic_inner.cpp are expected to include this header.

#include <vector>

#include "source_base/complexmatrix.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/unitcell.h"

#include <ostream>

namespace module_charge
{
namespace detail
{

/// Context bundling commonly-used data for rho_g3d fill helpers.
struct RhoG3dCtx
{
    ModuleBase::ComplexMatrix& rho_g3d;
    const ModuleBase::ComplexMatrix& strucFac;
    const std::vector<double>& rho_lgl;
    const ModulePW::PW_Basis* rhopw;
    int it;
};

/// Compute the radial atomic charge density rhoatm from pseudopotential.
std::vector<double> compute_rhoatm(const Atom& atom,
                                   int mesh,
                                   std::ostream& ofs_warning);

/// Compute the 1D charge in G space from rhoatm.
std::vector<double> compute_rho_lgl(const Atom& atom,
                                    const ModulePW::PW_Basis* rhopw,
                                    const UnitCell& ucell,
                                    const std::vector<double>& rhoatm,
                                    int test_charge,
                                    double omega,
                                    std::ostream& ofs_warning);

/// Fill rho_g3d for nspin==1 case.
void fill_rho_g3d_nspin1(RhoG3dCtx& ctx);

/// Fill rho_g3d for nspin==2 case (both startmag_type 1 and 2).
void fill_rho_g3d_nspin2(RhoG3dCtx& ctx,
                         int startmag_type,
                         double start_mag,
                         const Atom& atom);

/// Fill rho_g3d for nspin==4, startmag_type==1 case.
void fill_rho_g3d_nspin4_type1(RhoG3dCtx& ctx,
                               double start_mag,
                               const Atom& atom,
                               bool domag,
                               bool domag_z);

/// Fill rho_g3d for nspin==4, startmag_type==2 case.
void fill_rho_g3d_nspin4_type2(RhoG3dCtx& ctx,
                               const Atom& atom,
                               bool domag,
                               bool domag_z);

/// FFT rho_g3d to real space, check for negative/imaginary charge,
/// and normalize to target electron number.
void normalize_and_check(double** rho_in,
                         const ModuleBase::ComplexMatrix& rho_g3d,
                         const ModulePW::PW_Basis* rhopw,
                         int spin_number_need,
                         double omega,
                         std::ostream& ofs_warning,
                         double nelec);

} // namespace detail
} // namespace module_charge

#endif // CHG_ATOMIC_DETAIL_H
