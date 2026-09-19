#ifndef CHG_SYMM_H
#define CHG_SYMM_H

// TODO: make cal_rhog_symm / cal_rhog_symm_soc internal (detail or anonymous
// namespace) so that external callers only use symmetrize_rho.  Blocked by:
//   1. get_pchg_lcao/pw call the double** overload — need a symmetrize_rho
//      overload that accepts raw arrays (with nspin=4 branch).
//   2. write_mlkedf_desc symmetrizes a single component of a temporary array
//      — symmetrize_rho cannot express that yet.
//   3. setup_pot, ctrl_output_fp, read_wf2rho, update_state_rdmft already use
//      the Charge& overload and can be migrated directly.

#include <complex>

#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/module_symmetry/symmetry.h"

class Charge;

/**
 * @brief Charge-density symmetrization free functions.
 *
 * The functions are stateless: every input is passed explicitly. The
 * reciprocal-space helpers shared between translation units live in
 * module_charge::detail (see chg_symm_detail.h).
 */
namespace module_charge
{

/**
 * @brief Symmetrize charge density for all spin channels
 *
 * This is a helper function that symmetrizes the charge density
 * for all spin channels by calling cal_rhog_symm() for each spin.
 *
 * @param nspin Number of spin channels
 * @param chr Charge object containing the density
 * @param pw Plane wave basis
 * @param symm Symmetry object
 */
void symmetrize_rho(const int nspin,
                    const Charge& chr,
                    const ModulePW::PW_Basis* pw,
                    ModuleSymmetry::Symmetry& symm);

/**
 * @brief Symmetrize one spin channel of the charge density.
 *
 * The kinetic-energy density kin_r is symmetrized as well whenever its
 * buffer is allocated (meta-GGA functionals, or ELF output requested).
 */
void cal_rhog_symm(const int& spin_now,
                   const Charge& CHR,
                   const ModulePW::PW_Basis* pw,
                   ModuleSymmetry::Symmetry& symm);

/**
 * @brief Symmetrize one spin channel of raw density arrays.
 */
void cal_rhog_symm(const int& spin_now,
                   double** rho,
                   std::complex<double>** rhog,
                   int ngmc,
                   double** kin_r,
                   const ModulePW::PW_Basis* pw,
                   ModuleSymmetry::Symmetry& symm);

/**
 * @brief Symmetrize raw nspin=4 spin-density arrays with coupled spin rotations.
 *
 * @param rho Real-space density components ordered as rho0, mx, my, mz.
 * @param rhog Reciprocal-space work arrays with the same component ordering.
 * @param pw Plane-wave basis used for the Fourier transforms.
 * @param symm Symmetry operations and spin rotations.
 */
void cal_rhog_symm_soc(double** rho,
                       std::complex<double>** rhog,
                       const ModulePW::PW_Basis* pw,
                       ModuleSymmetry::Symmetry& symm);

} // namespace module_charge

#endif
