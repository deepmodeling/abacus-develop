#pragma once
#include "module_basis/module_nao/two_center_bundle.h"
#include "module_elecstate/module_dm/density_matrix.h"
#include "module_hamilt_lcao/module_gint/gint_gamma.h"
#include "module_hamilt_lcao/module_gint/gint_k.h"
#include "module_elecstate/potentials/potential_new.h"
#include "module_cell/unitcell.h"
#include "module_hamilt_lcao/hamilt_lcaodft/stress_tools.h"
#ifndef TGINT_H
#define TGINT_H
template <typename T>
struct TGint;
template <> struct TGint<double> { using type = Gint_Gamma; };
template <> struct TGint<std::complex<double>> { using type = Gint_k; };
#endif

/// calculate the abstract formulas: 
/// $Tr[D*dH/dx]$ (force) and $1/V Tr[D*(dH/dx_a*x_b)]$ (stress)
/// where D can be any (energy) density matrix
/// and H can be any operator
namespace PulayForceStress
{
    /// for 2-center-integration terms, provided force and stress derivatives
    template<typename TK, typename TR>
    void cal_pulay_fs(
        ModuleBase::matrix& f,  ///< [out] force
        ModuleBase::matrix& s,  ///< [out] stress
        const elecstate::DensityMatrix<TK, TR>& dm,  ///< [in] density matrix or energy density matrix
        const UnitCell& ucell,  ///< [in] unit cell
        const Parallel_Orbitals& pv,  ///< [in] parallel orbitals
        const double* (&dHSx)[3],  ///< [in] dHSx x, y, z, for force
        const double* (&dHSxy)[6],  ///< [in] dHSxy 11, 12, 13, 22, 23, 33, for stress
        const bool& isforce,
        const bool& isstress,
        Record_adj* ra = nullptr,
        const double& factor_force = 1.0,
        const double& factor_stress = 1.0);

    /// for 2-center-integration terms, provided force derivatives and coordinate difference
    template<typename TK, typename TR>
    void cal_pulay_fs(
        ModuleBase::matrix& f,  ///< [out] force
        ModuleBase::matrix& s,  ///< [out] stress
        const elecstate::DensityMatrix<TK, TR>& dm,  ///< [in] density matrix or energy density matrix
        const UnitCell& ucell,  ///< [in] unit cell
        const Parallel_Orbitals& pv,  ///< [in] parallel orbitals
        const double* (&dHSx)[3],  ///< [in] dHSx x, y, z, for force and stress
        const double* dtau,  ///< [in] dr x, y, z, for stress
        const bool& isforce,
        const bool& isstress,
        Record_adj* ra = nullptr,
        const double& factor_force = 1.0,
        const double& factor_stress = 1.0);

    /// for grid-integration terms
    template<typename TK, typename TR>
    void cal_pulay_fs(
        ModuleBase::matrix& f,  ///< [out] force
        ModuleBase::matrix& s,  ///< [out] stress
        const elecstate::DensityMatrix<TK, TR>& dm,  ///< [in] density matrix or energy density matrix
        const UnitCell& ucell,  ///< [in] unit cell
        const elecstate::Potential* pot, ///< [in] potential on grid
        typename TGint<TK>::type& gint, ///< [in] Gint object
        const bool& isforce,
        const bool& isstress,
        const bool& set_dmr_gint = true);
}
#include "pulay_force_stress_center2_template.hpp"
#include "pulay_force_stress_gint.hpp"