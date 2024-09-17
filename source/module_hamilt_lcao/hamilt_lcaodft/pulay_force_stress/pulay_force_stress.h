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

/// calculate $Tr[D*dH/dx]$ and $1/V Tr[D*(dH/dx_a*x_b)]
/// where D can be either density matrix or energy density matrix
namespace PulayForceStress
{
    /// for 2-center terms
    // template<typename TK, typename TR>
    // void cal_pulay_fs(
    //     ModuleBase::matrix& f,  ///< [out] force
    //     ModuleBase::matrix& s,  ///< [out] stress
    //     const elecstate::DensityMatrix<TK, TR>& dm,  ///< [in] density matrix or energy density matrix
    //     const UnitCell& ucell,  ///< [in] unit cell
    //     const ForceStressArrays& fsr,
    //     const bool& isstress
    // );
    /// for grid terms
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
        const bool& set_dmr_gint = true
    );
}

#include "pulay_force_stress.hpp"