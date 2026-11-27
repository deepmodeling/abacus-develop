#ifndef CHG_SYMM_DETAIL_H
#define CHG_SYMM_DETAIL_H

#include <complex>

#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/module_symmetry/symmetry.h"

/**
 * @brief Internal reciprocal-space helpers for charge-density symmetrization.
 *
 * Not part of the public module_charge API: only chg_symm.cpp and
 * chg_symm_detail.cpp are expected to include this header.
 */
namespace module_charge
{
namespace detail
{

/**
 * @brief Symmetrize one reciprocal-space density component.
 */
void psymmg(std::complex<double>* rhog_part,
            const ModulePW::PW_Basis* rho_basis,
            ModuleSymmetry::Symmetry& symm);

/**
 * @brief Symmetrize the three coupled reciprocal-space spin components for nspin=4.
 */
void psymmg_soc(std::complex<double>* rhog_x,
                std::complex<double>* rhog_y,
                std::complex<double>* rhog_z,
                const ModulePW::PW_Basis* rho_basis,
                ModuleSymmetry::Symmetry& symm);

} // namespace detail
} // namespace module_charge

#endif
