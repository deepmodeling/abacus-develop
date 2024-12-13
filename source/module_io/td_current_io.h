#ifndef W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_IO_TD_CURRENT_IO_H
#define W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_IO_TD_CURRENT_IO_H

#include "module_basis/module_nao/two_center_bundle.h"
#include "module_elecstate/elecstate_lcao.h"
#include "module_elecstate/module_dm/density_matrix.h"
#include "module_psi/psi.h"

namespace ModuleIO
{
#ifdef __LCAO
/// @brief func to output current, only used in tddft
void write_current(const UnitCell& ucell,
                   const Grid_Driver& gd,
                   const int istep,
                   const psi::Psi<std::complex<double>>* psi,
                   const elecstate::ElecState* pelec,
                   const K_Vectors& kv,
                   const TwoCenterIntegrator* intor,
                   const Parallel_Orbitals* pv,
                   const LCAO_Orbitals& orb,
                   Record_adj& ra);

/// @brief calculate sum_n[𝜌_(𝑛𝑘,𝜇𝜈)] for current calculation
void cal_tmp_DM(elecstate::DensityMatrix<std::complex<double>, double>& DM_real,
                elecstate::DensityMatrix<std::complex<double>, double>& DM_imag,
                const int ik,
                const int nspin,
                const int is);

#endif // __LCAO
} // namespace ModuleIO
#endif // W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_IO_TD_CURRENT_IO_H
