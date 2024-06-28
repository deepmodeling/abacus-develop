#ifndef W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_IO_TD_CURRENT_IO_H
#define W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_IO_TD_CURRENT_IO_H

#include "module_basis/module_nao/two_center_bundle.h"
#include "module_elecstate/elecstate_lcao.h"
#include "module_elecstate/module_dm/density_matrix.h"
#include "module_hamilt_lcao/hamilt_lcaodft/force_stress_arrays.h"
#include "module_psi/psi.h"

namespace ModuleIO
{
#ifdef __LCAO
/// @brief func to output current, only used in tddft
void write_current(const int istep,
                   const psi::Psi<std::complex<double>>* psi,
                   const elecstate::ElecState* pelec,
                   const K_Vectors& kv,
                   const TwoCenterBundle& two_center_bundle,
                   const Parallel_Orbitals* pv,
                   Record_adj& ra,
                   LCAO_Matrix& lm); // mohan add 2024-04-02

/// @brief calculate sum_n[𝜌_(𝑛𝑘,𝜇𝜈)] for current calculation
void cal_tmp_DM(elecstate::DensityMatrix<std::complex<double>, double>& DM, const int ik, const int nspin);

/// @brief Init DS_locR for currrent calculation
void Init_DS_tmp(const Parallel_Orbitals& pv,
                 LCAO_Matrix& lm,
                 ForceStressArrays& fsr,
                 const TwoCenterBundle& two_center_bundle);

/// @brief DS_locR will be initialized again in force calculation, so it must be destoryed here.
void destory_DS_tmp(ForceStressArrays& fsr);

#endif // __LCAO
} // namespace ModuleIO
#endif // W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_IO_TD_CURRENT_IO_H
