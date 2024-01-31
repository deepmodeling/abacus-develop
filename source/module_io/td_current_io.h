#ifndef TD_CURRENT_H
#define TD_CURRENT_H

#include "module_elecstate/module_dm/density_matrix.h"
#include "module_elecstate/elecstate_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_hamilt.h"
#include "module_psi/psi.h"

namespace ModuleIO
{
void write_current(const int istep,
                    const psi::Psi<std::complex<double>>* psi,
                    const elecstate::ElecState* pelec,
                    const K_Vectors& kv,
                    const Parallel_Orbitals* pv,
                    Record_adj& ra,
                    LCAO_Hamilt& UHM);
void cal_tmp_DM(elecstate::DensityMatrix<std::complex<double>, double>& DM, const int ik, const int nspin);
void Init_DS_tmp(const Parallel_Orbitals& pv,LCAO_Hamilt& UHM);
void destory_DS_tmp(LCAO_Hamilt& UHM);
}

#endif // TD_CURRENT_H