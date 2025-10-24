#ifndef TD_CURRENT_IO_H
#define TD_CURRENT_IO_H

#include "source_basis/module_nao/two_center_bundle.h"
#include "source_estate/elecstate_lcao.h"
#include "source_estate/module_dm/density_matrix.h"
#include "source_psi/psi.h"
#include "source_lcao/module_rt/velocity_op.h"

namespace ModuleIO
{
#ifdef __LCAO
/// @brief func to output current, only used in tddft
template <typename TR>
void write_current_eachk(const UnitCell& ucell,
                        const int istep,
                        const psi::Psi<std::complex<double>>* psi,
                        const elecstate::ElecState* pelec,
                        const K_Vectors& kv,
                        const TwoCenterIntegrator* intor,
                        const Parallel_Orbitals* pv,
                        const LCAO_Orbitals& orb,
                        const Velocity_op<TR>* cal_current,
                        Record_adj& ra);
template <typename TR>
void write_current(const UnitCell& ucell,
                const int istep,
                const psi::Psi<std::complex<double>>* psi,
                const elecstate::ElecState* pelec,
                const K_Vectors& kv,
                const TwoCenterIntegrator* intor,
                const Parallel_Orbitals* pv,
                const LCAO_Orbitals& orb,
                const Velocity_op<TR>* cal_current,
                Record_adj& ra);

/// @brief calculate sum_n[𝜌_(𝑛𝑘,𝜇𝜈)] for current calculation
void cal_tmp_DM_k(const UnitCell& ucell,
                elecstate::DensityMatrix<std::complex<double>, double>& DM_real,
                elecstate::DensityMatrix<std::complex<double>, double>& DM_imag,
                const int ik,
                const int nspin,
                const int is,
                const bool reset = true);

void cal_tmp_DM(const UnitCell& ucell,
                elecstate::DensityMatrix<std::complex<double>, double>& DM_real,
                elecstate::DensityMatrix<std::complex<double>, double>& DM_imag,
                const int nspin);

#endif // __LCAO
} // namespace ModuleIO
#endif 
