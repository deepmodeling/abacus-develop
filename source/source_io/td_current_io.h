#ifndef TD_CURRENT_IO_H
#define TD_CURRENT_IO_H

#include "source_basis/module_nao/two_center_bundle.h"
#include "source_estate/elecstate_lcao.h"
#include "source_estate/module_dm/density_matrix.h"
#include "source_psi/psi.h"
#include "source_lcao/module_rt/velocity_op.h"
#include "source_lcao/setup_exx.h"
#ifdef __EXX
#include <RI/global/Tensor.h>
#endif

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
/// @brief func to output current calculated using i[r,H] directly
template <typename TR>
void write_current(
    const UnitCell& ucell,
    const Grid_Driver& GridD,
    const int istep,
    const psi::Psi<std::complex<double>>* psi,
    const elecstate::ElecState* pelec,
    const K_Vectors& kv,
    const Parallel_Orbitals* pv,
    const LCAO_Orbitals& orb,
    cal_r_overlap_R& r_calculator,
    const hamilt::HContainer<TR>* sR,
    const hamilt::HContainer<TR>* hR,
    const Exx_NAO<std::complex<double>>& exx_nao
);
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
void set_rR_from_hR(const UnitCell& ucell,
                    const Grid_Driver& GridD,
                    const LCAO_Orbitals& orb,
                    const Parallel_Orbitals* pv,
                    cal_r_overlap_R& r_calculator,
                    const hamilt::HContainer<std::complex<double>>* hR,
                    ModuleBase::Vector3<hamilt::HContainer<double>*>& rR);
template <typename TR>
void sum_HR(
    const UnitCell& ucell,
    const Parallel_Orbitals& pv,
    const K_Vectors& kv,
    const hamilt::HContainer<TR>* hR,
    hamilt::HContainer<std::complex<double>>* full_hR,
    const Exx_NAO<std::complex<double>>& exx_nao
);

template <typename Tadd, typename Tfull>
void add_HR(const hamilt::HContainer<Tadd>* hR, hamilt::HContainer<Tfull>* full_hR);

void init_from_adj(const UnitCell& ucell,
                   const Grid_Driver& GridD,
                   const LCAO_Orbitals& orb,
                   const Parallel_Orbitals* pv,
                   std::vector<AdjacentAtomInfo>& adjs_all,
                   ModuleBase::Vector3<hamilt::HContainer<double>*>& rR);
template <typename TR, typename TA>
void init_from_hR(const hamilt::HContainer<TR>* hR, hamilt::HContainer<TA>* aimR);
template <typename TR>
void cal_velocity_basis_k(const UnitCell& ucell,
                          const LCAO_Orbitals& orb,
                          const Parallel_Orbitals* pv,
                          const K_Vectors& kv,
                          const ModuleBase::Vector3<hamilt::HContainer<double>*>& rR,
                          const hamilt::HContainer<TR>& sR,
                          const hamilt::HContainer<std::complex<double>>& hR,
                          std::vector<ModuleBase::Vector3<std::complex<double>*>>& velocity_basis_k);

void cal_velocity_matrix(const psi::Psi<std::complex<double>>* psi,
                         const Parallel_Orbitals* pv,
                         const K_Vectors& kv,
                         const std::vector<ModuleBase::Vector3<std::complex<double>*>>& velocity_basis_k,
                         std::vector<std::array<ModuleBase::ComplexMatrix, 3>>& velocity_k);
template <typename TR>
void cal_current_comm_k(const UnitCell& ucell,
                        const Grid_Driver& GridD,
                        const LCAO_Orbitals& orb,
                        const Parallel_Orbitals* pv,
                        const K_Vectors& kv,
                        cal_r_overlap_R& r_calculator,
                        const hamilt::HContainer<TR>& sR,
                        const hamilt::HContainer<std::complex<double>>& hR,
                        const psi::Psi<std::complex<double>>* psi,
                        const elecstate::ElecState* pelec,
                        std::vector<ModuleBase::Vector3<double>>& current_k);
#endif // __LCAO
} // namespace ModuleIO
#endif 
