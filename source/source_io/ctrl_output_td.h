#ifndef CTRL_OUTPUT_TD_H
#define CTRL_OUTPUT_TD_H

#include "source_basis/module_ao/ORB_read.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_basis/module_nao/two_center_bundle.h"
#include "source_cell/unitcell.h"
#include "source_estate/elecstate_lcao.h"
#include "source_estate/module_pot/H_TDDFT_pw.h"
#include "source_lcao/module_rt/td_info.h"
#include "source_lcao/module_rt/velocity_op.h"
#include "source_lcao/record_adj.h"
#include "source_psi/psi.h"
#include "source_lcao/hamilt_lcao.cpp"
#ifdef __EXX
#include <RI/global/Tensor.h>
#endif

namespace ModuleIO
{

template <typename TR>
void ctrl_output_td(const UnitCell& ucell,
                    double** rho_save,
                    const ModulePW::PW_Basis* rhopw,
                    const int istep,
                    const psi::Psi<std::complex<double>>* psi,
                    const elecstate::ElecState* pelec,
                    const K_Vectors& kv,
                    const TwoCenterIntegrator* intor,
                    const Parallel_Orbitals* pv,
                    const LCAO_Orbitals& orb,
                    const Velocity_op<TR>* velocity_mat,
                    const Grid_Driver& grid,
                    hamilt::HamiltLCAO<std::complex<double>, TR>* p_hamilt,
                    Record_adj& RA,
                    TD_info* td_p
#ifdef __EXX
                    ,
                    const std::vector<std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<std::complex<double>>>>>* hexxs = nullptr
#endif
                    );

} // namespace ModuleIO

#endif // CTRL_OUTPUT_TD_H