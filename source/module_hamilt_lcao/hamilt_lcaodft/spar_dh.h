#ifndef W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_HAMILT_LCAO_HAMILT_LCAODFT_SPAR_DH_H
#define W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_HAMILT_LCAO_HAMILT_LCAODFT_SPAR_DH_H

#include "module_cell/module_neighbor/sltk_atom_arrange.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_HS_arrays.hpp"
#include "module_hamilt_lcao/hamilt_lcaodft/force_stress_arrays.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include <vector>

namespace sparse_format
{
void cal_dH(const UnitCell& ucell,
            const Parallel_Orbitals& pv,
            LCAO_HS_Arrays& HS_Arrays,
            Grid_Driver& grid,
            const TwoCenterBundle& two_center_bundle,
            const LCAO_Orbitals& orb,
            const int& current_spin,
            const double& sparse_thr,
            Gint_k& gint_k);

// be called by 'cal_dH_sparse'
void set_R_range(std::set<Abfs::Vector3_Order<int>>& all_R_coor, Grid_Driver& grid);

// be called by 'cal_dH_sparse'
void cal_dSTN_R(const UnitCell& ucell,
                const Parallel_Orbitals& pv,
                LCAO_HS_Arrays& HS_Arrays,
                ForceStressArrays& fsr, // mohan add 2024-06-16
                Grid_Driver& grid,
                const std::vector<double>& orb_cutoff,
                const int& current_spin,
                const double& sparse_thr);

void destroy_dH_R_sparse(LCAO_HS_Arrays& HS_Arrays);

} // namespace sparse_format

#endif
