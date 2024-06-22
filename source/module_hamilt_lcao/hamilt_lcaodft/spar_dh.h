#ifndef SPARSE_FORMAT_DH_H
#define SPARSE_FORMAT_DH_H

#include "module_cell/module_neighbor/sltk_atom_arrange.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_lcao/hamilt_lcaodft/force_stress_arrays.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

namespace sparse_format
{
void cal_dH(LCAO_Matrix& lm,
            Grid_Driver& grid,
            const TwoCenterBundle& two_center_bundle,
            const int& current_spin,
            const double& sparse_thr,
            Gint_k& gint_k);

// be called by 'cal_dH_sparse'
void set_R_range(std::set<Abfs::Vector3_Order<int>>& all_R_coor, Grid_Driver& grid);

// be called by 'cal_dH_sparse'
void cal_dSTN_R(LCAO_Matrix& lm,
                ForceStressArrays& fsr, // mohan add 2024-06-16
                Grid_Driver& grid,
                const int& current_spin,
                const double& sparse_thr);
} // namespace sparse_format

#endif
