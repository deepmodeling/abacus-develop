#ifndef SPAR_DH_H
#define SPAR_DH_H

#include "source_cell/module_neighbor/sltk_atom_arrange.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_lcao/LCAO_HS_arrays.hpp"
#include "source_lcao/force_stress_arrays.h"
#include "source_lcao/hamilt_lcao.h"
#include "source_pw/module_pwdft/global.h"
#include <vector>

namespace sparse_format
{
void cal_dH(const UnitCell& ucell,
            const Parallel_Orbitals& pv,
            LCAO_HS_Arrays& HS_Arrays,
            const Grid_Driver& grid,
            const TwoCenterBundle& two_center_bundle,
            const LCAO_Orbitals& orb,
            const int& current_spin,
            const double& sparse_thr,
            const ModuleBase::matrix& v_eff,
            Gint_k& gint_k);

// calculated the derivative of the overlap matrix: <phi|dphi>
void cal_dS(const UnitCell& ucell,
            const Parallel_Orbitals& pv,
            LCAO_HS_Arrays& HS_Arrays,
            const Grid_Driver& grid,
            const TwoCenterBundle& two_center_bundle,
            const LCAO_Orbitals& orb,
            const double& sparse_thr);

// be called by 'cal_dH_sparse'
void set_R_range(std::set<Abfs::Vector3_Order<int>>& all_R_coor, const Grid_Driver& grid);

// be called by 'cal_dH_sparse'
void cal_dSTN_R(const UnitCell& ucell,
                const Parallel_Orbitals& pv,
                LCAO_HS_Arrays& HS_Arrays,
                ForceStressArrays& fsr, // mohan add 2024-06-16
                const Grid_Driver& grid,
                const std::vector<double>& orb_cutoff,
                const int& current_spin,
                const double& sparse_thr);

void destroy_dH_R_sparse(LCAO_HS_Arrays& HS_Arrays);

} // namespace sparse_format

#endif
