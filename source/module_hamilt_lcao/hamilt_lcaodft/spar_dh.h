#ifndef SPARSE_FORMAT_DH_H 
#define SPARSE_FORMAT_DH_H

#include "module_cell/module_neighbor/sltk_atom_arrange.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"

namespace sparse_format
{
	void cal_dH(
			LCAO_Matrix &lm,
			Grid_Driver &grid,
			LCAO_gen_fixedH &gen_h, 
			const int &current_spin, 
			const double &sparse_thr,
			Gint_k &gint_k);

	// be called by 'cal_dH_sparse'
	void set_R_range(
			std::set<Abfs::Vector3_Order<int>> &all_R_coor,
			Grid_Driver &grid);

	// be called by 'cal_dH_sparse' 
	void cal_dSTN_R(
			LCAO_Matrix &lm,
		    Grid_Driver &grid,
			const int &current_spin, 
			const double &sparse_thr);
}

#endif
