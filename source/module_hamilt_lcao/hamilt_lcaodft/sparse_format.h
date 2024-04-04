#ifndef SPARSE_FORMAT_H 
#define SPARSE_FORMAT_H

#include "module_cell/module_neighbor/sltk_atom_arrange.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"


namespace sparse_format
{

void cal_dH(
		LCAO_Matrix &lm,
		LCAO_gen_fixedH &gen_h, 
		const int &current_spin, 
		const double &sparse_threshold,
		Gint_k &gint_k);

// be called by 'cal_dH_sparse'
void set_R_range(LCAO_Matrix &lm);

// be called by 'cal_dH_sparse' 
void cal_dSTN_R(
		LCAO_Matrix &lm,
		const int &current_spin, 
		const double &sparse_threshold);


}

#endif
