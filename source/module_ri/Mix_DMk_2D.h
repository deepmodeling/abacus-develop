//=======================
// AUTHOR : Peize Lin
// DATE :   2023-05-09
//=======================

#ifndef MIX_DMK_2D_H
#define MIX_DMK_2D_H

#include "Mix_Data.h"
#include "module_base/matrix.h"
#include "module_base/complexmatrix.h"

#include <vector>

class Mix_DMk_2D
{
public:
	Mix_DMk_2D &set_nks(const int nks);
	Mix_DMk_2D &set_mixing_mode(const Mixing_Mode mixing_mode);
	Mix_DMk_2D &set_mixing_beta(const double mixing_beta);
	Mix_DMk_2D &set_coef_pulay(const int iter, const Charge_Mixing &chr_mix);

	void mix(const std::vector<ModuleBase::matrix> &dm, const bool flag_restart);
	void mix(const std::vector<ModuleBase::ComplexMatrix> &dm, const bool flag_restart);

	std::vector<const ModuleBase::matrix*> get_DMk_gamma_out() const;
	std::vector<const ModuleBase::ComplexMatrix*> get_DMk_k_out() const;

private:
	std::vector<Mix_Data<ModuleBase::matrix>> mix_DMk_gamma;
	std::vector<Mix_Data<ModuleBase::ComplexMatrix>> mix_DMk_k;
};

#endif