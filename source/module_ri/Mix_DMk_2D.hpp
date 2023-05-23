//=======================
// AUTHOR : Peize Lin
// DATE :   2023-05-09
//=======================

#include "Mix_DMk_2D.h"
#include "module_base/tool_title.h"

#ifndef MIX_DMK_2D_HPP
#define MIX_DMK_2D_HPP

template<typename ChgMix>
Mix_DMk_2D& Mix_DMk_2D::set_coef_pulay(const int iter, const ChgMix& chr_mix)
{
	ModuleBase::TITLE("Mix_DMk_2D","set_coef_pulay");
	if(this->gamma_only)
		for(Mix_Data<ModuleBase::matrix> &mix_one : this->mix_DMk_gamma)
			mix_one.set_coef_pulay(iter, chr_mix);
	else
		for(Mix_Data<ModuleBase::ComplexMatrix> &mix_one : this->mix_DMk_k)
			mix_one.set_coef_pulay(iter, chr_mix);
	return *this;
}

#endif