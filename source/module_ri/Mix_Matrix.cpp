//=======================
// AUTHOR : Peize Lin
// DATE :   2023-05-09
//=======================

#ifndef MIX_DATA_HPP
#define MIX_DATA_HPP

#include "Mix_Matrix.h"
#include "module_base/matrix.h"
#include "module_base/complexmatrix.h"
#include "module_base/tool_title.h"

template<typename Tdata>
void Mix_Matrix<Tdata>::mix(const Tdata &data_in, const bool flag_restart)
{
	ModuleBase::TITLE("Mix_Matrix","mix");
	if(separate_loop)
	{
			this->mixing = new Base_Mixing::Plain_Mixing(this->mixing_beta);
	}

	if(flag_restart)
	{
		this->data_out = data_in;
		this->mixing->init_mixing_data(this->matrix_data, data_in.nc*data_in.nr, sizeof(*data_in.c));
	}
	else
	{
		this->mixing->push_data(this->matrix_data, data_out.c, data_in.c, nullptr, false);
		this->mixing->mix_data(this->matrix_data, data_in.c);
	}
	
	if(separate_loop)
	{
		delete this->mixing;
		this->mixing = nullptr;
	}
}

template class Mix_Matrix<ModuleBase::matrix>;
template class Mix_Matrix<ModuleBase::ComplexMatrix>;

#endif