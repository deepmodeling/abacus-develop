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

template<>
void Mix_Matrix<ModuleBase::matrix>::mix(const ModuleBase::matrix& data_in, const bool flag_restart)
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
		this->mixing->mix_data(this->matrix_data, data_out.c);
	}
	
	if(separate_loop)
	{
		delete this->mixing;
		this->mixing = nullptr;
	}
}
template<>
void Mix_Matrix<ModuleBase::ComplexMatrix>::mix(const ModuleBase::ComplexMatrix& data_in, const bool flag_restart)
{
    ModuleBase::TITLE("Mix_Matrix", "mix");
    if (separate_loop)
    {
        this->mixing = new Base_Mixing::Plain_Mixing(this->mixing_beta);
    }

    if (flag_restart)
    {
        this->data_out = data_in;
        this->mixing->init_mixing_data(this->matrix_data, data_in.nc * data_in.nr, sizeof(*data_in.c));
    }
    else
    {
        this->mixing->push_data(this->matrix_data, data_out.c, data_in.c, nullptr, false);
        this->mixing->mix_data(this->matrix_data, data_out.c);
    }

    if (separate_loop)
    {
        delete this->mixing;
        this->mixing = nullptr;
    }
}

template<>
void Mix_Matrix<std::vector<double>>::mix(const std::vector<double>& data_in, const bool flag_restart)
{
    ModuleBase::TITLE("Mix_Matrix", "mix");
    if (separate_loop)
    {
        this->mixing = new Base_Mixing::Plain_Mixing(this->mixing_beta);
    }

    if (flag_restart)
    {
        this->data_out = data_in;
        this->mixing->init_mixing_data(this->matrix_data, data_in.size(), sizeof(*data_in.data()));
    }
    else
    {
        this->mixing->push_data(this->matrix_data, data_out.data(), data_in.data(), nullptr, false);
        this->mixing->mix_data(this->matrix_data, data_out.data());
    }

    if (separate_loop)
    {
        delete this->mixing;
        this->mixing = nullptr;
    }
}
template<>
void Mix_Matrix<std::vector<std::complex<double>>>::mix(const std::vector<std::complex<double>>& data_in, const bool flag_restart)
{
    ModuleBase::TITLE("Mix_Matrix", "mix");
    if (separate_loop)
    {
        this->mixing = new Base_Mixing::Plain_Mixing(this->mixing_beta);
    }

    if (flag_restart)
    {
        this->data_out = data_in;
        this->mixing->init_mixing_data(this->matrix_data, data_in.size(), sizeof(*data_in.data()));
    }
    else
    {
        this->mixing->push_data(this->matrix_data, data_out.data(), data_in.data(), nullptr, false);
        this->mixing->mix_data(this->matrix_data, data_out.data());
    }

    if (separate_loop)
    {
        delete this->mixing;
        this->mixing = nullptr;
    }
}

/*
// mix for ct::Tensor if useful in the future. comment out for UT-coverage
template<>
void Mix_Matrix<ct::Tensor>::mix(const ct::Tensor& data_in, const bool flag_restart)
{
    ModuleBase::TITLE("Mix_Matrix", "mix");
    if (separate_loop)
    {
        this->mixing = new Base_Mixing::Plain_Mixing(this->mixing_beta);
    }

    if (flag_restart)
    {
        this->data_out = data_in;
        this->mixing->init_mixing_data(this->matrix_data, data_in.NumElements(), sizeof(data_in.SizeOfType(data_in.data_type())));
    }
    else
    {
        switch (data_in.data_type())
        {
        case ct::DataType::DT_DOUBLE:
            this->mixing->push_data(this->matrix_data, data_out.data<double>(), data_in.data<double>(), nullptr, false);
            this->mixing->mix_data(this->matrix_data, data_out.data<double>());
            break;
        case ct::DataType::DT_COMPLEX_DOUBLE:
            this->mixing->push_data(this->matrix_data, data_out.data<std::complex<double>>(), data_in.data<std::complex<double>>(), nullptr, false);
            this->mixing->mix_data(this->matrix_data, data_out.data<std::complex<double>>());
            break;
        default:
            throw std::invalid_argument("Mix_Matrix: data type not supported");
            break;
        }
    }

    if (separate_loop)
    {
        delete this->mixing;
        this->mixing = nullptr;
    }
}
*/

template class Mix_Matrix<ModuleBase::matrix>;
template class Mix_Matrix<ModuleBase::ComplexMatrix>;
template class Mix_Matrix<std::vector<double>>;
template class Mix_Matrix<std::vector<std::complex<double>>>;
// template class Mix_Matrix<ct::Tensor>;

#endif