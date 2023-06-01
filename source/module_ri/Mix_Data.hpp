//=======================
// AUTHOR : Peize Lin
// DATE :   2023-05-09
//=======================

#ifndef MIX_DATA_HPP
#define MIX_DATA_HPP

#include "Mix_Data.h"
#include "module_base/tool_title.h"
#include "module_elecstate/module_charge/charge_mixing.h"

template<typename Tdata>
void Mix_Data<Tdata>::mix(const Tdata &data_in, const bool flag_restart)
{
	ModuleBase::TITLE("Mix_Data","mix");
	switch(this->mixing_mode)
	{
		case Mixing_Mode::No:
			this->data_out = data_in;
			break;
		case Mixing_Mode::Plain:
			if(flag_restart)
				this->data_out = data_in;
			else
				this->data_out =     this->mixing_beta  * data_in
				                + (1-this->mixing_beta) * this->data_out;
			break;
		case Mixing_Mode::Pulay:
			this->pulay_mixing(data_in, flag_restart);
			break;
		default:
			throw std::domain_error(std::string(__FILE__)+" line "+std::to_string(__LINE__));	break;
	}	
}

template<typename Tdata>
void Mix_Data<Tdata>::pulay_mixing(const Tdata &data_in, const bool flag_restart)
{
	ModuleBase::TITLE("Mix_Data","pulay_mixing");
	if(flag_restart)
	{
		this->data_pulay_list.clear();
		this->data_out = data_in;
	}
	else
	{
		this->data_pulay_list.push_back(
		         this->mixing_beta  * data_in
			+ (1-this->mixing_beta) * this->data_out);
		//while(this->data_pulay_list.size() > this->coef_pulay_list.size())
		if(this->data_pulay_list.size() == 1+this->coef_pulay_list.size())
			this->data_pulay_list.pop_front();
		assert(this->data_pulay_list.size() == this->coef_pulay_list.size());
		assert(this->coef_pulay_list.size() > 0);
		
		this->data_out = this->coef_pulay_list[0] * this->data_pulay_list[0];
		for(std::size_t i=1; i<this->data_pulay_list.size(); ++i)
			this->data_out += this->coef_pulay_list[i] * this->data_pulay_list[i];
	}
}

template<typename Tdata>
template<typename ChgMix>
void Mix_Data<Tdata>::set_coef_pulay(const int iter, const ChgMix& chr_mix)
{
	ModuleBase::TITLE("Mix_Data","set_coef_pulay");
	const std::size_t coef_size = std::min(iter-1, chr_mix.get_mixing_ndim());
	assert(coef_size>=0);
	this->coef_pulay_list.resize(coef_size);
	if(coef_size==1)
	{
		this->coef_pulay_list[0] = 1;
	}
	else if(coef_size>1)
	{
		auto pos_mod = [](const int i, const int N) ->int { return (i%N+N)%N; };
		auto index = [&](const int i) -> int
		{
			const int alpha_size = coef_size-1;
			const int alpha_begin = chr_mix.get_idstep() - alpha_size;
			return pos_mod( alpha_begin+i, chr_mix.get_dstep() );
		};		
		this->coef_pulay_list[0] = -chr_mix.get_alpha()[index(0)];
		for(std::size_t i=1; i<=coef_size-2; ++i)
			this->coef_pulay_list[i] = ( chr_mix.get_alpha()[index(i-1)] - chr_mix.get_alpha()[index(i)] );
		this->coef_pulay_list[coef_size-1] = 1+chr_mix.get_alpha()[index(coef_size-2)];
	}
}

#endif