#ifndef EXX_GLOBAL_H
#define EXX_GLOBAL_H

#include "xc_functional.h"

struct Exx_Global
{
	enum class Hybrid_Type {No,HF,PBE0,HSE,Generate_Matrix};
	struct Exx_Info
	{
		Exx_Global::Hybrid_Type hybrid_type;

		double hybrid_alpha = 0.25;
		double hse_omega = 0.11;
		
		bool separate_loop = true;
		size_t hybrid_step = 1;
	};	
	Exx_Info info;
};

#endif
