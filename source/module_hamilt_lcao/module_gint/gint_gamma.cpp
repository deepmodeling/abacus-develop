#include "gint_gamma.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_base/ylm.h"
#include "module_cell/module_neighbor/sltk_atom_arrange.h"
#include "module_base/timer.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __MKL
#include <mkl_service.h>
#endif

Gint_Gamma::Gint_Gamma()
{
   
    sender_index_size = 1;
	sender_local_index = nullptr;
    sender_size_process = nullptr;
    sender_displacement_process = nullptr;
    sender_size=1;
    sender_buffer=nullptr;

    receiver_index_size=1;
    receiver_global_index = nullptr;
    receiver_size_process = nullptr;
    receiver_displacement_process = nullptr;
    receiver_size=1;
    receiver_buffer=nullptr;
}

Gint_Gamma::~Gint_Gamma()
{
    // mohan add if 2024-04-09
	if(sender_local_index != nullptr)
	{
		delete[] sender_local_index;
	}

	if(sender_size_process != nullptr)
	{
		delete[] sender_size_process;
	}

	if(sender_displacement_process != nullptr)
	{
		delete[] sender_displacement_process;
	}

	if(sender_buffer != nullptr)
	{
		delete[] sender_buffer;
	}

    if(receiver_global_index != nullptr)
	{
		delete[] receiver_global_index;
	}

    if(receiver_size_process != nullptr)
	{
		delete[] receiver_size_process;
	}

    if(receiver_displacement_process != nullptr)
	{
		delete[] receiver_displacement_process;
	}

    if(receiver_buffer != nullptr)
	{
		delete[] receiver_buffer;
	}
}
