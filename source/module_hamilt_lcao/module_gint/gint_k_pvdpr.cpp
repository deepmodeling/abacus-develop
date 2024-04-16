#include "gint_k.h"
#include "grid_technique.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/memory.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_base/tool_threading.h"
#include "module_base/ylm.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

void Gint_k::allocate_pvdpR(void)
{
    ModuleBase::TITLE("Gint_k","allocate_pvpR");

    const int nspin = GlobalV::NSPIN;
    assert(nspin>0);

    //xiaohui modify 2015-05-30
    // the number of matrix element <phi_0 | V | dphi_R> is this->gridt->nnrg.
    this->pvdpRx_reduced = new double*[nspin];
    this->pvdpRy_reduced = new double*[nspin];
    this->pvdpRz_reduced = new double*[nspin];
    for(int is =0;is<nspin;is++)
    {
        this->pvdpRx_reduced[is] = new double[this->gridt->nnrg];	
        ModuleBase::GlobalFunc::ZEROS( pvdpRx_reduced[is], this->gridt->nnrg);
        this->pvdpRy_reduced[is] = new double[this->gridt->nnrg];	
        ModuleBase::GlobalFunc::ZEROS( pvdpRy_reduced[is], this->gridt->nnrg);
        this->pvdpRz_reduced[is] = new double[this->gridt->nnrg];	
        ModuleBase::GlobalFunc::ZEROS( pvdpRz_reduced[is], this->gridt->nnrg);
    }

    ModuleBase::Memory::record("pvdpR_reduced", 3 * sizeof(double) * this->gridt->nnrg * nspin);

    return;
}

void Gint_k::destroy_pvdpR(void)
{
    ModuleBase::TITLE("Gint_k","destroy_pvpR");

    const int nspin = GlobalV::NSPIN;
    assert(nspin>0);
    
	for(int is =0;is<nspin;is++) 
	{
		delete[] pvdpRx_reduced[is];
	}
    delete[] pvdpRx_reduced;

	for(int is =0;is<nspin;is++) 
	{
		delete[] pvdpRy_reduced[is];
	}
    delete[] pvdpRy_reduced;

	for(int is =0;is<nspin;is++) 
	{
		delete[] pvdpRz_reduced[is];
	}
    delete[] pvdpRz_reduced;

    return;
}
