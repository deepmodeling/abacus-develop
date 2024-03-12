#include "local_orbital_charge.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/read_wfc_nao.h"
#include "module_base/parallel_reduce.h"
#include "module_base/parallel_common.h"
#include "module_base/memory.h"
#include "module_base/timer.h"


// allocate density kernel may change once the ion
// positions change
void Local_Orbital_Charge::allocate_gamma(
    const int& lgd,
    psi::Psi<double>* psid,
    elecstate::ElecState* pelec,
    const int& nks,
    const int& istep)
{
     ModuleBase::TITLE("Local_Orbital_Charge","allocate_gamma");

    // mohan fix serious bug 2010-09-06
    this->lgd_now = lgd;
    //xiaohui add 'GlobalV::OUT_LEVEL' line, 2015-09-16
    if(GlobalV::OUT_LEVEL != "m") ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"lgd_last",lgd_last);
    if(GlobalV::OUT_LEVEL != "m") ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"lgd_now",lgd_now);

    // mohan add 2010-07-01
    if(this->init_DM)
    {
		assert(lgd_last > 0);
		for (int is=0; is<GlobalV::NSPIN; is++)
		{
			delete[] DM[is];
			delete[] DM_pool[is];
		}
		delete[] DM;
		delete[] DM_pool;
		init_DM = false;
    }

    assert(lgd_now <= GlobalV::NLOCAL);

    // mohan update 2010-09-06
    if(lgd_now > 0)
    {
		const std::size_t size_lgd_squre = static_cast<std::size_t>(lgd_now) * lgd_now;
		this->DM = new double**[GlobalV::NSPIN];
		this->DM_pool = new double *[GlobalV::NSPIN];
		for(int is=0; is<GlobalV::NSPIN; is++)
		{
			this->DM_pool[is]=new double [size_lgd_squre];
			ModuleBase::GlobalFunc::ZEROS(DM_pool[is], size_lgd_squre);
			this->DM[is] = new double*[lgd_now];

			for (int i=0; i<lgd_now; i++)
			{
				DM[is][i] = &DM_pool[is][i*lgd_now];
			}
		}
		this->init_DM = true;
        this->lgd_last = lgd_now;
        ModuleBase::Memory::record("LOC::DM", sizeof(double) * GlobalV::NSPIN*size_lgd_squre);
        //xiaohui add 'GlobalV::OUT_LEVEL', 2015-09-16
        if(GlobalV::OUT_LEVEL != "m") GlobalV::ofs_running << " allocate DM , the dimension is " << lgd_now << std::endl;
    }
    else if(lgd_now == 0)
    {
        this->init_DM = false;
    }
    else
    {
        ModuleBase::WARNING_QUIT("Local_Orbital_Charge::allocate","lgd<0!Something Wrong!");
    }
    
#ifdef __MPI
    this->dm2g.setAlltoallvParameter(this->ParaV->comm_2D, GlobalV::NLOCAL, this->ParaV->blacs_ctxt, this->ParaV->nb, this->lgd_now, this->LOWF->gridt->trace_lo);
#endif

	// Peize Lin test 2019-01-16
    this->init_dm_2d(nks);

    if (istep == 0 && INPUT.init_wfc == "file")
    {
        this->LOWF->gamma_file(psid, pelec);
    }
    return;
}

void Local_Orbital_Charge::cal_dk_gamma_from_2D_pub(void)
{
    ModuleBase::TITLE("Local_Orbital_Charge","cal_dk_gamma_from_2D_pub");

    this->dm2g.cal_dk_gamma_from_2D(this->dm_gamma, this->DM, GlobalV::NSPIN, GlobalV::NLOCAL, this->lgd_now, GlobalV::ofs_running);
}

