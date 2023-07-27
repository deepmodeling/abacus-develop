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
                const int& nks)
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
		this->DM = new double**[GlobalV::NSPIN];
		this->DM_pool = new double *[GlobalV::NSPIN];
		for(int is=0; is<GlobalV::NSPIN; is++)
		{
			this->DM_pool[is]=new double [lgd_now*lgd_now];
			ModuleBase::GlobalFunc::ZEROS(DM_pool[is], lgd_now*lgd_now);
			this->DM[is] = new double*[lgd_now];

			for (int i=0; i<lgd_now; i++)
			{
				DM[is][i] = &DM_pool[is][i*lgd_now];
			}
		}
		this->init_DM = true;
        this->lgd_last = lgd_now;
        ModuleBase::Memory::record("LOC::DM", sizeof(double) * GlobalV::NSPIN*lgd_now*lgd_now);
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

    if (INPUT.init_wfc == "file")
    {
        this->gamma_file(psid, this->LOWF[0], pelec);
    }
    return;
}

void Local_Orbital_Charge::gamma_file(psi::Psi<double>* psid, Local_Orbital_wfc &lowf, elecstate::ElecState* pelec)
{
	ModuleBase::TITLE("Local_Orbital_Charge","gamma_file");

	int error;
	std::cout << " Read in gamma point wave function files " << std::endl;

	double **ctot;

    //allocate psi
    int ncol = this->ParaV->ncol_bands;
    if(GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="lapack_gvx" || GlobalV::KS_SOLVER == "scalapack_gvx"
#ifdef __CUSOLVER_LCAO
    ||GlobalV::KS_SOLVER=="cusolver"
#endif
    )
    {
        ncol = this->ParaV->ncol;
    }
    if(psid == nullptr)
    {
        ModuleBase::WARNING_QUIT("gamma_file", "psid should be allocated first!");
    }
    else
    {
        psid->resize(GlobalV::NSPIN, ncol, this->ParaV->nrow);
    }
    ModuleBase::GlobalFunc::ZEROS( psid->get_pointer(), psid->size() );

	for(int is=0; is<GlobalV::NSPIN; ++is)
	{

		GlobalV::ofs_running << " Read in wave functions " << is << std::endl;
		error = ModuleIO::read_wfc_nao( ctot , is, this->ParaV, psid, pelec);
#ifdef __MPI
		Parallel_Common::bcast_int(error);
#endif
		GlobalV::ofs_running << " Error=" << error << std::endl;
		if(error==1)
		{
			ModuleBase::WARNING_QUIT("Local_Orbital_wfc","Can't find the wave function file: LOWF.dat");
		}
		else if(error==2)
		{
			ModuleBase::WARNING_QUIT("Local_Orbital_wfc","In wave function file, band number doesn't match");
		}
		else if(error==3)
		{
			ModuleBase::WARNING_QUIT("Local_Orbital_wfc","In wave function file, nlocal doesn't match");
		}
		else if(error==4)
		{
			ModuleBase::WARNING_QUIT("Local_Orbital_wfc","In k-dependent wave function file, k point is not correct");
		}

	}//loop ispin
}

void Local_Orbital_Charge::cal_dk_gamma_from_2D_pub(void)
{
    ModuleBase::TITLE("Local_Orbital_Charge","cal_dk_gamma_from_2D_pub");

    this->dm2g.cal_dk_gamma_from_2D(this->dm_gamma, this->DM, GlobalV::NSPIN, GlobalV::NLOCAL, this->lgd_now, GlobalV::ofs_running);
}

