#include "local_orbital_charge.h"
#include "../src_pw/global.h"
#include "../module_base/blas_connector.h"
#include "../module_base/timer.h"

//#include "../src_onscaling/on_tests.h"
//#include "../src_siao/selinv.h"
// 2014.10.29 add memory pool for DM and by yshen

// Shen Yu add 2019/5/9
extern "C"
{
    void Cblacs_gridinfo(int icontxt, int* nprow, int *npcol, int *myprow, int *mypcol);
    void Cblacs_pinfo(int *myid, int *nprocs);
    void Cblacs_pcoord(int icontxt, int pnum, int *prow, int *pcol);
    int Cblacs_pnum(int icontxt, int prow, int pcol);
}

int Local_Orbital_Charge::out_dm = 0;

Local_Orbital_Charge::Local_Orbital_Charge()
{
    // for gamma algorithms.
    this->init_DM = false;  
    this->lgd_now = 0;
    this->lgd_last = 0;

    // for k-dependent algorithms.
    this->init_DM_R = false;

    //xiaohui add 2014-06-19
    //band_local = nullptr;
    //Z_wg = nullptr;
    //Z_LOC = nullptr;
    sender_2D_index = nullptr;
    sender_size_process = nullptr;
    sender_displacement_process = nullptr;
    sender_buffer = nullptr;

    receiver_local_index = nullptr;
    receiver_size_process = nullptr;
    receiver_displacement_process = nullptr;
    receiver_buffer = nullptr;
}

Local_Orbital_Charge::~Local_Orbital_Charge()
{
    // with gamma point only
     if (this->init_DM)
	 {
		 for (int is=0; is<GlobalV::NSPIN; is++)
		 {
			 delete[] DM[is];
			 delete[] DM_pool[is];
		 }
		 delete[] DM;
		 delete[] DM_pool;
		 delete[] sender_2D_index;
		 delete[] sender_size_process;
		 delete[] sender_displacement_process;

		 delete[] receiver_local_index;
		 delete[] receiver_size_process;
		 delete[] receiver_displacement_process;
	 }

    // with k points
    if (this->init_DM_R)
    {
        for(int is=0; is<GlobalV::NSPIN; is++)
        {
            delete[] DM_R[is];
        }
        delete[] DM_R;
    }
}



void Local_Orbital_Charge::allocate_dm_wfc(const int& lgd,
    Local_Orbital_wfc &lowf, 
    psi::Psi<double>* psid,
    psi::Psi<std::complex<double>>* psi)
{
    ModuleBase::TITLE("Local_Orbital_Charge", "allocate_dm_wfc");

    this->LOWF = &lowf;
    if (GlobalV::GAMMA_ONLY_LOCAL)
	{
		// here we reset the density matrix dimension.
		this->allocate_gamma(lgd, psid);
	}
	else
    {
		lowf.allocate_k(lgd, psi);
		this->allocate_DM_k();
	}
    
    return;
}


void Local_Orbital_Charge::sum_bands(LCAO_Hamilt &uhm)
{
    ModuleBase::TITLE("Local_Orbital_Charge","sum_bands");
    ModuleBase::timer::tick("Local_Orbital_Cha","sum_bands");

    GlobalC::en.eband = 0.0;

    //xiaohui add 2013-09-02
    for(int ik=0; ik<GlobalC::kv.nks; ik++)
    {
        for (int ib=0; ib<GlobalV::NBANDS; ib++)
        {
            GlobalC::en.eband += GlobalC::wf.ekb[ik][ib] * GlobalC::wf.wg(ik, ib);
        }
    } 

    //xiaohui add 2013-09-02
    if(GlobalV::GAMMA_ONLY_LOCAL)
    {
        if(GlobalV::KS_SOLVER=="selinv")
        {
            //density matrix has already been calcualted.
        }
        else if(GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="scalapack_gvx" || GlobalV::KS_SOLVER=="lapack" || GlobalV::KS_SOLVER=="cusolver")
        {
            //LiuXh modify 2021-09-06, clear memory, cal_dk_gamma() not used for genelpa solver.
            //density matrix has already been calculated.
            ModuleBase::timer::tick("LCAO_Charge","cal_dm_2d");

            //caution:wfc and dm 
            this->cal_dm(GlobalC::wf.wg,
                this->LOWF->wfc_gamma, this->dm_gamma);        // Peize Lin test 2019-01-16

            ModuleBase::timer::tick("LCAO_Charge","cal_dm_2d");

            this->cal_dk_gamma_from_2D(); // transform dm_gamma[is].c to this->DM[is]
        }
        else if(GlobalV::KS_SOLVER=="hpseps") //LiuXh add 2021-09-06, used for hpseps solver
        {
            this->cal_dk_gamma();//calculate the density matrix.
        }
    }
    else
    {
        ModuleBase::GlobalFunc::NOTE("Calculate the density matrix.");
        this->cal_dk_k( GlobalC::GridT, GlobalC::wf.wg );
        if(GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="scalapack_gvx" || GlobalV::KS_SOLVER=="lapack" || GlobalV::KS_SOLVER=="cusolver")        // Peize Lin test 2019-05-15
		{
            this->cal_dm(GlobalC::wf.wg,
                this->LOWF->wfc_k, this->dm_k);
        }
    }


    for(int is=0; is<GlobalV::NSPIN; is++)
    {
        ModuleBase::GlobalFunc::ZEROS( GlobalC::CHR.rho[is], GlobalC::rhopw->nrxx ); // mohan 2009-11-10
    }

    //------------------------------------------------------------
    //calculate the charge density on real space grid.
    //------------------------------------------------------------
     time_t start = time(NULL);

    if(GlobalV::GAMMA_ONLY_LOCAL)
    {
        Gint_inout inout(this->DM, (Charge*)(&GlobalC::CHR), Gint_Tools::job_type::rho);
        uhm.GG.cal_gint(&inout);
        if (XC_Functional::get_func_type() == 3)
        {
            for(int is=0; is<GlobalV::NSPIN; is++)
            {
                ModuleBase::GlobalFunc::ZEROS(GlobalC::CHR.kin_r[0], GlobalC::rhopw->nrxx);
            }
            Gint_inout inout1(this->DM, (Charge*)(&GlobalC::CHR), Gint_Tools::job_type::tau);
            uhm.GG.cal_gint(&inout1);
        }
    }
    else
    {
        ModuleBase::GlobalFunc::NOTE("Calculate the charge on real space grid!");
        Gint_inout inout(this->DM_R, (Charge*)(&GlobalC::CHR), Gint_Tools::job_type::rho);
        uhm.GK.cal_gint(&inout);
        if (XC_Functional::get_func_type() == 3)
        {
            for(int is=0; is<GlobalV::NSPIN; is++)
            {
                ModuleBase::GlobalFunc::ZEROS(GlobalC::CHR.kin_r[0], GlobalC::rhopw->nrxx);
            }
            Gint_inout inout1(this->DM_R, (Charge*)(&GlobalC::CHR), Gint_Tools::job_type::tau);
            uhm.GK.cal_gint(&inout1);
        }
    }

     time_t end = time(NULL);

     //GlobalV::ofs_running << " START_Charge Time : " << ctime(&time_charge_start);
     //GlobalV::ofs_running << " END_Charge  Time : " << ctime(&time_charge_end);
     //GlobalV::ofs_running << " FINAL_Charge Time : " << difftime(time_charge_end, time_charge_start) << " (Seconds)" << std::endl;

    ModuleBase::GlobalFunc::OUT_TIME("charge grid integration", start, end);

	//BLOCK_HERE("sum_bands::before renormalize rho");  

	GlobalC::CHR.renormalize_rho();

	ModuleBase::timer::tick("Local_Orbital_Cha","sum_bands");
	return;
}
