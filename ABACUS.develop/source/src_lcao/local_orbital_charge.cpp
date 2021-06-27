#include "local_orbital_charge.h"
#include "../src_pw/global.h"
#include "module_base/blas_connector.h"

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

Local_Orbital_Charge::Local_Orbital_Charge()
{
    // for gamma algorithms.
    this->init_DM = false;  
    this->lgd_now = 0;
    this->lgd_last = 0;

    // for k-dependent algorithms.
    this->init_DM_R = false;

	// whether to printout density matrix
    this->out_dm = 0;

    //xiaohui add 2014-06-19
    //band_local = nullptr;
    //Z_wg = nullptr;
    //Z_LOC = nullptr;
    sender_2D_index = nullptr;
    sender_size_process = nullptr;
    sender_displacement_process = nullptr;

    receiver_local_index = nullptr;
    receiver_size_process = nullptr;
    receiver_displacement_process = nullptr;
}

Local_Orbital_Charge::~Local_Orbital_Charge()
{
    // with gamma point only
     if (this->init_DM)
	 {
		 for (int is=0; is<NSPIN; is++)
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
        for(int is=0; is<NSPIN; is++)
        {
            delete[] DM_R[is];
        }
        delete[] DM_R;
    }
}



void Local_Orbital_Charge::allocate_dm_wfc(const Grid_Technique &gt)
{
    TITLE("Local_Orbital_Charge","allocate_dm_wfc");

	if(GAMMA_ONLY_LOCAL)
	{
		// here we reset the density matrix dimension.
		this->allocate_gamma(gt);
	}
	else
	{
		LOWF.allocate_k(gt);
		this->allocate_DM_k();
	}

	return;
}


void Local_Orbital_Charge::sum_bands(void)
{
    TITLE("Local_Orbital_Charge","sum_bands");
    timer::tick("Local_Orbital_Cha","sum_bands",'E');

    en.eband = 0.0;

    //xiaohui add 2013-09-02
    for(int ik=0; ik<kv.nks; ik++)
    {
        for (int ib=0; ib<NBANDS; ib++)
        {
            en.eband += wf.ekb[ik][ib] * wf.wg(ik, ib);
        }
    } 

    //xiaohui add 2013-09-02
    if(GAMMA_ONLY_LOCAL)
    {
        if(KS_SOLVER=="selinv")
        {
            //density matrix has already been calcualted.
        }
        else if(KS_SOLVER=="genelpa" || KS_SOLVER=="scalapack_gvx")
        {
            if(INPUT.new_dm>0)
            {
                //density matrix has already been calcualted.
                timer::tick("LCAO_Charge","cal_dm_2d",'F');
                wfc_dm_2d.cal_dm(wf.wg);        // Peize Lin test 2019-01-16
                timer::tick("LCAO_Charge","cal_dm_2d",'F');
                this->cal_dk_gamma_from_2D(); // transform dm_gamma[is].c to this->DM[is]
            }
            else
            {
                //xiaohui modify 2014-06-18
                this->cal_dk_gamma();//calculate the density matrix.
            }
        }
    }
    else
    {
        NOTE("Calculate the density matrix!");
        this->cal_dk_k( GridT );
        if(KS_SOLVER=="genelpa" || KS_SOLVER=="scalapack_gvx")        // Peize Lin test 2019-05-15
		{
            wfc_dm_2d.cal_dm(wf.wg);
		}
    }
            
    for(int is=0; is<NSPIN; is++)
    {
        ZEROS( CHR.rho[is], pw.nrxx ); // mohan 2009-11-10
    }

    //------------------------------------------------------------
    //calculate the density in real space grid.
    //------------------------------------------------------------
     time_t start = time(NULL);

    if(GAMMA_ONLY_LOCAL)
    {
        UHM.GG.cal_rho();
    }
    else
    {
        NOTE("Calculate the charge on real space grid!");
        UHM.GK.calculate_charge();
    }

     time_t end = time(NULL);

     //ofs_running << " START_Charge Time : " << ctime(&time_charge_start);
     //ofs_running << " END_Charge  Time : " << ctime(&time_charge_end);
     //ofs_running << " FINAL_Charge Time : " << difftime(time_charge_end, time_charge_start) << " (Seconds)" << endl;

    OUT_TIME("charge grid integration", start, end);

	//BLOCK_HERE("sum_bands::before renormalize rho");  

	CHR.renormalize_rho();

	timer::tick("Local_Orbital_Cha","sum_bands",'E');
	return;
}
