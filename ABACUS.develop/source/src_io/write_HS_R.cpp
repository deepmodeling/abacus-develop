#include "src_lcao/LOOP_ions.h"
#include "src_io/cal_r_overlap_R.h"
#include "src_pw/global.h"
#include "src_io/write_HS.h"


void LOOP_ions::output_HS_R(void)
{
    TITLE("LOOP_ions","output_HS_R"); 
    timer::tick("LOOP_ions","output_HS_R",'D'); 
	
	// add by jingan for out r_R matrix 2019.8.14
	if(INPUT.out_r_matrix)
	{
		cal_r_overlap_R r_matrix;
		r_matrix.init();
		r_matrix.out_r_overlap_R(NSPIN);
	}

    if(NSPIN==1||NSPIN==4)
    {
        UHM.calculate_STN_R();
        UHM.GK.cal_vlocal_R(0);
        UHM.GK.distribute_pvpR_tr();
        HS_Matrix::save_HSR_tr(0);
    }
    ///*
    else if(NSPIN==2)
    {
        UHM.calculate_STN_R();
        for(int ik=0; ik<kv.nks; ik++)
        {
            if(ik==0 || ik==kv.nks/2)
            {
                if(NSPIN==2)CURRENT_SPIN = kv.isk[ik];
                for(int ir=0; ir<pw.nrxx; ir++)
                {
                    pot.vr_eff1[ir] = pot.vr_eff( CURRENT_SPIN, ir);
                }
        	    	
                if(!GAMMA_ONLY_LOCAL)
                {
                    if(VL_IN_H)
                    {
						//UHM.GK.cal_vlocal_k(pot.vrs1,GridT);
						UHM.GK.cal_vlocal_k(pot.vr_eff1, GridT, CURRENT_SPIN);
                    }
                }
                UHM.GK.cal_vlocal_R(CURRENT_SPIN);
                UHM.GK.distribute_pvpR_tr();
                HS_Matrix::save_HSR_tr(CURRENT_SPIN);
            }
        }
    }

    if(!GAMMA_ONLY_LOCAL) //LiuXh 20181011
    {
        UHM.GK.destroy_pvpR();
    } //LiuXh 20181011

    timer::tick("LOOP_ions","output_HS_R",'D'); 
    return;
}
