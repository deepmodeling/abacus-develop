#include "FORCE_gamma.h"
#include "../src_pw/global.h"
#include "dftu.h"  //Quxin add for DFT+U on 20201029

Force_LCAO_gamma::Force_LCAO_gamma ()
{}

Force_LCAO_gamma::~Force_LCAO_gamma ()
{}

// be called in force_lo.cpp
void Force_LCAO_gamma::ftable_gamma (
	const bool isforce,
	const bool isstress,
	matrix& foverlap,
	matrix& ftvnl_dphi,
	matrix& fvnl_dbeta,	
	matrix& fvl_dphi,
	matrix& soverlap,
	matrix& stvnl_dphi,
	matrix& svnl_dbeta,
	matrix& svl_dphi)
{
    TITLE("Force_LCAO_gamma", "ftable");
    timer::tick("Force_LCAO_gamma","ftable_gamma",'F');
    
    // allocate DSloc_x, DSloc_y, DSloc_z
    // allocate DHloc_fixed_x, DHloc_fixed_y, DHloc_fixed_z
    this->allocate_gamma();

    // calculate the 'energy density matrix' here.
    this->cal_foverlap(isforce, isstress, foverlap, soverlap);

    if(INPUT.new_dm>0)
    {
        this->cal_ftvnl_dphi(LOC.wfc_dm_2d.dm_gamma, isforce, isstress, ftvnl_dphi, stvnl_dphi);
        this->cal_fvnl_dbeta(LOC.wfc_dm_2d.dm_gamma, isforce, isstress, fvnl_dbeta, svnl_dbeta);
        this->cal_fvl_dphi(LOC.wfc_dm_2d.dm_gamma, isforce, isstress, fvl_dphi, svl_dphi);

        //quxin added for DFT+U
        if(INPUT.dft_plus_u) 
	    {
	    	dftu.force_stress();
	    	dftu.erase_force_stress();
	    }
    }
    else
    {
        timer::tick("Force_LCAO_gamma","cal_dm_grid",'G');
        // calculate the 'density matrix' here.
        matrix dm2d;
		dm2d.create(NSPIN, ParaO.nloc);
        Memory::record ("Force_LCAO_gamma", "dm2d", ParaO.nloc*NSPIN, "double");    

        bool with_energy = false;

		this->set_EDM_gamma(dm2d, with_energy);

        timer::tick("Force_LCAO_gamma","cal_dm_grid",'G');

        this->cal_ftvnl_dphi(dm2d, isforce, isstress, ftvnl_dphi, stvnl_dphi);
        this->cal_fvnl_dbeta(dm2d, isforce, isstress, fvnl_dbeta, svnl_dbeta);

        //quxin added for DFT+U
        if(INPUT.dft_plus_u) 
	    {
	    	dftu.force_stress();
	    	dftu.erase_force_stress();
	    }

        // calculate < dphi | V | phi > on real space grid.
        this->cal_fvl_dphi(dm2d, isforce, isstress, fvl_dphi, svl_dphi);

    }
	if(isforce)
	{
        Parallel_Reduce::reduce_double_pool( foverlap.c, foverlap.nr * foverlap.nc);
        Parallel_Reduce::reduce_double_pool( ftvnl_dphi.c, ftvnl_dphi.nr * ftvnl_dphi.nc);
        Parallel_Reduce::reduce_double_pool( fvnl_dbeta.c, fvnl_dbeta.nr * fvnl_dbeta.nc);
        Parallel_Reduce::reduce_double_pool( fvl_dphi.c, fvl_dphi.nr * fvl_dphi.nc);
	}
    if(isstress)
    {
        Parallel_Reduce::reduce_double_pool( soverlap.c, soverlap.nr * soverlap.nc);
        Parallel_Reduce::reduce_double_pool( stvnl_dphi.c, stvnl_dphi.nr * stvnl_dphi.nc);
        Parallel_Reduce::reduce_double_pool( svnl_dbeta.c, svnl_dbeta.nr * svnl_dbeta.nc);
        Parallel_Reduce::reduce_double_pool( svl_dphi.c, svl_dphi.nr * svl_dphi.nc);
    }

    // delete DSloc_x, DSloc_y, DSloc_z
    // delete DHloc_fixed_x, DHloc_fixed_y, DHloc_fixed_z
    this->finish_ftable_gamma();

    timer::tick("Force_LCAO_gamma","ftable_gamma",'F');
    return;
}

void Force_LCAO_gamma::allocate_gamma(void)
{
    TITLE("Force_LCAO_gamma","allocate_gamma");
    timer::tick("Force_LCAO_gamma","allocate_gamma",'G');

    // need to calculate the derivative in build_ST_new
    bool cal_deri = true;

    //calculate dS in LCAO
    //liaochen add on 2010/7/12
    //save the results in dense matrix by now.
    //ParaO.nloc: number of H elements in this proc.
    LM.DSloc_x = new double [ParaO.nloc];
    LM.DSloc_y = new double [ParaO.nloc];
    LM.DSloc_z = new double [ParaO.nloc];
    ZEROS(LM.DSloc_x, ParaO.nloc);
    ZEROS(LM.DSloc_y, ParaO.nloc);
    ZEROS(LM.DSloc_z, ParaO.nloc);
    //allocate stress part in gamma_only-line, added by zhengdy-stress
    if(STRESS)
    {
        LM.DSloc_11 = new double [ParaO.nloc];
        LM.DSloc_12 = new double [ParaO.nloc];
        LM.DSloc_13 = new double [ParaO.nloc];
        LM.DSloc_22 = new double [ParaO.nloc];
        LM.DSloc_23 = new double [ParaO.nloc];
        LM.DSloc_33 = new double [ParaO.nloc];
        ZEROS(LM.DSloc_11, ParaO.nloc);
        ZEROS(LM.DSloc_12, ParaO.nloc);
        ZEROS(LM.DSloc_13, ParaO.nloc);
        ZEROS(LM.DSloc_22, ParaO.nloc);
        ZEROS(LM.DSloc_23, ParaO.nloc);
        ZEROS(LM.DSloc_33, ParaO.nloc);
        LM.DHloc_fixed_11 = new double [ParaO.nloc];
        LM.DHloc_fixed_12 = new double [ParaO.nloc];
        LM.DHloc_fixed_13 = new double [ParaO.nloc];
        LM.DHloc_fixed_22 = new double [ParaO.nloc];
        LM.DHloc_fixed_23 = new double [ParaO.nloc];
        LM.DHloc_fixed_33 = new double [ParaO.nloc];
        ZEROS (LM.DHloc_fixed_11, ParaO.nloc);
        ZEROS (LM.DHloc_fixed_12, ParaO.nloc);
        ZEROS (LM.DHloc_fixed_13, ParaO.nloc);
        ZEROS (LM.DHloc_fixed_22, ParaO.nloc);
        ZEROS (LM.DHloc_fixed_23, ParaO.nloc);
        ZEROS (LM.DHloc_fixed_33, ParaO.nloc);
    }
    //calculate dS in LCAO basis
    // tips: build_ST_new --> ParaO.set_force 
    //timer::tick("Force_LCAO_gamma","build_S_new",'H');
    UHM.genH.build_ST_new ('S', cal_deri);
    //timer::tick("Force_LCAO_gamma","build_S_new",'H');

    Memory::record("force_lo", "dS", ParaO.nloc*3, "double");

    //calculate dT in LCAP
    //allocation dt
    //liaochen add on 2010/7/12
    LM.DHloc_fixed_x = new double [ParaO.nloc];
    LM.DHloc_fixed_y = new double [ParaO.nloc];
    LM.DHloc_fixed_z = new double [ParaO.nloc];
    ZEROS (LM.DHloc_fixed_x, ParaO.nloc);
    ZEROS (LM.DHloc_fixed_y, ParaO.nloc);
    ZEROS (LM.DHloc_fixed_z, ParaO.nloc);
    
    //calculate dT
    //calculate T + VNL(P1) in LCAO basis
    //timer::tick("Force_LCAO_gamma","build_T_new",'H');
    UHM.genH.build_ST_new ('T', cal_deri);
    //timer::tick("Force_LCAO_gamma","build_T_new",'H');
    //test_gamma(LM.DHloc_fixed_x, "dHloc_fixed_x T part");
    
    //UHM.genH.build_Nonlocal_beta (cal_deri);
    //timer::tick("Force_LCAO_gamma","build_Nonlocal_mu",'H');
    UHM.genH.build_Nonlocal_mu (cal_deri);
    //timer::tick("Force_LCAO_gamma","build_Nonlocal_mu",'H');
    //test_gamma(LM.DHloc_fixed_x, "dHloc_fixed_x Vnl part");

    Memory::record("force_lo", "dTVNL", ParaO.nloc*3, "double");

    timer::tick("Force_LCAO_gamma","allocate_gamma",'G');
    return;
}

void Force_LCAO_gamma::finish_ftable_gamma(void)
{
    delete [] LM.DSloc_x;
    delete [] LM.DSloc_y;
    delete [] LM.DSloc_z;
    delete [] LM.DHloc_fixed_x;
    delete [] LM.DHloc_fixed_y;
    delete [] LM.DHloc_fixed_z;
    if(STRESS)//added by zhengdy-stress
    {
        delete [] LM.DSloc_11;
        delete [] LM.DSloc_12;
        delete [] LM.DSloc_13;
        delete [] LM.DHloc_fixed_11;
        delete [] LM.DHloc_fixed_12;
        delete [] LM.DHloc_fixed_13;
        delete [] LM.DSloc_22;
        delete [] LM.DSloc_23;
        delete [] LM.DSloc_33;
        delete [] LM.DHloc_fixed_22;
        delete [] LM.DHloc_fixed_23;
        delete [] LM.DHloc_fixed_33;
    }
    return;
}


void Force_LCAO_gamma::test_gamma(double* mm, const string &name)
{
    cout << "\n PRINT " << name << endl;
    cout << setprecision(6) << endl;
    for(int i=0; i<NLOCAL; i++)
    {
        for(int j=0; j<NLOCAL; j++)
        {
            if( abs(mm[i*NLOCAL+j])>1.0e-5)
            {
                cout << setw(12) << mm[i*NLOCAL+j];
            }
            else
            {
                cout << setw(12) << "0";
            }
        }
        cout << endl;
    }
    return;
}
