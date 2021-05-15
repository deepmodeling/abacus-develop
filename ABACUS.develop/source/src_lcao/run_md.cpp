#include "run_md.h"
#include "src_pw/global.h"
#include "src_parallel/parallel_orbitals.h"
#include "src_pdiag/pdiag_double.h"
#include "LCAO_nnr.h"
#include "FORCE_STRESS.h"
#include "src_global/global_function.h"
#include "src_io/write_HS.h"
#include "src_io/cal_r_overlap_R.h"
#include "src_ions/variable_cell.h" // mohan add 2021-02-01
#include "src_ri/exx_abfs.h"
#include "src_ri/exx_opt_orb.h"
#include "ELEC_scf.h"
#include "src_global/sltk_atom_arrange.h"
#include "src_pw/vdwd2.h"

Run_MD::Run_MD()
{}

Run_MD::~Run_MD() 
{}


void Run_MD::opt_cell(void)
{
	TITLE("Run_MD","opt_cell");


	

    // Initialize the local wave functions.
    // npwx, eigenvalues, and weights
    // npwx may change according to cell change
    // this function belongs to cell LOOP
    wf.allocate_ekb_wg(kv.nks);

    // Initialize the FFT.
    // this function belongs to cell LOOP
    UFFT.allocate();

    // output is ppcell.vloc 3D local pseudopotentials
    // this function belongs to cell LOOP
    ppcell.init_vloc(pw.nggm, ppcell.vloc);

    // Initialize the sum of all local potentials.
    // if ion_step==0, read in/initialize the potentials
    // this function belongs to ions LOOP
    int ion_step=0;
    pot.init_pot(ion_step, pw.strucFac);

	
	opt_ions();
	return;
}


void Run_MD::opt_ions(void)
{
    TITLE("Run_MD","opt_ions"); 
    timer::tick("Run_MD","opt_ions",'B'); 
		
    if(OUT_LEVEL=="i")
    {
        cout << setprecision(12);
        cout<< " " << setw(7)<< "ISTEP"
        <<setw(5)<< "NE"
        <<setw(18)<< "ETOT(eV)"
        <<setw(10)<< "dE(meV)"
        <<setw(10)<< "F(eV/A)"
        <<setw(10)<< "T(MIN)"
        <<endl;
    }

    // Geometry optimization algorithm setup.
    if(FORCE)
    {
        //Ions_Move_Methods 
        IMM.allocate();
        //Charge_Extrapolation
        CE.allocate_ions();
    }

	// pengfei Li 2018-05-14
    if(STRESS)
    {
		// allocate arrays related to changes of lattice vectors
        LCM.allocate();
    }

    MD_basic mdb(INPUT.mdp, ucell);
    int mdtype = INPUT.mdp.mdtype;

    this->istep = 1;
    int force_step = 1;
    int stress_step = 1;
    bool stop = false;
    while(istep <= NSTEP && !stop)
    {
        time_t estart = time(NULL);
	
		// xiaohui add "m" option, 2015-09-16
        if(OUT_LEVEL=="ie" || OUT_LEVEL=="m")
        {
			cout << " ---------------------------------------------------------" << endl;
			cout<<" Molecular Dynamics STEP "<< mdb.getRealStep()<<endl;
			cout << " ---------------------------------------------------------" << endl;
        }

		// solve electronic structures in terms of LCAO
		// mohan add 2021-02-09
		LOE.solve_elec_stru(this->istep);

		time_t eend = time(NULL);

        //xiaohui add 2014-07-07, for second-order extrapolation
		CE.update_all_pos(ucell);

		if(mdtype==1||mdtype==2)   
		{
			mdb.runNVT(istep);
		}
		else if(mdtype==0)  
		{
			mdb.runNVE(istep);
		}
        else if(mdtype==-1)
        {
            stop = mdb.runFIRE(istep);
        }
        else
        {
            WARNING_QUIT("opt_ions", "mdtype should be -1~2!");
        }

        if(pot.out_potential == 2)
        {
            stringstream ssp;
            stringstream ssp_ave;
            ssp << global_out_dir << "ElecStaticPot";
            ssp_ave << global_out_dir << "ElecStaticPot_AVE";
            pot.write_elecstat_pot(ssp.str(), ssp_ave.str()); //output 'Hartree + local pseudopot'
        }

        time_t fstart = time(NULL);
        time_t fend = time(NULL);

        //xiaohui add 2014-07-07, for second-order extrapolation
		CE.save_pos_next(ucell);

		//xiaohui add CE.istep = istep 2014-07-07
		CE.update_istep(istep);

		// charge extrapolation if istep>0.
		CE.extrapolate_charge();

		if(pot.extra_pot=="dm")//xiaohui modify 2015-02-01
		{
			// done after grid technique.
		}
		else
		{
			pot.init_pot( istep, pw.strucFac );
		}


        if(OUT_LEVEL=="i")
        {
            double etime_min = difftime(eend, estart)/60.0;
            double ftime_min = difftime(fend, fstart)/60.0;
            stringstream ss;
            ss << MOVE_IONS << istep;

            cout << setiosflags(ios::scientific)
            << " " << setw(7) << ss.str()
            << setw(5) << ELEC_scf::iter
            << setw(18) << setprecision(6) << en.etot * Ry_to_eV;

            cout << setprecision(2) << setiosflags(ios::scientific)
            << setw(10) << IMM.get_ediff() * Ry_to_eV * 1000
            << setw(10) << IMM.get_largest_grad() * Ry_to_eV / BOHR_TO_A;

            cout << resetiosflags(ios::scientific)
            << setprecision(2) << setw(10) << etime_min + ftime_min;
            cout << endl;
        }

        ++istep;
    }

    if(istep>1) final_scf();

	// mohan update 2021-02-10
    hm.orb_con.clear_after_ions(UOT, ORB, INPUT.out_descriptor);

    timer::tick("Run_MD","opt_ions",'B'); 
    return;
}

void Run_MD::final_scf(void)
{
    TITLE("Run_MD","final_scf"); 

    FINAL_SCF = true;
    Variable_Cell::final_calculation_after_vc();
    atom_arrange::set_sr_NL();
    atom_arrange::search( SEARCH_RADIUS );
    GridT.set_pbc_grid(
        pw.ncx, pw.ncy, pw.ncz,
        pw.bx, pw.by, pw.bz,
        pw.nbx, pw.nby, pw.nbz,
        pw.nbxx, pw.nbzp_start, pw.nbzp);

    // (2) If k point is used here, allocate HlocR after atom_arrange.
    if(!GAMMA_ONLY_LOCAL)
    {
        // For each atom, calculate the adjacent atoms in different cells 
        // and allocate the space for H(R) and S(R).
        LNNR.cal_nnr();
        LM.allocate_HS_R(LNNR.nnr);
        
		// need to first calculae lgd.
        // using GridT.init.
        LNNR.cal_nnrg(GridT);
    }

    // (4) set the augmented orbitals index.
    // after ParaO and GridT, 
    // this information is used to calculate
    // the force.
    LOWF.set_trace_aug(GridT);
		
    // (5) init density kernel
    // (6) init wave functions.
    if(GAMMA_ONLY_LOCAL)
    {
        // here we reset the density matrix dimension.
        LOC.allocate_gamma(GridT);
    }
    else
    {
        LOWF.allocate_k(GridT);
        LOC.allocate_DM_k();
    }

    UHM.set_lcao_matrices();

    if(vdwd2_para.flag_vdwd2)							//Peize Lin add 2014-04-04, update 2021-03-09
    {
        Vdwd2 vdwd2(ucell,vdwd2_para);
        vdwd2.cal_energy();
        en.evdw = vdwd2.get_energy();
    }
	else if(vdwd3_para.flag_vdwd3)							//jiyy add 2019-05-18, update 2021-05-02
    {
        Vdwd3 vdwd3(ucell,vdwd3_para);
        vdwd3.cal_energy();
        en.evdw = vdwd3.get_energy();
    }												  
    
	ELEC_scf es;
	es.scf(0);

    ofs_running << "\n\n --------------------------------------------" << endl;
    ofs_running << setprecision(16);
    ofs_running << " !FINAL_ETOT_IS " << en.etot * Ry_to_eV << " eV" << endl; 
    ofs_running << " --------------------------------------------\n\n" << endl;

    return;
}
