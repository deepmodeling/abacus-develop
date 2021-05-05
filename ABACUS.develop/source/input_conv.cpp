#include "src_pw/global.h"
#include "src_pw/tools.h"
#include "input.h"
#include "input_conv.h"
#include "src_ions/ions_move_basic.h"
#include "src_io/optical.h"
#include "src_lcao/FORCE_STRESS.h"
#include "src_lcao/local_orbital_charge.h"
#include "src_lcao/ORB_read.h"
#include "src_pw/efield.h"
#include "src_lcao/global_fp.h"
#include "src_pw/vdwd2_parameters.h"
#include "src_pw/vdwd3.h"
#include "src_io/chi0_hilbert.h"
#include "src_io/chi0_standard.h"
#include "src_io/epsilon0_pwscf.h"
#include "src_io/epsilon0_vasp.h"
#include "src_pw/unitcell.h"
#include "src_ri/exx_abfs-jle.h"
#include "src_pw/occupy.h"
#include "src_io/berryphase.h"
#include "src_pw/symmetry.h"

void Input_Conv::Convert(void)
{
    TITLE("Input_Conv","Convert");
	timer::tick("Input_Conv","Convert",'B');
//----------------------------------------------------------
// main parameters / electrons / spin ( 10/16 )
//----------------------------------------------------------
//  suffix
    if(INPUT.atom_file!="") global_atom_card = INPUT.atom_file;
	global_wannier_card = INPUT.wannier_card;
    if(INPUT.kpoint_file!= "") global_kpoint_card = INPUT.kpoint_file;
    if(INPUT.pseudo_dir != "") global_pseudo_dir = INPUT.pseudo_dir + "/";
    global_pseudo_type = INPUT.pseudo_type;
	ucell.latName = INPUT.latname; 
	ucell.ntype = INPUT.ntype;
	ucell.nelec = INPUT.nelec;
//  ucell.lmaxmax = INPUT.lmaxmax;

    NBANDS = INPUT.nbands;
    NBANDS_ISTATE = INPUT.nbands_istate;
	NPOOL = INPUT.npool;
	CALCULATION = INPUT.calculation;

	// qianrui add 2021-2-5
	STO_WF.nchi = INPUT.nbands_sto;
	STO_WF.nche_sto = INPUT.nche_sto;

	// Electrical Field
	EFIELD = INPUT.efield;
	Efield::edir = INPUT.edir;
	Efield::emaxpos = INPUT.emaxpos;
	Efield::eopreg = INPUT.eopreg;
	Efield::eamp = INPUT.eamp;

	// optical
	Optical::opt_epsilon2=INPUT.opt_epsilon2;	// mohan add 2010-03-24	
	Optical::opt_nbands=INPUT.opt_nbands;		// number of bands for optical transition.	

	DFT_FUNCTIONAL = INPUT.dft_functional;
    NSPIN = INPUT.nspin;
    CURRENT_SPIN = 0;

    FORCE = INPUT.force;
    FORCE_THR = INPUT.force_thr;

    STRESS_THR = INPUT.stress_thr;
    PRESS1 = INPUT.press1;
    PRESS2 = INPUT.press2;
    PRESS3 = INPUT.press3;
	Force_Stress_LCAO::force_invalid_threshold_ev = INPUT.force_thr_ev2;

	BFGS_Basic::w1 = INPUT.bfgs_w1;
	BFGS_Basic::w2 = INPUT.bfgs_w2;

	Ions_Move_Basic::trust_radius_max = INPUT.trust_radius_max;
	Ions_Move_Basic::trust_radius_min = INPUT.trust_radius_min;
	Ions_Move_Basic::trust_radius_ini = INPUT.trust_radius_ini;
	Ions_Move_Basic::out_stru = INPUT.out_stru; //mohan add 2012-03-23

	STRESS = INPUT.stress;


	// pengfei Li add 2018-11-11
	if(INPUT.fixed_axes == "None")
	{
		ucell.lc[0] = 1; ucell.lc[1] = 1; ucell.lc[2] = 1;
	}
	else if(INPUT.fixed_axes == "volume")
	{
		ucell.lc[0] = 1; ucell.lc[1] = 1; ucell.lc[2] = 1;
	}
	else if(INPUT.fixed_axes == "a")
	{
		ucell.lc[0] = 0; ucell.lc[1] = 1; ucell.lc[2] = 1;
	}
	else if(INPUT.fixed_axes == "b")
	{
		ucell.lc[0] = 1; ucell.lc[1] = 0; ucell.lc[2] = 1;
	}
	else if(INPUT.fixed_axes == "c")
	{
		ucell.lc[0] = 1; ucell.lc[1] = 1; ucell.lc[2] = 0;
	}
	else if(INPUT.fixed_axes == "ab")
	{
		ucell.lc[0] = 0; ucell.lc[1] = 0; ucell.lc[2] = 1;
	}
	else if(INPUT.fixed_axes == "ac")
	{
		ucell.lc[0] = 0; ucell.lc[1] = 1; ucell.lc[2] = 0;
	}
	else if(INPUT.fixed_axes == "bc")
	{
		ucell.lc[0] = 1; ucell.lc[1] = 0; ucell.lc[2] = 0;
	}
	else if(INPUT.fixed_axes == "abc")
	{
		ucell.lc[0] = 0; ucell.lc[1] = 0; ucell.lc[2] = 0;
	}
	else
	{
		WARNING_QUIT("Input", "fixed_axes should be None,a,b,c,ab,ac,bc or abc!");
	}

	MOVE_IONS = INPUT.ion_dynamics;
	OUT_LEVEL = INPUT.out_level;
	Ions_Move_CG::CG_THRESHOLD = INPUT.cg_threshold; // pengfei add 2013-09-09

    Symmetry::symm_flag = INPUT.symmetry;						// 9
	BASIS_TYPE = INPUT.basis_type;
	KS_SOLVER = INPUT.ks_solver;
	SEARCH_RADIUS = INPUT.search_radius;
	SEARCH_PBC = INPUT.search_pbc;

//----------------------------------------------------------
// planewave (8/8)
//----------------------------------------------------------
    pw.set(
        INPUT.gamma_only,
        INPUT.ecutwfc,
        INPUT.ecutrho,
        INPUT.nx,
        INPUT.ny,
        INPUT.nz,
        INPUT.ncx,
        INPUT.ncy,
        INPUT.ncz,
		INPUT.bx,
		INPUT.by,
		INPUT.bz
    );
	GAMMA_ONLY_LOCAL = INPUT.gamma_only_local;

//----------------------------------------------------------
// diagonalization  (5/5)
//----------------------------------------------------------
	DIAGO_PROC = INPUT.diago_proc;
    DIAGO_CG_MAXITER = INPUT.diago_cg_maxiter;
	DIAGO_CG_PREC = INPUT.diago_cg_prec;
    DIAGO_DAVID_NDIM = INPUT.diago_david_ndim;
    ETHR = INPUT.ethr;
	NB2D = INPUT.nb2d;
	NURSE = INPUT.nurse;
	COLOUR = INPUT.colour;
	T_IN_H = INPUT.t_in_h;
	VL_IN_H = INPUT.vl_in_h;
	VNL_IN_H = INPUT.vnl_in_h;
	TEST_FORCE = INPUT.test_force;
	TEST_STRESS = INPUT.test_stress;

//----------------------------------------------------------
// iteration (1/3)
//----------------------------------------------------------
    DRHO2 = INPUT.dr2;

//----------------------------------------------------------
// wavefunction / charge / potential / (2/4)
//----------------------------------------------------------
    RESTART_MODE = INPUT.restart_mode;
    wf.start_wfc = INPUT.start_wfc;
	wf.mem_saver = INPUT.mem_saver; //mohan add 2010-09-07
	en.printe    = INPUT.printe; // mohan add 2011-03-16

//----------------------------------------------------------
// about vdw, jiyy add 2019-08-04
//----------------------------------------------------------	
	if(INPUT.vdw_method=="d2")
	{
		vdwd2_para.flag_vdwd2 = true;
		vdwd2_para.scaling = std::stod(INPUT.vdw_s6);
		vdwd2_para.damping = INPUT.vdw_d;
		vdwd2_para.C6_input(INPUT.vdw_C6_file, INPUT.vdw_C6_unit);
		vdwd2_para.R0_input(INPUT.vdw_R0_file, INPUT.vdw_R0_unit);
		vdwd2_para.model = INPUT.vdw_model;
		if(INPUT.vdw_model=="radius")
		{
			if(INPUT.vdw_radius_unit=="Bohr")
			{
				vdwd2_para.radius = std::stod(INPUT.vdw_radius);
			}
			else
			{
				vdwd2_para.radius = std::stod(INPUT.vdw_radius) * BOHR_TO_A;
			}
		}
		else if(INPUT.vdw_model=="period")
		{
			vdwd2_para.period = INPUT.vdw_period;
		}
	}
    if(INPUT.vdw_method=="d3_0" || INPUT.vdw_method=="d3_bj")
    {
		vdwd3.vdwD3 = true;
		vdwd3.s6 = std::stod(INPUT.vdw_s6);
		vdwd3.s18 = std::stod(INPUT.vdw_s8);
		vdwd3.rs6 = std::stod(INPUT.vdw_a1);
		vdwd3.rs18 = std::stod(INPUT.vdw_a2);					
		vdwd3.abc = INPUT.vdw_abc;
		if(INPUT.vdw_method=="d3_0")
		{
			vdwd3.version = 1;
		}
		if(INPUT.vdw_method=="d3_bj")
		{
			vdwd3.version = 2;
		}
		vdwd3.model = INPUT.vdw_model;
		if(INPUT.vdw_model=="radius")
	    {
			if(INPUT.vdw_radius_unit=="Bohr")
			{
			    vdwd3.rthr2 = pow(std::stod(INPUT.vdw_radius),2);
			}
			else
			{
				vdwd3.rthr2 = pow((std::stod(INPUT.vdw_radius) * BOHR_TO_A),2);       
			}
		    if(INPUT.vdw_cn_thr_unit=="Bohr")
			{
			    vdwd3.cn_thr2 = pow(INPUT.vdw_cn_thr,2);
			}
			else
			{  
				vdwd3.cn_thr2 = pow((INPUT.vdw_cn_thr * BOHR_TO_A),2);			
			}
		}
		else if(INPUT.vdw_model=="period")
		{
			vdwd3.rep_vdw[0] = INPUT.vdw_period.x;
			vdwd3.rep_vdw[1] = INPUT.vdw_period.y;
			vdwd3.rep_vdw[2] = INPUT.vdw_period.z;
		}
	}
    
//----------------------------------------------------------
// about spectrum, pengfei 2016-12-14
//----------------------------------------------------------

	if( (INPUT.spectral_type == "eels" && INPUT.eels_method == 0) 
		|| (INPUT.spectral_type == "None" 
		&& INPUT.eels_method == 0 
		&& INPUT.kmesh_interpolation) )
	{
		if(INPUT.spectral_type == "eels")
		{
			chi0_hilbert.epsilon = true;
		}
		else if(INPUT.spectral_type == "None")
		{
			chi0_hilbert.epsilon = false;
		}
		//chi0_hilbert.epsilon = INPUT.epsilon;
		chi0_hilbert.kernel_type = INPUT.kernel_type;
		chi0_hilbert.system = INPUT.system;
		chi0_hilbert.eta = INPUT.eta;
		chi0_hilbert.domega = INPUT.domega;
		chi0_hilbert.nomega = INPUT.nomega;
		chi0_hilbert.dim = INPUT.ecut_chi; 
		//chi0_hilbert.oband = INPUT.oband;

		chi0_hilbert.q_start[0] = INPUT.q_start[0];  
		chi0_hilbert.q_start[1] = INPUT.q_start[1]; 
		chi0_hilbert.q_start[2] = INPUT.q_start[2];

		chi0_hilbert.direct[0] = INPUT.q_direct[0]; 
		chi0_hilbert.direct[1] = INPUT.q_direct[1]; 
		chi0_hilbert.direct[2] = INPUT.q_direct[2];

		//chi0_hilbert.start_q = INPUT.start_q;
		//chi0_hilbert.interval_q = INPUT.interval_q;
		chi0_hilbert.nq = INPUT.nq;
		chi0_hilbert.out_epsilon = INPUT.out_epsilon;
		chi0_hilbert.out_chi = INPUT.out_chi;
		chi0_hilbert.out_chi0 = INPUT.out_chi0;
		chi0_hilbert.fermi_level = INPUT.fermi_level;
		chi0_hilbert.coulomb_cutoff = INPUT.coulomb_cutoff;
		chi0_hilbert.kmesh_interpolation = INPUT.kmesh_interpolation;
		for(int i=0; i<100; i++)
		{
			chi0_hilbert.qcar[i][0] = INPUT.qcar[i][0]; 
			chi0_hilbert.qcar[i][1] = INPUT.qcar[i][1]; 
			chi0_hilbert.qcar[i][2] = INPUT.qcar[i][2]; 
		}
		chi0_hilbert.lcao_box[0] = INPUT.lcao_box[0]; 
		chi0_hilbert.lcao_box[1] = INPUT.lcao_box[1]; 
		chi0_hilbert.lcao_box[2] = INPUT.lcao_box[2];
	}
	
	//if( INPUT.epsilon && (INPUT.epsilon_choice == 1))
	if( INPUT.spectral_type == "eels" && INPUT.eels_method == 1)
	{
		//chi0_standard.epsilon = INPUT.epsilon;
		chi0_standard.epsilon = true;
		chi0_standard.system = INPUT.system;
		chi0_standard.eta = INPUT.eta;
		chi0_standard.domega = INPUT.domega;
		chi0_standard.nomega = INPUT.nomega;
		chi0_standard.dim = INPUT.ecut_chi;
		//chi0_standard.oband = INPUT.oband;
		chi0_standard.q_start[0] = INPUT.q_start[0]; 
	 	chi0_standard.q_start[1] = INPUT.q_start[1]; 
		chi0_standard.q_start[2] = INPUT.q_start[2];
		chi0_standard.direct[0] = INPUT.q_direct[0];  
		chi0_standard.direct[1] = INPUT.q_direct[1]; 
		chi0_standard.direct[2] = INPUT.q_direct[2];
		//chi0_standard.start_q = INPUT.start_q;
		//chi0_standard.interval_q = INPUT.interval_q;
		chi0_standard.nq = INPUT.nq;
		chi0_standard.out_epsilon = INPUT.out_epsilon;		
	}
	
	//if( INPUT.epsilon0 && (INPUT.epsilon0_choice == 1) )
	if( INPUT.spectral_type == "absorption" && INPUT.absorption_method == 1)
	{
		//epsilon0_pwscf.epsilon = INPUT.epsilon0;
		epsilon0_pwscf.epsilon = true;
		epsilon0_pwscf.intersmear = INPUT.eta;
		epsilon0_pwscf.intrasmear = INPUT.intrasmear;
		epsilon0_pwscf.domega = INPUT.domega;
		epsilon0_pwscf.nomega = INPUT.nomega;
		epsilon0_pwscf.shift = INPUT.shift;
		epsilon0_pwscf.metalcalc = INPUT.metalcalc;
		epsilon0_pwscf.degauss = INPUT.eps_degauss;
	}
	
	//if( INPUT.epsilon0 && (INPUT.epsilon0_choice == 0))
	if( INPUT.spectral_type == "absorption" && INPUT.absorption_method == 0)
	{
		//epsilon0_vasp.epsilon = INPUT.epsilon0;
		epsilon0_vasp.epsilon = true;
		epsilon0_vasp.domega = INPUT.domega;
		epsilon0_vasp.nomega = INPUT.nomega;
		epsilon0_vasp.eta = INPUT.eta;
	}

//--------------------------------------------
// added by zhengdy-soc
//--------------------------------------------
	if(INPUT.noncolin||INPUT.lspinorb) 
	{
		NSPIN = 4;
	}
	if(NSPIN == 4)
	{
		NONCOLIN = INPUT.noncolin;
		//wavefunctions are spinors with 2 components
		NPOL = 2;
		//set the domag variable to make a spin-orbit calculation with zero magnetization
		if(NONCOLIN)
		{
			DOMAG = true;
			DOMAG_Z = false;
		}
		else{
			DOMAG = false;
			DOMAG_Z = true;
		}
		LSPINORB = INPUT.lspinorb;

		delete[] soc.m_loc;
		delete[] soc.angle1;
		delete[] soc.angle2;
		soc.m_loc = new Vector3<double> [INPUT.ntype];
		soc.angle1 = new double[INPUT.ntype];
		soc.angle2 = new double[INPUT.ntype];
		for(int i = 0;i<INPUT.ntype;i++)
		{
			soc.angle1[i] = INPUT.angle1[i]/180*PI;
			soc.angle2[i] = INPUT.angle2[i]/180*PI;
		}
#ifdef __MPI
//			Parallel_Common::bcast_double(soc.angle1[i]);
//			Parallel_Common::bcast_double(soc.angle2[i]);
#endif
	}
	else{
		LSPINORB = false;
		NONCOLIN = false;
		DOMAG = false;
		DOMAG_Z = false;
		NPOL = 1;
	}
	
//----------------------------------------------------------
// Fuxiang He add 2016-10-26
//----------------------------------------------------------
	tddft = INPUT.tddft;
	td_dr2 = INPUT.td_dr2;
	td_dt = INPUT.td_dt;
	td_force_dt = INPUT.td_force_dt;
	val_elec_01 = INPUT.val_elec_01;
	val_elec_02 = INPUT.val_elec_02;
	val_elec_03 = INPUT.val_elec_03;
	vext = INPUT.vext;
	vext_dire = INPUT.vext_dire;	


	// jiyy add 2020.10.11	
	ocp = INPUT.ocp;
     //ocp_n = INPUT.ocp_n;
    ocp_set = INPUT.ocp_set;
    if(ocp == 1)
	{
		int count = 0;
		string pattern("([0-9]+\\*[0-9.]+|[0-9,.]+)");
		vector<string> str;
		string::size_type pos1, pos2;
		string c = " ";
		pos2 = ocp_set.find(c);
		pos1 = 0;
		while(string::npos != pos2)
		{
			str.push_back(ocp_set.substr(pos1, pos2-pos1));
 
			pos1 = pos2 + c.size();
			pos2 = ocp_set.find(c, pos1);
		}
		if(pos1 != ocp_set.length())
		{
			str.push_back(ocp_set.substr(pos1));
		}

		regex_t reg;
		regcomp(&reg, pattern.c_str(), REG_EXTENDED);
		regmatch_t pmatch[1];
		const size_t nmatch=1;
		for(int i=0; i<str.size(); ++i)
		{
			if(str[i] == "") 
			{
				continue;
			}
			int status = regexec(&reg, str[i].c_str(), nmatch, pmatch, 0);
			string sub_str = "";
			for(int j=pmatch[0].rm_so; j!=pmatch[0].rm_eo; ++j)
			{
				sub_str += str[i][j];
			}
			string sub_pattern("\\*");
			regex_t sub_reg;
			regcomp(&sub_reg, sub_pattern.c_str(), REG_EXTENDED);
			regmatch_t sub_pmatch[1];
			const size_t sub_nmatch=1;
			if(regexec(&sub_reg, sub_str.c_str(), sub_nmatch, sub_pmatch, 0) == 0)
			{
				int pos = sub_str.find("*");
				int num = stoi(sub_str.substr(0, pos));
				double occ = stof(sub_str.substr(pos+1, sub_str.size()));
				vector<double> ocp_temp(num, occ);
				const vector<double>::iterator dest = ocp_kb.begin()+count;
				copy(ocp_temp.begin(), ocp_temp.end(), dest);
				count += num;
			}
			else
			{
				ocp_kb[count] = stof(sub_str);
				count += 1;
			}
		}
	}
		
    mulliken = INPUT. mulliken;//qifeng add 2019/9/10

//----------------------------------------------------------
// about restart, // Peize Lin add 2020-04-04
//----------------------------------------------------------	
	if(INPUT.restart_save)
	{
		restart.folder = global_out_dir + "restart/";
		const string command0 =  "test -d " + restart.folder + " || mkdir " + restart.folder;
		if(MY_RANK==0)
			system( command0.c_str() );
		if(INPUT.exx_hybrid_type=="no")
		{
			restart.info_save.save_charge = true;
		}
		else
		{
			restart.info_save.save_charge = true;
			restart.info_save.save_H = true;
		}
	}
	if(INPUT.restart_load)
	{
		restart.folder = global_out_dir + "restart/";
		if(INPUT.exx_hybrid_type=="no")
		{
			restart.info_load.load_charge = true;
		}
		else
		{
			restart.info_load.load_charge = true;
			restart.info_load.load_H = true;
		}
	}

//----------------------------------------------------------
// about exx, Peize Lin add 2018-06-20
//----------------------------------------------------------	
	if(INPUT.exx_hybrid_type=="no")
	{
		exx_global.info.hybrid_type = Exx_Global::Hybrid_Type::No;
	}
	else
	{
		if(INPUT.exx_hybrid_type=="hf")
		{
			exx_global.info.hybrid_type = Exx_Global::Hybrid_Type::HF;
		}
		else if(INPUT.exx_hybrid_type=="pbe0")
		{
			exx_global.info.hybrid_type = Exx_Global::Hybrid_Type::PBE0;
		}
		else if(INPUT.exx_hybrid_type=="hse")
		{
			exx_global.info.hybrid_type = Exx_Global::Hybrid_Type::HSE;
		}
		else if(INPUT.exx_hybrid_type=="opt_orb")
		{
			exx_global.info.hybrid_type = Exx_Global::Hybrid_Type::Generate_Matrix;
		}
		exx_global.info.hybrid_alpha    = INPUT.exx_hybrid_alpha      ;
		exx_global.info.hse_omega       = INPUT.exx_hse_omega         ;
		exx_global.info.separate_loop   = INPUT.exx_separate_loop     ;
		exx_global.info.hybrid_step     = INPUT.exx_hybrid_step       ;
		exx_lip.info.lambda             = INPUT.exx_lambda            ;
		exx_lcao.info.pca_threshold     = INPUT.exx_pca_threshold     ;
		exx_lcao.info.c_threshold       = INPUT.exx_c_threshold       ;
		exx_lcao.info.v_threshold       = INPUT.exx_v_threshold       ;
		exx_lcao.info.dm_threshold      = INPUT.exx_dm_threshold      ;
		exx_lcao.info.schwarz_threshold = INPUT.exx_schwarz_threshold ;
		exx_lcao.info.cauchy_threshold  = INPUT.exx_cauchy_threshold  ;
		exx_lcao.info.ccp_threshold     = INPUT.exx_ccp_threshold     ;
		exx_lcao.info.ccp_rmesh_times   = INPUT.exx_ccp_rmesh_times   ;
		if(INPUT.exx_distribute_type=="htime")
		{
			exx_lcao.info.distribute_type = Exx_Lcao::Distribute_Type::Htime;
		}
		else if(INPUT.exx_distribute_type=="kmeans2")
		{
			exx_lcao.info.distribute_type = Exx_Lcao::Distribute_Type::Kmeans2;
		}
		else if(INPUT.exx_distribute_type=="kmeans1")
		{
			exx_lcao.info.distribute_type = Exx_Lcao::Distribute_Type::Kmeans1;
		}
		else if(INPUT.exx_distribute_type=="order")
		{
			exx_lcao.info.distribute_type = Exx_Lcao::Distribute_Type::Order;
		}
		Exx_Abfs::Jle::Lmax      = INPUT.exx_opt_orb_lmax;
		Exx_Abfs::Jle::Ecut_exx  = INPUT.exx_opt_orb_ecut;
		Exx_Abfs::Jle::tolerence = INPUT.exx_opt_orb_tolerence;
	}

    ppcell.cell_factor = INPUT.cell_factor; //LiuXh add 20180619

    NEW_DM=INPUT.newDM;  // Shen Yu add 2019/5/9

//----------------------------------------------------------
// main parameters / electrons / spin ( 2/16 )
//----------------------------------------------------------
//	electrons::nelup = INPUT.nelup;
//	electrons::neldw = INPUT.neldw;

//----------------------------------------------------------
// occupation (3/3)
//----------------------------------------------------------
    Occupy::decision(INPUT.occupations,INPUT.smearing,INPUT.degauss);
//----------------------------------------------------------
// charge mixing(3/3)
//----------------------------------------------------------
    CHR.set_mixing(INPUT.mixing_mode, INPUT.mixing_beta, 
	INPUT.mixing_ndim, INPUT.mixing_gg0); //mohan modify 2014-09-27, add mixing_gg0

//----------------------------------------------------------
// iteration 
//----------------------------------------------------------
    NITER = INPUT.niter;
    NSTEP = INPUT.nstep;

//----------------------------------------------------------
// wavefunction / charge / potential / (2/4)
//----------------------------------------------------------
    pot.start_pot = INPUT.start_pot;
	pot.extra_pot = INPUT.charge_extrap;//xiaohui modify 2015-02-01
    CHR.out_charge = INPUT.out_charge;
	pot.out_potential = INPUT.out_potential;
    wf.out_wf = INPUT.out_wf;
	en.out_dos = INPUT.out_dos;
    en.out_band = INPUT.out_band;
	LOC.out_dm = INPUT.out_dm;
	ParaO.out_hs = INPUT.out_hs;
	ParaO.out_hsR = INPUT.out_hs2; //LiuXh add 2019-07-16
	ParaO.out_lowf = INPUT.out_lowf;

	en.dos_emin_ev = INPUT.dos_emin_ev;
	en.dos_emax_ev = INPUT.dos_emax_ev;
	en.dos_edelta_ev = INPUT.dos_edelta_ev;
    en.bcoeff = INPUT.b_coef;

//----------------------------------------------------------
// About LCAO
//----------------------------------------------------------
	ORB.ecutwfc = INPUT.lcao_ecut;
	ORB.dk = INPUT.lcao_dk;
	ORB.dR = INPUT.lcao_dr;
	ORB.Rmax = INPUT.lcao_rmax; 

	// mohan add 2021-02-16
	berryphase::berry_phase_flag = INPUT.berry_phase;

	timer::tick("Input_Conv","Convert",'B');
    return;
}

