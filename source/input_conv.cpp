#include "input.h"
#include "input_conv.h"

#include "src_io/optical.h"
#include "src_io/chi0_hilbert.h"
#include "src_io/chi0_standard.h"
#include "src_io/epsilon0_pwscf.h"
#include "src_io/epsilon0_vasp.h"
#include "src_io/berryphase.h"
#include "src_ions/ions_move_basic.h"
#include "src_pw/global.h"
#include "src_pw/tools.h"
#include "module_symmetry/symmetry.h"
#include "src_pw/efield.h"
#include "src_pw/occupy.h"
#include "module_cell/unitcell.h"
#include "src_ri/exx_abfs-jle.h"
#ifdef __LCAO
#include "module_orbital/ORB_read.h"
#include "src_lcao/global_fp.h"
#include "src_lcao/FORCE_STRESS.h"
#include "src_lcao/local_orbital_charge.h"
#include "src_lcao/ELEC_evolve.h"
#endif 

void Input_Conv::Convert(void)
{
    TITLE("Input_Conv","Convert");
	timer::tick("Input_Conv","Convert");
//----------------------------------------------------------
// main parameters / electrons / spin ( 10/16 )
//----------------------------------------------------------
//  suffix
    if(INPUT.atom_file!="") GlobalV::global_atom_card = INPUT.atom_file;
	GlobalV::global_wannier_card = INPUT.wannier_card;
    if(INPUT.kpoint_file!= "") GlobalV::global_kpoint_card = INPUT.kpoint_file;
    if(INPUT.pseudo_dir != "") GlobalV::global_pseudo_dir = INPUT.pseudo_dir + "/";
    GlobalV::global_pseudo_type = INPUT.pseudo_type;
	GlobalC::ucell.latName = INPUT.latname; 
	GlobalC::ucell.ntype = INPUT.ntype;
	GlobalC::ucell.lmaxmax = INPUT.lmaxmax;
	GlobalC::ucell.set_vel = INPUT.set_vel;

    GlobalV::NBANDS = INPUT.nbands;
    GlobalV::NBANDS_ISTATE = INPUT.nbands_istate;
	GlobalV::NPOOL = INPUT.npool;
	GlobalV::CALCULATION = INPUT.calculation;

	GlobalV::PSEUDORCUT = INPUT.pseudo_rcut; 
    GlobalV::RENORMWITHMESH = INPUT.renormwithmesh;

	// qianrui add 2021-2-5
	GlobalC::sto_wf.nchi = INPUT.nbands_sto;
	GlobalC::sto_wf.nche_sto = INPUT.nche_sto;
	GlobalC::sto_wf.seed_sto = INPUT.seed_sto;
	GlobalC::sto_wf.emax_sto = INPUT.emax_sto;
	GlobalC::sto_wf.emin_sto = INPUT.emin_sto;
	GlobalC::sto_wf.stotype = INPUT.stotype;

	// Electrical Field
	GlobalV::EFIELD = INPUT.efield;
	Efield::edir = INPUT.edir;
	Efield::emaxpos = INPUT.emaxpos;
	Efield::eopreg = INPUT.eopreg;
	Efield::eamp = INPUT.eamp;

	// optical
	Optical::opt_epsilon2=INPUT.opt_epsilon2;	// mohan add 2010-03-24	
	Optical::opt_nbands=INPUT.opt_nbands;		// number of bands for optical transition.	

	GlobalV::DFT_FUNCTIONAL = INPUT.dft_functional;
    GlobalV::NSPIN = INPUT.nspin;
    GlobalV::CURRENT_SPIN = 0;

    GlobalV::FORCE = INPUT.force;
    GlobalV::FORCE_THR = INPUT.force_thr;

    GlobalV::STRESS_THR = INPUT.stress_thr;
    GlobalV::PRESS1 = INPUT.press1;
    GlobalV::PRESS2 = INPUT.press2;
    GlobalV::PRESS3 = INPUT.press3;
#ifdef __LCAO
	Force_Stress_LCAO::force_invalid_threshold_ev = INPUT.force_thr_ev2;
#endif

	BFGS_Basic::w1 = INPUT.bfgs_w1;
	BFGS_Basic::w2 = INPUT.bfgs_w2;

	Ions_Move_Basic::trust_radius_max = INPUT.trust_radius_max;
	Ions_Move_Basic::trust_radius_min = INPUT.trust_radius_min;
	Ions_Move_Basic::trust_radius_ini = INPUT.trust_radius_ini;
	Ions_Move_Basic::out_stru = INPUT.out_stru; //mohan add 2012-03-23

	GlobalV::STRESS = INPUT.stress;


	// pengfei Li add 2018-11-11
	if(INPUT.fixed_axes == "None")
	{
		GlobalC::ucell.lc[0] = 1; GlobalC::ucell.lc[1] = 1; GlobalC::ucell.lc[2] = 1;
	}
	else if(INPUT.fixed_axes == "volume")
	{
		GlobalC::ucell.lc[0] = 1; GlobalC::ucell.lc[1] = 1; GlobalC::ucell.lc[2] = 1;
	}
	else if(INPUT.fixed_axes == "a")
	{
		GlobalC::ucell.lc[0] = 0; GlobalC::ucell.lc[1] = 1; GlobalC::ucell.lc[2] = 1;
	}
	else if(INPUT.fixed_axes == "b")
	{
		GlobalC::ucell.lc[0] = 1; GlobalC::ucell.lc[1] = 0; GlobalC::ucell.lc[2] = 1;
	}
	else if(INPUT.fixed_axes == "c")
	{
		GlobalC::ucell.lc[0] = 1; GlobalC::ucell.lc[1] = 1; GlobalC::ucell.lc[2] = 0;
	}
	else if(INPUT.fixed_axes == "ab")
	{
		GlobalC::ucell.lc[0] = 0; GlobalC::ucell.lc[1] = 0; GlobalC::ucell.lc[2] = 1;
	}
	else if(INPUT.fixed_axes == "ac")
	{
		GlobalC::ucell.lc[0] = 0; GlobalC::ucell.lc[1] = 1; GlobalC::ucell.lc[2] = 0;
	}
	else if(INPUT.fixed_axes == "bc")
	{
		GlobalC::ucell.lc[0] = 1; GlobalC::ucell.lc[1] = 0; GlobalC::ucell.lc[2] = 0;
	}
	else if(INPUT.fixed_axes == "abc")
	{
		GlobalC::ucell.lc[0] = 0; GlobalC::ucell.lc[1] = 0; GlobalC::ucell.lc[2] = 0;
	}
	else
	{
		WARNING_QUIT("Input", "fixed_axes should be None,a,b,c,ab,ac,bc or abc!");
	}

	GlobalV::MOVE_IONS = INPUT.ion_dynamics;
	GlobalV::OUT_LEVEL = INPUT.out_level;
	Ions_Move_CG::CG_THRESHOLD = INPUT.cg_threshold; // pengfei add 2013-09-09

    ModuleSymmetry::Symmetry::symm_flag = INPUT.symmetry;						// 9
	GlobalV::BASIS_TYPE = INPUT.basis_type;
	GlobalV::KS_SOLVER = INPUT.ks_solver;
	GlobalV::SEARCH_RADIUS = INPUT.search_radius;
	GlobalV::SEARCH_PBC = INPUT.search_pbc;

//----------------------------------------------------------
// planewave (8/8)
//----------------------------------------------------------
    GlobalC::pw.set(
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
	GlobalV::GAMMA_ONLY_LOCAL = INPUT.gamma_only_local;

//----------------------------------------------------------
// diagonalization  (5/5)
//----------------------------------------------------------
	GlobalV::DIAGO_PROC = INPUT.diago_proc;
    GlobalV::DIAGO_CG_MAXITER = INPUT.diago_cg_maxiter;
	GlobalV::DIAGO_CG_PREC = INPUT.diago_cg_prec;
    GlobalV::DIAGO_DAVID_NDIM = INPUT.diago_david_ndim;
    GlobalV::ETHR = INPUT.ethr;
	GlobalV::NB2D = INPUT.nb2d;
	GlobalV::NURSE = INPUT.nurse;
	GlobalV::COLOUR = INPUT.colour;
	GlobalV::T_IN_H = INPUT.t_in_h;
	GlobalV::VL_IN_H = INPUT.vl_in_h;
	GlobalV::VNL_IN_H = INPUT.vnl_in_h;
	GlobalV::VH_IN_H = INPUT.vh_in_h;
	GlobalV::VXC_IN_H = INPUT.vxc_in_h;
	GlobalV::VION_IN_H = INPUT.vion_in_h;
	GlobalV::TEST_FORCE = INPUT.test_force;
	GlobalV::TEST_STRESS = INPUT.test_stress;

//----------------------------------------------------------
// iteration (1/3)
//----------------------------------------------------------
    GlobalV::DRHO2 = INPUT.dr2;

//----------------------------------------------------------
// wavefunction / charge / potential / (2/4)
//----------------------------------------------------------
    GlobalV::RESTART_MODE = INPUT.restart_mode;
    GlobalC::wf.start_wfc = INPUT.start_wfc;
	GlobalC::wf.mem_saver = INPUT.mem_saver; //mohan add 2010-09-07
	GlobalC::en.printe    = INPUT.printe; // mohan add 2011-03-16

    
//----------------------------------------------------------
// about spectrum, pengfei 2016-12-14
//----------------------------------------------------------
#ifdef __LCAO
	if( (INPUT.spectral_type == "eels" && INPUT.eels_method == 0) 
		|| (INPUT.spectral_type == "None" 
		&& INPUT.eels_method == 0 
		&& INPUT.kmesh_interpolation) )
	{
		if(INPUT.spectral_type == "eels")
		{
			GlobalC::chi0_hilbert.epsilon = true;
		}
		else if(INPUT.spectral_type == "None")
		{
			GlobalC::chi0_hilbert.epsilon = false;
		}
		//GlobalC::chi0_hilbert.epsilon = INPUT.epsilon;
		GlobalC::chi0_hilbert.kernel_type = INPUT.kernel_type;
		GlobalC::chi0_hilbert.system = INPUT.system_type;
		GlobalC::chi0_hilbert.eta = INPUT.eta;
		GlobalC::chi0_hilbert.domega = INPUT.domega;
		GlobalC::chi0_hilbert.nomega = INPUT.nomega;
		GlobalC::chi0_hilbert.dim = INPUT.ecut_chi; 
		//GlobalC::chi0_hilbert.oband = INPUT.oband;

		GlobalC::chi0_hilbert.q_start[0] = INPUT.q_start[0];  
		GlobalC::chi0_hilbert.q_start[1] = INPUT.q_start[1]; 
		GlobalC::chi0_hilbert.q_start[2] = INPUT.q_start[2];

		GlobalC::chi0_hilbert.direct[0] = INPUT.q_direct[0]; 
		GlobalC::chi0_hilbert.direct[1] = INPUT.q_direct[1]; 
		GlobalC::chi0_hilbert.direct[2] = INPUT.q_direct[2];

		//GlobalC::chi0_hilbert.start_q = INPUT.start_q;
		//GlobalC::chi0_hilbert.interval_q = INPUT.interval_q;
		GlobalC::chi0_hilbert.nq = INPUT.nq;
		GlobalC::chi0_hilbert.out_epsilon = INPUT.out_epsilon;
		GlobalC::chi0_hilbert.out_chi = INPUT.out_chi;
		GlobalC::chi0_hilbert.out_chi0 = INPUT.out_chi0;
		GlobalC::chi0_hilbert.fermi_level = INPUT.fermi_level;
		GlobalC::chi0_hilbert.coulomb_cutoff = INPUT.coulomb_cutoff;
		GlobalC::chi0_hilbert.kmesh_interpolation = INPUT.kmesh_interpolation;
		for(int i=0; i<100; i++)
		{
			GlobalC::chi0_hilbert.qcar[i][0] = INPUT.qcar[i][0]; 
			GlobalC::chi0_hilbert.qcar[i][1] = INPUT.qcar[i][1]; 
			GlobalC::chi0_hilbert.qcar[i][2] = INPUT.qcar[i][2]; 
		}
		GlobalC::chi0_hilbert.lcao_box[0] = INPUT.lcao_box[0]; 
		GlobalC::chi0_hilbert.lcao_box[1] = INPUT.lcao_box[1]; 
		GlobalC::chi0_hilbert.lcao_box[2] = INPUT.lcao_box[2];
	}
#endif
	
	//if( INPUT.epsilon && (INPUT.epsilon_choice == 1))
	if( INPUT.spectral_type == "eels" && INPUT.eels_method == 1)
	{
		//GlobalC::chi0_standard.epsilon = INPUT.epsilon;
		GlobalC::chi0_standard.epsilon = true;
		GlobalC::chi0_standard.system = INPUT.system_type;
		GlobalC::chi0_standard.eta = INPUT.eta;
		GlobalC::chi0_standard.domega = INPUT.domega;
		GlobalC::chi0_standard.nomega = INPUT.nomega;
		GlobalC::chi0_standard.dim = INPUT.ecut_chi;
		//GlobalC::chi0_standard.oband = INPUT.oband;
		GlobalC::chi0_standard.q_start[0] = INPUT.q_start[0]; 
	 	GlobalC::chi0_standard.q_start[1] = INPUT.q_start[1]; 
		GlobalC::chi0_standard.q_start[2] = INPUT.q_start[2];
		GlobalC::chi0_standard.direct[0] = INPUT.q_direct[0];  
		GlobalC::chi0_standard.direct[1] = INPUT.q_direct[1]; 
		GlobalC::chi0_standard.direct[2] = INPUT.q_direct[2];
		//GlobalC::chi0_standard.start_q = INPUT.start_q;
		//GlobalC::chi0_standard.interval_q = INPUT.interval_q;
		GlobalC::chi0_standard.nq = INPUT.nq;
		GlobalC::chi0_standard.out_epsilon = INPUT.out_epsilon;		
	}
	
	//if( INPUT.epsilon0 && (INPUT.epsilon0_choice == 1) )
	if( INPUT.spectral_type == "absorption" && INPUT.absorption_method == 1)
	{
		//GlobalC::epsilon0_pwscf.epsilon = INPUT.epsilon0;
		GlobalC::epsilon0_pwscf.epsilon = true;
		GlobalC::epsilon0_pwscf.intersmear = INPUT.eta;
		GlobalC::epsilon0_pwscf.intrasmear = INPUT.intrasmear;
		GlobalC::epsilon0_pwscf.domega = INPUT.domega;
		GlobalC::epsilon0_pwscf.nomega = INPUT.nomega;
		GlobalC::epsilon0_pwscf.shift = INPUT.shift;
		GlobalC::epsilon0_pwscf.metalcalc = INPUT.metalcalc;
		GlobalC::epsilon0_pwscf.degauss = INPUT.eps_degauss;
	}
	
	//if( INPUT.epsilon0 && (INPUT.epsilon0_choice == 0))
	if( INPUT.spectral_type == "absorption" && INPUT.absorption_method == 0)
	{
		//GlobalC::epsilon0_vasp.epsilon = INPUT.epsilon0;
		GlobalC::epsilon0_vasp.epsilon = true;
		GlobalC::epsilon0_vasp.domega = INPUT.domega;
		GlobalC::epsilon0_vasp.nomega = INPUT.nomega;
		GlobalC::epsilon0_vasp.eta = INPUT.eta;
	}

//--------------------------------------------
// added by zhengdy-soc
//--------------------------------------------
	if(INPUT.noncolin||INPUT.lspinorb) 
	{
		GlobalV::NSPIN = 4;
	}
	if(GlobalV::NSPIN == 4)
	{
		GlobalV::NONCOLIN = INPUT.noncolin;
		//wavefunctions are spinors with 2 components
		GlobalV::NPOL = 2;
		//set the domag variable to make a spin-orbit calculation with zero magnetization
		if(GlobalV::NONCOLIN)
		{
			GlobalV::DOMAG = true;
			GlobalV::DOMAG_Z = false;
		}
		else{
			GlobalV::DOMAG = false;
			GlobalV::DOMAG_Z = true;
		}
		GlobalV::LSPINORB = INPUT.lspinorb;

		delete[] GlobalC::ucell.magnet.m_loc_;
		delete[] GlobalC::ucell.magnet.angle1_;
		delete[] GlobalC::ucell.magnet.angle2_;
		GlobalC::ucell.magnet.m_loc_ = new Vector3<double> [INPUT.ntype];
		GlobalC::ucell.magnet.angle1_ = new double[INPUT.ntype];
		GlobalC::ucell.magnet.angle2_ = new double[INPUT.ntype];
		for(int i = 0;i<INPUT.ntype;i++)
		{
			GlobalC::ucell.magnet.angle1_[i] = INPUT.angle1[i]/180*PI;
			GlobalC::ucell.magnet.angle2_[i] = INPUT.angle2[i]/180*PI;
		}
#ifdef __MPI
//			Parallel_Common::bcast_double(GlobalC::ucell.magnet.angle1_[i]);
//			Parallel_Common::bcast_double(GlobalC::ucell.magnet.angle2_[i]);
#endif
	}
	else{
		delete[] GlobalC::ucell.magnet.m_loc_;
		GlobalC::ucell.magnet.m_loc_ = new Vector3<double> [INPUT.ntype];
		GlobalV::LSPINORB = false;
		GlobalV::NONCOLIN = false;
		GlobalV::DOMAG = false;
		GlobalV::DOMAG_Z = false;
		GlobalV::NPOL = 1;
	}
	
//----------------------------------------------------------
// Fuxiang He add 2016-10-26
//----------------------------------------------------------
#ifdef __LCAO
	ELEC_evolve::tddft = INPUT.tddft;
	ELEC_evolve::td_dr2 = INPUT.td_dr2;
	ELEC_evolve::td_dt = INPUT.td_dt;
	ELEC_evolve::td_force_dt = INPUT.td_force_dt;
	ELEC_evolve::td_val_elec_01 = INPUT.td_val_elec_01;
	ELEC_evolve::td_val_elec_02 = INPUT.td_val_elec_02;
	ELEC_evolve::td_val_elec_03 = INPUT.td_val_elec_03;
	ELEC_evolve::td_vext = INPUT.td_vext;
	ELEC_evolve::td_vext_dire = INPUT.td_vext_dire;	
	ELEC_evolve::td_timescale = INPUT.td_timescale;
	ELEC_evolve::td_vexttype = INPUT.td_vexttype;
	ELEC_evolve::td_vextout = INPUT.td_vextout;
	ELEC_evolve::td_dipoleout = INPUT.td_dipoleout;
#endif



	// jiyy add 2020.10.11	
	GlobalV::ocp = INPUT.ocp;
     //ocp_n = INPUT.ocp_n;
    GlobalV::ocp_set = INPUT.ocp_set;
    if(GlobalV::ocp == 1)
	{
		int count = 0;
		string pattern("([0-9]+\\*[0-9.]+|[0-9,.]+)");
		vector<string> str;
		string::size_type pos1, pos2;
		string c = " ";
		pos2 = GlobalV::ocp_set.find(c);
		pos1 = 0;
		while(string::npos != pos2)
		{
			str.push_back(GlobalV::ocp_set.substr(pos1, pos2-pos1));
 
			pos1 = pos2 + c.size();
			pos2 = GlobalV::ocp_set.find(c, pos1);
		}
		if(pos1 != GlobalV::ocp_set.length())
		{
			str.push_back(GlobalV::ocp_set.substr(pos1));
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
				const vector<double>::iterator dest = GlobalV::ocp_kb.begin()+count;
				copy(ocp_temp.begin(), ocp_temp.end(), dest);
				count += num;
			}
			else
			{
				GlobalV::ocp_kb[count] = stof(sub_str);
				count += 1;
			}
		}
	}
		
    GlobalV::mulliken = INPUT.mulliken;//qifeng add 2019/9/10

//----------------------------------------------------------
// about restart, // Peize Lin add 2020-04-04
//----------------------------------------------------------	
	if(INPUT.restart_save)
	{
		GlobalC::restart.folder = GlobalV::global_out_dir + "restart/";
		const string command0 =  "test -d " + GlobalC::restart.folder + " || mkdir " + GlobalC::restart.folder;
		if(GlobalV::MY_RANK==0)
			system( command0.c_str() );
		if(INPUT.exx_hybrid_type=="no")
		{
			GlobalC::restart.info_save.save_charge = true;
		}
		else
		{
			GlobalC::restart.info_save.save_charge = true;
			GlobalC::restart.info_save.save_H = true;
		}
	}
	if(INPUT.restart_load)
	{
		GlobalC::restart.folder = GlobalV::global_out_dir + "restart/";
		if(INPUT.exx_hybrid_type=="no")
		{
			GlobalC::restart.info_load.load_charge = true;
		}
		else
		{
			GlobalC::restart.info_load.load_charge = true;
			GlobalC::restart.info_load.load_H = true;
		}
	}

//----------------------------------------------------------
// about exx, Peize Lin add 2018-06-20
//----------------------------------------------------------	
#ifdef __LCAO
	if(INPUT.exx_hybrid_type=="no")
	{
		GlobalC::exx_global.info.hybrid_type = Exx_Global::Hybrid_Type::No;
	}
	else
	{
		if(INPUT.exx_hybrid_type=="hf")
		{
			GlobalC::exx_global.info.hybrid_type = Exx_Global::Hybrid_Type::HF;
		}
		else if(INPUT.exx_hybrid_type=="pbe0")
		{
			GlobalC::exx_global.info.hybrid_type = Exx_Global::Hybrid_Type::PBE0;
		}
		else if(INPUT.exx_hybrid_type=="hse")
		{
			GlobalC::exx_global.info.hybrid_type = Exx_Global::Hybrid_Type::HSE;
		}
		else if(INPUT.exx_hybrid_type=="opt_orb")
		{
			GlobalC::exx_global.info.hybrid_type = Exx_Global::Hybrid_Type::Generate_Matrix;
		}
		GlobalC::exx_global.info.hybrid_alpha    = INPUT.exx_hybrid_alpha      ;
		GlobalC::exx_global.info.hse_omega       = INPUT.exx_hse_omega         ;
		GlobalC::exx_global.info.separate_loop   = INPUT.exx_separate_loop     ;
		GlobalC::exx_global.info.hybrid_step     = INPUT.exx_hybrid_step       ;
		GlobalC::exx_lip.info.lambda             = INPUT.exx_lambda            ;
		GlobalC::exx_lcao.info.pca_threshold     = INPUT.exx_pca_threshold     ;
		GlobalC::exx_lcao.info.c_threshold       = INPUT.exx_c_threshold       ;
		GlobalC::exx_lcao.info.v_threshold       = INPUT.exx_v_threshold       ;
		GlobalC::exx_lcao.info.dm_threshold      = INPUT.exx_dm_threshold      ;
		GlobalC::exx_lcao.info.schwarz_threshold = INPUT.exx_schwarz_threshold ;
		GlobalC::exx_lcao.info.cauchy_threshold  = INPUT.exx_cauchy_threshold  ;
		GlobalC::exx_lcao.info.ccp_threshold     = INPUT.exx_ccp_threshold     ;
		GlobalC::exx_lcao.info.ccp_rmesh_times   = INPUT.exx_ccp_rmesh_times   ;
		if(INPUT.exx_distribute_type=="htime")
		{
			GlobalC::exx_lcao.info.distribute_type = Exx_Lcao::Distribute_Type::Htime;
		}
		else if(INPUT.exx_distribute_type=="kmeans2")
		{
			GlobalC::exx_lcao.info.distribute_type = Exx_Lcao::Distribute_Type::Kmeans2;
		}
		else if(INPUT.exx_distribute_type=="kmeans1")
		{
			GlobalC::exx_lcao.info.distribute_type = Exx_Lcao::Distribute_Type::Kmeans1;
		}
		else if(INPUT.exx_distribute_type=="order")
		{
			GlobalC::exx_lcao.info.distribute_type = Exx_Lcao::Distribute_Type::Order;
		}
		Exx_Abfs::Jle::Lmax      = INPUT.exx_opt_orb_lmax;
		Exx_Abfs::Jle::Ecut_exx  = INPUT.exx_opt_orb_ecut;
		Exx_Abfs::Jle::tolerence = INPUT.exx_opt_orb_tolerence;
	}
#endif

    GlobalC::ppcell.cell_factor = INPUT.cell_factor; //LiuXh add 20180619

//    NEW_DM=INPUT.new_dm;  // Shen Yu add 2019/5/9

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
    GlobalC::CHR.set_mixing(INPUT.mixing_mode, INPUT.mixing_beta, 
	INPUT.mixing_ndim, INPUT.mixing_gg0); //mohan modify 2014-09-27, add mixing_gg0

//----------------------------------------------------------
// iteration 
//----------------------------------------------------------
    GlobalV::NITER = INPUT.niter;
    GlobalV::NSTEP = INPUT.nstep;

//----------------------------------------------------------
// wavefunction / charge / potential / (2/4)
//----------------------------------------------------------
    GlobalC::pot.start_pot = INPUT.start_pot;
	GlobalC::pot.extra_pot = INPUT.charge_extrap;//xiaohui modify 2015-02-01
    GlobalC::CHR.out_charge = INPUT.out_charge;
	GlobalC::CHR.nelec = INPUT.nelec;
	GlobalC::pot.out_potential = INPUT.out_potential;
    GlobalC::wf.out_wf = INPUT.out_wf;
	GlobalC::en.out_dos = INPUT.out_dos;
    GlobalC::en.out_band = INPUT.out_band;
#ifdef __LCAO
	GlobalC::LOC.out_dm = INPUT.out_dm;
	GlobalC::ParaO.out_hs = INPUT.out_hs;
	GlobalC::ParaO.out_hsR = INPUT.out_hs2; //LiuXh add 2019-07-16
	GlobalC::ParaO.out_lowf = INPUT.out_lowf;
#endif

	GlobalC::en.dos_emin_ev = INPUT.dos_emin_ev;
	GlobalC::en.dos_emax_ev = INPUT.dos_emax_ev;
	GlobalC::en.dos_edelta_ev = INPUT.dos_edelta_ev;
	GlobalC::en.dos_scale= INPUT.dos_scale;
    GlobalC::en.bcoeff = INPUT.b_coef;

//----------------------------------------------------------
// About LCAO
//----------------------------------------------------------
// mohan add 2021-04-16
//	ORB.ecutwfc = INPUT.lcao_ecut;
//	ORB.dk = INPUT.lcao_dk;
//	ORB.dR = INPUT.lcao_dr;
//	ORB.Rmax = INPUT.lcao_rmax; 

	// mohan add 2021-02-16
	berryphase::berry_phase_flag = INPUT.berry_phase;

	// wenfei 2021-7-28
	if(GlobalV::DFT_FUNCTIONAL == "scan")
	{
		if(GlobalV::BASIS_TYPE != "pw")
		{
			 WARNING_QUIT("Input_conv","add metaGGA for pw first");
		}
		GlobalV::DFT_META = 1;
	}

	timer::tick("Input_Conv","Convert");
    return;
}

