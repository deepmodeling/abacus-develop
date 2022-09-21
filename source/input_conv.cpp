#include "input_conv.h"

#include "input.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_cell/unitcell.h"
#include "module_surchem/surchem.h"
#include "module_symmetry/symmetry.h"
#include "src_io/berryphase.h"
#include "src_io/chi0_hilbert.h"
#include "src_io/chi0_standard.h"
#include "src_io/epsilon0_pwscf.h"
#include "src_io/epsilon0_vasp.h"
#include "src_io/optical.h"
#include "module_relaxation/ions_move_basic.h"
#include "src_pw/global.h"
#include "src_pw/occupy.h"
#ifdef __EXX
#include "src_ri/exx_abfs-jle.h"
#endif
#ifdef __LCAO
#include "module_orbital/ORB_read.h"
#include "src_lcao/ELEC_evolve.h"
#include "src_lcao/FORCE_STRESS.h"
#include "src_lcao/dftu.h"
#include "src_lcao/global_fp.h"
#include "src_lcao/local_orbital_charge.h"
#endif
#include "module_base/timer.h"
#include "module_elecstate/elecstate_lcao.h"
#include "module_hsolver/hsolver_lcao.h"
#include "module_surchem/efield.h"
#include "module_surchem/gatefield.h"

void Input_Conv::Convert(void)
{
    ModuleBase::TITLE("Input_Conv", "Convert");
    ModuleBase::timer::tick("Input_Conv", "Convert");
    //----------------------------------------------------------
    // main parameters / electrons / spin ( 10/16 )
    //----------------------------------------------------------
    //  suffix
    if (INPUT.stru_file != "")
        GlobalV::stru_file = INPUT.stru_file;
    GlobalV::global_wannier_card = INPUT.wannier_card;
    if (INPUT.kpoint_file != "")
        GlobalV::global_kpoint_card = INPUT.kpoint_file;
    if (INPUT.pseudo_dir != "")
        GlobalV::global_pseudo_dir = INPUT.pseudo_dir + "/";
    if (INPUT.orbital_dir != "")
        GlobalV::global_orbital_dir = INPUT.orbital_dir + "/";
    // GlobalV::global_pseudo_type = INPUT.pseudo_type;
    GlobalC::ucell.setup(INPUT.latname, INPUT.ntype, INPUT.lmaxmax, INPUT.init_vel, INPUT.fixed_axes);

    GlobalV::KSPACING = INPUT.kspacing;
    GlobalV::MIN_DIST_COEF = INPUT.min_dist_coef;
    GlobalV::NBANDS = INPUT.nbands;
    GlobalC::wf.pw_seed = INPUT.pw_seed;
    GlobalV::NBANDS_ISTATE = INPUT.nbands_istate;
#if ((defined __CUDA) || (defined __ROCM))
    int temp_nproc;
    MPI_Comm_size(MPI_COMM_WORLD, &temp_nproc);
    if (temp_nproc != INPUT.kpar)
    {
        ModuleBase::WARNING("Input_conv", "None kpar set in INPUT file, auto set kpar value.");
    }
    GlobalV::KPAR = temp_nproc;
#else
    GlobalV::KPAR = INPUT.kpar;
    GlobalV::NSTOGROUP = INPUT.bndpar;
#endif
    GlobalV::CALCULATION = INPUT.calculation;

    GlobalV::PSEUDORCUT = INPUT.pseudo_rcut;
    GlobalV::PSEUDO_MESH = INPUT.pseudo_mesh;

    GlobalV::DFT_FUNCTIONAL = INPUT.dft_functional;
    GlobalV::XC_TEMPERATURE = INPUT.xc_temperature;
    GlobalV::NSPIN = INPUT.nspin;
    GlobalV::CURRENT_SPIN = 0;

    GlobalV::CAL_FORCE = INPUT.cal_force;
    GlobalV::FORCE_THR = INPUT.force_thr;

    GlobalV::STRESS_THR = INPUT.stress_thr;
    GlobalV::PRESS1 = INPUT.press1;
    GlobalV::PRESS2 = INPUT.press2;
    GlobalV::PRESS3 = INPUT.press3;
    GlobalV::out_element_info = INPUT.out_element_info;
#ifdef __LCAO
    Force_Stress_LCAO::force_invalid_threshold_ev = INPUT.force_thr_ev2;
#endif

    if((INPUT.calculation=="relax" || INPUT.calculation=="cell-relax") && INPUT.chg_extrap!="atomic")
    {
        std::cout << " For relaxation, charge extrapolation is set to atomic." << std::endl;
        INPUT.chg_extrap="atomic";
    }

    BFGS_Basic::relax_bfgs_w1 = INPUT.relax_bfgs_w1;
    BFGS_Basic::relax_bfgs_w2 = INPUT.relax_bfgs_w2;

    Ions_Move_Basic::relax_bfgs_rmax = INPUT.relax_bfgs_rmax;
    Ions_Move_Basic::relax_bfgs_rmin = INPUT.relax_bfgs_rmin;
    Ions_Move_Basic::relax_bfgs_init = INPUT.relax_bfgs_init;
    Ions_Move_Basic::out_stru = INPUT.out_stru; // mohan add 2012-03-23
    Lattice_Change_Basic::out_stru = INPUT.out_stru;

    GlobalV::CAL_STRESS = INPUT.cal_stress;

    GlobalV::RELAX_METHOD = INPUT.relax_method;
    GlobalV::OUT_LEVEL = INPUT.out_level;
    Ions_Move_CG::RELAX_CG_THR = INPUT.relax_cg_thr; // pengfei add 2013-09-09

    ModuleSymmetry::Symmetry::symm_flag = INPUT.symmetry; // 9
    GlobalC::symm.epsilon = INPUT.symmetry_prec; // LiuXh add 2021-08-12, accuracy for symmetry
    GlobalV::BASIS_TYPE = INPUT.basis_type;
    GlobalV::KS_SOLVER = INPUT.ks_solver;
    GlobalV::SEARCH_RADIUS = INPUT.search_radius;
    GlobalV::SEARCH_PBC = INPUT.search_pbc;

    //----------------------------------------------------------
    // planewave (8/8)
    //----------------------------------------------------------
    GlobalC::sf.set(INPUT.nbspline);
    GlobalV::GAMMA_ONLY_LOCAL = INPUT.gamma_only_local;

    //----------------------------------------------------------
    // diagonalization  (5/5)
    //----------------------------------------------------------
    GlobalV::DIAGO_PROC = INPUT.diago_proc;
    GlobalV::PW_DIAG_NMAX = INPUT.pw_diag_nmax;
    GlobalV::DIAGO_CG_PREC = INPUT.diago_cg_prec;
    GlobalV::PW_DIAG_NDIM = INPUT.pw_diag_ndim;
    GlobalV::PW_DIAG_THR = INPUT.pw_diag_thr;
    GlobalV::NB2D = INPUT.nb2d;
    GlobalV::NURSE = INPUT.nurse;
    GlobalV::COLOUR = INPUT.colour;
    GlobalV::T_IN_H = INPUT.t_in_h;
    GlobalV::VL_IN_H = INPUT.vl_in_h;
    GlobalV::VNL_IN_H = INPUT.vnl_in_h;
    GlobalV::VH_IN_H = INPUT.vh_in_h;
    GlobalV::VION_IN_H = INPUT.vion_in_h;
    GlobalV::TEST_FORCE = INPUT.test_force;
    GlobalV::TEST_STRESS = INPUT.test_stress;
    GlobalV::test_skip_ewald = INPUT.test_skip_ewald;

    //----------------------------------------------------------
    // iteration (1/3)
    //----------------------------------------------------------
    GlobalV::SCF_THR = INPUT.scf_thr;

    //----------------------------------------------------------
    // wavefunction / charge / potential / (2/4)
    //----------------------------------------------------------
    GlobalC::wf.init_wfc = INPUT.init_wfc;
    GlobalC::wf.mem_saver = INPUT.mem_saver; // mohan add 2010-09-07
    GlobalC::en.printe = INPUT.printe; // mohan add 2011-03-16

    if (INPUT.dft_plus_u)
    {
        GlobalC::dftu.dftu_type = INPUT.dftu_type;
        GlobalC::dftu.double_counting = INPUT.double_counting;
        GlobalC::dftu.Yukawa = INPUT.yukawa_potential;
        GlobalC::dftu.omc = INPUT.omc;
        GlobalC::dftu.orbital_corr = INPUT.orbital_corr;
        if (!INPUT.yukawa_potential)
        {
            GlobalC::dftu.U = INPUT.hubbard_u; // Hubbard Coulomb interaction parameter U(ev)
            GlobalC::dftu.J = INPUT.hund_j; // Hund exchange parameter J(ev)
        }
    }
    /*
#ifndef __CMD
    GlobalC::ucell.n_mag_at=INPUT.n_mag_at;
    GlobalC::ucell.atom_mag=INPUT.atom_mag;
#endif*/
    //--------------------------------------------
    // added by zhengdy-soc
    //--------------------------------------------
    if (INPUT.noncolin || INPUT.lspinorb)
    {
        GlobalV::NSPIN = 4;
    }

    if (GlobalV::NSPIN == 4)
    {
        GlobalV::NONCOLIN = INPUT.noncolin;
        // wavefunctions are spinors with 2 components
        GlobalV::NPOL = 2;
        // set the domag variable to make a spin-orbit calculation with zero magnetization
        if (GlobalV::NONCOLIN)
        {
            GlobalV::DOMAG = true;
            GlobalV::DOMAG_Z = false;
        }
        else
        {
            GlobalV::DOMAG = false;
            GlobalV::DOMAG_Z = true;
        }
        GlobalV::LSPINORB = INPUT.lspinorb;
        GlobalV::soc_lambda = INPUT.soc_lambda;
    }
    else
    {
        GlobalV::LSPINORB = false;
        GlobalV::NONCOLIN = false;
        GlobalV::DOMAG = false;
        GlobalV::DOMAG_Z = false;
        GlobalV::NPOL = 1;
    }

    //----------------------------------------------------------
    // Yu Liu add 2022-05-18
    //----------------------------------------------------------
    GlobalV::EFIELD_FLAG = INPUT.efield_flag;
    GlobalV::DIP_COR_FLAG = INPUT.dip_cor_flag;
    Efield::efield_dir = INPUT.efield_dir;
    Efield::efield_pos_max = INPUT.efield_pos_max;
    Efield::efield_pos_dec = INPUT.efield_pos_dec;
    Efield::efield_amp = INPUT.efield_amp;

    //----------------------------------------------------------
    // Yu Liu add 2022-09-13
    //----------------------------------------------------------
    GlobalV::GATE_FLAG = INPUT.gate_flag;
    GlobalV::NELEC = INPUT.nelec;
    Gatefield::zgate = INPUT.zgate;
    Gatefield::relax = INPUT.relax;
    Gatefield::block = INPUT.block;
    Gatefield::block_down = INPUT.block_down;
    Gatefield::block_up = INPUT.block_up;
    Gatefield::block_height = INPUT.block_height;

//----------------------------------------------------------
// Fuxiang He add 2016-10-26
//----------------------------------------------------------
#ifdef __LCAO
    ELEC_evolve::tddft = INPUT.tddft;
    ELEC_evolve::td_scf_thr = INPUT.td_scf_thr;
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

    // setting for constrained DFT, jiyy add 2020.10.11
    // For example, when we studying nitrogen-vacancy center,
    // it requires an additional excitation of an electron conduction band to simulate the excited state,
    // used for TDDFT only.
    if (GlobalV::ocp == 1)
    {
        int count = 0;
        std::string pattern("([0-9]+\\*[0-9.]+|[0-9,.]+)");
        std::vector<std::string> str;
        std::string::size_type pos1, pos2;
        std::string c = " ";
        pos2 = GlobalV::ocp_set.find(c);
        pos1 = 0;
        while (std::string::npos != pos2)
        {
            str.push_back(GlobalV::ocp_set.substr(pos1, pos2 - pos1));

            pos1 = pos2 + c.size();
            pos2 = GlobalV::ocp_set.find(c, pos1);
        }
        if (pos1 != GlobalV::ocp_set.length())
        {
            str.push_back(GlobalV::ocp_set.substr(pos1));
        }

        regex_t reg;
        regcomp(&reg, pattern.c_str(), REG_EXTENDED);
        regmatch_t pmatch[1];
        const size_t nmatch = 1;
        for (int i = 0; i < str.size(); ++i)
        {
            if (str[i] == "")
            {
                continue;
            }
            int status = regexec(&reg, str[i].c_str(), nmatch, pmatch, 0);
            std::string sub_str = "";
            for (int j = pmatch[0].rm_so; j != pmatch[0].rm_eo; ++j)
            {
                sub_str += str[i][j];
            }
            std::string sub_pattern("\\*");
            regex_t sub_reg;
            regcomp(&sub_reg, sub_pattern.c_str(), REG_EXTENDED);
            regmatch_t sub_pmatch[1];
            const size_t sub_nmatch = 1;
            if (regexec(&sub_reg, sub_str.c_str(), sub_nmatch, sub_pmatch, 0) == 0)
            {
                int pos = sub_str.find("*");
                int num = stoi(sub_str.substr(0, pos));
                double occ = stof(sub_str.substr(pos + 1, sub_str.size()));
                std::vector<double> ocp_temp(num, occ);
                const std::vector<double>::iterator dest = GlobalV::ocp_kb.begin() + count;
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

    GlobalV::out_mul = INPUT.out_mul; // qifeng add 2019/9/10

    //----------------------------------------------------------
    // about restart, // Peize Lin add 2020-04-04
    //----------------------------------------------------------
    if (INPUT.restart_save)
    {
        GlobalC::restart.folder = GlobalV::global_out_dir + "restart/";
        const std::string command0 = "test -d " + GlobalC::restart.folder + " || mkdir " + GlobalC::restart.folder;
        if (GlobalV::MY_RANK == 0)
            system(command0.c_str());
        if (INPUT.dft_functional == "hf" || INPUT.dft_functional == "pbe0" || INPUT.dft_functional == "hse"
            || INPUT.dft_functional == "opt_orb" || INPUT.dft_functional == "scan0")
        {
            GlobalC::restart.info_save.save_charge = true;
            GlobalC::restart.info_save.save_H = true;
        }
        else
        {
            GlobalC::restart.info_save.save_charge = true;
        }
    }
    if (INPUT.restart_load)
    {
        GlobalC::restart.folder = GlobalV::global_out_dir + "restart/";
        if (INPUT.dft_functional == "hf" || INPUT.dft_functional == "pbe0" || INPUT.dft_functional == "hse"
            || INPUT.dft_functional == "opt_orb" || INPUT.dft_functional == "scan0")
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
#ifdef __MPI // liyuanbo 2022/2/23
#ifdef __LCAO

    if (INPUT.dft_functional == "hf")
    {
        GlobalC::exx_global.info.hybrid_type = Exx_Global::Hybrid_Type::HF;
    }
    else if (INPUT.dft_functional == "pbe0")
    {
        GlobalC::exx_global.info.hybrid_type = Exx_Global::Hybrid_Type::PBE0;
    }
    else if (INPUT.dft_functional == "scan0")
    {
        GlobalC::exx_global.info.hybrid_type = Exx_Global::Hybrid_Type::SCAN0;
    }
    else if (INPUT.dft_functional == "hse")
    {
        GlobalC::exx_global.info.hybrid_type = Exx_Global::Hybrid_Type::HSE;
    }
    else if (INPUT.dft_functional == "opt_orb")
    {
        GlobalC::exx_global.info.hybrid_type = Exx_Global::Hybrid_Type::Generate_Matrix;
    }
    else
    {
        GlobalC::exx_global.info.hybrid_type = Exx_Global::Hybrid_Type::No;
    }

    if (GlobalC::exx_global.info.hybrid_type != Exx_Global::Hybrid_Type::No)
    {
        //EXX case, convert all EXX related variables 
        GlobalC::exx_global.info.hybrid_alpha = INPUT.exx_hybrid_alpha;
        XC_Functional::get_hybrid_alpha(INPUT.exx_hybrid_alpha);
        GlobalC::exx_global.info.hse_omega = INPUT.exx_hse_omega;
        GlobalC::exx_global.info.separate_loop = INPUT.exx_separate_loop;
        GlobalC::exx_global.info.hybrid_step = INPUT.exx_hybrid_step;
        GlobalC::exx_lip.info.lambda = INPUT.exx_lambda;
        GlobalC::exx_lcao.info.pca_threshold = INPUT.exx_pca_threshold;
        GlobalC::exx_lcao.info.c_threshold = INPUT.exx_c_threshold;
        GlobalC::exx_lcao.info.v_threshold = INPUT.exx_v_threshold;
        GlobalC::exx_lcao.info.dm_threshold = INPUT.exx_dm_threshold;
        GlobalC::exx_lcao.info.schwarz_threshold = INPUT.exx_schwarz_threshold;
        GlobalC::exx_lcao.info.cauchy_threshold = INPUT.exx_cauchy_threshold;
        GlobalC::exx_lcao.info.ccp_threshold = INPUT.exx_ccp_threshold;
        GlobalC::exx_lcao.info.ccp_rmesh_times = INPUT.exx_ccp_rmesh_times;
        if (INPUT.exx_distribute_type == "htime")
        {
            GlobalC::exx_lcao.info.distribute_type = Exx_Lcao::Distribute_Type::Htime;
        }
        else if (INPUT.exx_distribute_type == "kmeans2")
        {
            GlobalC::exx_lcao.info.distribute_type = Exx_Lcao::Distribute_Type::Kmeans2;
        }
        else if (INPUT.exx_distribute_type == "kmeans1")
        {
            GlobalC::exx_lcao.info.distribute_type = Exx_Lcao::Distribute_Type::Kmeans1;
        }
        else if (INPUT.exx_distribute_type == "order")
        {
            GlobalC::exx_lcao.info.distribute_type = Exx_Lcao::Distribute_Type::Order;
        }
        Exx_Abfs::Jle::Lmax = INPUT.exx_opt_orb_lmax;
        Exx_Abfs::Jle::Ecut_exx = INPUT.exx_opt_orb_ecut;
        Exx_Abfs::Jle::tolerence = INPUT.exx_opt_orb_tolerence;

        //EXX does not support any symmetry analyse, force symmetry setting to -1
        ModuleSymmetry::Symmetry::symm_flag = -1;
    }
#endif
#endif
    GlobalC::ppcell.cell_factor = INPUT.cell_factor; // LiuXh add 20180619

    //----------------------------------------------------------
    // main parameters / electrons / spin ( 2/16 )
    //----------------------------------------------------------
    //	electrons::nelup = INPUT.nelup;
    //	electrons::neldw = INPUT.neldw;

    //----------------------------------------------------------
    // occupation (3/3)
    //----------------------------------------------------------
    Occupy::decision(INPUT.occupations, INPUT.smearing_method, INPUT.smearing_sigma);
    //----------------------------------------------------------
    // charge mixing(3/3)
    //----------------------------------------------------------
    GlobalC::CHR.set_mixing(INPUT.mixing_mode,
                            INPUT.mixing_beta,
                            INPUT.mixing_ndim,
                            INPUT.mixing_gg0); // mohan modify 2014-09-27, add mixing_gg0

    //----------------------------------------------------------
    // iteration
    //----------------------------------------------------------
    GlobalV::SCF_NMAX = INPUT.scf_nmax;
    GlobalV::RELAX_NMAX = INPUT.relax_nmax;
    GlobalV::MD_NSTEP = INPUT.mdp.md_nstep;

    //----------------------------------------------------------
    // wavefunction / charge / potential / (2/4)
    //----------------------------------------------------------
    GlobalV::OUT_FREQ_ELEC = INPUT.out_freq_elec;
    GlobalV::OUT_FREQ_ION = INPUT.out_freq_ion;
    GlobalC::pot.init_chg = INPUT.init_chg;
    GlobalC::pot.chg_extrap = INPUT.chg_extrap; // xiaohui modify 2015-02-01
    GlobalC::CHR.out_chg = INPUT.out_chg;
    GlobalC::CHR.nelec = INPUT.nelec;
    GlobalC::pot.out_pot = INPUT.out_pot;
    GlobalC::wf.out_wfc_pw = INPUT.out_wfc_pw;
    GlobalC::wf.out_wfc_r = INPUT.out_wfc_r;
    GlobalC::en.out_dos = INPUT.out_dos;
    GlobalC::en.out_band = INPUT.out_band;
    GlobalC::en.out_proj_band = INPUT.out_proj_band;
#ifdef __LCAO
    Local_Orbital_Charge::out_dm = INPUT.out_dm;
    hsolver::HSolverLCAO::out_mat_hs = INPUT.out_mat_hs;
    hsolver::HSolverLCAO::out_mat_hsR = INPUT.out_mat_hs2; // LiuXh add 2019-07-16
    elecstate::ElecStateLCAO::out_wfc_lcao = INPUT.out_wfc_lcao;
    if (INPUT.calculation == "nscf" && !INPUT.towannier90 && !INPUT.berry_phase)
    {
        elecstate::ElecStateLCAO::need_psi_grid = false;
    }
    if(INPUT.calculation == "test_neighbour" && GlobalV::NPROC>1)
    {
        ModuleBase::WARNING_QUIT("Input_conv", "test_neighbour must be done with 1 processor");
    }
#endif

    GlobalC::en.dos_emin_ev = INPUT.dos_emin_ev;
    GlobalC::en.dos_emax_ev = INPUT.dos_emax_ev;
    GlobalC::en.dos_edelta_ev = INPUT.dos_edelta_ev;
    GlobalC::en.dos_scale = INPUT.dos_scale;
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

//-----------------------------------------------
// caoyu add for DeePKS
//-----------------------------------------------
#ifdef __DEEPKS
    GlobalV::deepks_scf = INPUT.deepks_scf;
    GlobalV::deepks_bandgap = INPUT.deepks_bandgap; // QO added for bandgap label 2021-12-15
    GlobalV::deepks_out_unittest = INPUT.deepks_out_unittest;
    GlobalV::deepks_out_labels = INPUT.deepks_out_labels;
    if (GlobalV::deepks_out_unittest)
    {
        GlobalV::deepks_out_labels = 1;
        GlobalV::deepks_scf = 1;
        if (GlobalV::NPROC > 1)
            ModuleBase::WARNING_QUIT("Input_conv", "generate deepks unittest with only 1 processor");
        if (GlobalV::CAL_FORCE != 1)
            ModuleBase::WARNING_QUIT("Input_conv", "force is required in generating deepks unittest");
        if (GlobalV::CAL_STRESS != 1)
            ModuleBase::WARNING_QUIT("Input_conv", "stress is required in generating deepks unittest");
    }
    if (GlobalV::deepks_scf || GlobalV::deepks_out_labels)
        GlobalV::deepks_setorb = 1;
#endif
    //-----------------------------------------------
    // sunml add for implicit solvation model
    //-----------------------------------------------
    GlobalV::imp_sol = INPUT.imp_sol;
    GlobalV::eb_k = INPUT.eb_k;
    GlobalV::tau = INPUT.tau;
    GlobalV::sigma_k = INPUT.sigma_k;
    GlobalV::nc_k = INPUT.nc_k;

    ModuleBase::timer::tick("Input_Conv", "Convert");
    return;
}
