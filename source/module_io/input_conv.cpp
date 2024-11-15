#include "module_io/input_conv.h"

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_cell/module_symmetry/symmetry.h"
#include "module_cell/unitcell.h"
#include "module_elecstate/occupy.h"
#include "module_hamilt_general/module_surchem/surchem.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/berryphase.h"
#include "module_parameter/parameter.h"
#include "module_relax/relax_old/ions_move_basic.h"
#include "module_relax/relax_old/lattice_change_basic.h"

#include <algorithm>

#ifdef __EXX
#include "module_ri/exx_abfs-jle.h"
#endif

#ifdef __LCAO
#include "module_basis/module_ao/ORB_read.h"
#include "module_elecstate/potentials/H_TDDFT_pw.h"
#include "module_hamilt_lcao/hamilt_lcaodft/FORCE_STRESS.h"
#include "module_hamilt_lcao/module_dftu/dftu.h"
#include "module_hamilt_lcao/module_tddft/evolve_elec.h"
#include "module_hamilt_lcao/module_tddft/td_velocity.h"
#endif
#ifdef __PEXSI
#include "module_hsolver/module_pexsi/pexsi_solver.h"
#endif
#ifdef __MPI
#include "module_hsolver/diago_elpa.h"
#include "module_hsolver/diago_elpa_native.h"
#endif

#include "module_base/module_device/device.h"
#include "module_base/timer.h"
#include "module_elecstate/elecstate_lcao.h"
#include "module_elecstate/potentials/efield.h"
#include "module_elecstate/potentials/gatefield.h"
#include "module_hsolver/hsolver_lcao.h"
#include "module_hsolver/hsolver_pw.h"
#include "module_md/md_func.h"

#ifdef __LCAO
std::vector<double> Input_Conv::convert_units(std::string params, double c) {
    std::vector<double> params_ori;
    std::vector<double> params_out;
    parse_expression(params, params_ori);
    for (auto param: params_ori)
        params_out.emplace_back(param * c);

    return params_out;
}

void Input_Conv::read_td_efield()
{
    elecstate::H_TDDFT_pw::stype = PARAM.inp.td_stype;
    if (PARAM.inp.esolver_type == "tddft" && elecstate::H_TDDFT_pw::stype == 1)
    {
        TD_Velocity::tddft_velocity = true;
    } else {
        TD_Velocity::tddft_velocity = false;
    }
    if (PARAM.inp.out_mat_hs2 == 1)
    {
        TD_Velocity::out_mat_R = true;
    } else {
        TD_Velocity::out_mat_R = false;
    }
    parse_expression(PARAM.inp.td_ttype, elecstate::H_TDDFT_pw::ttype);

    elecstate::H_TDDFT_pw::tstart = PARAM.inp.td_tstart;
    elecstate::H_TDDFT_pw::tend = PARAM.inp.td_tend;

    elecstate::H_TDDFT_pw::dt = PARAM.mdp.md_dt / ModuleBase::AU_to_FS;
    elecstate::H_TDDFT_pw::dt_int = elecstate::H_TDDFT_pw::dt;

    // space domain parameters

    // length gauge
    elecstate::H_TDDFT_pw::lcut1 = PARAM.inp.td_lcut1;
    elecstate::H_TDDFT_pw::lcut2 = PARAM.inp.td_lcut2;

    // time domain parameters

    // Gauss
    elecstate::H_TDDFT_pw::gauss_omega = convert_units(PARAM.inp.td_gauss_freq,
                                                       2 * ModuleBase::PI * ModuleBase::AU_to_FS); // time(a.u.)^-1
    elecstate::H_TDDFT_pw::gauss_phase = convert_units(PARAM.inp.td_gauss_phase, 1.0);
    elecstate::H_TDDFT_pw::gauss_sigma = convert_units(PARAM.inp.td_gauss_sigma, 1 / ModuleBase::AU_to_FS);
    elecstate::H_TDDFT_pw::gauss_t0 = convert_units(PARAM.inp.td_gauss_t0, 1.0);
    elecstate::H_TDDFT_pw::gauss_amp = convert_units(PARAM.inp.td_gauss_amp,
                                                     ModuleBase::BOHR_TO_A / ModuleBase::Ry_to_eV); // Ry/bohr
    // init ncut for velocity gauge integral
    for (auto omega: elecstate::H_TDDFT_pw::gauss_omega) {
        int ncut
            = int(100.0 * omega * elecstate::H_TDDFT_pw::dt / ModuleBase::PI);
        if (ncut % 2 == 0) {
            ncut += 2;
        } else {
            ncut += 1;
        }
        if (elecstate::H_TDDFT_pw::stype == 0)
            ncut = 1;
        elecstate::H_TDDFT_pw::gauss_ncut.push_back(ncut);
    }
    // trapezoid
    elecstate::H_TDDFT_pw::trape_omega = convert_units(PARAM.inp.td_trape_freq,
                                                       2 * ModuleBase::PI * ModuleBase::AU_to_FS); // time(a.u.)^-1
    elecstate::H_TDDFT_pw::trape_phase = convert_units(PARAM.inp.td_trape_phase, 1.0);
    elecstate::H_TDDFT_pw::trape_t1 = convert_units(PARAM.inp.td_trape_t1, 1.0);
    elecstate::H_TDDFT_pw::trape_t2 = convert_units(PARAM.inp.td_trape_t2, 1.0);
    elecstate::H_TDDFT_pw::trape_t3 = convert_units(PARAM.inp.td_trape_t3, 1.0);
    elecstate::H_TDDFT_pw::trape_amp = convert_units(PARAM.inp.td_trape_amp,
                                                     ModuleBase::BOHR_TO_A / ModuleBase::Ry_to_eV); // Ry/bohr
    // init ncut for velocity gauge integral
    for (auto omega: elecstate::H_TDDFT_pw::trape_omega) {
        int ncut
            = int(100.0 * omega * elecstate::H_TDDFT_pw::dt / ModuleBase::PI);
        if (ncut % 2 == 0) {
            ncut += 2;
        } else {
            ncut += 1;
        }
        if (elecstate::H_TDDFT_pw::stype == 0)
            ncut = 1;
        elecstate::H_TDDFT_pw::trape_ncut.push_back(ncut);
    }
    // Trigonometric
    elecstate::H_TDDFT_pw::trigo_omega1 = convert_units(PARAM.inp.td_trigo_freq1,
                                                        2 * ModuleBase::PI * ModuleBase::AU_to_FS); // time(a.u.)^-1
    elecstate::H_TDDFT_pw::trigo_omega2 = convert_units(PARAM.inp.td_trigo_freq2,
                                                        2 * ModuleBase::PI * ModuleBase::AU_to_FS); // time(a.u.)^-1
    elecstate::H_TDDFT_pw::trigo_phase1 = convert_units(PARAM.inp.td_trigo_phase1, 1.0);
    elecstate::H_TDDFT_pw::trigo_phase2 = convert_units(PARAM.inp.td_trigo_phase2, 1.0);
    elecstate::H_TDDFT_pw::trigo_amp = convert_units(PARAM.inp.td_trigo_amp,
                                                     ModuleBase::BOHR_TO_A / ModuleBase::Ry_to_eV); // Ry/bohr
    // init ncut for velocity gauge integral
    for (auto omega: elecstate::H_TDDFT_pw::trigo_omega1) {
        int ncut
            = int(100.0 * omega * elecstate::H_TDDFT_pw::dt / ModuleBase::PI);
        if (ncut % 2 == 0) {
            ncut += 2;
        } else {
            ncut += 1;
        }
        if (elecstate::H_TDDFT_pw::stype == 0)
            ncut = 1;
        elecstate::H_TDDFT_pw::trigo_ncut.push_back(ncut);
    }
    // Heaviside
    elecstate::H_TDDFT_pw::heavi_t0 = convert_units(PARAM.inp.td_heavi_t0, 1.0);
    elecstate::H_TDDFT_pw::heavi_amp = convert_units(PARAM.inp.td_heavi_amp,
                                                     ModuleBase::BOHR_TO_A / ModuleBase::Ry_to_eV); // Ry/bohr

    return;
}
#endif

void Input_Conv::Convert()
{
    ModuleBase::TITLE("Input_Conv", "Convert");
    ModuleBase::timer::tick("Input_Conv", "Convert");
    //----------------------------------------------------------
    // main parameters / electrons / spin ( 10/16 )
    //----------------------------------------------------------

    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "pseudo_dir", PARAM.inp.pseudo_dir);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "orbital_dir", PARAM.inp.orbital_dir);
    // GlobalV::global_pseudo_type = PARAM.inp.pseudo_type;
    GlobalC::ucell.setup(PARAM.inp.latname,
                         PARAM.inp.ntype,
                         PARAM.inp.lmaxmax,
                         PARAM.inp.init_vel,
                         PARAM.inp.fixed_axes);

    if (PARAM.inp.calculation == "relax" || PARAM.inp.calculation == "cell-relax")
    {
    }


    if (PARAM.inp.device  == "gpu" && PARAM.inp.basis_type == "pw")
    {
        GlobalV::KPAR = base_device::information::get_device_kpar(PARAM.inp.kpar, PARAM.inp.bndpar);
    }
#ifdef __LCAO
    else if (PARAM.inp.basis_type == "lcao") {
        /// GlobalV::KPAR_LCAO is used in LCAO diagonalization only
        GlobalV::KPAR_LCAO = PARAM.inp.kpar;
        /// all other parts of the code use GlobalV::KPAR = 1
        GlobalV::KPAR = 1;
    }
#endif
    else
    {
        GlobalV::KPAR = PARAM.inp.kpar;
    }
    if (PARAM.inp.device  == "cpu" and PARAM.inp.precision == "single")
    {
// cpu single precision is not supported while float_fftw lib is not available
#ifndef __ENABLE_FLOAT_FFTW
        ModuleBase::WARNING_QUIT(
            "Input_Conv",
            "Single precision with cpu is not supported while float_fftw lib is not available; \
            \n Please recompile with cmake flag \"-DENABLE_FLOAT_FFTW=ON\".\n");
#endif // __ENABLE_FLOAT_FFTW
    }


#ifdef __LCAO
    Force_Stress_LCAO<double>::force_invalid_threshold_ev = PARAM.inp.force_thr_ev2;
    Force_Stress_LCAO<std::complex<double>>::force_invalid_threshold_ev = PARAM.inp.force_thr_ev2;
#endif

    BFGS_Basic::relax_bfgs_w1 = PARAM.inp.relax_bfgs_w1;
    BFGS_Basic::relax_bfgs_w2 = PARAM.inp.relax_bfgs_w2;

    Ions_Move_Basic::relax_bfgs_rmax = PARAM.inp.relax_bfgs_rmax;
    Ions_Move_Basic::relax_bfgs_rmin = PARAM.inp.relax_bfgs_rmin;
    Ions_Move_Basic::relax_bfgs_init = PARAM.inp.relax_bfgs_init;
    Ions_Move_Basic::out_stru = PARAM.inp.out_stru; // mohan add 2012-03-23
    Ions_Move_Basic::relax_method = PARAM.inp.relax_method;
    Lattice_Change_Basic::fixed_axes = PARAM.inp.fixed_axes;


    Ions_Move_CG::RELAX_CG_THR = PARAM.inp.relax_cg_thr; // pengfei add 2013-09-09

    ModuleSymmetry::Symmetry::symm_flag = std::stoi(PARAM.inp.symmetry);
    ModuleSymmetry::Symmetry::symm_autoclose = PARAM.inp.symmetry_autoclose;

    //----------------------------------------------------------
    // planewave (8/8)
    //----------------------------------------------------------

    //----------------------------------------------------------
    // diagonalization  (5/5)
    //----------------------------------------------------------


    //----------------------------------------------------------
    // iteration (1/3)
    //----------------------------------------------------------

#ifdef __LCAO
    if (PARAM.inp.dft_plus_u)
    {
        GlobalC::dftu.Yukawa = PARAM.inp.yukawa_potential;
        GlobalC::dftu.omc = PARAM.inp.omc;
        GlobalC::dftu.orbital_corr = PARAM.inp.orbital_corr;
        GlobalC::dftu.uramping = PARAM.globalv.uramping;
        GlobalC::dftu.mixing_dftu = PARAM.inp.mixing_dftu;
        GlobalC::dftu.U = PARAM.globalv.hubbard_u;
        GlobalC::dftu.U0 = PARAM.globalv.hubbard_u;
        if (PARAM.globalv.uramping > 0.01)
        {
            ModuleBase::GlobalFunc::ZEROS(GlobalC::dftu.U.data(), GlobalC::ucell.ntype);
        }
    }
#endif

    //----------------------------------------------------------
    // Yu Liu add 2022-05-18
    //----------------------------------------------------------
    elecstate::Efield::efield_dir = PARAM.inp.efield_dir;
    elecstate::Efield::efield_pos_max = PARAM.inp.efield_pos_max;
    elecstate::Efield::efield_pos_dec = PARAM.inp.efield_pos_dec;
    elecstate::Efield::efield_amp = PARAM.inp.efield_amp;

    //----------------------------------------------------------
    // Yu Liu add 2022-09-13
    //----------------------------------------------------------
    elecstate::Gatefield::zgate = PARAM.inp.zgate;
    elecstate::Gatefield::relax = PARAM.inp.relax;
    elecstate::Gatefield::block = PARAM.inp.block;
    elecstate::Gatefield::block_down = PARAM.inp.block_down;
    elecstate::Gatefield::block_up = PARAM.inp.block_up;
    elecstate::Gatefield::block_height = PARAM.inp.block_height;

//----------------------------------------------------------
// Fuxiang He add 2016-10-26
//----------------------------------------------------------
#ifdef __LCAO
    module_tddft::Evolve_elec::td_force_dt = PARAM.inp.td_force_dt;
    module_tddft::Evolve_elec::td_vext = PARAM.inp.td_vext;
    if (module_tddft::Evolve_elec::td_vext)
    {
        parse_expression(PARAM.inp.td_vext_dire, module_tddft::Evolve_elec::td_vext_dire_case);
    }
    module_tddft::Evolve_elec::out_dipole = PARAM.inp.out_dipole;
    module_tddft::Evolve_elec::out_efield = PARAM.inp.out_efield;
    module_tddft::Evolve_elec::td_print_eij = PARAM.inp.td_print_eij;
    module_tddft::Evolve_elec::td_edm = PARAM.inp.td_edm;
    TD_Velocity::out_current = PARAM.inp.out_current;
    TD_Velocity::out_current_k = PARAM.inp.out_current_k;
    TD_Velocity::out_vecpot = PARAM.inp.out_vecpot;
    TD_Velocity::init_vecpot_file = PARAM.inp.init_vecpot_file;
    read_td_efield();
#endif

   

    //----------------------------------------------------------
    // about restart, // Peize Lin add 2020-04-04
    //----------------------------------------------------------
    if (PARAM.inp.restart_save)
    {
        std::string dft_functional_lower = PARAM.inp.dft_functional;
        std::transform(PARAM.inp.dft_functional.begin(),
                       PARAM.inp.dft_functional.end(),
                       dft_functional_lower.begin(),
                       tolower);
        GlobalC::restart.folder = PARAM.globalv.global_readin_dir + "restart/";
        ModuleBase::GlobalFunc::MAKE_DIR(GlobalC::restart.folder);
        if (dft_functional_lower == "hf" || dft_functional_lower == "pbe0"
            || dft_functional_lower == "hse"
            || dft_functional_lower == "opt_orb"
            || dft_functional_lower == "scan0") {
            GlobalC::restart.info_save.save_charge = true;
            GlobalC::restart.info_save.save_H = true;
        } else {
            GlobalC::restart.info_save.save_charge = true;
        }
    }
    if (PARAM.inp.restart_load)
    {
        std::string dft_functional_lower = PARAM.inp.dft_functional;
        std::transform(PARAM.inp.dft_functional.begin(),
                       PARAM.inp.dft_functional.end(),
                       dft_functional_lower.begin(),
                       tolower);
        GlobalC::restart.folder = PARAM.globalv.global_readin_dir + "restart/";
        if (dft_functional_lower == "hf" || dft_functional_lower == "pbe0"
            || dft_functional_lower == "hse"
            || dft_functional_lower == "opt_orb"
            || dft_functional_lower == "scan0") {
            GlobalC::restart.info_load.load_charge = true;
            GlobalC::restart.info_load.load_H = true;
        } else {
            GlobalC::restart.info_load.load_charge = true;
        }
    }

//----------------------------------------------------------
// about exx, Peize Lin add 2018-06-20
//----------------------------------------------------------
#ifdef __EXX
#ifdef __LCAO

    std::string dft_functional_lower = PARAM.inp.dft_functional;
    std::transform(PARAM.inp.dft_functional.begin(),
                   PARAM.inp.dft_functional.end(),
                   dft_functional_lower.begin(),
                   tolower);
    if (dft_functional_lower == "hf" || dft_functional_lower == "pbe0" || dft_functional_lower == "scan0")
    {
        GlobalC::exx_info.info_global.cal_exx = true;
        GlobalC::exx_info.info_global.ccp_type
            = Conv_Coulomb_Pot_K::Ccp_Type::Hf;
    } else if (dft_functional_lower == "hse") {
        GlobalC::exx_info.info_global.cal_exx = true;
        GlobalC::exx_info.info_global.ccp_type
            = Conv_Coulomb_Pot_K::Ccp_Type::Hse;
    } else if (dft_functional_lower == "opt_orb") {
        GlobalC::exx_info.info_global.cal_exx = false;
        Exx_Abfs::Jle::generate_matrix = true;
    } else {
        GlobalC::exx_info.info_global.cal_exx = false;
    }

    if (GlobalC::exx_info.info_global.cal_exx || Exx_Abfs::Jle::generate_matrix || PARAM.inp.rpa)
    {
        // EXX case, convert all EXX related variables
        // GlobalC::exx_info.info_global.cal_exx = true;
        GlobalC::exx_info.info_global.hybrid_alpha = std::stod(PARAM.inp.exx_hybrid_alpha);
        XC_Functional::set_hybrid_alpha(std::stod(PARAM.inp.exx_hybrid_alpha));
        GlobalC::exx_info.info_global.hse_omega = PARAM.inp.exx_hse_omega;
        GlobalC::exx_info.info_global.separate_loop = PARAM.inp.exx_separate_loop;
        GlobalC::exx_info.info_global.hybrid_step = PARAM.inp.exx_hybrid_step;
        GlobalC::exx_info.info_global.mixing_beta_for_loop1 = PARAM.inp.exx_mixing_beta;
        GlobalC::exx_info.info_global.exx_symmetry_realspace = PARAM.inp.exx_symmetry_realspace;
        GlobalC::exx_info.info_lip.lambda = PARAM.inp.exx_lambda;

        GlobalC::exx_info.info_ri.real_number = std::stoi(PARAM.inp.exx_real_number);
        GlobalC::exx_info.info_ri.pca_threshold = PARAM.inp.exx_pca_threshold;
        GlobalC::exx_info.info_ri.C_threshold = PARAM.inp.exx_c_threshold;
        GlobalC::exx_info.info_ri.V_threshold = PARAM.inp.exx_v_threshold;
        GlobalC::exx_info.info_ri.dm_threshold = PARAM.inp.exx_dm_threshold;
        GlobalC::exx_info.info_ri.cauchy_threshold = PARAM.inp.exx_cauchy_threshold;
        GlobalC::exx_info.info_ri.C_grad_threshold = PARAM.inp.exx_c_grad_threshold;
        GlobalC::exx_info.info_ri.V_grad_threshold = PARAM.inp.exx_v_grad_threshold;
        GlobalC::exx_info.info_ri.C_grad_R_threshold = PARAM.inp.exx_c_grad_r_threshold;
        GlobalC::exx_info.info_ri.V_grad_R_threshold = PARAM.inp.exx_v_grad_r_threshold;
        GlobalC::exx_info.info_ri.cauchy_force_threshold = PARAM.inp.exx_cauchy_force_threshold;
        GlobalC::exx_info.info_ri.cauchy_stress_threshold = PARAM.inp.exx_cauchy_stress_threshold;
        GlobalC::exx_info.info_ri.ccp_rmesh_times = std::stod(PARAM.inp.exx_ccp_rmesh_times);

        Exx_Abfs::Jle::Lmax = PARAM.inp.exx_opt_orb_lmax;
        Exx_Abfs::Jle::Ecut_exx = PARAM.inp.exx_opt_orb_ecut;
        Exx_Abfs::Jle::tolerence = PARAM.inp.exx_opt_orb_tolerence;

        // EXX does not support symmetry for nspin==4
        if (PARAM.inp.calculation != "nscf" && PARAM.inp.symmetry == "1" && PARAM.inp.nspin == 4)
        {
            ModuleSymmetry::Symmetry::symm_flag = -1;
        }
    }
#endif                                                   // __LCAO
#endif                                                   // __EXX
    GlobalC::ppcell.cell_factor = PARAM.inp.cell_factor; // LiuXh add 20180619

    //----------------------------------------------------------
    // reset symmetry flag to avoid error
    //----------------------------------------------------------
    // In these case, symmetry should be reset to 0
    // efield does not support symmetry=1
    if (PARAM.inp.efield_flag && ModuleSymmetry::Symmetry::symm_flag == 1)
    {
        ModuleSymmetry::Symmetry::symm_flag = 0;
    }
    // In these case, inversion symmetry is also not allowed, symmetry should be
    // reset to -1
    if (PARAM.inp.lspinorb)
    {
        ModuleSymmetry::Symmetry::symm_flag = -1;
    }
    // end of symmetry reset

    //----------------------------------------------------------
    // main parameters / electrons / spin ( 2/16 )
    //----------------------------------------------------------
    //	electrons::nelup = PARAM.inp.nelup;
    //	electrons::neldw = PARAM.inp.neldw;

    //----------------------------------------------------------
    // occupation (3/3)
    //----------------------------------------------------------
    std::string occupations = "smearing";
    Occupy::decision(occupations, PARAM.inp.smearing_method, PARAM.inp.smearing_sigma);

    //----------------------------------------------------------
    // iteration
    //----------------------------------------------------------

    //----------------------------------------------------------
    // wavefunction / charge / potential / (2/4)
    //----------------------------------------------------------

#ifdef __LCAO

    if (PARAM.globalv.gamma_only_local)
    {
        elecstate::ElecStateLCAO<double>::out_wfc_lcao = PARAM.inp.out_wfc_lcao;
    }
    else if (!PARAM.globalv.gamma_only_local)
    {
        elecstate::ElecStateLCAO<std::complex<double>>::out_wfc_lcao = PARAM.inp.out_wfc_lcao;
    }
    if (PARAM.inp.calculation == "nscf" && !PARAM.inp.towannier90 && !PARAM.inp.berry_phase)
    {
        if (PARAM.globalv.gamma_only_local)
        {
            elecstate::ElecStateLCAO<double>::need_psi_grid = false;
        } else if (!PARAM.globalv.gamma_only_local) {
            elecstate::ElecStateLCAO<std::complex<double>>::need_psi_grid
                = false;
        }
    }
    if (PARAM.inp.calculation == "test_neighbour" && GlobalV::NPROC > 1)
    {
        ModuleBase::WARNING_QUIT("Input_conv", "test_neighbour must be done with 1 processor");
    }
#endif

    //----------------------------------------------------------
    // About LCAO
    //----------------------------------------------------------
    // mohan add 2021-04-16
    //	ORB.ecutwfc = PARAM.inp.lcao_ecut;
    //	ORB.dk = PARAM.inp.lcao_dk;
    //	ORB.dR = PARAM.inp.lcao_dr;
    //	ORB.Rmax = PARAM.inp.lcao_rmax;

    // mohan add 2021-02-16
    berryphase::berry_phase_flag = PARAM.inp.berry_phase;

//-----------------------------------------------
// caoyu add for DeePKS
//-----------------------------------------------
    //-----------------------------------------------
    // sunml add for implicit solvation model
    //-----------------------------------------------

    //-----------------------------------------------
    // Deltaspin related parameters
    //-----------------------------------------------

    // mixing parameters

    //-----------------------------------------------
    // Quasiatomic Orbital analysis
    //-----------------------------------------------

    //-----------------------------------------------
    // PEXSI related parameters
    //-----------------------------------------------
#ifdef __PEXSI
    pexsi::PEXSI_Solver::pexsi_npole = PARAM.inp.pexsi_npole;
    pexsi::PEXSI_Solver::pexsi_inertia = PARAM.inp.pexsi_inertia;
    pexsi::PEXSI_Solver::pexsi_nmax = PARAM.inp.pexsi_nmax;
    // pexsi::PEXSI_Solver::pexsi_symbolic = PARAM.inp.pexsi_symbolic;
    pexsi::PEXSI_Solver::pexsi_comm = PARAM.inp.pexsi_comm;
    pexsi::PEXSI_Solver::pexsi_storage = PARAM.inp.pexsi_storage;
    pexsi::PEXSI_Solver::pexsi_ordering = PARAM.inp.pexsi_ordering;
    pexsi::PEXSI_Solver::pexsi_row_ordering = PARAM.inp.pexsi_row_ordering;
    pexsi::PEXSI_Solver::pexsi_nproc = PARAM.inp.pexsi_nproc;
    pexsi::PEXSI_Solver::pexsi_symm = PARAM.inp.pexsi_symm;
    pexsi::PEXSI_Solver::pexsi_trans = PARAM.inp.pexsi_trans;
    pexsi::PEXSI_Solver::pexsi_method = PARAM.inp.pexsi_method;
    pexsi::PEXSI_Solver::pexsi_nproc_pole = PARAM.inp.pexsi_nproc_pole;
    // pexsi::PEXSI_Solver::pexsi_spin = PARAM.inp.pexsi_spin;
    pexsi::PEXSI_Solver::pexsi_temp = PARAM.inp.pexsi_temp;
    pexsi::PEXSI_Solver::pexsi_gap = PARAM.inp.pexsi_gap;
    pexsi::PEXSI_Solver::pexsi_delta_e = PARAM.inp.pexsi_delta_e;
    pexsi::PEXSI_Solver::pexsi_mu_lower = PARAM.inp.pexsi_mu_lower;
    pexsi::PEXSI_Solver::pexsi_mu_upper = PARAM.inp.pexsi_mu_upper;
    pexsi::PEXSI_Solver::pexsi_mu = PARAM.inp.pexsi_mu;
    pexsi::PEXSI_Solver::pexsi_mu_thr = PARAM.inp.pexsi_mu_thr;
    pexsi::PEXSI_Solver::pexsi_mu_expand = PARAM.inp.pexsi_mu_expand;
    pexsi::PEXSI_Solver::pexsi_mu_guard = PARAM.inp.pexsi_mu_guard;
    pexsi::PEXSI_Solver::pexsi_elec_thr = PARAM.inp.pexsi_elec_thr;
    pexsi::PEXSI_Solver::pexsi_zero_thr = PARAM.inp.pexsi_zero_thr;
#endif

    // elpa related
#ifdef __MPI
    hsolver::DiagoElpa<std::complex<double>>::elpa_num_thread = PARAM.inp.elpa_num_thread;
    ;
    hsolver::DiagoElpa<double>::elpa_num_thread = PARAM.inp.elpa_num_thread;
    ;
    hsolver::DiagoElpaNative<std::complex<double>>::elpa_num_thread = PARAM.inp.elpa_num_thread;
    ;
    hsolver::DiagoElpaNative<double>::elpa_num_thread = PARAM.inp.elpa_num_thread;
    ;
#endif
    ModuleBase::timer::tick("Input_Conv", "Convert");
    return;
}
