#include "input_conv.h"

#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_cell/module_symmetry/symmetry.h"
#include "source_cell/unitcell.h"
#include "source_estate/occupy.h"
#include "source_hamilt/module_surchem/surchem.h"
#include "source_hamilt/module_xc/exx_info.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "../module_unk/berryphase.h"
#include "source_io/module_parameter/parameter.h"
#include "source_io/module_restart/restart.h"
#include "source_relax/ions_move_basic.h"
#include "source_relax/lattice_change_basic.h"

#include <algorithm>

#ifdef __EXX
#include "source_lcao/module_ri/exx_abfs_jle.h"
#endif

#include "source_lcao/module_dftu/dftu.h"
#ifdef __LCAO
#include "source_basis/module_ao/orb_read.h"
#include "source_lcao/force_stress_lcao.h"
#include "source_lcao/module_rt/td_info.h"
#endif
#ifdef __PEXSI
#include "source_hsolver/module_pexsi/pexsi_solver.h"
#endif
#ifdef __MPI
#include "source_hsolver/diago_elpa.h"
#include "source_hsolver/diago_elpa_native.h"
#endif

#include "source_base/module_device/device.h"
#include "source_base/timer.h"
#include "source_estate/elecstate_lcao.h"
#include "source_estate/module_pot/efield.h"
#include "source_estate/module_pot/gatefield.h"
#include "source_hsolver/hsolver_lcao.h"
#include "source_hsolver/hsolver_pw.h"
#include "source_relax/bfgs_basic.h"
#include "source_relax/ions_move_cg.h"

void Input_Conv::Convert()
{
    ModuleBase::TITLE("Input_Conv", "Convert");
    ModuleBase::timer::start("Input_Conv", "Convert");
    //----------------------------------------------------------
    // main parameters / electrons / spin ( 10/16 )
    //----------------------------------------------------------

    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "pseudo_dir", PARAM.inp.pseudo_dir);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "orbital_dir", PARAM.inp.orbital_dir);
    // GlobalV::global_pseudo_type = PARAM.inp.pseudo_type;


    GlobalV::KPAR = PARAM.inp.kpar;



#ifdef __LCAO
    Force_Stress_LCAO<double>::force_invalid_threshold_ev = PARAM.inp.force_zero_out;
    Force_Stress_LCAO<std::complex<double>>::force_invalid_threshold_ev = PARAM.inp.force_zero_out;
#endif

    BFGS_Basic::relax_bfgs_w1 = PARAM.inp.relax_bfgs_w1;
    BFGS_Basic::relax_bfgs_w2 = PARAM.inp.relax_bfgs_w2;

    Ions_Move_Basic::relax_bfgs_rmax = PARAM.inp.relax_bfgs_rmax;
    Ions_Move_Basic::relax_bfgs_rmin = PARAM.inp.relax_bfgs_rmin;
    Ions_Move_Basic::relax_bfgs_init = PARAM.inp.relax_bfgs_init;
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

    if (PARAM.inp.dft_plus_u)
    {
        Plus_U::omc = PARAM.inp.omc;
        Plus_U::orbital_corr = PARAM.inp.orbital_corr;
        Plus_U::uramping = PARAM.globalv.uramping;
        Plus_U::mixing_dftu = PARAM.inp.mixing_dftu;
        Plus_U::U = PARAM.globalv.hubbard_u;
        Plus_U::U0 = PARAM.globalv.hubbard_u;
        if (PARAM.globalv.uramping > 0.01)
        {
            ModuleBase::GlobalFunc::ZEROS(Plus_U::U.data(), PARAM.inp.ntype);
        }
    }

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
    TD_info::out_current = PARAM.inp.out_current;
    TD_info::out_current_k = PARAM.inp.out_current_k;
    TD_info::out_vecpot = PARAM.inp.out_vecpot;
    TD_info::init_vecpot_file = PARAM.inp.init_vecpot_file;
    TD_info::out_mat_R = PARAM.inp.out_hsr[0] == 1 || PARAM.inp.out_hsr[0] == 3 || PARAM.inp.out_hsr_npz_compat;
#endif // __LCAO



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
        }
        else if ( dft_functional_lower == "muller" || dft_functional_lower == "power"
            || dft_functional_lower == "wp22"
            || dft_functional_lower == "cwp22" ) // added by jghan, 2024-07-07
        {
            GlobalC::restart.info_save.save_charge = true;
            GlobalC::restart.info_save.save_H = true;
        }
        else {
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
            || dft_functional_lower == "scan0"
            || dft_functional_lower == "lc_pbe"
            || dft_functional_lower == "lc_wpbe"
            || dft_functional_lower == "lrc_wpbe"
            || dft_functional_lower == "lrc_wpbeh"
            || dft_functional_lower == "cam_pbeh") {
            GlobalC::restart.info_load.load_charge = true;
            GlobalC::restart.info_load.load_H = true;
        }
        else if ( dft_functional_lower == "muller" || dft_functional_lower == "power"
            || dft_functional_lower == "wp22"
            || dft_functional_lower == "cwp22" ) // added by jghan, 2024-07-07
        {
            GlobalC::restart.info_load.load_charge = true;
            GlobalC::restart.info_load.load_H = true;
        }
        else {
            GlobalC::restart.info_load.load_charge = true;
        }
    }

//----------------------------------------------------------
// about exx, Peize Lin add 2018-06-20
//----------------------------------------------------------
    // Initialize a local Exx_Info from input parameters. The ESolver will
    // own its own Exx_Info copy (initialized the same way); the global
    // GlobalC::exx_info has been removed.
    Exx_Info local_exx_info;
    const bool generate_opt_orb = init_exx_info(local_exx_info, PARAM.inp);

    // Set XC_Functional global parameters from the local Exx_Info.
    if (local_exx_info.info_global.cal_exx
#ifdef __EXX
        || generate_opt_orb
        || PARAM.inp.rpa
#endif
        )
    {
        XC_Functional::set_hybrid_alpha(local_exx_info.info_global.hybrid_alpha);
        XC_Functional::set_hse_omega(local_exx_info.info_global.hse_omega);
    }

    // Local aliases: keep this PR's global-state reference budget non-increasing.
    const auto& inp = PARAM.inp;
    const bool cal_exx = local_exx_info.info_global.cal_exx;

    if (cal_exx && inp.basis_type == "pw")
    {
        if (ModuleSymmetry::Symmetry::symm_flag != -1)
        {
            ModuleBase::WARNING("Input_Conv", "EXX PW works only with symmetry=-1");
            ModuleSymmetry::Symmetry::symm_flag = -1;
        }

        if (inp.nspin != 1 && inp.nspin != 2)
        {
            ModuleBase::WARNING_QUIT("Input_Conv", "EXX PW works only with nspin=1 and 2");
        }

        if (inp.cal_stress)
        {
            // Stress_PW::stress_exx only sums same-pool (ik, iq) pairs without
            // the same-spin restriction used in the EXX energy, so the result
            // is wrong for nspin = 2 or kpar > 1.
            if (inp.nspin != 1)
            {
                ModuleBase::WARNING_QUIT("Input_Conv", "EXX PW stress supports only nspin = 1");
            }
            if (inp.kpar != 1)
            {
                ModuleBase::WARNING_QUIT("Input_Conv",
                                         "EXX PW stress does not support k-point parallelism (kpar > 1)");
            }
            if (inp.device == "gpu")
            {
                ModuleBase::WARNING_QUIT("Input_Conv", "EXX PW stress is not supported on GPU");
            }
        }
    }

    if (cal_exx && inp.basis_type == "lcao_in_pw" && inp.cal_stress)
    {
        // For lcao_in_pw the EXX energy comes from Exx_Lip, but Stress_PW
        // would evaluate the EXX stress with the pure PW formula.
        ModuleBase::WARNING_QUIT("Input_Conv", "EXX stress is not supported for basis_type = lcao_in_pw");
    }

    //----------------------------------------------------------
    // reset symmetry flag to avoid error
    //----------------------------------------------------------
    // In these case, symmetry should be reset to 0
    // efield does not support symmetry=1
    if (PARAM.inp.efield_flag && ModuleSymmetry::Symmetry::symm_flag == 1)
    {
        ModuleSymmetry::Symmetry::symm_flag = 0;
    }
    // end of symmetry reset

    //----------------------------------------------------------
    // main parameters / electrons / spin ( 2/16 )
    //----------------------------------------------------------
    //    electrons::nelup = PARAM.inp.nelup;
    //    electrons::neldw = PARAM.inp.neldw;

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
    //    ORB.ecutwfc = PARAM.inp.lcao_ecut;
    //    ORB.dk = PARAM.inp.lcao_dk;
    //    ORB.dR = PARAM.inp.lcao_dr;
    //    ORB.Rmax = PARAM.inp.lcao_rmax;

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
    ModuleBase::timer::end("Input_Conv", "Convert");
    return;
}
