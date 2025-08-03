#include "source_io/input_conv.h"

#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_cell/module_symmetry/symmetry.h"
#include "source_cell/unitcell.h"
#include "source_estate/occupy.h"
#include "source_hamilt/module_surchem/surchem.h"
#include "source_pw/module_pwdft/global.h"
#include "source_io/berryphase.h"
#include "source_io/module_parameter/parameter.h"
#include "source_relax/ions_move_basic.h"
#include "source_relax/lattice_change_basic.h"

#include <algorithm>

#ifdef __EXX
#include "source_lcao/module_ri/exx_abfs-jle.h"
#endif

#include "source_lcao/module_dftu/dftu.h"
#ifdef __LCAO
#include "source_basis/module_ao/ORB_read.h"
#include "source_estate/module_pot/H_TDDFT_pw.h"
#include "source_lcao/FORCE_STRESS.h"
#include "source_lcao/module_rt/evolve_elec.h"
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
#include "source_md/md_func.h"

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
    if (PARAM.inp.out_mat_hs2 == 1)
    {
        TD_info::out_mat_R = true;
    } else {
        TD_info::out_mat_R = false;
    }
    parse_expression(PARAM.inp.td_ttype, elecstate::H_TDDFT_pw::ttype);

    elecstate::H_TDDFT_pw::tstart = PARAM.inp.td_tstart;
    elecstate::H_TDDFT_pw::tend = PARAM.inp.td_tend;
    if(PARAM.inp.td_dt!=-1.0)
    {
        elecstate::H_TDDFT_pw::dt = PARAM.inp.td_dt / ModuleBase::AU_to_FS;
    }
    else
    {
        elecstate::H_TDDFT_pw::dt = PARAM.mdp.md_dt / PARAM.inp.estep_per_md / ModuleBase::AU_to_FS;
    }

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
            ModuleBase::GlobalFunc::ZEROS(GlobalC::dftu.U.data(), PARAM.inp.ntype);
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
    read_td_efield();
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
    std::string dft_functional_lower = PARAM.inp.dft_functional;
    std::transform(PARAM.inp.dft_functional.begin(),
                   PARAM.inp.dft_functional.end(),
                   dft_functional_lower.begin(),
                   tolower);
    if (dft_functional_lower == "hf"
    || dft_functional_lower == "pbe0" || dft_functional_lower == "b3lyp" || dft_functional_lower == "hse"
    || dft_functional_lower == "scan0"
    || dft_functional_lower == "muller" || dft_functional_lower == "power"
    || dft_functional_lower == "cwp22" || dft_functional_lower == "wp22" 
    || dft_functional_lower == "lc_pbe"
    || dft_functional_lower == "lc_wpbe" 
    || dft_functional_lower == "lrc_wpbe"
    || dft_functional_lower == "lrc_wpbeh"
    || dft_functional_lower == "cam_pbeh")
    {
        GlobalC::exx_info.info_global.cal_exx = true;

        GlobalC::exx_info.info_global.hybrid_alpha = 0;
        std::vector<double> fock_alpha(PARAM.inp.exx_fock_alpha.size());
        for(std::size_t i=0; i<fock_alpha.size(); ++i)
        {
            fock_alpha[i] = std::stod(PARAM.inp.exx_fock_alpha[i]);
            GlobalC::exx_info.info_global.hybrid_alpha = std::max(std::abs(fock_alpha[i]), GlobalC::exx_info.info_global.hybrid_alpha);
        }
        std::vector<double> erfc_alpha(PARAM.inp.exx_erfc_alpha.size());
        for(std::size_t i=0; i<erfc_alpha.size(); ++i)
        {
            erfc_alpha[i] = std::stod(PARAM.inp.exx_erfc_alpha[i]);
            GlobalC::exx_info.info_global.hybrid_alpha = std::max(std::abs(erfc_alpha[i]), GlobalC::exx_info.info_global.hybrid_alpha);
        }
        assert(GlobalC::exx_info.info_global.hybrid_alpha>0);
        for(std::size_t i=0; i<fock_alpha.size(); ++i)
            { fock_alpha[i] /= GlobalC::exx_info.info_global.hybrid_alpha; }
        for(std::size_t i=0; i<erfc_alpha.size(); ++i)
            { erfc_alpha[i] /= GlobalC::exx_info.info_global.hybrid_alpha; }        

        if(!fock_alpha.empty())
        {
            if(PARAM.inp.basis_type == "lcao")
            {
                GlobalC::exx_info.info_global.coulomb_param[Conv_Coulomb_Pot_K::Coulomb_Type::Fock].resize(fock_alpha.size());
                for(std::size_t i=0; i<fock_alpha.size(); ++i)
                {
                    GlobalC::exx_info.info_global.coulomb_param[Conv_Coulomb_Pot_K::Coulomb_Type::Fock] = {{
                        {"alpha", ModuleBase::GlobalFunc::TO_STRING(fock_alpha[i])},
                        {"singularity_correction", PARAM.inp.exx_singularity_correction} }};
                }
            }
            else if(PARAM.inp.basis_type == "lcao_in_pw")
            {
                assert(fock_alpha.size() == PARAM.inp.exx_fock_lambda.size());
                GlobalC::exx_info.info_global.coulomb_param[Conv_Coulomb_Pot_K::Coulomb_Type::Fock].resize(fock_alpha.size());
                for(std::size_t i=0; i<fock_alpha.size(); ++i)
                {
                    GlobalC::exx_info.info_global.coulomb_param[Conv_Coulomb_Pot_K::Coulomb_Type::Fock] = {{
                        {"alpha", ModuleBase::GlobalFunc::TO_STRING(fock_alpha[i])},
                        {"lambda", PARAM.inp.exx_fock_lambda[i]} }};
                }
            }
            else if(PARAM.inp.basis_type == "pw")
            {
                GlobalC::exx_info.info_global.coulomb_param[Conv_Coulomb_Pot_K::Coulomb_Type::Fock].resize(fock_alpha.size());
                for(std::size_t i=0; i<fock_alpha.size(); ++i)
                {
                    GlobalC::exx_info.info_global.coulomb_param[Conv_Coulomb_Pot_K::Coulomb_Type::Fock] = {{
                        {"alpha", ModuleBase::GlobalFunc::TO_STRING(fock_alpha[i])} }};
                }
            }
            else
            {
                throw std::invalid_argument(std::string(__FILE__)+" line "+std::to_string(__LINE__));
            }
        }
        if(!erfc_alpha.empty())
        {
            assert(erfc_alpha.size() == PARAM.inp.exx_erfc_omega.size());
            if(PARAM.inp.basis_type == "lcao")
            {
                GlobalC::exx_info.info_global.coulomb_param[Conv_Coulomb_Pot_K::Coulomb_Type::Erfc].resize(erfc_alpha.size());
                for(std::size_t i=0; i<erfc_alpha.size(); ++i)
                {
                    GlobalC::exx_info.info_global.coulomb_param[Conv_Coulomb_Pot_K::Coulomb_Type::Erfc] = {{
                        {"alpha", ModuleBase::GlobalFunc::TO_STRING(erfc_alpha[i])},
                        {"omega", ModuleBase::GlobalFunc::TO_STRING(PARAM.inp.exx_erfc_omega[i])},
                        {"singularity_correction", PARAM.inp.exx_singularity_correction} }};
                }
            }
            else if(PARAM.inp.basis_type == "pw")
            {
                GlobalC::exx_info.info_global.coulomb_param[Conv_Coulomb_Pot_K::Coulomb_Type::Erfc].resize(erfc_alpha.size());
                for(std::size_t i=0; i<erfc_alpha.size(); ++i)
                {
                    GlobalC::exx_info.info_global.coulomb_param[Conv_Coulomb_Pot_K::Coulomb_Type::Erfc] = {{
                        {"alpha", ModuleBase::GlobalFunc::TO_STRING(erfc_alpha[i])},
                        {"omega", ModuleBase::GlobalFunc::TO_STRING(PARAM.inp.exx_erfc_omega[i])} }};
                }
            }
        }
    }
    else
    {
        GlobalC::exx_info.info_global.cal_exx = false;
    }

    // info_global.ccp_type will be removed in the future. these codes for pw and lcao_in_pw temporarily
    if (dft_functional_lower == "hf"
     || dft_functional_lower == "pbe0" || dft_functional_lower == "b3lyp"
     || dft_functional_lower == "scan0"
     || dft_functional_lower == "muller" || dft_functional_lower == "power")
    {
        GlobalC::exx_info.info_global.ccp_type = Conv_Coulomb_Pot_K::Ccp_Type::Hf;
    }
    // use the error function erf(w|r-r'|), exx just has the short-range part
    else if (dft_functional_lower == "hse"
          || dft_functional_lower == "cwp22")
    {
        GlobalC::exx_info.info_global.ccp_type = Conv_Coulomb_Pot_K::Ccp_Type::Erfc;
    }
    // use the error function erf(w|r-r'|), exx just has the long-range part
    else if ( dft_functional_lower == "wp22" )
    {
        GlobalC::exx_info.info_global.ccp_type = Conv_Coulomb_Pot_K::Ccp_Type::Erf;
    }

    if (GlobalC::exx_info.info_global.cal_exx
#ifdef __EXX
        || PARAM.inp.rpa
#endif
        )
    {
        // EXX case, convert all EXX related variables
        XC_Functional::set_hybrid_alpha(GlobalC::exx_info.info_global.hybrid_alpha);
        if(!PARAM.inp.exx_erfc_omega.empty())
            { GlobalC::exx_info.info_global.hse_omega = std::stod(PARAM.inp.exx_erfc_omega[0]); }
        if(!PARAM.inp.exx_fock_lambda.empty())
            { GlobalC::exx_info.info_lip.lambda = std::stod(PARAM.inp.exx_fock_lambda[0]); }
        GlobalC::exx_info.info_global.separate_loop = PARAM.inp.exx_separate_loop;
        GlobalC::exx_info.info_global.hybrid_step = PARAM.inp.exx_hybrid_step;
        GlobalC::exx_info.info_global.mixing_beta_for_loop1 = PARAM.inp.exx_mixing_beta;

        GlobalC::exx_info.info_ri.real_number = std::stoi(PARAM.inp.exx_real_number);
        GlobalC::exx_info.info_ri.pca_threshold = PARAM.inp.exx_pca_threshold;
        GlobalC::exx_info.info_ri.C_threshold = PARAM.inp.exx_c_threshold;
        GlobalC::exx_info.info_ri.V_threshold = PARAM.inp.exx_v_threshold;
        GlobalC::exx_info.info_ri.dm_threshold = PARAM.inp.exx_dm_threshold;
        GlobalC::exx_info.info_ri.C_grad_threshold = PARAM.inp.exx_c_grad_threshold;
        GlobalC::exx_info.info_ri.V_grad_threshold = PARAM.inp.exx_v_grad_threshold;
        GlobalC::exx_info.info_ri.C_grad_R_threshold = PARAM.inp.exx_c_grad_r_threshold;
        GlobalC::exx_info.info_ri.V_grad_R_threshold = PARAM.inp.exx_v_grad_r_threshold;
        GlobalC::exx_info.info_ri.ccp_rmesh_times = std::stod(PARAM.inp.exx_ccp_rmesh_times);
        GlobalC::exx_info.info_ri.exx_symmetry_realspace = PARAM.inp.exx_symmetry_realspace;

        GlobalC::exx_info.info_opt_abfs.abfs_Lmax  = PARAM.inp.exx_opt_orb_lmax;
        GlobalC::exx_info.info_opt_abfs.ecut_exx = PARAM.inp.exx_opt_orb_ecut;
        GlobalC::exx_info.info_opt_abfs.tolerence = PARAM.inp.exx_opt_orb_tolerence;

        // EXX does not support symmetry for nspin==4
        if (PARAM.inp.calculation != "nscf" && PARAM.inp.symmetry == "1" && PARAM.inp.nspin == 4 && PARAM.inp.basis_type == "lcao")
        {
            ModuleSymmetry::Symmetry::symm_flag = -1;
        }
    }

    if (GlobalC::exx_info.info_global.cal_exx && PARAM.inp.basis_type == "pw")
    {
        if (ModuleSymmetry::Symmetry::symm_flag != -1)
        {
            ModuleBase::WARNING_QUIT("Input_Conv", "EXX PW works only with symmetry=-1");
            // ModuleSymmetry::Symmetry::symm_flag = -1;
        }

        if (PARAM.inp.nspin != 1 && PARAM.inp.nspin != 2)
        {
            ModuleBase::WARNING_QUIT("Input_Conv", "EXX PW works only with nspin=1 and 2");
        }

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
    ModuleBase::timer::tick("Input_Conv", "Convert");
    return;
}
