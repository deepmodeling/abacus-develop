#include "module_io/input_conv.h"

#include <algorithm>

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_cell/module_symmetry/symmetry.h"
#include "module_cell/unitcell.h"
#include "module_elecstate/occupy.h"
#include "module_hamilt_general/module_surchem/surchem.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/berryphase.h"
#include "module_io/input.h"
#include "module_relax/relax_old/ions_move_basic.h"
#include "module_relax/relax_old/lattice_change_basic.h"
#ifdef __EXX
#include "module_ri/exx_abfs-jle.h"
#endif
#ifdef __LCAO
#include "module_basis/module_ao/ORB_read.h"
#include "module_elecstate/potentials/H_TDDFT_pw.h"
#include "module_hamilt_lcao/hamilt_lcaodft/FORCE_STRESS.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_charge.h"
#include "module_hamilt_lcao/module_dftu/dftu.h"
#include "module_hamilt_lcao/module_tddft/evolve_elec.h"
#endif

#include "module_base/timer.h"
#include "module_elecstate/elecstate_lcao.h"
#include "module_elecstate/potentials/efield.h"
#include "module_elecstate/potentials/gatefield.h"
#include "module_hsolver/hsolver_lcao.h"
#include "module_md/md_func.h"
#include "module_psi/kernels/device.h"

template <typename T>
void Input_Conv::parse_expression(const std::string &fn, std::vector<T> &vec)
{
    ModuleBase::TITLE("Input_Conv", "parse_expression");
    int count = 0;
    std::string pattern("([0-9]+\\*[0-9.]+|[0-9,.]+)");
    std::vector<std::string> str;
    std::stringstream ss(fn);
    std::string section;
    while (ss >> section) {
        int index = 0;
        if (str.empty()) {
            while (index < section.size() && std::isspace(section[index])) {
                index++;
            }
        }
        section.erase(0, index);
        str.push_back(section);
    }
    // std::string::size_type pos1, pos2;
    // std::string c = " ";
    // pos2 = fn.find(c);
    // pos1 = 0;
    // while (std::string::npos != pos2)
    // {
    //     str.push_back(fn.substr(pos1, pos2 - pos1));
    //     pos1 = pos2 + c.size();
    //     pos2 = fn.find(c, pos1);
    // }
    // if (pos1 != fn.length())
    // {
    //     str.push_back(fn.substr(pos1));
    // }
    regex_t reg;
    regcomp(&reg, pattern.c_str(), REG_EXTENDED);
    regmatch_t pmatch[1];
    const size_t nmatch = 1;
    for (size_t i = 0; i < str.size(); ++i)
    {
        if (str[i] == "")
        {
            continue;
        }
        int status = regexec(&reg, str[i].c_str(), nmatch, pmatch, 0);
        std::string sub_str = "";
        for (size_t j = pmatch[0].rm_so; j != pmatch[0].rm_eo; ++j)
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
            T occ = stof(sub_str.substr(pos + 1, sub_str.size()));
            // std::vector<double> ocp_temp(num, occ);
            // const std::vector<double>::iterator dest = vec.begin() + count;
            // copy(ocp_temp.begin(), ocp_temp.end(), dest);
            // count += num;
            for (size_t k = 0; k != num; k++)
                vec.emplace_back(occ);
        }
        else
        {
            // vec[count] = stof(sub_str);
            // count += 1;
            std::stringstream convert;
            convert << sub_str;
            T occ;
            convert >> occ;
            vec.emplace_back(occ);
        }
        regfree(&sub_reg);
    }
    regfree(&reg);
}

#ifdef __LCAO
std::vector<double> Input_Conv::convert_units(std::string params, double c)
{
    std::vector<double> params_ori;
    std::vector<double> params_out;
    parse_expression(params, params_ori);
    for (auto param: params_ori)
        params_out.emplace_back(param * c);

    return params_out;
}

void Input_Conv::read_td_efield()
{
    elecstate::H_TDDFT_pw::stype = INPUT.td_stype;

    parse_expression(INPUT.td_ttype, elecstate::H_TDDFT_pw::ttype);

    elecstate::H_TDDFT_pw::tstart = INPUT.td_tstart;
    elecstate::H_TDDFT_pw::tend = INPUT.td_tend;

    elecstate::H_TDDFT_pw::dt = INPUT.mdp.md_dt / ModuleBase::AU_to_FS;

    // space domain parameters

    // length gauge
    elecstate::H_TDDFT_pw::lcut1 = INPUT.td_lcut1;
    elecstate::H_TDDFT_pw::lcut2 = INPUT.td_lcut2;

    // time domain parameters

    // Gauss
    elecstate::H_TDDFT_pw::gauss_omega
        = convert_units(INPUT.td_gauss_freq, 2 * ModuleBase::PI * ModuleBase::AU_to_FS); // time(a.u.)^-1
    elecstate::H_TDDFT_pw::gauss_phase = convert_units(INPUT.td_gauss_phase, 1.0);
    elecstate::H_TDDFT_pw::gauss_sigma = convert_units(INPUT.td_gauss_sigma, 1 / ModuleBase::AU_to_FS);
    elecstate::H_TDDFT_pw::gauss_t0 = convert_units(INPUT.td_gauss_t0, 1.0);
    elecstate::H_TDDFT_pw::gauss_amp
        = convert_units(INPUT.td_gauss_amp, ModuleBase::BOHR_TO_A / ModuleBase::Ry_to_eV); // Ry/bohr

    // trapezoid
    elecstate::H_TDDFT_pw::trape_omega
        = convert_units(INPUT.td_trape_freq, 2 * ModuleBase::PI * ModuleBase::AU_to_FS); // time(a.u.)^-1
    elecstate::H_TDDFT_pw::trape_phase = convert_units(INPUT.td_trape_phase, 1.0);
    elecstate::H_TDDFT_pw::trape_t1 = convert_units(INPUT.td_trape_t1, 1.0);
    elecstate::H_TDDFT_pw::trape_t2 = convert_units(INPUT.td_trape_t2, 1.0);
    elecstate::H_TDDFT_pw::trape_t3 = convert_units(INPUT.td_trape_t3, 1.0);
    elecstate::H_TDDFT_pw::trape_amp
        = convert_units(INPUT.td_trape_amp, ModuleBase::BOHR_TO_A / ModuleBase::Ry_to_eV); // Ry/bohr

    // Trigonometric
    elecstate::H_TDDFT_pw::trigo_omega1
        = convert_units(INPUT.td_trigo_freq1, 2 * ModuleBase::PI * ModuleBase::AU_to_FS); // time(a.u.)^-1
    elecstate::H_TDDFT_pw::trigo_omega2
        = convert_units(INPUT.td_trigo_freq2, 2 * ModuleBase::PI * ModuleBase::AU_to_FS); // time(a.u.)^-1
    elecstate::H_TDDFT_pw::trigo_phase1 = convert_units(INPUT.td_trigo_phase1, 1.0);
    elecstate::H_TDDFT_pw::trigo_phase2 = convert_units(INPUT.td_trigo_phase2, 1.0);
    elecstate::H_TDDFT_pw::trigo_amp
        = convert_units(INPUT.td_trigo_amp, ModuleBase::BOHR_TO_A / ModuleBase::Ry_to_eV); // Ry/bohr

    // Heaviside
    elecstate::H_TDDFT_pw::heavi_t0 = convert_units(INPUT.td_heavi_t0, 1.0);
    elecstate::H_TDDFT_pw::heavi_amp
        = convert_units(INPUT.td_heavi_amp, ModuleBase::BOHR_TO_A / ModuleBase::Ry_to_eV); // Ry/bohr

    return;
}
#endif

void Input_Conv::Convert(void)
{
    ModuleBase::TITLE("Input_Conv", "Convert");
    ModuleBase::timer::tick("Input_Conv", "Convert");
    //-----------------------------------------------
    // set read_file_dir
    //-----------------------------------------------
    if (INPUT.read_file_dir == "auto")
    {
        GlobalV::global_readin_dir = GlobalV::global_out_dir;
    }
    else
    {
        GlobalV::global_readin_dir = INPUT.read_file_dir + '/';
    }
    //----------------------------------------------------------
    // main parameters / electrons / spin ( 10/16 )
    //----------------------------------------------------------
    //  suffix
    if (INPUT.calculation == "md" && INPUT.mdp.md_restart) // md restart  liuyu add 2023-04-12
    {
        int istep = 0;
        MD_func::current_md_info(GlobalV::MY_RANK, GlobalV::global_readin_dir, istep, INPUT.mdp.md_tfirst);
        INPUT.mdp.md_tfirst *= ModuleBase::Hartree_to_K;
        if (INPUT.read_file_dir == "auto")
        {
            GlobalV::stru_file = INPUT.stru_file = GlobalV::global_stru_dir + "STRU_MD_" + std::to_string(istep);
        }
        else
        {
            GlobalV::stru_file = INPUT.stru_file = GlobalV::global_readin_dir + "STRU_MD_" + std::to_string(istep);
        }
    }
    else if (INPUT.stru_file != "")
    {
        GlobalV::stru_file = INPUT.stru_file;
    }
    GlobalV::global_wannier_card = INPUT.wannier_card;
    if (INPUT.kpoint_file != "")
        GlobalV::global_kpoint_card = INPUT.kpoint_file;
    if (INPUT.pseudo_dir != "")
        GlobalV::global_pseudo_dir = INPUT.pseudo_dir + "/";
    if (INPUT.orbital_dir != "")
        GlobalV::global_orbital_dir = INPUT.orbital_dir + "/";
    // GlobalV::global_pseudo_type = INPUT.pseudo_type;
    GlobalC::ucell.setup(INPUT.latname, INPUT.ntype, INPUT.lmaxmax, INPUT.init_vel, INPUT.fixed_axes);

    if (INPUT.calculation == "relax" || INPUT.calculation == "cell-relax")
    {
        if (INPUT.fixed_ibrav && !INPUT.relax_new)
        {
            ModuleBase::WARNING_QUIT("Input_Conv", "fixed_ibrav only available for relax_new = 1");
        }
        if (INPUT.latname == "none" && INPUT.fixed_ibrav)
        {
            ModuleBase::WARNING_QUIT("Input_Conv", "to use fixed_ibrav, latname must be provided");
        }
        if (INPUT.calculation == "relax" && INPUT.fixed_atoms)
        {
            ModuleBase::WARNING_QUIT("Input_Conv", "fixed_atoms is not meant to be used for calculation = relax");
        }
        if (INPUT.relax_new && INPUT.relax_method != "cg")
        {
            INPUT.relax_new = false;
        }
        if (!INPUT.relax_new && (INPUT.fixed_axes == "shape" || INPUT.fixed_axes == "volume"))
        {
            ModuleBase::WARNING_QUIT("Input_Conv", "fixed shape and fixed volume only supported for relax_new = 1");
        }
        GlobalV::fixed_atoms = INPUT.fixed_atoms;
    }

    for(int i=0;i<3;i++)
    {
        GlobalV::KSPACING[i] = INPUT.kspacing[i];
    }
    GlobalV::MIN_DIST_COEF = INPUT.min_dist_coef;
    GlobalV::NBANDS = INPUT.nbands;
    GlobalV::NBANDS_ISTATE = INPUT.nbands_istate;
    GlobalV::device_flag = psi::device::get_device_flag(INPUT.device, INPUT.ks_solver, INPUT.basis_type);

    if (GlobalV::device_flag == "gpu")
    {
        GlobalV::KPAR = psi::device::get_device_kpar(INPUT.kpar);
    }
    else
    {
        GlobalV::KPAR = INPUT.kpar;
        GlobalV::NSTOGROUP = INPUT.bndpar;
    }
    GlobalV::precision_flag = INPUT.precision;
    GlobalV::CALCULATION = INPUT.calculation;
    GlobalV::ESOLVER_TYPE = INPUT.esolver_type;

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

    BFGS_Basic::relax_bfgs_w1 = INPUT.relax_bfgs_w1;
    BFGS_Basic::relax_bfgs_w2 = INPUT.relax_bfgs_w2;

    Ions_Move_Basic::relax_bfgs_rmax = INPUT.relax_bfgs_rmax;
    Ions_Move_Basic::relax_bfgs_rmin = INPUT.relax_bfgs_rmin;
    Ions_Move_Basic::relax_bfgs_init = INPUT.relax_bfgs_init;
    Ions_Move_Basic::out_stru = INPUT.out_stru; // mohan add 2012-03-23
    Lattice_Change_Basic::fixed_axes = INPUT.fixed_axes;

    GlobalV::CAL_STRESS = INPUT.cal_stress;

    GlobalV::RELAX_METHOD = INPUT.relax_method;
    GlobalV::relax_scale_force = INPUT.relax_scale_force;
    GlobalV::relax_new = INPUT.relax_new;

    GlobalV::use_paw = INPUT.use_paw;

    GlobalV::OUT_LEVEL = INPUT.out_level;
    Ions_Move_CG::RELAX_CG_THR = INPUT.relax_cg_thr; // pengfei add 2013-09-09

    ModuleSymmetry::Symmetry::symm_flag = std::stoi(INPUT.symmetry);
    ModuleSymmetry::Symmetry::symm_autoclose = INPUT.symmetry_autoclose;
    GlobalV::BASIS_TYPE = INPUT.basis_type;
    GlobalV::KS_SOLVER = INPUT.ks_solver;
    GlobalV::SEARCH_RADIUS = INPUT.search_radius;
    GlobalV::SEARCH_PBC = INPUT.search_pbc;

    //----------------------------------------------------------
    // planewave (8/8)
    //----------------------------------------------------------
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
    GlobalV::SCF_THR_TYPE = INPUT.scf_thr_type;

#ifdef __LCAO
    if (INPUT.dft_plus_u)
    {
        GlobalV::dft_plus_u = INPUT.dft_plus_u;
        GlobalC::dftu.Yukawa = INPUT.yukawa_potential;
        GlobalC::dftu.omc = INPUT.omc;
        GlobalC::dftu.orbital_corr = INPUT.orbital_corr;
        GlobalC::dftu.mixing_dftu = INPUT.mixing_dftu;
        if (!INPUT.yukawa_potential)
        {
            // Duradev's rotational invariant formulation is implemented
            // where only an effective U given by U-J is used
            // unit is in eV
            GlobalC::dftu.U = INPUT.hubbard_u;
        }
    }
#endif
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
        GlobalV::DOMAG = false;
        GlobalV::DOMAG_Z = true;
        GlobalV::LSPINORB = INPUT.lspinorb;
        GlobalV::soc_lambda = INPUT.soc_lambda;

        if (INPUT.cal_force || INPUT.cal_stress)
        {
            ModuleBase::WARNING_QUIT("input_conv", "force & stress not ready for soc yet!");
        }

        if(INPUT.gamma_only_local)
        {
            ModuleBase::WARNING_QUIT("input_conv", "soc does not support gamma only calculation");
        }
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
    elecstate::Efield::efield_dir = INPUT.efield_dir;
    elecstate::Efield::efield_pos_max = INPUT.efield_pos_max;
    elecstate::Efield::efield_pos_dec = INPUT.efield_pos_dec;
    elecstate::Efield::efield_amp = INPUT.efield_amp;

    //----------------------------------------------------------
    // Yu Liu add 2022-09-13
    //----------------------------------------------------------
    GlobalV::GATE_FLAG = INPUT.gate_flag;
    GlobalV::nelec = INPUT.nelec;
    if (std::abs(INPUT.nupdown) > 1e-6)
    {
        GlobalV::TWO_EFERMI = true;
        GlobalV::nupdown = INPUT.nupdown;
    }
    elecstate::Gatefield::zgate = INPUT.zgate;
    elecstate::Gatefield::relax = INPUT.relax;
    elecstate::Gatefield::block = INPUT.block;
    elecstate::Gatefield::block_down = INPUT.block_down;
    elecstate::Gatefield::block_up = INPUT.block_up;
    elecstate::Gatefield::block_height = INPUT.block_height;

    //----------------------------------------------------------
    // Yu Liu add 2023-05-09
    //----------------------------------------------------------
    INPUT.mdp.force_thr = GlobalV::FORCE_THR;
    INPUT.mdp.my_rank = GlobalV::MY_RANK;
    INPUT.mdp.cal_stress = GlobalV::CAL_STRESS;

//----------------------------------------------------------
// Fuxiang He add 2016-10-26
//----------------------------------------------------------
#ifdef __LCAO
    module_tddft::Evolve_elec::td_force_dt = INPUT.td_force_dt;
    module_tddft::Evolve_elec::td_vext = INPUT.td_vext;
    if (module_tddft::Evolve_elec::td_vext)
    {
        parse_expression(INPUT.td_vext_dire, module_tddft::Evolve_elec::td_vext_dire_case);
    }
    module_tddft::Evolve_elec::out_dipole = INPUT.out_dipole;
    module_tddft::Evolve_elec::out_efield = INPUT.out_efield;
    module_tddft::Evolve_elec::td_print_eij = INPUT.td_print_eij;
    module_tddft::Evolve_elec::td_edm = INPUT.td_edm;
    read_td_efield();
#endif

    // setting for constrained DFT, jiyy add 2020.10.11
    // For example, when we studying nitrogen-vacancy center,
    // it requires an additional excitation of an electron conduction band to simulate the excited state,
    // used for TDDFT only.
    GlobalV::ocp = INPUT.ocp;
    GlobalV::ocp_set = INPUT.ocp_set;
    if (GlobalV::ocp == 1)
    {
        parse_expression(GlobalV::ocp_set, GlobalV::ocp_kb);
    }

    GlobalV::out_mul = INPUT.out_mul; // qifeng add 2019/9/10

    //----------------------------------------------------------
    // about restart, // Peize Lin add 2020-04-04
    //----------------------------------------------------------
    if (INPUT.restart_save)
    {
        std::string dft_functional_lower = INPUT.dft_functional;
        std::transform(INPUT.dft_functional.begin(), INPUT.dft_functional.end(), dft_functional_lower.begin(), tolower);
        GlobalC::restart.folder = GlobalV::global_readin_dir + "restart/";
        ModuleBase::GlobalFunc::MAKE_DIR(GlobalC::restart.folder);
        if (dft_functional_lower == "hf" || dft_functional_lower == "pbe0" || dft_functional_lower == "hse"
            || dft_functional_lower == "opt_orb" || dft_functional_lower == "scan0")
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
        std::string dft_functional_lower = INPUT.dft_functional;
        std::transform(INPUT.dft_functional.begin(), INPUT.dft_functional.end(), dft_functional_lower.begin(), tolower);
        GlobalC::restart.folder = GlobalV::global_readin_dir + "restart/";
        if (dft_functional_lower == "hf" || dft_functional_lower == "pbe0" || dft_functional_lower == "hse"
            || dft_functional_lower == "opt_orb" || dft_functional_lower == "scan0")
        {
            GlobalC::restart.info_load.load_charge = true;
        }
        else
        {
            GlobalC::restart.info_load.load_charge = true;
            GlobalC::restart.info_load.load_H = true;
        }
    }

    if (GlobalV::CALCULATION == "cell-relax" && INPUT.cell_factor < 2.0)
    {
        INPUT.cell_factor = 2.0; // follows QE
    }

//----------------------------------------------------------
// about exx, Peize Lin add 2018-06-20
//----------------------------------------------------------
#ifdef __EXX
#ifdef __LCAO

    std::string dft_functional_lower = INPUT.dft_functional;
    std::transform(INPUT.dft_functional.begin(), INPUT.dft_functional.end(), dft_functional_lower.begin(), tolower);
    if (dft_functional_lower == "hf" || dft_functional_lower == "pbe0" || dft_functional_lower == "scan0")
    {
        GlobalC::exx_info.info_global.cal_exx = true;
        GlobalC::exx_info.info_global.ccp_type = Conv_Coulomb_Pot_K::Ccp_Type::Hf;
    }
    else if (dft_functional_lower == "hse")
    {
        GlobalC::exx_info.info_global.cal_exx = true;
        GlobalC::exx_info.info_global.ccp_type = Conv_Coulomb_Pot_K::Ccp_Type::Hse;
    }
    else if (dft_functional_lower == "opt_orb")
    {
        GlobalC::exx_info.info_global.cal_exx = false;
        Exx_Abfs::Jle::generate_matrix = true;
    }
    else if (INPUT.rpa)
    {
        GlobalC::exx_info.info_global.ccp_type = Conv_Coulomb_Pot_K::Ccp_Type::Hf;
    }
    else
    {
        GlobalC::exx_info.info_global.cal_exx = false;
    }

    if (GlobalC::exx_info.info_global.cal_exx || Exx_Abfs::Jle::generate_matrix || INPUT.rpa)
    {
        // EXX case, convert all EXX related variables
        // GlobalC::exx_info.info_global.cal_exx = true;
        GlobalC::exx_info.info_global.hybrid_alpha = std::stod(INPUT.exx_hybrid_alpha);
        XC_Functional::get_hybrid_alpha(std::stod(INPUT.exx_hybrid_alpha));
        GlobalC::exx_info.info_global.hse_omega = INPUT.exx_hse_omega;
        GlobalC::exx_info.info_global.separate_loop = INPUT.exx_separate_loop;
        GlobalC::exx_info.info_global.hybrid_step = INPUT.exx_hybrid_step;
        GlobalC::exx_info.info_global.mixing_beta_for_loop1 = INPUT.exx_mixing_beta;
        GlobalC::exx_info.info_lip.lambda = INPUT.exx_lambda;

        GlobalC::exx_info.info_ri.real_number = std::stoi(INPUT.exx_real_number);
        GlobalC::exx_info.info_ri.pca_threshold = INPUT.exx_pca_threshold;
        GlobalC::exx_info.info_ri.C_threshold = INPUT.exx_c_threshold;
        GlobalC::exx_info.info_ri.V_threshold = INPUT.exx_v_threshold;
        GlobalC::exx_info.info_ri.dm_threshold = INPUT.exx_dm_threshold;
        GlobalC::exx_info.info_ri.cauchy_threshold = INPUT.exx_cauchy_threshold;
        GlobalC::exx_info.info_ri.C_grad_threshold = INPUT.exx_c_grad_threshold;
        GlobalC::exx_info.info_ri.V_grad_threshold = INPUT.exx_v_grad_threshold;
        GlobalC::exx_info.info_ri.cauchy_force_threshold = INPUT.exx_cauchy_force_threshold;
        GlobalC::exx_info.info_ri.cauchy_stress_threshold = INPUT.exx_cauchy_stress_threshold;
        GlobalC::exx_info.info_ri.ccp_threshold = INPUT.exx_ccp_threshold;
        GlobalC::exx_info.info_ri.ccp_rmesh_times = std::stod(INPUT.exx_ccp_rmesh_times);

        Exx_Abfs::Jle::Lmax = INPUT.exx_opt_orb_lmax;
        Exx_Abfs::Jle::Ecut_exx = INPUT.exx_opt_orb_ecut;
        Exx_Abfs::Jle::tolerence = INPUT.exx_opt_orb_tolerence;

        // EXX does not support symmetry=1
        if (INPUT.calculation != "nscf" && INPUT.symmetry == "1")
            ModuleSymmetry::Symmetry::symm_flag = 0;
    }
#endif // __LCAO
#endif // __EXX
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
    // iteration
    //----------------------------------------------------------
    GlobalV::SCF_NMAX = INPUT.scf_nmax;
    GlobalV::RELAX_NMAX = INPUT.relax_nmax;
    GlobalV::md_prec_level = INPUT.mdp.md_prec_level;

    //----------------------------------------------------------
    // wavefunction / charge / potential / (2/4)
    //----------------------------------------------------------
    GlobalV::OUT_FREQ_ELEC = INPUT.out_freq_elec;
    GlobalV::OUT_FREQ_ION = INPUT.out_freq_ion;
    GlobalV::init_chg = INPUT.init_chg;
    GlobalV::init_wfc = INPUT.init_wfc;
    GlobalV::psi_initializer = INPUT.psi_initializer;
    GlobalV::chg_extrap = INPUT.chg_extrap; // xiaohui modify 2015-02-01
    GlobalV::out_chg = INPUT.out_chg;
    GlobalV::nelec = INPUT.nelec;
    GlobalV::out_pot = INPUT.out_pot;
    GlobalV::out_app_flag = INPUT.out_app_flag;

    GlobalV::out_bandgap = INPUT.out_bandgap; // QO added for bandgap printing
    GlobalV::out_interval = INPUT.out_interval;
#ifdef __LCAO
    Local_Orbital_Charge::out_dm = INPUT.out_dm;
    Local_Orbital_Charge::out_dm1 = INPUT.out_dm1;
    hsolver::HSolverLCAO::out_mat_hs = INPUT.out_mat_hs;
    hsolver::HSolverLCAO::out_mat_hsR = INPUT.out_mat_hs2; // LiuXh add 2019-07-16
    hsolver::HSolverLCAO::out_mat_t = INPUT.out_mat_t;
    hsolver::HSolverLCAO::out_mat_dh = INPUT.out_mat_dh;
    if (GlobalV::GAMMA_ONLY_LOCAL)
    {
        elecstate::ElecStateLCAO<double>::out_wfc_lcao = INPUT.out_wfc_lcao;
    }
    else if (!GlobalV::GAMMA_ONLY_LOCAL)
    {
        elecstate::ElecStateLCAO<std::complex<double>>::out_wfc_lcao = INPUT.out_wfc_lcao;
    }
    if (INPUT.calculation == "nscf" && !INPUT.towannier90 && !INPUT.berry_phase)
    {
        if (GlobalV::GAMMA_ONLY_LOCAL)
        {
            elecstate::ElecStateLCAO<double>::need_psi_grid = false;
        }
        else if (!GlobalV::GAMMA_ONLY_LOCAL)
        {
            elecstate::ElecStateLCAO<std::complex<double>>::need_psi_grid = false;
        }
    }
    if (INPUT.calculation == "test_neighbour" && GlobalV::NPROC > 1)
    {
        ModuleBase::WARNING_QUIT("Input_conv", "test_neighbour must be done with 1 processor");
    }
#endif

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
#else
    if (INPUT.deepks_scf || INPUT.deepks_out_labels || INPUT.deepks_bandgap)
    {
        ModuleBase::WARNING_QUIT("Input_conv", "please compile with DeePKS");
    }
#endif
    //-----------------------------------------------
    // sunml add for implicit solvation model
    //-----------------------------------------------
    GlobalV::imp_sol = INPUT.imp_sol;
    GlobalV::eb_k = INPUT.eb_k;
    GlobalV::tau = INPUT.tau;
    GlobalV::sigma_k = INPUT.sigma_k;
    GlobalV::nc_k = INPUT.nc_k;

    //-----------------------------------------------
    // sunliang add for ofdft 2022-05-11
    //-----------------------------------------------
    GlobalV::of_kinetic = INPUT.of_kinetic;
    GlobalV::of_method = INPUT.of_method;
    GlobalV::of_conv = INPUT.of_conv;
    GlobalV::of_tole = INPUT.of_tole;
    GlobalV::of_tolp = INPUT.of_tolp;
    GlobalV::of_tf_weight = INPUT.of_tf_weight;
    GlobalV::of_vw_weight = INPUT.of_vw_weight;
    GlobalV::of_wt_alpha = INPUT.of_wt_alpha;
    GlobalV::of_wt_beta = INPUT.of_wt_beta;
    GlobalV::of_wt_rho0 = INPUT.of_wt_rho0;
    GlobalV::of_hold_rho0 = INPUT.of_hold_rho0;
    GlobalV::of_lkt_a = INPUT.of_lkt_a;
    GlobalV::of_full_pw = INPUT.of_full_pw;
    GlobalV::of_full_pw_dim = INPUT.of_full_pw_dim;
    GlobalV::of_read_kernel = INPUT.of_read_kernel;
    GlobalV::of_kernel_file = INPUT.of_kernel_file;

    //-----------------------------------------------
    // NCSCD related parameters
    //-----------------------------------------------
    GlobalV::sc_mag_switch = INPUT.sc_mag_switch;
    GlobalV::sc_thr = INPUT.sc_thr;
    GlobalV::nsc = INPUT.nsc;
    GlobalV::nsc_min = INPUT.nsc_min;
    GlobalV::sc_file = INPUT.sc_file;

    ModuleBase::timer::tick("Input_Conv", "Convert");
    return;
}
