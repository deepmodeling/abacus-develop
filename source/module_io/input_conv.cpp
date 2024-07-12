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
#include "module_parameter/parameter.h"
#include "module_relax/relax_old/ions_move_basic.h"
#include "module_relax/relax_old/lattice_change_basic.h"

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
#endif

#include "module_base/module_device/device.h"
#include "module_base/timer.h"
#include "module_elecstate/elecstate_lcao.h"
#include "module_elecstate/potentials/efield.h"
#include "module_elecstate/potentials/gatefield.h"
#include "module_hsolver/hsolver_lcao.h"
#include "module_hsolver/hsolver_pw.h"
#include "module_md/md_func.h"

template <typename T>
void Input_Conv::parse_expression(const std::string& fn, std::vector<T>& vec)
{
    ModuleBase::TITLE("Input_Conv", "parse_expression");
    int count = 0;
    std::string pattern("([0-9]+\\*[0-9.]+|[0-9,.]+)");
    std::vector<std::string> str;
    std::stringstream ss(fn);
    std::string section;
    while (ss >> section)
    {
        int index = 0;
        if (str.empty())
        {
            while (index < section.size() && std::isspace(section[index]))
            {
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
            for (size_t k = 0; k != num; k++) {
                vec.emplace_back(occ);
}
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
    elecstate::H_TDDFT_pw::stype = PARAM.inp.td_stype;
    if (INPUT.esolver_type == "tddft" && elecstate::H_TDDFT_pw::stype == 1)
    {
        TD_Velocity::tddft_velocity = true;
    }
    else
    {
        TD_Velocity::tddft_velocity = false;
    }
    if (PARAM.inp.out_mat_hs2 == 1)
    {
        TD_Velocity::out_mat_R = true;
    }
    else
    {
        TD_Velocity::out_mat_R = false;
    }
    parse_expression(PARAM.inp.td_ttype, elecstate::H_TDDFT_pw::ttype);

    elecstate::H_TDDFT_pw::tstart = PARAM.inp.td_tstart;
    elecstate::H_TDDFT_pw::tend = PARAM.inp.td_tend;

    elecstate::H_TDDFT_pw::dt = INPUT.mdp.md_dt / ModuleBase::AU_to_FS;
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
    for (auto omega: elecstate::H_TDDFT_pw::gauss_omega)
    {
        int ncut = int(100.0 * omega * elecstate::H_TDDFT_pw::dt / ModuleBase::PI);
        if (ncut % 2 == 0)
        {
            ncut += 2;
        }
        else
        {
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
    for (auto omega: elecstate::H_TDDFT_pw::trape_omega)
    {
        int ncut = int(100.0 * omega * elecstate::H_TDDFT_pw::dt / ModuleBase::PI);
        if (ncut % 2 == 0)
        {
            ncut += 2;
        }
        else
        {
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
    for (auto omega: elecstate::H_TDDFT_pw::trigo_omega1)
    {
        int ncut = int(100.0 * omega * elecstate::H_TDDFT_pw::dt / ModuleBase::PI);
        if (ncut % 2 == 0)
        {
            ncut += 2;
        }
        else
        {
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
    GlobalV::CALCULATION = PARAM.inp.sup.global_calculation;
    GlobalV::double_grid = PARAM.inp.sup.double_grid;
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
    {
        GlobalV::global_kpoint_card = INPUT.kpoint_file;
    }
    if (INPUT.pseudo_dir != "")
    {
        GlobalV::global_pseudo_dir = INPUT.pseudo_dir + "/";
    }
    if (INPUT.orbital_dir != "")
    {
        GlobalV::global_orbital_dir = INPUT.orbital_dir + "/";
    }
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "pseudo_dir", GlobalV::global_pseudo_dir);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "orbital_dir", GlobalV::global_orbital_dir);
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
            ModuleBase::WARNING_QUIT("Input_Conv",
                                     "fixed shape and fixed volume only "
                                     "supported for relax_new = 1");
        }
        GlobalV::fixed_atoms = INPUT.fixed_atoms;
    }

    for (int i = 0; i < 3; i++)
    {
        GlobalV::KSPACING[i] = INPUT.kspacing[i];
    }
    GlobalV::MIN_DIST_COEF = INPUT.min_dist_coef;
    GlobalV::NBANDS = INPUT.nbands;
    GlobalV::NBANDS_ISTATE = INPUT.nbands_istate;

    GlobalV::device_flag = base_device::information::get_device_flag(INPUT.device,
                                                                     INPUT.ks_solver,
                                                                     INPUT.basis_type,
                                                                     INPUT.gamma_only_local);

    if (GlobalV::device_flag == "gpu" && INPUT.basis_type == "pw")
    {
        GlobalV::KPAR = base_device::information::get_device_kpar(INPUT.kpar);
    }
    else
    {
        GlobalV::KPAR = INPUT.kpar;
        GlobalV::NSTOGROUP = INPUT.bndpar;
    }
    GlobalV::precision_flag = INPUT.precision;
    if (GlobalV::device_flag == "cpu" and GlobalV::precision_flag == "single")
    {
// cpu single precision is not supported while float_fftw lib is not available
#ifndef __ENABLE_FLOAT_FFTW
        ModuleBase::WARNING_QUIT("Input_Conv",
                                 "Single precision with cpu is not supported while float_fftw lib is not available; \
            \n Please recompile with cmake flag \"-DENABLE_FLOAT_FFTW=ON\".\n");
#endif // __ENABLE_FLOAT_FFTW
    }
    GlobalV::CALCULATION = INPUT.calculation;
    GlobalV::ESOLVER_TYPE = INPUT.esolver_type;

    GlobalV::PSEUDORCUT = INPUT.pseudo_rcut;
    GlobalV::PSEUDO_MESH = INPUT.pseudo_mesh;

    GlobalV::DFT_FUNCTIONAL = INPUT.dft_functional;
    GlobalV::XC_TEMPERATURE = INPUT.xc_temperature;
    GlobalV::NSPIN = INPUT.nspin;

    GlobalV::CAL_FORCE = INPUT.cal_force;
    GlobalV::FORCE_THR = INPUT.force_thr;

    GlobalV::STRESS_THR = INPUT.stress_thr;
    GlobalV::PRESS1 = INPUT.press1;
    GlobalV::PRESS2 = INPUT.press2;
    GlobalV::PRESS3 = INPUT.press3;
    GlobalV::out_element_info = INPUT.out_element_info;
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
    Lattice_Change_Basic::fixed_axes = INPUT.fixed_axes;

    GlobalV::CAL_STRESS = INPUT.cal_stress;

    GlobalV::NUM_STREAM = INPUT.nstream;

    GlobalV::RELAX_METHOD = INPUT.relax_method;
    GlobalV::relax_scale_force = INPUT.relax_scale_force;
    GlobalV::relax_new = INPUT.relax_new;

    GlobalV::use_paw = INPUT.use_paw;

    GlobalV::OUT_LEVEL = INPUT.out_level;
    Ions_Move_CG::RELAX_CG_THR = INPUT.relax_cg_thr; // pengfei add 2013-09-09

    ModuleSymmetry::Symmetry::symm_flag = std::stoi(PARAM.inp.symmetry);
    ModuleSymmetry::Symmetry::symm_autoclose = PARAM.inp.symmetry_autoclose;
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

    hsolver::HSolverPW<std::complex<float>, base_device::DEVICE_CPU>::diago_full_acc = PARAM.inp.diago_full_acc;
    hsolver::HSolverPW<std::complex<double>, base_device::DEVICE_CPU>::diago_full_acc = PARAM.inp.diago_full_acc;

#if ((defined __CUDA) || (defined __ROCM))
    hsolver::HSolverPW<std::complex<float>, base_device::DEVICE_GPU>::diago_full_acc = PARAM.inp.diago_full_acc;
    hsolver::HSolverPW<std::complex<double>, base_device::DEVICE_GPU>::diago_full_acc = PARAM.inp.diago_full_acc;
#endif

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
        GlobalC::dftu.uramping = INPUT.uramping;
        GlobalC::dftu.mixing_dftu = PARAM.inp.mixing_dftu;
        if (INPUT.yukawa_potential && INPUT.hubbard_u == nullptr)
        {
            // Duradev's rotational invariant formulation is implemented
            // where only an effective U given by U-J is used
            // unit is in eV
            INPUT.hubbard_u = new double[GlobalC::ucell.ntype];
        }
        GlobalC::dftu.U = INPUT.hubbard_u;
        GlobalC::dftu.U0 = std::vector<double>(INPUT.hubbard_u, INPUT.hubbard_u + GlobalC::ucell.ntype);
        if (INPUT.uramping > 0.01)
        {
            ModuleBase::GlobalFunc::ZEROS(GlobalC::dftu.U, GlobalC::ucell.ntype);
        }
    }
    GlobalV::onsite_radius = INPUT.onsite_radius;
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
        // set the domag variable to make a spin-orbit calculation with zero
        // magnetization
        GlobalV::DOMAG = false;
        GlobalV::DOMAG_Z = true;
        GlobalV::LSPINORB = INPUT.lspinorb;
        GlobalV::soc_lambda = INPUT.soc_lambda;
        if (INPUT.gamma_only_local)
        {
            ModuleBase::WARNING_QUIT("input_conv",
                                     "nspin=4(soc or noncollinear-spin) does "
                                     "not support gamma only calculation");
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
    GlobalV::EFIELD_FLAG = PARAM.inp.efield_flag;
    GlobalV::DIP_COR_FLAG = PARAM.inp.dip_cor_flag;
    elecstate::Efield::efield_dir = PARAM.inp.efield_dir;
    elecstate::Efield::efield_pos_max = PARAM.inp.efield_pos_max;
    elecstate::Efield::efield_pos_dec = PARAM.inp.efield_pos_dec;
    elecstate::Efield::efield_amp = PARAM.inp.efield_amp;

    //----------------------------------------------------------
    // Yu Liu add 2022-09-13
    //----------------------------------------------------------
    GlobalV::GATE_FLAG = PARAM.inp.gate_flag;
    GlobalV::nelec = PARAM.inp.nelec;
    if (PARAM.inp.sup.two_fermi)
    {
        GlobalV::TWO_EFERMI = true;
        GlobalV::nupdown = PARAM.inp.nupdown;
    }
    elecstate::Gatefield::zgate = PARAM.inp.zgate;
    elecstate::Gatefield::relax = PARAM.inp.relax;
    elecstate::Gatefield::block = PARAM.inp.block;
    elecstate::Gatefield::block_down = PARAM.inp.block_down;
    elecstate::Gatefield::block_up = PARAM.inp.block_up;
    elecstate::Gatefield::block_height = PARAM.inp.block_height;

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

    // setting for constrained DFT, jiyy add 2020.10.11
    // For example, when we studying nitrogen-vacancy center,
    // it requires an additional excitation of an electron conduction band to
    // simulate the excited state, used for TDDFT only.
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
            GlobalC::restart.info_load.load_H = true;
        }
        else
        {
            GlobalC::restart.info_load.load_charge = true;
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
    else
    {
        GlobalC::exx_info.info_global.cal_exx = false;
    }

    if (GlobalC::exx_info.info_global.cal_exx || Exx_Abfs::Jle::generate_matrix || INPUT.rpa)
    {
        // EXX case, convert all EXX related variables
        // GlobalC::exx_info.info_global.cal_exx = true;
        GlobalC::exx_info.info_global.hybrid_alpha = std::stod(PARAM.inp.exx_hybrid_alpha);
        XC_Functional::get_hybrid_alpha(std::stod(PARAM.inp.exx_hybrid_alpha));
        GlobalC::exx_info.info_global.hse_omega = PARAM.inp.exx_hse_omega;
        GlobalC::exx_info.info_global.separate_loop = PARAM.inp.exx_separate_loop;
        GlobalC::exx_info.info_global.hybrid_step = PARAM.inp.exx_hybrid_step;
        GlobalC::exx_info.info_global.mixing_beta_for_loop1 = PARAM.inp.exx_mixing_beta;
        GlobalC::exx_info.info_lip.lambda = PARAM.inp.exx_lambda;

        GlobalC::exx_info.info_ri.real_number = std::stoi(PARAM.inp.exx_real_number);
        GlobalC::exx_info.info_ri.pca_threshold = PARAM.inp.exx_pca_threshold;
        GlobalC::exx_info.info_ri.C_threshold = PARAM.inp.exx_c_threshold;
        GlobalC::exx_info.info_ri.V_threshold = PARAM.inp.exx_v_threshold;
        GlobalC::exx_info.info_ri.dm_threshold = PARAM.inp.exx_dm_threshold;
        GlobalC::exx_info.info_ri.cauchy_threshold = PARAM.inp.exx_cauchy_threshold;
        GlobalC::exx_info.info_ri.C_grad_threshold = PARAM.inp.exx_c_grad_threshold;
        GlobalC::exx_info.info_ri.V_grad_threshold = PARAM.inp.exx_v_grad_threshold;
        GlobalC::exx_info.info_ri.cauchy_force_threshold = PARAM.inp.exx_cauchy_force_threshold;
        GlobalC::exx_info.info_ri.cauchy_stress_threshold = PARAM.inp.exx_cauchy_stress_threshold;
        GlobalC::exx_info.info_ri.ccp_rmesh_times = std::stod(PARAM.inp.exx_ccp_rmesh_times);

        Exx_Abfs::Jle::Lmax = PARAM.inp.exx_opt_orb_lmax;
        Exx_Abfs::Jle::Ecut_exx = PARAM.inp.exx_opt_orb_ecut;
        Exx_Abfs::Jle::tolerence = PARAM.inp.exx_opt_orb_tolerence;

        // EXX does not support symmetry=1
        if (INPUT.calculation != "nscf" && PARAM.inp.symmetry == "1") {
            ModuleSymmetry::Symmetry::symm_flag = 0;
}
    }
#endif                                               // __LCAO
#endif                                               // __EXX
    GlobalC::ppcell.cell_factor = INPUT.cell_factor; // LiuXh add 20180619

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
    if (GlobalV::LSPINORB)
    {
        ModuleSymmetry::Symmetry::symm_flag = -1;
    }
    // end of symmetry reset

    //----------------------------------------------------------
    // main parameters / electrons / spin ( 2/16 )
    //----------------------------------------------------------
    //	electrons::nelup = INPUT.nelup;
    //	electrons::neldw = INPUT.neldw;

    //----------------------------------------------------------
    // occupation (3/3)
    //----------------------------------------------------------
    std::string occupations = "smearing";
    Occupy::decision(occupations, INPUT.smearing_method, INPUT.smearing_sigma);

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
    GlobalV::nelec = PARAM.inp.nelec;
    GlobalV::nelec_delta = PARAM.inp.nelec_delta;
    GlobalV::out_pot = INPUT.out_pot;
    GlobalV::out_app_flag = INPUT.out_app_flag;
    GlobalV::out_ndigits = INPUT.out_ndigits;

    GlobalV::out_bandgap = INPUT.out_bandgap; // QO added for bandgap printing
    GlobalV::out_interval = INPUT.out_interval;
#ifdef __LCAO
    hsolver::HSolverLCAO<double>::out_mat_hs = PARAM.inp.out_mat_hs;
    hsolver::HSolverLCAO<double>::out_mat_hsR = PARAM.inp.out_mat_hs2; // LiuXh add 2019-07-16
    hsolver::HSolverLCAO<double>::out_mat_t = PARAM.inp.out_mat_t;
    hsolver::HSolverLCAO<double>::out_mat_dh = PARAM.inp.out_mat_dh;
    hsolver::HSolverLCAO<std::complex<double>>::out_mat_hs = PARAM.inp.out_mat_hs;
    hsolver::HSolverLCAO<std::complex<double>>::out_mat_hsR = PARAM.inp.out_mat_hs2; // LiuXh add 2019-07-16
    hsolver::HSolverLCAO<std::complex<double>>::out_mat_t = PARAM.inp.out_mat_t;
    hsolver::HSolverLCAO<std::complex<double>>::out_mat_dh = PARAM.inp.out_mat_dh;
    if (GlobalV::GAMMA_ONLY_LOCAL)
    {
        elecstate::ElecStateLCAO<double>::out_wfc_lcao = PARAM.inp.out_wfc_lcao;
    }
    else if (!GlobalV::GAMMA_ONLY_LOCAL)
    {
        elecstate::ElecStateLCAO<std::complex<double>>::out_wfc_lcao = PARAM.inp.out_wfc_lcao;
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
    GlobalV::deepks_scf = PARAM.inp.deepks_scf;
    GlobalV::deepks_bandgap = PARAM.inp.deepks_bandgap; // QO added for bandgap label 2021-12-15
    GlobalV::deepks_out_unittest = PARAM.inp.deepks_out_unittest;
    GlobalV::deepks_out_labels = PARAM.inp.deepks_out_labels;
    GlobalV::deepks_equiv = PARAM.inp.deepks_equiv;

    if (GlobalV::deepks_equiv && GlobalV::deepks_bandgap)
    {
        ModuleBase::WARNING_QUIT("Input_conv", "deepks_equiv and deepks_bandgap cannot be used together");
    }
    if (GlobalV::deepks_out_unittest)
    {
        GlobalV::deepks_out_labels = true;
        GlobalV::deepks_scf = true;
        if (GlobalV::NPROC > 1) {
            ModuleBase::WARNING_QUIT("Input_conv", "generate deepks unittest with only 1 processor");
}
        if (GlobalV::CAL_FORCE != 1) {
            ModuleBase::WARNING_QUIT("Input_conv", "force is required in generating deepks unittest");
}
        if (GlobalV::CAL_STRESS != 1) {
            ModuleBase::WARNING_QUIT("Input_conv", "stress is required in generating deepks unittest");
}
    }
    if (GlobalV::deepks_scf || GlobalV::deepks_out_labels) {
        GlobalV::deepks_setorb = true;
}
#else
    if (PARAM.inp.deepks_scf || PARAM.inp.deepks_out_labels || PARAM.inp.deepks_bandgap)
    {
        ModuleBase::WARNING_QUIT("Input_conv", "please compile with DeePKS");
    }
#endif
    //-----------------------------------------------
    // sunml add for implicit solvation model
    //-----------------------------------------------
    GlobalV::imp_sol = PARAM.inp.imp_sol;
    GlobalV::eb_k = PARAM.inp.eb_k;
    GlobalV::tau = PARAM.inp.tau;
    GlobalV::sigma_k = PARAM.inp.sigma_k;
    GlobalV::nc_k = PARAM.inp.nc_k;

    //-----------------------------------------------
    // Deltaspin related parameters
    //-----------------------------------------------
    GlobalV::sc_mag_switch = PARAM.inp.sc_mag_switch;
    GlobalV::decay_grad_switch = PARAM.inp.decay_grad_switch;
    GlobalV::sc_thr = PARAM.inp.sc_thr;
    GlobalV::nsc = PARAM.inp.nsc;
    GlobalV::nsc_min = PARAM.inp.nsc_min;
    GlobalV::sc_scf_nmin = PARAM.inp.sc_scf_nmin;
    GlobalV::alpha_trial = PARAM.inp.alpha_trial;
    GlobalV::sccut = PARAM.inp.sccut;
    GlobalV::sc_file = PARAM.inp.sc_file;

    // mixing parameters
    GlobalV::MIXING_MODE = PARAM.inp.mixing_mode;
    GlobalV::MIXING_BETA = PARAM.inp.mixing_beta;
    GlobalV::MIXING_NDIM = PARAM.inp.mixing_ndim;
    GlobalV::MIXING_RESTART = PARAM.inp.mixing_restart;
    GlobalV::MIXING_GG0 = PARAM.inp.mixing_gg0;
    GlobalV::MIXING_BETA_MAG = PARAM.inp.mixing_beta_mag;
    GlobalV::MIXING_GG0_MAG = PARAM.inp.mixing_gg0_mag;
    GlobalV::MIXING_GG0_MIN = PARAM.inp.mixing_gg0_min;
    GlobalV::MIXING_ANGLE = PARAM.inp.mixing_angle;
    GlobalV::MIXING_TAU = PARAM.inp.mixing_tau;
    GlobalV::MIXING_DMR = PARAM.inp.mixing_dmr;

    //-----------------------------------------------
    // Quasiatomic Orbital analysis
    //-----------------------------------------------
    GlobalV::qo_switch = PARAM.inp.qo_switch;
    GlobalV::qo_basis = PARAM.inp.qo_basis;
    GlobalV::qo_strategy = PARAM.inp.qo_strategy;
    GlobalV::qo_thr = PARAM.inp.qo_thr;
    GlobalV::qo_screening_coeff = PARAM.inp.qo_screening_coeff;

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
#endif
    ModuleBase::timer::tick("Input_Conv", "Convert");
    return;
}
