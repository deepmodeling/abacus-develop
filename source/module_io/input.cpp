//==========================================================
// Author: Lixin He,mohan
// DATE : 2008-11-6
//==========================================================
// #include "global.h"
#include "module_io/input.h"

#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include "module_base/constants.h"
#include "module_base/global_file.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/parallel_common.h"
#include "module_base/timer.h"
#include "version.h"
Input INPUT;

void Input::Init(const std::string& fn)
{
    ModuleBase::timer::tick("Input", "Init");
    this->Default();

    // only rank 0 reads the input file, check the code in this->Read()
    bool success = this->Read(fn);

    this->Default_2();

// xiaohui add 2015-09-16
#ifdef __MPI
    Parallel_Common::bcast_bool(input_error);
#endif
    if (input_error == 1)
    {
        ModuleBase::WARNING_QUIT("Input", "Bad parameter, please check the input parameters in file INPUT", 1);
    }

#ifdef __MPI
    Parallel_Common::bcast_bool(success);
#endif
    if (!success)
    {
        ModuleBase::WARNING_QUIT("Input::Init", "Error during readin parameters.", 1);
    }
#ifdef __MPI
    this->Bcast();
#endif

    // mohan move forward 2011-02-26
    //----------------------------------------------------------
    // OTHRE CLASS MEMBER FUNCTION :
    // NAME : Run::make_dir( dir name : OUT.suffix)
    //----------------------------------------------------------
    bool out_dir = false;
    if (!out_app_flag && (out_mat_hs2 || out_mat_r || out_mat_t || out_mat_dh))
        out_dir = true;
    ModuleBase::Global_File::make_dir_out(this->suffix,
                                          this->calculation,
                                          out_dir,
                                          GlobalV::MY_RANK,
                                          this->mdp.md_restart,
                                          this->out_alllog); // xiaohui add 2013-09-01
    Check();
#ifdef VERSION
    const char* version = VERSION;
#else
    const char* version = "unknown";
#endif
#ifdef COMMIT
    const char* commit = COMMIT;
#else
    const char* commit = "unknown";
#endif
    time_t time_now = time(NULL);
    GlobalV::ofs_running << "                                                                                     "
                         << std::endl;
    GlobalV::ofs_running << "                              ABACUS " << version << std::endl << std::endl;
    GlobalV::ofs_running << "               Atomic-orbital Based Ab-initio Computation at UStc                    "
                         << std::endl
                         << std::endl;
    GlobalV::ofs_running << "                     Website: http://abacus.ustc.edu.cn/                             "
                         << std::endl;
    GlobalV::ofs_running << "               Documentation: https://abacus.deepmodeling.com/                       "
                         << std::endl;
    GlobalV::ofs_running << "                  Repository: https://github.com/abacusmodeling/abacus-develop       "
                         << std::endl;
    GlobalV::ofs_running << "                              https://github.com/deepmodeling/abacus-develop         "
                         << std::endl;
    GlobalV::ofs_running << "                      Commit: " << commit << std::endl << std::endl;
    GlobalV::ofs_running << std::setiosflags(std::ios::right);

#ifdef __MPI
    // GlobalV::ofs_running << "    Version: Parallel, under ALPHA test" << std::endl;
    // GlobalV::ofs_running << "    Version: Parallel, in development" << std::endl;
    // GlobalV::ofs_running << "    Processor Number is " << GlobalV::NPROC << std::endl;
    ModuleBase::TITLE("Input", "init");
    ModuleBase::TITLE("Input", "Bcast");
#else
    GlobalV::ofs_running << "    This is SERIES version." << std::endl;
    ModuleBase::TITLE("Input", "init");
#endif
    GlobalV::ofs_running << "    Start Time is " << ctime(&time_now);
    GlobalV::ofs_running << "                                                                                     "
                         << std::endl;
    GlobalV::ofs_running << " ------------------------------------------------------------------------------------"
                         << std::endl;

    GlobalV::ofs_running << std::setiosflags(std::ios::left);
    std::cout << std::setiosflags(std::ios::left);

    GlobalV::ofs_running << "\n READING GENERAL INFORMATION" << std::endl;
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "global_out_dir", GlobalV::global_out_dir);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "global_in_card", GlobalV::global_in_card);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "pseudo_dir", GlobalV::global_pseudo_dir);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "orbital_dir", GlobalV::global_orbital_dir);

    // ModuleBase::GlobalFunc::OUT(
    //     GlobalV::ofs_running,
    //     "pseudo_type",
    //     pseudo_type); // mohan add 2013-05-20 (xiaohui add 2013-06-23, GlobalV::global_pseudo_type -> pseudo_type)

    ModuleBase::timer::tick("Input", "Init");
    return;
}

void Input::Default(void)
{
    ModuleBase::TITLE("Input", "Default");
    //----------------------------------------------------------
    // main parameters
    //----------------------------------------------------------
    // xiaohui modify 2015-03-25
    // suffix = "MESIA";
    suffix = "ABACUS";
    stru_file = "";   // xiaohui modify 2015-02-01
    kpoint_file = ""; // xiaohui modify 2015-02-01
    pseudo_dir = "";
    orbital_dir = ""; // liuyu add 2021-08-14
    read_file_dir = "auto";
    // pseudo_type = "auto"; // mohan add 2013-05-20 (xiaohui add 2013-06-23)
    wannier_card = "none";
    latname = "none";
    // xiaohui modify 2015-09-15, relax -> scf
    // calculation = "relax";
    calculation = "scf";
    esolver_type = "ksdft";
    pseudo_rcut = 15.0;  // qianrui add this parameter 2021-5
    pseudo_mesh = false; // qianrui add this pararmeter
    ntype = 0;
    nbands = 0;
    nbands_sto = 256;
    nbndsto_str = "256";
    nbands_istate = 5;
    pw_seed = 1;
    emin_sto = 0.0;
    emax_sto = 0.0;
    nche_sto = 100;
    seed_sto = 0;
    bndpar = 1;
    kpar = 1;
    initsto_freq = 0;
    method_sto = 2;
    npart_sto = 1;
    cal_cond = false;
    dos_nche = 100;
    cond_che_thr = 1e-8;
    cond_dw = 0.1;
    cond_wcut = 10;
    cond_dt = 0.02;
    cond_dtbatch = 0;
    cond_smear = 1;
    cond_fwhm = 0.4;
    cond_nonlocal = true;
    berry_phase = false;
    gdir = 3;
    towannier90 = false;
    nnkpfile = "seedname.nnkp";
    wannier_spin = "up";
    wannier_method = 1;
    out_wannier_amn = true;
    out_wannier_eig = true;
    out_wannier_mmn = true;
    out_wannier_unk = false;
    out_wannier_wvfn_formatted = true;
    for (int i = 0; i < 3; i++)
    {
        kspacing[i] = 0;
    }
    min_dist_coef = 0.2;
    //----------------------------------------------------------
    // electrons / spin
    //----------------------------------------------------------
    dft_functional = "default";
    xc_temperature = 0.0;
    nspin = 1;
    nelec = 0.0;
    lmaxmax = 2;
    //----------------------------------------------------------
    // new function
    //----------------------------------------------------------
    basis_type = "pw";     // xiaohui add 2013-09-01
    ks_solver = "default"; // xiaohui add 2013-09-01
    search_radius = -1.0;  // unit: a.u. -1.0 has no meaning.
    search_pbc = true;
    symmetry = "default";
    init_vel = false;
    ref_cell_factor = 1.0;
    symmetry_prec = 1.0e-6;    // LiuXh add 2021-08-12, accuracy for symmetry
    symmetry_autoclose = true; // whether to close symmetry automatically when error occurs in symmetry analysis
    cal_force = 0;
    force_thr = 1.0e-3;
    force_thr_ev2 = 0;
    stress_thr = 0.5; // LiuXh add 20180515 liuyu update 2023-05-10
    press1 = 0.0;
    press2 = 0.0;
    press3 = 0.0;
    cal_stress = false;
    fixed_axes = "None"; // pengfei 2018-11-9
    fixed_ibrav = false;
    fixed_atoms = false;
    relax_method = "cg"; // pengfei  2014-10-13
    relax_cg_thr = 0.5;  // pengfei add 2013-08-15
    out_level = "ie";
    out_md_control = false;
    relax_new = true;
    relax_bfgs_w1 = 0.01; // mohan add 2011-03-13
    relax_bfgs_w2 = 0.5;
    relax_bfgs_rmax = 0.8; // bohr
    relax_bfgs_rmin = 1e-5;
    relax_bfgs_init = 0.5; // bohr
    relax_scale_force = 0.5;
    nbspline = -1;

    use_paw = false;
    //----------------------------------------------------------
    // ecutwfc
    //----------------------------------------------------------
    // gamma_only = false;
    gamma_only = false;
    gamma_only_local = false;
    ecutwfc = 50.0;
    ecutrho = 0.0;
    erf_ecut = 0.0;
    erf_height = 0.0;
    erf_sigma = 0.1;
    ncx = 0;
    ncy = 0;
    ncz = 0;
    nx = 0;
    ny = 0;
    nz = 0;
    bx = 0;
    by = 0;
    bz = 0;
    ndx = 0;
    ndy = 0;
    ndz = 0;
    //----------------------------------------------------------
    // diagonalization
    //----------------------------------------------------------
    // diago_type = "default"; xiaohui modify 2013-09-01 //mohan update 2012-02-06
    diago_proc = 0; // if 0, then diago_proc = GlobalV::NPROC
    pw_diag_nmax = 50;
    diago_cg_prec = 1; // mohan add 2012-03-31
    pw_diag_ndim = 4;
    pw_diag_thr = 1.0e-2;
    nb2d = 0;
    nurse = 0;
    colour = 0;
    t_in_h = 1;
    vl_in_h = 1;
    vnl_in_h = 1;
    vh_in_h = 1;
    vion_in_h = 1;
    test_force = 0;
    test_stress = 0;
    //----------------------------------------------------------
    // iteration
    //----------------------------------------------------------
    scf_thr = -1.0;    // the default value (1e-9 for pw, and 1e-7 for lcao) will be set in Default_2
    scf_thr_type = -1; // the default value (1 for pw, and 2 for lcao) will be set in Default_2
    scf_nmax = 100;
    relax_nmax = 0;
    out_stru = 0;
    //----------------------------------------------------------
    // occupation
    //----------------------------------------------------------
    occupations = "smearing";
    smearing_method = "gauss"; // this setting is based on the report in Issue #2847
    smearing_sigma = 0.015;    // this setting is based on the report in Issue #2847
    //----------------------------------------------------------
    //  charge mixing
    //----------------------------------------------------------
    mixing_mode = "broyden";
    mixing_beta = -10;
    mixing_ndim = 8;
    mixing_gg0 = 1.00;       // use Kerker defaultly
    mixing_beta_mag = -10.0; // only set when nspin == 2 || nspin == 4
    mixing_gg0_mag = 0.0;    // defaultly exclude Kerker from mixing magnetic density
    mixing_gg0_min = 0.1;    // defaultly minimum kerker coefficient
    mixing_angle = -10.0;    // defaultly close for npsin = 4
    mixing_tau = false;
    mixing_dftu = false;
    //----------------------------------------------------------
    // potential / charge / wavefunction / energy
    //----------------------------------------------------------
    init_wfc = "atomic";
    psi_initializer = false;
    mem_saver = 0;
    printe = 100; // must > 0
    init_chg = "atomic";
    chg_extrap = "default"; // xiaohui modify 2015-02-01
    out_freq_elec = 0;
    out_freq_ion = 0;
    out_chg = 0;
    out_dm = 0;
    out_dm1 = 0;

    out_bandgap = 0; // QO added for bandgap printing

    band_print_num = 0;

    deepks_out_labels = 0; // caoyu added 2020-11-24, mohan added 2021-01-03
    deepks_scf = 0;
    deepks_bandgap = 0;
    deepks_out_unittest = 0;

    out_pot = 0;
    out_wfc_pw = 0;
    out_wfc_r = 0;
    out_dos = 0;
    out_band = {0, 8};
    out_proj_band = 0;
    out_mat_hs = {0, 8};
    out_mat_xc = 0;
    cal_syns = 0;
    dmax = 0.01;
    out_mat_hs2 = 0; // LiuXh add 2019-07-15
    out_mat_t = 0;
    out_interval = 1;
    out_app_flag = true;
    out_ndigits = 8;
    out_mat_r = 0; // jingan add 2019-8-14
    out_mat_dh = 0;
    out_wfc_lcao = 0;
    out_alllog = false;
    dos_emin_ev = -15;    //(ev)
    dos_emax_ev = 15;     //(ev)
    dos_edelta_ev = 0.01; //(ev)
    dos_scale = 0.01;
    dos_sigma = 0.07;
    out_element_info = false;
    //----------------------------------------------------------
    // LCAO
    //----------------------------------------------------------
    lcao_ecut = 0; // (Ry)
    lcao_dk = 0.01;
    lcao_dr = 0.01;
    lcao_rmax = 30; // (a.u.)
    //----------------------------------------------------------
    // efield and dipole correction     Yu Liu add 2022-05-18
    //----------------------------------------------------------
    efield_flag = false;
    dip_cor_flag = false;
    efield_dir = 2;
    efield_pos_max = -1.0;
    efield_pos_dec = -1.0;
    efield_amp = 0.0;
    //----------------------------------------------------------
    // gatefield                        Yu Liu add 2022-09-13
    //----------------------------------------------------------
    gate_flag = false;
    zgate = 0.5;
    relax = false;
    block = false;
    block_down = 0.45;
    block_up = 0.55;
    block_height = 0.1;
    //----------------------------------------------------------
    // vdw									//jiyy add 2019-08-04
    //----------------------------------------------------------
    vdw_method = "none";
    vdw_s6 = "default";
    vdw_s8 = "default";
    vdw_a1 = "default";
    vdw_a2 = "default";
    vdw_d = 20;
    vdw_abc = false;
    vdw_cutoff_radius = "default";
    vdw_radius_unit = "Bohr";
    vdw_cn_thr = 40.0;
    vdw_cn_thr_unit = "Bohr";
    vdw_C6_file = "default";
    vdw_C6_unit = "Jnm6/mol";
    vdw_R0_file = "default";
    vdw_R0_unit = "A";
    vdw_cutoff_type = "radius";
    vdw_cutoff_period = {3, 3, 3};

    //----------------------------------------------------------
    // exx										//Peize Lin add 2018-06-20
    //----------------------------------------------------------
    exx_hybrid_alpha = "default";
    exx_hse_omega = 0.11;

    exx_separate_loop = true;
    exx_hybrid_step = 100;
    exx_mixing_beta = 1.0;

    exx_lambda = 0.3;

    exx_real_number = "default";
    exx_pca_threshold = 1E-4;
    exx_c_threshold = 1E-4;
    exx_v_threshold = 1E-1;
    exx_dm_threshold = 1E-4;
    exx_schwarz_threshold = 0;
    exx_cauchy_threshold = 1E-7;
    exx_c_grad_threshold = 1E-4;
    exx_v_grad_threshold = 1E-1;
    exx_cauchy_force_threshold = 1E-7;
    exx_cauchy_stress_threshold = 1E-7;
    exx_ccp_threshold = 1E-8;
    exx_ccp_rmesh_times = "default";

    exx_distribute_type = "htime";

    exx_opt_orb_lmax = 0;
    exx_opt_orb_ecut = 0.0;
    exx_opt_orb_tolerence = 0.0;

    // added by zhengdy-soc
    noncolin = false;
    lspinorb = false;
    soc_lambda = 1.0;

    // xiaohui add 2015-09-16
    input_error = 0;

    //----------------------------------------------------------			//Fuxiang He add 2016-10-26
    // tddft
    //----------------------------------------------------------
    td_force_dt = 0.02;
    td_vext = false;
    td_vext_dire = "1";

    propagator = 0;

    out_dipole = false;
    out_efield = false;

    td_print_eij = -1.0;
    td_edm = 0;

    td_stype = 0; // 0 : length gauge  1: velocity gauge

    td_ttype = '0';
    //  0  Gauss type function.
    //  1  trapezoid type function.
    //  2  Trigonometric function, sin^2.
    //  3  heaviside function.
    //  4  HHG function.

    td_tstart = 1;
    td_tend = 1000;

    // space domain parameters

    // length gauge
    td_lcut1 = 0.05;
    td_lcut2 = 0.95;

    // time domain parameters

    // Gauss
    td_gauss_freq = "22.13"; // fs^-1
    td_gauss_phase = "0.0";
    td_gauss_sigma = "30.0"; // fs
    td_gauss_t0 = "100.0";
    td_gauss_amp = "0.25"; // V/A

    // Trapezoid
    td_trape_freq = "1.60"; // fs^-1
    td_trape_phase = "0.0";
    td_trape_t1 = "1875";
    td_trape_t2 = "5625";
    td_trape_t3 = "7500";
    td_trape_amp = "2.74"; // V/A

    // Trigonometric
    td_trigo_freq1 = "1.164656"; // time(fs)^-1
    td_trigo_freq2 = "0.029116"; // time(fs)^-1
    td_trigo_phase1 = "0.0";
    td_trigo_phase2 = "0.0";
    td_trigo_amp = "2.74"; // V/A

    // Heaviside
    td_heavi_t0 = "100";
    td_heavi_amp = "1.0"; // V/A

    // HHG
    // td_hhg_amp1 = "2.74"; // V/A
    // td_hhg_amp2 = "2.74"; // V/A
    // td_hhg_freq1 = "1.164656"; // time(fs)^-1
    // td_hhg_freq2 = "0.029116"; // time(fs)^-1
    // td_hhg_phase1 = "0.0";
    // td_hhg_phase2 = "0.0";
    // td_hhg_t0 = "700";
    // td_hhg_sigma = "30"; // fs

    //----------------------------------------------------------			//Fuxiang He add 2016-10-26
    // constrained DFT
    //----------------------------------------------------------
    // GlobalV::ocp = 0;
    // ocp_n = 0;
    // GlobalV::ocp_set = "none";
    // for(int i=0; i<10000; i++)
    // {
    // GlobalV::ocp_kb[i] = 0.0;
    // }

    cell_factor = 1.2; // LiuXh add 20180619

    out_mul = false; // qi feng add 2019/9/10

    //----------------------------------------------------------			//Peize Lin add 2020-04-04
    // restart
    //----------------------------------------------------------
    restart_save = false;
    restart_load = false;

    //==========================================================
    // test only
    //==========================================================
    test_skip_ewald = false;

    //==========================================================
    //    DFT+U     Xin Qu added on 2020-10-29
    //==========================================================
    dft_plus_u = false; // 1:DFT+U correction; 0: standard DFT calcullation
    yukawa_potential = false;
    yukawa_lambda = -1.0;
    omc = 0;

    //==========================================================
    //    DFT+DMFT     Xin Qu added on 2020-08
    //==========================================================
    dft_plus_dmft = false;

    //==========================================================
    //    RPA    Rong Shi added on 2022-04
    //==========================================================
    rpa = false;
    coulomb_type = "full";

    //==========================================================
    //    implicit solvation model       sunml added on 2022-04-04
    //==========================================================
    imp_sol = 0;
    eb_k = 80.0;
    tau = 1.0798 * 1e-5;
    sigma_k = 0.6;
    nc_k = 0.00037;

    //==========================================================
    //    OFDFT sunliang added on 2022-05-05
    //==========================================================
    of_kinetic = "wt";
    of_method = "tn";
    of_conv = "energy";
    of_tole = 1e-6;
    of_tolp = 1e-5;
    of_tf_weight = 1.;
    of_vw_weight = 1.;
    of_wt_alpha = 5. / 6.;
    of_wt_beta = 5. / 6.;
    of_wt_rho0 = 0.;
    of_hold_rho0 = false;
    of_lkt_a = 1.3;
    of_full_pw = true;
    of_full_pw_dim = 0;
    of_read_kernel = false;
    of_kernel_file = "WTkernel.txt";

    //==========================================================
    // spherical bessel  Peize Lin added on 2022-12-15
    //==========================================================
    bessel_nao_smooth = true;
    bessel_nao_sigma = 0.1;
    bessel_nao_ecut = "default";
    bessel_nao_rcut = 6.0; // -1.0 for forcing manual setting
    bessel_nao_rcuts = {};
    bessel_nao_tolerence = 1E-12;

    bessel_descriptor_lmax = 2; // -1 for forcing manual setting
    bessel_descriptor_smooth = true;
    bessel_descriptor_sigma = 0.1;
    bessel_descriptor_ecut = "default";
    bessel_descriptor_rcut = 6.0; // -1.0 for forcing manual setting
    bessel_descriptor_tolerence = 1E-12;

    //==========================================================
    //    device control denghui added on 2022-11-15
    //==========================================================
    device = "cpu";
    //==========================================================
    //    precision control denghui added on 2023-01-01
    //==========================================================
    precision = "double";
    //==========================================================
    // variables for deltaspin
    //==========================================================
    sc_mag_switch = 0;
    decay_grad_switch = false;
    sc_thr = 1e-6;
    nsc = 100;
    nsc_min = 2;
    sc_scf_nmin = 2;
    alpha_trial = 0.01;
    sccut = 3.0;
    sc_file = "none";
    //==========================================================
    // variables for Quasiatomic Orbital analysis
    //==========================================================
    qo_switch = false;
    qo_basis = "hydrogen";
    qo_strategy = {};
    qo_thr = 1e-6;
    qo_screening_coeff = {};

    return;
}

bool Input::Read(const std::string& fn)
{
    ModuleBase::TITLE("Input", "Read");

    if (GlobalV::MY_RANK != 0)
        return false;

    std::ifstream ifs(fn.c_str(), std::ios::in); // "in_datas/input_parameters"

    if (!ifs)
    {
        std::cout << " Can't find the INPUT file." << std::endl;
        return false;
    }

    ifs.clear();
    ifs.seekg(0);

    char word[80];
    char word1[80];
    int ierr = 0;

    // ifs >> std::setiosflags(ios::uppercase);
    ifs.rdstate();
    while (ifs.good())
    {
        ifs >> word;
        ifs.ignore(150, '\n');
        if (strcmp(word, "INPUT_PARAMETERS") == 0)
        {
            ierr = 1;
            break;
        }
        ifs.rdstate();
    }

    if (ierr == 0)
    {
        std::cout << " Error parameter list. "
                  << " The parameter list always starts with key word 'INPUT_PARAMETERS'. " << std::endl;
        return false; // return error : false
    }

    ifs.rdstate();
    while (ifs.good())
    {
        ifs >> word1;
        if (ifs.eof())
            break;
        strtolower(word1, word);

        //----------------------------------------------------------
        // main parameters
        //----------------------------------------------------------
        if (strcmp("suffix", word) == 0) // out dir
        {
            read_value(ifs, suffix);
        }
        else if (strcmp("stru_file", word) == 0) // xiaohui modify 2015-02-01
        {
            read_value(ifs, stru_file); // xiaohui modify 2015-02-01
        }
        else if (strcmp("pseudo_dir", word) == 0)
        {
            read_value(ifs, pseudo_dir);
        }
        // else if (strcmp("pseudo_type", word) == 0) // mohan add 2013-05-20 (xiaohui add 2013-06-23)
        // {
        //     read_value(ifs, pseudo_type);
        // }
        else if (strcmp("orbital_dir", word) == 0) // liuyu add 2021-08-14
        {
            read_value(ifs, orbital_dir);
        }
        else if (strcmp("kpoint_file", word) == 0) // xiaohui modify 2015-02-01
        {
            read_value(ifs, kpoint_file); // xiaohui modify 2015-02-01
        }
        else if (strcmp("wannier_card", word) == 0) // mohan add 2009-12-25
        {
            read_value(ifs, wannier_card);
        }
        else if (strcmp("latname", word) == 0) // which material
        {
            read_value(ifs, latname);
        }
        else if (strcmp("pseudo_rcut", word) == 0) //
        {
            read_value(ifs, pseudo_rcut);
        }
        else if (strcmp("pseudo_mesh", word) == 0) //
        {
            read_bool(ifs, pseudo_mesh);
        }
        else if (strcmp("calculation", word) == 0) // which type calculation
        {
            read_value(ifs, calculation);
        }
        else if (strcmp("esolver_type", word) == 0)
        {
            read_value(ifs, esolver_type);
        }
        else if (strcmp("ntype", word) == 0) // number of atom types
        {
            read_value(ifs, ntype);
        }
        else if (strcmp("nbands", word) == 0) // number of atom bands
        {
            read_value(ifs, nbands);
        }
        else if (strcmp("nbands_sto", word) == 0) // number of stochastic bands
        {
            std::string nbsto_str;
            read_value(ifs, nbndsto_str);
            if (nbndsto_str != "all")
            {
                nbands_sto = std::stoi(nbndsto_str);
            }
        }
        else if (strcmp("kspacing", word) == 0)
        {
            read_kspacing(ifs);
        }
        else if (strcmp("min_dist_coef", word) == 0)
        {
            read_value(ifs, min_dist_coef);
        }
        else if (strcmp("nbands_istate", word) == 0) // number of atom bands
        {
            read_value(ifs, nbands_istate);
            // Originally disabled in line 2401.
            // if (nbands_istate < 0)
            // 	ModuleBase::WARNING_QUIT("Input", "NBANDS_ISTATE must > 0");
        }
        else if (strcmp("nche_sto", word) == 0) // Chebyshev expansion order
        {
            read_value(ifs, nche_sto);
        }
        else if (strcmp("seed_sto", word) == 0)
        {
            read_value(ifs, seed_sto);
        }
        else if (strcmp("initsto_ecut", word) == 0)
        {
            read_value(ifs, initsto_ecut);
        }
        else if (strcmp("pw_seed", word) == 0)
        {
            read_value(ifs, pw_seed);
        }
        else if (strcmp("emax_sto", word) == 0)
        {
            read_value(ifs, emax_sto);
        }
        else if (strcmp("emin_sto", word) == 0)
        {
            read_value(ifs, emin_sto);
        }
        else if (strcmp("initsto_freq", word) == 0)
        {
            read_value(ifs, initsto_freq);
        }
        else if (strcmp("method_sto", word) == 0)
        {
            read_value(ifs, method_sto);
        }
        else if (strcmp("npart_sto", word) == 0)
        {
            read_value(ifs, npart_sto);
        }
        else if (strcmp("cal_cond", word) == 0)
        {
            read_bool(ifs, cal_cond);
        }
        else if (strcmp("cond_che_thr", word) == 0)
        {
            read_value(ifs, cond_che_thr);
        }
        else if (strcmp("cond_dw", word) == 0)
        {
            read_value(ifs, cond_dw);
        }
        else if (strcmp("cond_wcut", word) == 0)
        {
            read_value(ifs, cond_wcut);
        }
        else if (strcmp("cond_dt", word) == 0)
        {
            read_value(ifs, cond_dt);
        }
        else if (strcmp("cond_dtbatch", word) == 0)
        {
            read_value(ifs, cond_dtbatch);
        }
        else if (strcmp("cond_smear", word) == 0)
        {
            read_value(ifs, cond_smear);
        }
        else if (strcmp("cond_fwhm", word) == 0)
        {
            read_value(ifs, cond_fwhm);
        }
        else if (strcmp("cond_nonlocal", word) == 0)
        {
            read_bool(ifs, cond_nonlocal);
        }
        else if (strcmp("bndpar", word) == 0)
        {
            read_value(ifs, bndpar);
        }
        else if (strcmp("kpar", word) == 0) // number of pools
        {
            read_value(ifs, kpar);
        }
        else if (strcmp("berry_phase", word) == 0) // berry phase calculation
        {
            read_bool(ifs, berry_phase);
        }
        else if (strcmp("gdir", word) == 0) // berry phase calculation
        {
            read_value(ifs, gdir);
        }
        else if (strcmp("towannier90", word) == 0) // add by jingan for wannier90
        {
            read_bool(ifs, towannier90);
        }
        else if (strcmp("nnkpfile", word) == 0) // add by jingan for wannier90
        {
            read_value(ifs, nnkpfile);
        }
        else if (strcmp("wannier_spin", word) == 0) // add by jingan for wannier90
        {
            read_value(ifs, wannier_spin);
        }
        else if (strcmp("wannier_method", word) == 0) // add by jingan for wannier90
        {
            read_value(ifs, wannier_method);
        }
        else if (strcmp("out_wannier_mmn", word) == 0) // add by renxi for wannier90
        {
            read_bool(ifs, out_wannier_mmn);
        }
        else if (strcmp("out_wannier_amn", word) == 0) // add by renxi for wannier90
        {
            read_bool(ifs, out_wannier_amn);
        }
        else if (strcmp("out_wannier_unk", word) == 0) // add by renxi for wannier90
        {
            read_bool(ifs, out_wannier_unk);
        }
        else if (strcmp("out_wannier_eig", word) == 0) // add by renxi for wannier90
        {
            read_bool(ifs, out_wannier_eig);
        }
        else if (strcmp("out_wannier_wvfn_formatted", word) == 0)
        {
            read_bool(ifs, out_wannier_wvfn_formatted);
        }
        //----------------------------------------------------------
        // electrons / spin
        //----------------------------------------------------------
        else if (strcmp("dft_functional", word) == 0)
        {
            read_value(ifs, dft_functional);
        }
        else if (strcmp("xc_temperature", word) == 0)
        {
            read_value(ifs, xc_temperature);
        }
        else if (strcmp("nspin", word) == 0)
        {
            read_value(ifs, nspin);
        }
        else if (strcmp("nelec", word) == 0)
        {
            read_value(ifs, nelec);
        }
        else if (strcmp("nupdown", word) == 0)
        {
            read_value(ifs, nupdown);
        }
        else if (strcmp("lmaxmax", word) == 0)
        {
            read_value(ifs, lmaxmax);
        }
        //----------------------------------------------------------
        // new function
        //----------------------------------------------------------
        else if (strcmp("basis_type", word) == 0)
        {
            read_value(ifs, basis_type);
        } // xiaohui add 2013-09-01
        else if (strcmp("ks_solver", word) == 0)
        {
            read_value(ifs, ks_solver);
        } // xiaohui add 2013-09-01
        else if (strcmp("search_radius", word) == 0)
        {
            read_value(ifs, search_radius);
        }
        else if (strcmp("search_pbc", word) == 0)
        {
            read_bool(ifs, search_pbc);
        }
        else if (strcmp("symmetry", word) == 0)
        {
            read_value(ifs, symmetry);
        }
        else if (strcmp("init_vel", word) == 0)
        {
            read_bool(ifs, init_vel);
        }
        else if (strcmp("ref_cell_factor", word) == 0)
        {
            read_value(ifs, ref_cell_factor);
        }
        else if (strcmp("symmetry_prec", word) == 0) // LiuXh add 2021-08-12, accuracy for symmetry
        {
            read_value(ifs, symmetry_prec);
        }
        else if (strcmp("symmetry_autoclose", word) == 0)
        {
            read_value(ifs, symmetry_autoclose);
        }
        else if (strcmp("cal_force", word) == 0)
        {
            read_bool(ifs, cal_force);
        }
        else if (strcmp("force_thr", word) == 0)
        {
            read_value(ifs, force_thr);
        }
        else if (strcmp("force_thr_ev", word) == 0)
        {
            read_value(ifs, force_thr);
            force_thr = force_thr / 13.6058 * 0.529177;
        }
        else if (strcmp("force_thr_ev2", word) == 0)
        {
            read_value(ifs, force_thr_ev2);
        }
        else if (strcmp("stress_thr", word) == 0)
        {
            read_value(ifs, stress_thr);
        }
        else if (strcmp("press1", word) == 0)
        {
            read_value(ifs, press1);
        }
        else if (strcmp("press2", word) == 0)
        {
            read_value(ifs, press2);
        }
        else if (strcmp("press3", word) == 0)
        {
            read_value(ifs, press3);
        }
        else if (strcmp("cal_stress", word) == 0)
        {
            read_bool(ifs, cal_stress);
        }
        else if (strcmp("fixed_axes", word) == 0)
        {
            read_value(ifs, fixed_axes);
        }
        else if (strcmp("fixed_ibrav", word) == 0)
        {
            read_bool(ifs, fixed_ibrav);
        }
        else if (strcmp("fixed_atoms", word) == 0)
        {
            read_bool(ifs, fixed_atoms);
        }
        else if (strcmp("relax_method", word) == 0)
        {
            read_value(ifs, relax_method);
        }
        else if (strcmp("relax_cg_thr", word) == 0) // pengfei add 2013-08-15
        {
            read_value(ifs, relax_cg_thr);
        }
        else if (strcmp("out_level", word) == 0)
        {
            read_value(ifs, out_level);
            out_md_control = true;
        }
        else if (strcmp("relax_bfgs_w1", word) == 0)
        {
            read_value(ifs, relax_bfgs_w1);
        }
        else if (strcmp("relax_bfgs_w2", word) == 0)
        {
            read_value(ifs, relax_bfgs_w2);
        }
        else if (strcmp("relax_bfgs_rmax", word) == 0)
        {
            read_value(ifs, relax_bfgs_rmax);
        }
        else if (strcmp("relax_bfgs_rmin", word) == 0)
        {
            read_value(ifs, relax_bfgs_rmin);
        }
        else if (strcmp("relax_bfgs_init", word) == 0)
        {
            read_value(ifs, relax_bfgs_init);
        }
        else if (strcmp("relax_scale_force", word) == 0)
        {
            read_value(ifs, relax_scale_force);
        }
        else if (strcmp("relax_new", word) == 0)
        {
            read_bool(ifs, relax_new);
        }
        else if (strcmp("use_paw", word) == 0)
        {
            read_bool(ifs, use_paw);
        }
        //----------------------------------------------------------
        // plane waves
        //----------------------------------------------------------
        else if (strcmp("gamma_only", word) == 0)
        {
            read_bool(ifs, gamma_only);
        }
        else if (strcmp("fft_mode", word) == 0)
        {
            read_value(ifs, fft_mode);
        }
        else if (strcmp("ecutwfc", word) == 0)
        {
            read_value(ifs, ecutwfc);
        }
        else if (strcmp("ecutrho", word) == 0)
        {
            read_value(ifs, ecutrho);
        }
        else if (strcmp("nx", word) == 0)
        {
            read_value(ifs, nx);
            ncx = nx;
        }
        else if (strcmp("ny", word) == 0)
        {
            read_value(ifs, ny);
            ncy = ny;
        }
        else if (strcmp("nz", word) == 0)
        {
            read_value(ifs, nz);
            ncz = nz;
        }
        else if (strcmp("bx", word) == 0)
        {
            read_value(ifs, bx);
        }
        else if (strcmp("by", word) == 0)
        {
            read_value(ifs, by);
        }
        else if (strcmp("bz", word) == 0)
        {
            read_value(ifs, bz);
        }
        else if (strcmp("ndx", word) == 0)
        {
            read_value(ifs, ndx);
        }
        else if (strcmp("ndy", word) == 0)
        {
            read_value(ifs, ndy);
        }
        else if (strcmp("ndz", word) == 0)
        {
            read_value(ifs, ndz);
        }
        else if (strcmp("erf_ecut", word) == 0)
        {
            read_value(ifs, erf_ecut);
        }
        else if (strcmp("erf_height", word) == 0)
        {
            read_value(ifs, erf_height);
        }
        else if (strcmp("erf_sigma", word) == 0)
        {
            read_value(ifs, erf_sigma);
        }
        //----------------------------------------------------------
        // diagonalization
        //----------------------------------------------------------
        // else if (strcmp("diago_type", word) == 0)
        //{
        //    read_value(ifs, diago_type);
        //} xiaohui modify 2013-09-01
        else if (strcmp("diago_proc", word) == 0)
        {
            read_value(ifs, diago_proc);
        }
        else if (strcmp("pw_diag_nmax", word) == 0)
        {
            read_value(ifs, pw_diag_nmax);
        }
        else if (strcmp("diago_cg_prec", word) == 0) // mohan add 2012-03-31
        {
            read_value(ifs, diago_cg_prec);
        }
        else if (strcmp("pw_diag_ndim", word) == 0)
        {
            read_value(ifs, pw_diag_ndim);
        }
        else if (strcmp("pw_diag_thr", word) == 0)
        {
            read_value(ifs, pw_diag_thr);
        }
        else if (strcmp("nb2d", word) == 0)
        {
            read_value(ifs, nb2d);
        }
        else if (strcmp("nurse", word) == 0)
        {
            read_value(ifs, nurse);
        }
        else if (strcmp("colour", word) == 0)
        {
            read_bool(ifs, colour);
        }
        else if (strcmp("nbspline", word) == 0)
        {
            read_value(ifs, nbspline);
        }
        else if (strcmp("t_in_h", word) == 0)
        {
            read_bool(ifs, t_in_h);
        }
        else if (strcmp("vl_in_h", word) == 0)
        {
            read_bool(ifs, vl_in_h);
        }
        else if (strcmp("vnl_in_h", word) == 0)
        {
            read_bool(ifs, vnl_in_h);
        }
        else if (strcmp("vh_in_h", word) == 0)
        {
            read_bool(ifs, vh_in_h);
        }
        else if (strcmp("vion_in_h", word) == 0)
        {
            read_bool(ifs, vion_in_h);
        }
        else if (strcmp("test_force", word) == 0)
        {
            read_bool(ifs, test_force);
        }
        else if (strcmp("test_stress", word) == 0)
        {
            read_bool(ifs, test_stress);
        }
        //----------------------------------------------------------
        // iteration
        //----------------------------------------------------------
        else if (strcmp("scf_thr", word) == 0)
        {
            read_value(ifs, scf_thr);
        }
        else if (strcmp("scf_thr_type", word) == 0)
        {
            read_value(ifs, scf_thr_type);
        }
        else if (strcmp("scf_nmax", word) == 0)
        {
            read_value(ifs, scf_nmax);
        }
        else if (strcmp("relax_nmax", word) == 0)
        {
            read_value(ifs, this->relax_nmax);
        }
        else if (strcmp("out_stru", word) == 0)
        {
            read_bool(ifs, out_stru);
        }
        //----------------------------------------------------------
        // occupation
        //----------------------------------------------------------
        // else if (strcmp("occupations", word) == 0)
        //{
        //    read_value(ifs, occupations);
        //}
        else if (strcmp("smearing_method", word) == 0)
        {
            read_value(ifs, smearing_method);
        }
        else if (strcmp("smearing_sigma", word) == 0)
        {
            read_value(ifs, smearing_sigma);
        }
        else if (strcmp("smearing_sigma_temp", word) == 0)
        {
            double smearing_sigma_temp;
            read_value(ifs, smearing_sigma_temp);
            smearing_sigma = smearing_sigma_temp * 3.166815e-6;
        }
        //----------------------------------------------------------
        // charge mixing
        //----------------------------------------------------------
        else if (strcmp("mixing_type", word) == 0)
        {
            read_value(ifs, mixing_mode);
        }
        else if (strcmp("mixing_beta", word) == 0)
        {
            read_value(ifs, mixing_beta);
        }
        else if (strcmp("mixing_ndim", word) == 0)
        {
            read_value(ifs, mixing_ndim);
        }
        else if (strcmp("mixing_gg0", word) == 0) // mohan add 2014-09-27
        {
            read_value(ifs, mixing_gg0);
        }
        else if (strcmp("mixing_beta_mag", word) == 0)
        {
            read_value(ifs, mixing_beta_mag);
        }
        else if (strcmp("mixing_gg0_mag", word) == 0)
        {
            read_value(ifs, mixing_gg0_mag);
        }
        else if (strcmp("mixing_gg0_min", word) == 0)
        {
            read_value(ifs, mixing_gg0_min);
        }
        else if (strcmp("mixing_angle", word) == 0)
        {
            read_value(ifs, mixing_angle);
        }
        else if (strcmp("mixing_tau", word) == 0)
        {
            read_bool(ifs, mixing_tau);
        }
        else if (strcmp("mixing_dftu", word) == 0)
        {
            read_bool(ifs, mixing_dftu);
        }
        //----------------------------------------------------------
        // charge / potential / wavefunction
        //----------------------------------------------------------
        else if (strcmp("read_file_dir", word) == 0)
        {
            read_value(ifs, read_file_dir);
        }
        else if (strcmp("init_wfc", word) == 0)
        {
            read_value(ifs, init_wfc);
        }
        else if (strcmp("psi_initializer", word) == 0)
        {
            read_value(ifs, psi_initializer);
        }
        else if (strcmp("mem_saver", word) == 0)
        {
            read_value(ifs, mem_saver);
        }
        else if (strcmp("printe", word) == 0)
        {
            read_value(ifs, printe);
        }
        else if (strcmp("init_chg", word) == 0)
        {
            read_value(ifs, init_chg);
        }
        else if (strcmp("chg_extrap", word) == 0) // xiaohui modify 2015-02-01
        {
            read_value(ifs, chg_extrap); // xiaohui modify 2015-02-01
        }
        else if (strcmp("out_freq_elec", word) == 0)
        {
            read_value(ifs, out_freq_elec);
        }
        else if (strcmp("out_freq_ion", word) == 0)
        {
            read_value(ifs, out_freq_ion);
        }
        else if (strcmp("out_chg", word) == 0)
        {
            read_bool(ifs, out_chg);
        }
        else if (strcmp("band_print_num", word) == 0)
        {
            read_value(ifs, band_print_num);
        }
        else if (strcmp("bands_to_print", word) == 0)
        {
            ifs.ignore(150, '\n');
        }
        else if (strcmp("out_dm", word) == 0)
        {
            read_bool(ifs, out_dm);
        }
        else if (strcmp("out_dm1", word) == 0)
        {
            read_bool(ifs, out_dm1);
        }
        else if (strcmp("out_bandgap", word) == 0) // for bandgap printing
        {
            read_bool(ifs, out_bandgap);
        }
        else if (strcmp("deepks_out_labels", word) == 0) // caoyu added 2020-11-24, mohan modified 2021-01-03
        {
            read_bool(ifs, deepks_out_labels);
        }
        else if (strcmp("deepks_scf", word) == 0) // caoyu added 2020-11-24, mohan modified 2021-01-03
        {
            read_bool(ifs, deepks_scf);
        }
        else if (strcmp("deepks_bandgap", word) == 0) // caoyu added 2020-11-24, mohan modified 2021-01-03
        {
            read_bool(ifs, deepks_bandgap);
        }
        else if (strcmp("deepks_out_unittest", word) == 0) // mohan added 2021-01-03
        {
            read_bool(ifs, deepks_out_unittest);
        }
        else if (strcmp("deepks_model", word) == 0) // caoyu added 2021-06-03
        {
            read_value(ifs, deepks_model);
        }
        else if (strcmp("out_pot", word) == 0)
        {
            read_value(ifs, out_pot);
        }
        else if (strcmp("out_wfc_pw", word) == 0)
        {
            read_value(ifs, out_wfc_pw);
        }
        else if (strcmp("out_wfc_r", word) == 0)
        {
            read_bool(ifs, out_wfc_r);
        }
        // mohan add 20090909
        else if (strcmp("out_dos", word) == 0)
        {
            read_value(ifs, out_dos);
        }
        else if (strcmp("out_band", word) == 0)
        {
            read_value2stdvector(ifs, out_band);
            if(out_band.size() == 1) out_band.push_back(8);
        }
        else if (strcmp("out_proj_band", word) == 0)
        {
            read_bool(ifs, out_proj_band);
        }
        else if (strcmp("out_mat_hs", word) == 0)
        {
            read_value2stdvector(ifs, out_mat_hs);
            if(out_mat_hs.size() == 1) out_mat_hs.push_back(8);
        }
        // LiuXh add 2019-07-15
        else if (strcmp("out_mat_hs2", word) == 0)
        {
            read_bool(ifs, out_mat_hs2);
        }
        else if (strcmp("out_mat_t", word) == 0)
        {
            read_bool(ifs, out_mat_t);
        }
        else if (strcmp("out_mat_dh", word) == 0)
        {
            read_bool(ifs, out_mat_dh);
        }
        else if (strcmp("out_mat_xc", word) == 0)
        {
            read_bool(ifs, out_mat_xc);
        }
        else if (strcmp("out_interval", word) == 0)
        {
            read_value(ifs, out_interval);
        }
        else if (strcmp("out_app_flag", word) == 0)
        {
            read_bool(ifs, out_app_flag);
        }
        else if (strcmp("out_ndigits", word) == 0)
        {
            read_value(ifs, out_ndigits);
        }
        else if (strcmp("out_mat_r", word) == 0)
        {
            read_bool(ifs, out_mat_r);
        }
        else if (strcmp("out_wfc_lcao", word) == 0)
        {
            read_value(ifs, out_wfc_lcao);
        }
        else if (strcmp("out_alllog", word) == 0)
        {
            read_bool(ifs, out_alllog);
        }
        else if (strcmp("out_element_info", word) == 0)
        {
            read_bool(ifs, out_element_info);
        }
        else if (strcmp("dos_emin_ev", word) == 0)
        {
            read_value(ifs, dos_emin_ev);
            dos_setemin = true;
        }
        else if (strcmp("dos_emax_ev", word) == 0)
        {
            read_value(ifs, dos_emax_ev);
            dos_setemax = true;
        }
        else if (strcmp("dos_edelta_ev", word) == 0)
        {
            read_value(ifs, dos_edelta_ev);
        }
        else if (strcmp("dos_scale", word) == 0)
        {
            read_value(ifs, dos_scale);
        }
        else if (strcmp("dos_sigma", word) == 0)
        {
            read_value(ifs, dos_sigma);
        }
        else if (strcmp("dos_nche", word) == 0)
        {
            read_value(ifs, dos_nche);
        }

        //----------------------------------------------------------
        // Parameters about LCAO
        // mohan add 2009-11-11
        //----------------------------------------------------------
        else if (strcmp("lcao_ecut", word) == 0)
        {
            read_value(ifs, lcao_ecut);
        }
        else if (strcmp("lcao_dk", word) == 0)
        {
            read_value(ifs, lcao_dk);
        }
        else if (strcmp("lcao_dr", word) == 0)
        {
            read_value(ifs, lcao_dr);
        }
        else if (strcmp("lcao_rmax", word) == 0)
        {
            read_value(ifs, lcao_rmax);
        }
        //----------------------------------------------------------
        // Molecule Dynamics
        // Yu Liu add 2021-07-30
        //----------------------------------------------------------
        else if (strcmp("md_type", word) == 0)
        {
            read_value(ifs, mdp.md_type);
        }
        else if (strcmp("md_thermostat", word) == 0)
        {
            read_value(ifs, mdp.md_thermostat);
        }
        else if (strcmp("md_nraise", word) == 0)
        {
            read_value(ifs, mdp.md_nraise);
        }
        else if (strcmp("cal_syns", word) == 0)
        {
            read_value(ifs, cal_syns);
        }
        else if (strcmp("dmax", word) == 0)
        {
            read_value(ifs, dmax);
        }
        else if (strcmp("md_tolerance", word) == 0)
        {
            read_value(ifs, mdp.md_tolerance);
        }
        else if (strcmp("md_nstep", word) == 0)
        {
            read_value(ifs, mdp.md_nstep);
        }
        else if (strcmp("md_dt", word) == 0)
        {
            read_value(ifs, mdp.md_dt);
        }
        else if (strcmp("md_tchain", word) == 0)
        {
            read_value(ifs, mdp.md_tchain);
        }
        else if (strcmp("md_tfirst", word) == 0)
        {
            read_value(ifs, mdp.md_tfirst);
        }
        else if (strcmp("md_tlast", word) == 0)
        {
            read_value(ifs, mdp.md_tlast);
        }
        else if (strcmp("md_dumpfreq", word) == 0)
        {
            read_value(ifs, mdp.md_dumpfreq);
        }
        else if (strcmp("md_restartfreq", word) == 0)
        {
            read_value(ifs, mdp.md_restartfreq);
        }
        else if (strcmp("md_seed", word) == 0)
        {
            read_value(ifs, mdp.md_seed);
        }
        else if (strcmp("md_prec_level", word) == 0)
        {
            read_value(ifs, mdp.md_prec_level);
        }
        else if (strcmp("md_restart", word) == 0)
        {
            read_bool(ifs, mdp.md_restart);
        }
        else if (strcmp("md_pmode", word) == 0)
        {
            read_value(ifs, mdp.md_pmode);
        }
        else if (strcmp("md_pcouple", word) == 0)
        {
            read_value(ifs, mdp.md_pcouple);
        }
        else if (strcmp("md_pchain", word) == 0)
        {
            read_value(ifs, mdp.md_pchain);
        }
        else if (strcmp("md_pfirst", word) == 0)
        {
            read_value(ifs, mdp.md_pfirst);
        }
        else if (strcmp("md_plast", word) == 0)
        {
            read_value(ifs, mdp.md_plast);
        }
        else if (strcmp("md_pfreq", word) == 0)
        {
            read_value(ifs, mdp.md_pfreq);
        }
        else if (strcmp("lj_rcut", word) == 0)
        {
            read_value(ifs, mdp.lj_rcut);
        }
        else if (strcmp("lj_epsilon", word) == 0)
        {
            read_value(ifs, mdp.lj_epsilon);
        }
        else if (strcmp("lj_sigma", word) == 0)
        {
            read_value(ifs, mdp.lj_sigma);
        }
        else if (strcmp("msst_direction", word) == 0)
        {
            read_value(ifs, mdp.msst_direction);
        }
        else if (strcmp("msst_vel", word) == 0)
        {
            read_value(ifs, mdp.msst_vel);
        }
        else if (strcmp("msst_vis", word) == 0)
        {
            read_value(ifs, mdp.msst_vis);
        }
        else if (strcmp("msst_tscale", word) == 0)
        {
            read_value(ifs, mdp.msst_tscale);
        }
        else if (strcmp("msst_qmass", word) == 0)
        {
            read_value(ifs, mdp.msst_qmass);
        }
        else if (strcmp("md_tfreq", word) == 0)
        {
            read_value(ifs, mdp.md_tfreq);
        }
        else if (strcmp("md_damp", word) == 0)
        {
            read_value(ifs, mdp.md_damp);
        }
        else if (strcmp("pot_file", word) == 0)
        {
            read_value(ifs, mdp.pot_file);
        }
        else if (strcmp("dump_force", word) == 0)
        {
            read_bool(ifs, mdp.dump_force);
        }
        else if (strcmp("dump_vel", word) == 0)
        {
            read_bool(ifs, mdp.dump_vel);
        }
        else if (strcmp("dump_virial", word) == 0)
        {
            read_bool(ifs, mdp.dump_virial);
        }
        //----------------------------------------------------------
        // efield and dipole correction
        // Yu Liu add 2022-05-18
        //----------------------------------------------------------
        else if (strcmp("efield_flag", word) == 0)
        {
            read_bool(ifs, efield_flag);
        }
        else if (strcmp("dip_cor_flag", word) == 0)
        {
            read_bool(ifs, dip_cor_flag);
        }
        else if (strcmp("efield_dir", word) == 0)
        {
            read_value(ifs, efield_dir);
        }
        else if (strcmp("efield_pos_max", word) == 0)
        {
            read_value(ifs, efield_pos_max);
        }
        else if (strcmp("efield_pos_dec", word) == 0)
        {
            read_value(ifs, efield_pos_dec);
        }
        else if (strcmp("efield_amp", word) == 0)
        {
            read_value(ifs, efield_amp);
        }
        //----------------------------------------------------------
        // gatefield (compensating charge)
        // Yu Liu add 2022-09-13
        //----------------------------------------------------------
        else if (strcmp("gate_flag", word) == 0)
        {
            read_bool(ifs, gate_flag);
        }
        else if (strcmp("zgate", word) == 0)
        {
            read_value(ifs, zgate);
        }
        else if (strcmp("relax", word) == 0)
        {
            read_bool(ifs, relax);
        }
        else if (strcmp("block", word) == 0)
        {
            read_bool(ifs, block);
        }
        else if (strcmp("block_down", word) == 0)
        {
            read_value(ifs, block_down);
        }
        else if (strcmp("block_up", word) == 0)
        {
            read_value(ifs, block_up);
        }
        else if (strcmp("block_height", word) == 0)
        {
            read_value(ifs, block_height);
        }
        //----------------------------------------------------------
        // tddft
        // Fuxiang He add 2016-10-26
        //----------------------------------------------------------
        else if (strcmp("td_force_dt", word) == 0)
        {
            read_value(ifs, td_force_dt);
        }
        else if (strcmp("td_vext", word) == 0)
        {
            read_value(ifs, td_vext);
        }
        else if (strcmp("td_vext_dire", word) == 0)
        {
            getline(ifs, td_vext_dire);
        }
        else if (strcmp("out_dipole", word) == 0)
        {
            read_value(ifs, out_dipole);
        }
        else if (strcmp("out_efield", word) == 0)
        {
            read_value(ifs, out_efield);
        }
        else if (strcmp("td_print_eij", word) == 0)
        {
            read_value(ifs, td_print_eij);
        }
        else if (strcmp("td_edm", word) == 0)
        {
            read_value(ifs, td_edm);
        }
        else if (strcmp("td_propagator", word) == 0)
        {
            read_value(ifs, propagator);
        }
        else if (strcmp("td_stype", word) == 0)
        {
            read_value(ifs, td_stype);
        }
        else if (strcmp("td_ttype", word) == 0)
        {
            getline(ifs, td_ttype);
        }
        else if (strcmp("td_tstart", word) == 0)
        {
            read_value(ifs, td_tstart);
        }
        else if (strcmp("td_tend", word) == 0)
        {
            read_value(ifs, td_tend);
        }
        else if (strcmp("td_lcut1", word) == 0)
        {
            read_value(ifs, td_lcut1);
        }
        else if (strcmp("td_lcut2", word) == 0)
        {
            read_value(ifs, td_lcut2);
        }
        else if (strcmp("td_gauss_freq", word) == 0)
        {
            getline(ifs, td_gauss_freq);
        }
        else if (strcmp("td_gauss_phase", word) == 0)
        {
            getline(ifs, td_gauss_phase);
        }
        else if (strcmp("td_gauss_sigma", word) == 0)
        {
            getline(ifs, td_gauss_sigma);
        }
        else if (strcmp("td_gauss_t0", word) == 0)
        {
            getline(ifs, td_gauss_t0);
        }
        else if (strcmp("td_gauss_amp", word) == 0)
        {
            getline(ifs, td_gauss_amp);
        }
        else if (strcmp("td_trape_freq", word) == 0)
        {
            getline(ifs, td_trape_freq);
        }
        else if (strcmp("td_trape_phase", word) == 0)
        {
            getline(ifs, td_trape_phase);
        }
        else if (strcmp("td_trape_t1", word) == 0)
        {
            getline(ifs, td_trape_t1);
        }
        else if (strcmp("td_trape_t2", word) == 0)
        {
            getline(ifs, td_trape_t2);
        }
        else if (strcmp("td_trape_t3", word) == 0)
        {
            getline(ifs, td_trape_t3);
        }
        else if (strcmp("td_trape_amp", word) == 0)
        {
            getline(ifs, td_trape_amp);
        }
        else if (strcmp("td_trigo_freq1", word) == 0)
        {
            getline(ifs, td_trigo_freq1);
        }
        else if (strcmp("td_trigo_freq2", word) == 0)
        {
            getline(ifs, td_trigo_freq2);
        }
        else if (strcmp("td_trigo_phase1", word) == 0)
        {
            getline(ifs, td_trigo_phase1);
        }
        else if (strcmp("td_trigo_phase2", word) == 0)
        {
            getline(ifs, td_trigo_phase2);
        }
        else if (strcmp("td_trigo_amp", word) == 0)
        {
            getline(ifs, td_trigo_amp);
        }
        else if (strcmp("td_heavi_t0", word) == 0)
        {
            getline(ifs, td_heavi_t0);
        }
        else if (strcmp("td_heavi_amp", word) == 0)
        {
            getline(ifs, td_heavi_amp);
        }
        // else if (strcmp("td_hhg_amp1", word) == 0)
        // {
        //     getline(ifs, td_hhg_amp1);
        // }
        // else if (strcmp("td_hhg_amp2", word) == 0)
        // {
        //     getline(ifs, td_hhg_amp2);
        // }
        // else if (strcmp("td_hhg_freq1", word) == 0)
        // {
        //     getline(ifs, td_hhg_freq1);
        // }
        // else if (strcmp("td_hhg_freq2", word) == 0)
        // {
        //     getline(ifs, td_hhg_freq2);
        // }
        // else if (strcmp("td_hhg_phase1", word) == 0)
        // {
        //     getline(ifs, td_hhg_phase1);
        // }
        // else if (strcmp("td_hhg_phase2", word) == 0)
        // {
        //     getline(ifs, td_hhg_phase2);
        // }
        // else if (strcmp("td_hhg_t0", word) == 0)
        // {
        //     getline(ifs, td_hhg_t0);
        // }
        // else if (strcmp("td_hhg_sigma", word) == 0)
        // {
        //     getline(ifs, td_hhg_sigma);
        // }
        //----------------------------------------------------------
        // vdw
        // jiyy add 2019-08-04
        //----------------------------------------------------------
        else if (strcmp("vdw_method", word) == 0)
        {
            read_value(ifs, vdw_method);
        }
        else if (strcmp("vdw_s6", word) == 0)
        {
            read_value(ifs, vdw_s6);
        }
        else if (strcmp("vdw_s8", word) == 0)
        {
            read_value(ifs, vdw_s8);
        }
        else if (strcmp("vdw_a1", word) == 0)
        {
            read_value(ifs, vdw_a1);
        }
        else if (strcmp("vdw_a2", word) == 0)
        {
            read_value(ifs, vdw_a2);
        }
        else if (strcmp("vdw_d", word) == 0)
        {
            read_value(ifs, vdw_d);
        }
        else if (strcmp("vdw_abc", word) == 0)
        {
            read_bool(ifs, vdw_abc);
        }
        else if (strcmp("vdw_cutoff_radius", word) == 0)
        {
            read_value(ifs, vdw_cutoff_radius);
        }
        else if (strcmp("vdw_radius_unit", word) == 0)
        {
            read_value(ifs, vdw_radius_unit);
        }
        else if (strcmp("vdw_cn_thr", word) == 0)
        {
            read_value(ifs, vdw_cn_thr);
        }
        else if (strcmp("vdw_cn_thr_unit", word) == 0)
        {
            read_value(ifs, vdw_cn_thr_unit);
        }
        else if (strcmp("vdw_c6_file", word) == 0)
        {
            read_value(ifs, vdw_C6_file);
        }
        else if (strcmp("vdw_c6_unit", word) == 0)
        {
            read_value(ifs, vdw_C6_unit);
        }
        else if (strcmp("vdw_r0_file", word) == 0)
        {
            read_value(ifs, vdw_R0_file);
        }
        else if (strcmp("vdw_r0_unit", word) == 0)
        {
            read_value(ifs, vdw_R0_unit);
        }
        else if (strcmp("vdw_cutoff_type", word) == 0)
        {
            read_value(ifs, vdw_cutoff_type);
        }
        else if (strcmp("vdw_cutoff_period", word) == 0)
        {
            ifs >> vdw_cutoff_period.x >> vdw_cutoff_period.y;
            read_value(ifs, vdw_cutoff_period.z);
        }
        //--------------------------------------------------------
        // restart           Peize Lin 2020-04-04
        //--------------------------------------------------------
        else if (strcmp("restart_save", word) == 0)
        {
            read_bool(ifs, restart_save);
        }
        else if (strcmp("restart_load", word) == 0)
        {
            read_bool(ifs, restart_load);
        }
        else if (strcmp("ocp", word) == 0)
        {
            read_bool(ifs, ocp);
        }
        else if (strcmp("ocp_set", word) == 0)
        {
            getline(ifs, ocp_set);
            //			ifs.ignore(150, '\n');
        }
        else if (strcmp("out_mul", word) == 0)
        {
            read_bool(ifs, out_mul);
        } // qifeng add 2019/9/10
        //----------------------------------------------------------
        // exx
        // Peize Lin add 2018-06-20
        //----------------------------------------------------------
        else if (strcmp("exx_hybrid_alpha", word) == 0)
        {
            read_value(ifs, exx_hybrid_alpha);
        }
        else if (strcmp("exx_hse_omega", word) == 0)
        {
            read_value(ifs, exx_hse_omega);
        }
        else if (strcmp("exx_separate_loop", word) == 0)
        {
            read_bool(ifs, exx_separate_loop);
        }
        else if (strcmp("exx_hybrid_step", word) == 0)
        {
            read_value(ifs, exx_hybrid_step);
        }
        else if (strcmp("exx_mixing_beta", word) == 0)
        {
            read_value(ifs, exx_mixing_beta);
        }
        else if (strcmp("exx_lambda", word) == 0)
        {
            read_value(ifs, exx_lambda);
        }
        else if (strcmp("exx_real_number", word) == 0)
        {
            read_value(ifs, exx_real_number);
        }
        else if (strcmp("exx_pca_threshold", word) == 0)
        {
            read_value(ifs, exx_pca_threshold);
        }
        else if (strcmp("exx_c_threshold", word) == 0)
        {
            read_value(ifs, exx_c_threshold);
        }
        else if (strcmp("exx_v_threshold", word) == 0)
        {
            read_value(ifs, exx_v_threshold);
        }
        else if (strcmp("exx_dm_threshold", word) == 0)
        {
            read_value(ifs, exx_dm_threshold);
        }
        else if (strcmp("exx_schwarz_threshold", word) == 0)
        {
            read_value(ifs, exx_schwarz_threshold);
        }
        else if (strcmp("exx_cauchy_threshold", word) == 0)
        {
            read_value(ifs, exx_cauchy_threshold);
        }
        else if (strcmp("exx_c_grad_threshold", word) == 0)
        {
            read_value(ifs, exx_c_grad_threshold);
        }
        else if (strcmp("exx_v_grad_threshold", word) == 0)
        {
            read_value(ifs, exx_v_grad_threshold);
        }
        else if (strcmp("exx_cauchy_force_threshold", word) == 0)
        {
            read_value(ifs, exx_cauchy_force_threshold);
        }
        else if (strcmp("exx_cauchy_stress_threshold", word) == 0)
        {
            read_value(ifs, exx_cauchy_stress_threshold);
        }
        else if (strcmp("exx_ccp_threshold", word) == 0)
        {
            read_value(ifs, exx_ccp_threshold);
        }
        else if (strcmp("exx_ccp_rmesh_times", word) == 0)
        {
            read_value(ifs, exx_ccp_rmesh_times);
        }
        else if (strcmp("exx_distribute_type", word) == 0)
        {
            read_value(ifs, exx_distribute_type);
        }
        else if (strcmp("exx_opt_orb_lmax", word) == 0)
        {
            read_value(ifs, exx_opt_orb_lmax);
        }
        else if (strcmp("exx_opt_orb_ecut", word) == 0)
        {
            read_value(ifs, exx_opt_orb_ecut);
        }
        else if (strcmp("exx_opt_orb_tolerence", word) == 0)
        {
            read_value(ifs, exx_opt_orb_tolerence);
        }
        else if (strcmp("noncolin", word) == 0)
        {
            read_bool(ifs, noncolin);
        }
        else if (strcmp("lspinorb", word) == 0)
        {
            read_bool(ifs, lspinorb);
        }
        else if (strcmp("soc_lambda", word) == 0)
        {
            read_value(ifs, soc_lambda);
        }
        else if (strcmp("cell_factor", word) == 0)
        {
            read_value(ifs, cell_factor);
        }
        else if (strcmp("test_skip_ewald", word) == 0)
        {
            read_bool(ifs, test_skip_ewald);
        }
        //--------------
        //----------------------------------------------------------------------------------
        //         Xin Qu added on 2020-10-29 for DFT+U
        //----------------------------------------------------------------------------------
        else if (strcmp("dft_plus_u", word) == 0)
        {
            read_bool(ifs, dft_plus_u);
        }
        else if (strcmp("yukawa_potential", word) == 0)
            ifs.ignore(150, '\n');
        else if (strcmp("hubbard_u", word) == 0)
            ifs.ignore(150, '\n');
        else if (strcmp("orbital_corr", word) == 0)
            ifs.ignore(150, '\n');
        else if (strcmp("omc", word) == 0)
            ifs.ignore(150, '\n');
        else if (strcmp("yukawa_lambda", word) == 0)
            ifs.ignore(150, '\n');
        //----------------------------------------------------------------------------------
        //         Xin Qu added on 2020-08 for DFT+DMFT
        //----------------------------------------------------------------------------------
        else if (strcmp("dft_plus_dmft", word) == 0)
        {
            read_bool(ifs, dft_plus_dmft);
        }
        //----------------------------------------------------------------------------------
        //         Rong Shi added for RPA
        //----------------------------------------------------------------------------------
        else if (strcmp("rpa", word) == 0)
        {
            read_bool(ifs, rpa);
            if (rpa)
                GlobalV::rpa_setorb = true;
        }
        //----------------------------------------------------------------------------------
        //    implicit solvation model       sunml added on 2022-04-04
        //----------------------------------------------------------------------------------
        else if (strcmp("imp_sol", word) == 0)
        {
            read_bool(ifs, imp_sol);
        }
        else if (strcmp("eb_k", word) == 0)
        {
            read_value(ifs, eb_k);
        }
        else if (strcmp("tau", word) == 0)
        {
            read_value(ifs, tau);
        }
        else if (strcmp("sigma_k", word) == 0)
        {
            read_value(ifs, sigma_k);
        }
        else if (strcmp("nc_k", word) == 0)
        {
            read_value(ifs, nc_k);
        }
        //----------------------------------------------------------------------------------
        //    OFDFT sunliang added on 2022-05-05
        //----------------------------------------------------------------------------------
        else if (strcmp("of_kinetic", word) == 0)
        {
            read_value(ifs, of_kinetic);
        }
        else if (strcmp("of_method", word) == 0)
        {
            read_value(ifs, of_method);
        }
        else if (strcmp("of_conv", word) == 0)
        {
            read_value(ifs, of_conv);
        }
        else if (strcmp("of_tole", word) == 0)
        {
            read_value(ifs, of_tole);
        }
        else if (strcmp("of_tolp", word) == 0)
        {
            read_value(ifs, of_tolp);
        }
        else if (strcmp("of_tf_weight", word) == 0)
        {
            read_value(ifs, of_tf_weight);
        }
        else if (strcmp("of_vw_weight", word) == 0)
        {
            read_value(ifs, of_vw_weight);
        }
        else if (strcmp("of_wt_alpha", word) == 0)
        {
            read_value(ifs, of_wt_alpha);
        }
        else if (strcmp("of_wt_beta", word) == 0)
        {
            read_value(ifs, of_wt_beta);
        }
        else if (strcmp("of_wt_rho0", word) == 0)
        {
            read_value(ifs, of_wt_rho0);
        }
        else if (strcmp("of_hold_rho0", word) == 0)
        {
            read_bool(ifs, of_hold_rho0);
        }
        else if (strcmp("of_lkt_a", word) == 0)
        {
            read_value(ifs, of_lkt_a);
        }
        else if (strcmp("of_full_pw", word) == 0)
        {
            read_bool(ifs, of_full_pw);
        }
        else if (strcmp("of_full_pw_dim", word) == 0)
        {
            read_value(ifs, of_full_pw_dim);
        }
        else if (strcmp("of_read_kernel", word) == 0)
        {
            read_bool(ifs, of_read_kernel);
        }
        else if (strcmp("of_kernel_file", word) == 0)
        {
            read_value(ifs, of_kernel_file);
        }
        else if (strcmp("bessel_nao_smooth", word) == 0)
        {
            read_value(ifs, bessel_nao_smooth);
        }
        else if (strcmp("bessel_nao_sigma", word) == 0)
        {
            read_value(ifs, bessel_nao_sigma);
        }
        else if (strcmp("bessel_nao_ecut", word) == 0)
        {
            read_value(ifs, bessel_nao_ecut);
        }
        else if (strcmp("bessel_nao_rcut", word) == 0)
        {
            //read_value(ifs, bessel_nao_rcut);
            read_value2stdvector(ifs, bessel_nao_rcuts);
            bessel_nao_rcut = bessel_nao_rcuts[0]; // also compatible with old input file
        }
        else if (strcmp("bessel_nao_tolerence", word) == 0)
        {
            read_value(ifs, bessel_nao_tolerence);
        }
        else if (strcmp("bessel_descriptor_lmax", word) == 0)
        {
            read_value(ifs, bessel_descriptor_lmax);
        }
        else if (strcmp("bessel_descriptor_smooth", word) == 0)
        {
            read_value(ifs, bessel_descriptor_smooth);
        }
        else if (strcmp("bessel_descriptor_sigma", word) == 0)
        {
            read_value(ifs, bessel_descriptor_sigma);
        }
        else if (strcmp("bessel_descriptor_ecut", word) == 0)
        {
            read_value(ifs, bessel_descriptor_ecut);
        }
        else if (strcmp("bessel_descriptor_rcut", word) == 0)
        {
            read_value(ifs, bessel_descriptor_rcut);
        }
        else if (strcmp("bessel_descriptor_tolerence", word) == 0)
        {
            read_value(ifs, bessel_descriptor_tolerence);
        }
        //----------------------------------------------------------------------------------
        //    device control denghui added on 2022-11-05
        //----------------------------------------------------------------------------------
        else if (strcmp("device", word) == 0)
        {
            read_value(ifs, device);
        }
        //----------------------------------------------------------------------------------
        //    precision control denghui added on 2023-01-01
        //----------------------------------------------------------------------------------
        else if (strcmp("precision", word) == 0)
        {
            read_value(ifs, precision);
        }
        //----------------------------------------------------------------------------------
        //    Deltaspin
        //----------------------------------------------------------------------------------
        else if (strcmp("sc_mag_switch", word) == 0)
        {
            read_bool(ifs, sc_mag_switch);
        }
        else if (strcmp("decay_grad_switch", word) == 0)
        {
            read_bool(ifs, decay_grad_switch);
        }
        else if (strcmp("sc_thr", word) == 0)
        {
            read_value(ifs, sc_thr);
        }
        else if (strcmp("nsc", word) == 0)
        {
            read_value(ifs, nsc);
        }
        else if (strcmp("nsc_min", word) == 0)
        {
            read_value(ifs, nsc_min);
        }
        else if (strcmp("sc_scf_nmin", word) == 0)
        {
            read_value(ifs, sc_scf_nmin);
        }
        else if (strcmp("alpha_trial", word) == 0)
        {
            read_value(ifs, alpha_trial);
        }
        else if (strcmp("sccut", word) == 0)
        {
            read_value(ifs, sccut);
        }
        else if (strcmp("sc_file", word) == 0)
        {
            read_value(ifs, sc_file);
        }
        else if (strcmp("qo_switch", word) == 0){
            read_bool(ifs, qo_switch);
        }
        else if (strcmp("qo_basis", word) == 0){
            read_value(ifs, qo_basis);
        }
        else if (strcmp("qo_thr", word) == 0){
            read_value(ifs, qo_thr);
        }
        else if (strcmp("qo_strategy", word) == 0){
            read_value2stdvector(ifs, qo_strategy);
        }
        else if (strcmp("qo_screening_coeff", word) == 0){
            read_value2stdvector(ifs, qo_screening_coeff);
        }
        else
        {
            // xiaohui add 2015-09-15
            if (word[0] != '#' && word[0] != '/')
            {
                input_error = 1;
                std::cout << " THE PARAMETER NAME '" << word << "' IS NOT USED!" << std::endl;
            }
            // mohan screen this 2012-06-30
            //			std::cout << " THE PARAMETER NAME '" << word
            //			<< "' IS NOT USED!" << std::endl;
            ifs.ignore(150, '\n');
        }

        ifs.rdstate();

        /*if(gamma_only == 1)
        {
           gamma_only_local = 1;      //pengfei 2014-10-15
           gamma_only = 0;
           std::cout << "gamma_only_local = " << gamma_only_local <<std::endl;
        }*/

        if (ifs.eof() != 0)
        {
            break;
        }
        else if (ifs.bad() != 0)
        {
            std::cout << " Bad input parameters. " << std::endl;
            return false;
        }
        else if (ifs.fail() != 0)
        {
            std::cout << " word = " << word << std::endl;
            std::cout << " Fail to read parameters. " << std::endl;
            ifs.clear();
            return false;
        }
        else if (ifs.good() == 0)
        {
            break;
        }
    }

    // sunliang added on 2022-12-06
    // To check if ntype in INPUT is equal to the atom species in STRU, if ntype is not set in INPUT, we will set it
    // according to STRU.
    if (this->stru_file == "")
    {
        this->stru_file = "STRU";
    }
    double ntype_stru = this->count_ntype(this->stru_file);
    if (this->ntype == 0)
    {
        this->ntype = ntype_stru;
        GlobalV::ofs_running << "ntype in INPUT is 0, and it is automatically set to " << this->ntype
                             << " according to STRU" << std::endl;
    }
    else if (this->ntype != ntype_stru)
    {
        ModuleBase::WARNING_QUIT("Input", "The ntype in INPUT is not equal to the ntype counted in STRU, check it.");
    }

    if(band_print_num > 0)
    {
        bands_to_print.resize(band_print_num);
        ifs.clear();
        ifs.seekg(0); // move to the beginning of the file
        ifs.rdstate();
        while (ifs.good())
        {
            ifs >> word1;
            if (ifs.eof() != 0)
                break;
            strtolower(word1, word); // convert uppercase std::string to lower case; word1 --> word

            if (strcmp("bands_to_print", word) == 0)
            {
                for(int i = 0; i < band_print_num; i ++)
                {
                    ifs >> bands_to_print[i];
                }
            }
        }
    }

    //----------------------------------------------------------
    //       DFT+U    Xin Qu  added on 2020-10-29
    //----------------------------------------------------------
    hubbard_u = new double[ntype];
    for (int i = 0; i < ntype; i++)
    {
        hubbard_u[i] = 0.0;
    }

    orbital_corr = new int[ntype];
    for (int i = 0; i < ntype; i++)
    {
        orbital_corr[i] = -1;
    }

    if (dft_plus_u)
    {
        ifs.clear();
        ifs.seekg(0); // move to the beginning of the file
        ifs.rdstate();
        while (ifs.good())
        {
            ifs >> word1;
            if (ifs.eof() != 0)
                break;
            strtolower(word1, word); // convert uppercase std::string to lower case; word1 --> word

            if (strcmp("yukawa_potential", word) == 0)
            {
                read_bool(ifs, yukawa_potential);
            }
            else if (strcmp("yukawa_lambda", word) == 0)
            {
                ifs >> yukawa_lambda;
            }
            else if (strcmp("hubbard_u", word) == 0)
            {
                for (int i = 0; i < ntype; i++)
                {
                    ifs >> hubbard_u[i];
                    hubbard_u[i] /= ModuleBase::Ry_to_eV;
                }
            }
            else if (strcmp("orbital_corr", word) == 0)
            {
                for (int i = 0; i < ntype; i++)
                {
                    ifs >> orbital_corr[i];
                }
            }
            else if (strcmp("omc", word) == 0)
            {
                read_value(ifs, omc);
            }
            else
            {
                ifs.ignore(150, '\n');
            }

            if (ifs.eof() != 0)
            {
                break;
            }
        }

        for (int i = 0; i < ntype; i++)
        {

            if (hubbard_u[i] < -1.0e-3)
            {
                std::cout << " WRONG ARGUMENTS OF hubbard_u " << std::endl;
                exit(0);
            }

            if ((orbital_corr[i] != -1) && (orbital_corr[i] != 0) && (orbital_corr[i] != 1) && (orbital_corr[i] != 2)
                && (orbital_corr[i] != 3))
            {
                std::cout << " WRONG ARGUMENTS OF orbital_corr " << std::endl;
                exit(0);
            }
        }

        dft_plus_u = 0;
        for (int i = 0; i < ntype; i++)
        {
            if (orbital_corr[i] != -1)
                dft_plus_u = 1;
        }

        if (strcmp("lcao", basis_type.c_str()) != 0)
        {
            std::cout << " WRONG ARGUMENTS OF basis_type, only lcao is support " << std::endl;
            exit(0);
        }

        if (strcmp("genelpa", ks_solver.c_str()) != 0 && strcmp(ks_solver.c_str(), "scalapack_gvx") != 0
            && strcmp(ks_solver.c_str(), "default") != 0)
        {
            std::cout << " WRONG ARGUMENTS OF ks_solver in DFT+U routine, only genelpa and scalapack_gvx are supported "
                      << std::endl;
            std::cout << " You'are using " << ks_solver.c_str() << std::endl;
            exit(0);
        }
    }

    //----------------------------------------------------------
    //       DFT+DMFT    Xin Qu  added on 2020-08
    //----------------------------------------------------------
    if (dft_plus_dmft)
    {
        ifs.clear();
        ifs.seekg(0); // move to the beginning of the file
        ifs.rdstate();
        while (ifs.good())
        {
            ifs >> word1;
            strtolower(word1, word); // convert uppercase std::string to lower case; word1 --> word

            if (strcmp("hubbard_u", word) == 0)
            {
                for (int i = 0; i < ntype; i++)
                {
                    ifs >> hubbard_u[i];
                    hubbard_u[i] /= ModuleBase::Ry_to_eV;
                }
            }
            else if (strcmp("orbital_corr", word) == 0)
            {
                for (int i = 0; i < ntype; i++)
                {
                    ifs >> orbital_corr[i];
                }
            }
            else
                ifs.ignore(150, '\n');

            if (ifs.eof() != 0)
                break;
        }

        for (int i = 0; i < ntype; i++)
        {

            if (hubbard_u[i] < -1.0e-3)
            {
                std::cout << " WRONG ARGUMENTS OF hubbard_u " << std::endl;
                exit(0);
            }

            if ((orbital_corr[i] != -1) && (orbital_corr[i] != 0) && (orbital_corr[i] != 1) && (orbital_corr[i] != 2)
                && (orbital_corr[i] != 3))
            {
                std::cout << " WRONG ARGUMENTS OF orbital_corr " << std::endl;
                exit(0);
            }
        }

        bool dmft_flag = false;
        for (int i = 0; i < ntype; i++)
        {
            if (orbital_corr[i] != -1)
            {
                dmft_flag = true;
                break;
            }
        }

        if (!dmft_flag)
        {
            std::cout << "No atoms are correlated!!!" << std::endl;
            exit(0);
        }

        if (strcmp("lcao", basis_type.c_str()) != 0)
        {
            std::cout << " WRONG ARGUMENTS OF basis_type, only lcao is support " << std::endl;
            exit(0);
        }

        /*
        if (strcmp("genelpa", ks_solver.c_str()) != 0)
        {
            std::cout << " WRONG ARGUMENTS OF ks_solver in DFT+DMFT routine, only genelpa is support " << std::endl;
            exit(0);
        }
         */
    }

    if (basis_type == "pw" && gamma_only != 0) // pengfei Li add 2015-1-31
    {
        gamma_only = 0;
        GlobalV::ofs_running << " WARNING : gamma_only has not been implemented for pw yet" << std::endl;
        GlobalV::ofs_running << " the INPUT parameter gamma_only has been reset to 0" << std::endl;
        GlobalV::ofs_running << " and a new KPT is generated with gamma point as the only k point" << std::endl;

        GlobalV::ofs_warning << " Auto generating k-points file: " << GlobalV::global_kpoint_card << std::endl;
        std::ofstream ofs(GlobalV::global_kpoint_card.c_str());
        ofs << "K_POINTS" << std::endl;
        ofs << "0" << std::endl;
        ofs << "Gamma" << std::endl;
        ofs << "1 1 1 0 0 0" << std::endl;
        ofs.close();

        // std::cout << "gamma_only =" << gamma_only << std::endl;
    }
    else if ((basis_type == "lcao" || basis_type == "lcao_in_pw") && (gamma_only == 1))
    {
        gamma_only_local = 1;
        // std::cout << "gamma_only_local =" << gamma_only_local << std::endl;
        if (esolver_type == "tddft")
        {
            GlobalV::ofs_running << " WARNING : gamma_only is not applicable for tddft" << std::endl;
            gamma_only_local = 0;
        }
    }
    if ((out_mat_r || out_mat_hs2 || out_mat_t || out_mat_dh) && gamma_only_local)
    {
        ModuleBase::WARNING_QUIT("Input",
                                 "printing of H(R)/S(R)/dH(R)/T(R) is not available for gamma only calculations");
    }
    if (out_mat_dh && nspin == 4)
    {
        ModuleBase::WARNING_QUIT("Input", "priting of dH not available for nspin = 4");
    }

    return true;
} // end read_parameters

void Input::Default_2(void) // jiyy add 2019-08-04
{
    if (GlobalV::MY_RANK != 0)
        return;
    //==========================================================
    // vdw
    // jiyy add 2019-08-04
    //==========================================================
    if (vdw_s6 == "default")
    {
        if (vdw_method == "d2")
        {
            vdw_s6 = "0.75";
        }
        else if (vdw_method == "d3_0" || vdw_method == "d3_bj")
        {
            vdw_s6 = "1.0";
        }
    }
    if (vdw_s8 == "default")
    {
        if (vdw_method == "d3_0")
        {
            vdw_s8 = "0.722";
        }
        else if (vdw_method == "d3_bj")
        {
            vdw_s8 = "0.7875";
        }
    }
    if (vdw_a1 == "default")
    {
        if (vdw_method == "d3_0")
        {
            vdw_a1 = "1.217";
        }
        else if (vdw_method == "d3_bj")
        {
            vdw_a1 = "0.4289";
        }
    }
    if (vdw_a2 == "default")
    {
        if (vdw_method == "d3_0")
        {
            vdw_a2 = "1.0";
        }
        else if (vdw_method == "d3_bj")
        {
            vdw_a2 = "4.4407";
        }
    }
    if (vdw_cutoff_radius == "default")
    {
        if (vdw_method == "d2")
        {
            vdw_cutoff_radius = "56.6918";
        }
        else if (vdw_method == "d3_0" || vdw_method == "d3_bj")
        {
            vdw_cutoff_radius = "95";
        }
    }

    if (nx * ny * nz && ndx * ndy * ndz == 0)
    {
        ndx = nx;
        ndy = ny;
        ndz = nz;
    }
    if (ndx * ndy * ndz && nx * ny * nz == 0)
    {
        nx = ndx;
        ny = ndy;
        nz = ndz;
    }
    if (ndx > nx || ndy > ny || ndz > nz)
    {
        GlobalV::double_grid = true;
    }

    if (ecutrho <= 0.0)
    {
        ecutrho = 4.0 * ecutwfc;
    }
    if (nx * ny * nz == 0 && ecutrho / ecutwfc > 4 + 1e-8)
    {
        GlobalV::double_grid = true;
    }

    if (esolver_type == "sdft" && psi_initializer)
    {
        GlobalV::ofs_warning << "psi_initializer is not available for sdft, it is automatically set to false"
                             << std::endl;
        psi_initializer = false;
    }

    if (nbndsto_str == "all")
    {
        nbands_sto = 0;
    }
    else if (nbndsto_str == "0" && esolver_type == "sdft")
    {
        esolver_type = "ksdft";
    }
    if (esolver_type != "sdft")
        bndpar = 1;
    if (bndpar > GlobalV::NPROC)
        bndpar = GlobalV::NPROC;
    if (method_sto != 1 && method_sto != 2)
    {
        method_sto = 2;
    }
    if (of_wt_rho0 != 0)
        of_hold_rho0 = true; // sunliang add 2022-06-17
    if (!of_full_pw)
        of_full_pw_dim = 0; // sunliang add 2022-08-31
    if (of_kinetic != "wt")
        of_read_kernel = false; // sunliang add 2022-09-12

    if (dft_functional == "default" && use_paw)
    {
        ModuleBase::WARNING_QUIT("Input", "dft_functional must be set when use_paw is true");
    }

    if (exx_hybrid_alpha == "default")
    {
        std::string dft_functional_lower = dft_functional;
        std::transform(dft_functional.begin(), dft_functional.end(), dft_functional_lower.begin(), tolower);
        if (dft_functional_lower == "hf" || rpa)
            exx_hybrid_alpha = "1";
        else if (dft_functional_lower == "pbe0" || dft_functional_lower == "hse" || dft_functional_lower == "scan0")
            exx_hybrid_alpha = "0.25";
    }
    if (exx_real_number == "default")
    {
        if (gamma_only)
            exx_real_number = "1";
        else
            exx_real_number = "0";
    }
    if (exx_ccp_rmesh_times == "default")
    {
        std::string dft_functional_lower = dft_functional;
        std::transform(dft_functional.begin(), dft_functional.end(), dft_functional_lower.begin(), tolower);
        if (dft_functional_lower == "hf" || dft_functional_lower == "pbe0" || dft_functional_lower == "scan0")
            exx_ccp_rmesh_times = "5";
        else if (dft_functional_lower == "hse")
            exx_ccp_rmesh_times = "1.5";
    }
    if (symmetry == "default")
    { // deal with no-forced default value
        if (gamma_only || calculation == "nscf" || calculation == "get_S" || calculation == "get_pchg"
            || calculation == "get_wf")
            symmetry = "0"; // if md or exx, symmetry will be force-set to 0 or -1 later
        else
            symmetry = "1";
    }
    if (diago_proc <= 0)
    {
        diago_proc = GlobalV::NPROC;
    }
    else if (diago_proc > GlobalV::NPROC)
    {
        diago_proc = GlobalV::NPROC;
    }
    if (calculation == "scf")
    {
        if (mem_saver == 1)
        {
            mem_saver = 0;
            ModuleBase::GlobalFunc::AUTO_SET("mem_saver", "0");
        }
        // xiaohui modify 2015-09-15, 0 -> 1
        // cal_force = 0;
        /*
                if(!noncolin)
                    cal_force = 1;
                else
                {
                    cal_force = 0;//modified by zhengdy-soc, can't calculate force now!
                    std::cout<<"sorry, can't calculate force with soc now, would be implement in next
           version!"<<std::endl;
                }
        */
        this->relax_nmax = 1;
    }
    else if (calculation == "relax") // pengfei 2014-10-13
    {
        if (mem_saver == 1)
        {
            mem_saver = 0;
            ModuleBase::GlobalFunc::AUTO_SET("mem_saver", "0");
        }
        cal_force = 1;
        if (!this->relax_nmax)
            this->relax_nmax = 50;
    }
    else if (calculation == "nscf" || calculation == "get_S")
    {
        GlobalV::CALCULATION = "nscf";
        this->relax_nmax = 1;
        out_stru = 0;

        if (basis_type == "pw"
            && calculation == "get_S") // xiaohui add 2013-09-01. Attention! maybe there is some problem
        {
            if (pw_diag_thr > 1.0e-3)
            {
                pw_diag_thr = 1.0e-5;
            }
        }
        if (cal_force) // mohan add 2010-09-07
        {
            cal_force = false;
            ModuleBase::GlobalFunc::AUTO_SET("cal_force", "false");
        }
        if (init_chg != "file")
        {
            init_chg = "file";
            ModuleBase::GlobalFunc::AUTO_SET("init_chg", init_chg);
        }
    }
    else if (calculation == "get_pchg")
    {
        GlobalV::CALCULATION = "get_pchg";
        this->relax_nmax = 1;
        out_stru = 0;
        out_dos = 0;
        out_band[0] = 0;
        out_proj_band = 0;
        cal_force = 0;
        init_wfc = "file";
        init_chg = "atomic";   // useless,
        chg_extrap = "atomic"; // xiaohui modify 2015-02-01
        out_chg = 1;           // this leads to the calculation of state charge.
        out_dm = 0;
        out_dm1 = 0;
        out_pot = 0;
    }
    else if (calculation == "get_wf")
    {
        GlobalV::CALCULATION = "get_wf"; // mohan fix 2011-11-04
        this->relax_nmax = 1;
        out_stru = 0;
        out_dos = 0;
        out_band[0] = 0;
        out_proj_band = 0;
        cal_force = 0;
        init_wfc = "file";
        init_chg = "atomic";
        chg_extrap = "atomic"; // xiaohui modify 2015-02-01
        out_chg = 1;
        out_dm = 0;
        out_dm1 = 0;
        out_pot = 0;
    }
    else if (calculation == "md") // mohan add 2011-11-04
    {
        GlobalV::CALCULATION = "md";
        symmetry = "0";
        cal_force = 1;
        if (mdp.md_nstep == 0)
        {
            GlobalV::ofs_running << "md_nstep should be set. Autoset md_nstep to 50!" << std::endl;
            mdp.md_nstep = 50;
        }
        if (!out_md_control)
            out_level = "m"; // zhengdy add 2019-04-07
        if (mdp.md_plast < 0.0)
            mdp.md_plast = mdp.md_pfirst;

        if (mdp.md_tfreq == 0)
        {
            mdp.md_tfreq = 1.0 / 40 / mdp.md_dt;
        }
        if (mdp.md_pfreq == 0)
        {
            mdp.md_pfreq = 1.0 / 400 / mdp.md_dt;
        }
        if (mdp.md_tfirst < 0 || mdp.md_restart)
        {
            init_vel = 1;
        }
        if (esolver_type == "lj" || esolver_type == "dp" || mdp.md_type == "msst" || mdp.md_type == "npt")
        {
            cal_stress = 1;
        }

        // md_prec_level only used in vc-md  liuyu 2023-03-27
        if (mdp.md_type != "msst" && mdp.md_type != "npt")
        {
            mdp.md_prec_level = 0;
        }
    }
    else if (calculation == "cell-relax") // mohan add 2011-11-04
    {
        cal_force = 1;
        cal_stress = 1;
        if (!this->relax_nmax)
            this->relax_nmax = 50;
    }
    else if (calculation == "test_memory")
    {
        this->relax_nmax = 1;
    }
    else if (calculation == "test_neighbour")
    {
        this->relax_nmax = 1;
    }
    else if (calculation == "gen_bessel")
    {
        this->relax_nmax = 1;
    }

    if (basis_type == "pw") // xiaohui add 2013-09-01
    {
        if (ks_solver == "default") // xiaohui add 2013-09-01
        {
            ks_solver = "cg";
            ModuleBase::GlobalFunc::AUTO_SET("ks_solver", "cg");
        }
        else if (ks_solver == "cg")
        {
            GlobalV::ofs_warning << " It's ok to use cg." << std::endl;
            // new rule, mohan add 2012-02-11
            // otherwise, there need wave functions transfers
            // if(diago_type=="cg") xiaohui modify 2013-09-01
            if (diago_proc != GlobalV::NPROC)
            {
                ModuleBase::WARNING("Input", "when CG is used for diago, diago_proc==GlobalV::NPROC");
                diago_proc = GlobalV::NPROC;
            }
        }
        else if (ks_solver == "dav")
        {
            GlobalV::ofs_warning << " It's ok to use dav." << std::endl;
        }
        //
        bx = 1;
        by = 1;
        bz = 1;
    }
    else if (basis_type == "lcao_in_pw")
    {
        if (ks_solver != "lapack")
        {
            ModuleBase::WARNING_QUIT("Input", "ks_solver must be lapack when basis_type is lcao_in_pw");
        }
        else
        {
            /*
                then psi initialization setting adjustment
            */
            if (!psi_initializer)
            {
                psi_initializer = true;
            }
            if (init_wfc != "nao")
            {
                init_wfc = "nao";
                GlobalV::ofs_warning << "init_wfc is set to nao when basis_type is lcao_in_pw" << std::endl;
            }
        }
        /*
            if bx, by and bz is not assigned here, undefined behavior will occur
        */
        bx = 1;
        by = 1;
        bz = 1;
    }
    else if (basis_type == "lcao")
    {
        if (ks_solver == "default")
        {
#ifdef __ELPA
            ks_solver = "genelpa";
            ModuleBase::GlobalFunc::AUTO_SET("ks_solver", "genelpa");
#else
            ks_solver = "scalapack_gvx";
            ModuleBase::GlobalFunc::AUTO_SET("ks_solver", "scalapack_gvx");
#endif
        }
        if (lcao_ecut == 0)
        {
            lcao_ecut = ecutwfc;
            ModuleBase::GlobalFunc::AUTO_SET("lcao_ecut", ecutwfc);
        }

        // if calculation is get_wf, function source/module_basis/module_pw/pw_basis_k_big.h/distribute_r()
        // will calculate nbx/nby/nbz by divide nx/ny/nz by bx/by/bz, so bx/by/bz should not be 0
        if (calculation == "get_wf")
        {
            if (!bx)
                bx = 1;
            if (!by)
                by = 1;
            if (!bz)
                bz = 1;
        }
    }

    if (basis_type == "pw" || basis_type == "lcao_in_pw")
    {
        if (gamma_only_local)
        {
            // means you can use > 1 number of k points.
            gamma_only_local = 0;
            ModuleBase::GlobalFunc::AUTO_SET("gamma_only_local", "0");
        }
    }
    // added by linpz 2023/02/13
    if (bessel_nao_ecut == "default")
    {
        bessel_nao_ecut = std::to_string(ecutwfc);
    }
    if (bessel_descriptor_ecut == "default")
    {
        bessel_descriptor_ecut = std::to_string(ecutwfc);
    }
    // charge extrapolation liuyu 2023/09/16
    if (chg_extrap == "default" && calculation == "md")
    {
        chg_extrap = "second-order";
    }
    else if (chg_extrap == "default" && (calculation == "relax" || calculation == "cell-relax"))
    {
        chg_extrap = "first-order";
    }
    else if (chg_extrap == "default")
    {
        chg_extrap = "atomic";
    }

    if (calculation != "md")
    {
        mdp.md_prec_level = 0;
    }

    if (scf_thr == -1.0)
    {
        if (basis_type == "lcao" || basis_type == "lcao_in_pw")
        {
            scf_thr = 1.0e-7;
        }
        else if (basis_type == "pw" and calculation != "nscf")
        {
            scf_thr = 1.0e-9;
        }
        else if (basis_type == "pw" and calculation == "nscf")
        {
            scf_thr = 1.0e-6;
            // In NSCF calculation, the diagonalization threshold is set to 0.1*scf/nelec.
            // In other words, the scf_thr is used to control diagonalization convergence
            // threthod in NSCF. In this case, the default 1.0e-9 is too strict.
            // renxi 20230908
        }
    }

    if (scf_thr_type == -1)
    {
        if (basis_type == "lcao" || basis_type == "lcao_in_pw")
        {
            scf_thr_type = 2;
        }
        else if (basis_type == "pw")
        {
            scf_thr_type = 1;
        }
    }

    if(qo_switch)
    {
        /* parameter logic of QO */
        out_mat_hs[0] = 1; // print H(k) and S(k)
        out_wfc_lcao = 1; // print wave function in lcao basis in kspace
        symmetry = "-1"; // disable kpoint reduce
    }
    if(qo_screening_coeff.size() != ntype)
    {
        double default_screening_coeff = (qo_screening_coeff.size() == 1)? qo_screening_coeff[0]: 0.1;
        qo_screening_coeff.resize(ntype, default_screening_coeff);
    }
    if(qo_strategy.size() != ntype)
    {
        if(qo_strategy.size() == 1)
        {
            qo_strategy.resize(ntype, qo_strategy[0]);
        }
        else
        {
            std::string default_strategy = (qo_basis == "hydrogen")? "minimal-valence": "all";
            qo_strategy.resize(ntype, default_strategy);
        }
    }

  
    // set nspin with noncolin
    if (noncolin || lspinorb)
    {
        nspin = 4;
    }

    // mixing parameters
    if (mixing_beta < 0.0)
    {
        if (nspin == 1)
        {
            mixing_beta = 0.8;
        }
        else if (nspin == 2)
        {
            mixing_beta = 0.4;
            mixing_beta_mag = 1.6;
            mixing_gg0_mag = 0.0;
        }
        else if (nspin == 4) // I will add this
        {
            mixing_beta = 0.4;
            mixing_beta_mag = 1.6;
            mixing_gg0_mag = 0.0;
        }
    }
    else
    {
        if ((nspin == 2 || nspin == 4) && mixing_beta_mag < 0.0)
        {
            if (mixing_beta <= 0.4)
            {
                mixing_beta_mag = 4 * mixing_beta;
            }
            else
            {
                mixing_beta_mag = 1.6; // 1.6 can be discussed
            }
        }
    }
}
#ifdef __MPI
void Input::Bcast()
{
    ModuleBase::TITLE("Input", "Bcast");

    //	std::cout << "\n Bcast()" << std::endl;
    //----------------------------------------------------------
    // main parameters
    //----------------------------------------------------------
    Parallel_Common::bcast_string(suffix);
    Parallel_Common::bcast_string(stru_file); // xiaohui modify 2015-02-01
    Parallel_Common::bcast_string(pseudo_dir);
    // Parallel_Common::bcast_string(pseudo_type); // mohan add 2013-05-20 (xiaohui add 2013-06-23)
    Parallel_Common::bcast_string(orbital_dir);
    Parallel_Common::bcast_string(kpoint_file); // xiaohui modify 2015-02-01
    Parallel_Common::bcast_string(wannier_card);
    Parallel_Common::bcast_string(latname);
    Parallel_Common::bcast_string(calculation);
    Parallel_Common::bcast_string(esolver_type);
    Parallel_Common::bcast_double(pseudo_rcut);
    Parallel_Common::bcast_bool(pseudo_mesh);
    Parallel_Common::bcast_int(ntype);
    Parallel_Common::bcast_int(nbands);
    Parallel_Common::bcast_int(nbands_sto);
    Parallel_Common::bcast_int(nbands_istate);
    for (int i = 0; i < 3; i++)
    {
        Parallel_Common::bcast_double(kspacing[i]);
    }
    Parallel_Common::bcast_double(min_dist_coef);
    Parallel_Common::bcast_int(nche_sto);
    Parallel_Common::bcast_int(seed_sto);
    Parallel_Common::bcast_double(initsto_ecut);
    Parallel_Common::bcast_int(pw_seed);
    Parallel_Common::bcast_double(emax_sto);
    Parallel_Common::bcast_double(emin_sto);
    Parallel_Common::bcast_int(initsto_freq);
    Parallel_Common::bcast_int(method_sto);
    Parallel_Common::bcast_int(npart_sto);
    Parallel_Common::bcast_bool(cal_cond);
    Parallel_Common::bcast_double(cond_che_thr);
    Parallel_Common::bcast_double(cond_dw);
    Parallel_Common::bcast_double(cond_wcut);
    Parallel_Common::bcast_double(cond_dt);
    Parallel_Common::bcast_int(cond_dtbatch);
    Parallel_Common::bcast_int(cond_smear);
    Parallel_Common::bcast_double(cond_fwhm);
    Parallel_Common::bcast_bool(cond_nonlocal);
    Parallel_Common::bcast_int(bndpar);
    Parallel_Common::bcast_int(kpar);
    Parallel_Common::bcast_bool(berry_phase);
    Parallel_Common::bcast_int(gdir);
    Parallel_Common::bcast_bool(towannier90);
    Parallel_Common::bcast_string(nnkpfile);
    Parallel_Common::bcast_string(wannier_spin);
    Parallel_Common::bcast_int(wannier_method);
    Parallel_Common::bcast_bool(out_wannier_mmn);
    Parallel_Common::bcast_bool(out_wannier_amn);
    Parallel_Common::bcast_bool(out_wannier_unk);
    Parallel_Common::bcast_bool(out_wannier_eig);
    Parallel_Common::bcast_bool(out_wannier_wvfn_formatted);

    Parallel_Common::bcast_string(dft_functional);
    Parallel_Common::bcast_double(xc_temperature);
    Parallel_Common::bcast_int(nspin);
    Parallel_Common::bcast_double(nelec);
    Parallel_Common::bcast_double(nupdown);
    Parallel_Common::bcast_int(lmaxmax);

    Parallel_Common::bcast_string(basis_type); // xiaohui add 2013-09-01
    Parallel_Common::bcast_string(ks_solver);  // xiaohui add 2013-09-01
    Parallel_Common::bcast_double(search_radius);
    Parallel_Common::bcast_bool(search_pbc);
    Parallel_Common::bcast_double(search_radius);
    Parallel_Common::bcast_string(symmetry);
    Parallel_Common::bcast_bool(init_vel); // liuyu 2021-07-14
    Parallel_Common::bcast_double(ref_cell_factor);
    Parallel_Common::bcast_double(symmetry_prec); // LiuXh add 2021-08-12, accuracy for symmetry
    Parallel_Common::bcast_bool(symmetry_autoclose);
    Parallel_Common::bcast_bool(cal_force);
    Parallel_Common::bcast_double(force_thr);
    Parallel_Common::bcast_double(force_thr_ev2);
    Parallel_Common::bcast_double(stress_thr); // LiuXh add 20180515
    Parallel_Common::bcast_double(press1);
    Parallel_Common::bcast_double(press2);
    Parallel_Common::bcast_double(press3);
    Parallel_Common::bcast_bool(cal_stress);
    Parallel_Common::bcast_string(fixed_axes);
    Parallel_Common::bcast_bool(fixed_ibrav);
    Parallel_Common::bcast_bool(fixed_atoms);
    Parallel_Common::bcast_string(relax_method);
    Parallel_Common::bcast_double(relax_cg_thr); // pengfei add 2013-08-15
    Parallel_Common::bcast_string(out_level);
    Parallel_Common::bcast_bool(out_md_control);
    Parallel_Common::bcast_double(relax_bfgs_w1);
    Parallel_Common::bcast_double(relax_bfgs_w2);
    Parallel_Common::bcast_double(relax_bfgs_rmax);
    Parallel_Common::bcast_double(relax_bfgs_rmin);
    Parallel_Common::bcast_double(relax_bfgs_init);
    Parallel_Common::bcast_double(relax_scale_force);
    Parallel_Common::bcast_bool(relax_new);

    Parallel_Common::bcast_bool(use_paw);

    Parallel_Common::bcast_bool(gamma_only);
    Parallel_Common::bcast_bool(gamma_only_local);
    Parallel_Common::bcast_int(fft_mode);
    Parallel_Common::bcast_double(ecutwfc);
    Parallel_Common::bcast_double(ecutrho);
    Parallel_Common::bcast_bool(GlobalV::double_grid);
    Parallel_Common::bcast_int(ncx);
    Parallel_Common::bcast_int(ncy);
    Parallel_Common::bcast_int(ncz);
    Parallel_Common::bcast_int(nx);
    Parallel_Common::bcast_int(ny);
    Parallel_Common::bcast_int(nz);
    Parallel_Common::bcast_int(bx);
    Parallel_Common::bcast_int(by);
    Parallel_Common::bcast_int(bz);
    Parallel_Common::bcast_int(ndx);
    Parallel_Common::bcast_int(ndy);
    Parallel_Common::bcast_int(ndz);
    Parallel_Common::bcast_double(erf_ecut);
    Parallel_Common::bcast_double(erf_height);
    Parallel_Common::bcast_double(erf_sigma);

    Parallel_Common::bcast_int(diago_proc); // mohan add 2012-01-03
    Parallel_Common::bcast_int(pw_diag_nmax);
    Parallel_Common::bcast_int(diago_cg_prec);
    Parallel_Common::bcast_int(pw_diag_ndim);
    Parallel_Common::bcast_double(pw_diag_thr);
    Parallel_Common::bcast_int(nb2d);
    Parallel_Common::bcast_int(nurse);
    Parallel_Common::bcast_bool(colour);
    Parallel_Common::bcast_int(nbspline);
    Parallel_Common::bcast_bool(t_in_h);
    Parallel_Common::bcast_bool(vl_in_h);
    Parallel_Common::bcast_bool(vnl_in_h);
    Parallel_Common::bcast_bool(vh_in_h);
    Parallel_Common::bcast_bool(vion_in_h);

    Parallel_Common::bcast_bool(test_force);
    Parallel_Common::bcast_bool(test_stress);

    Parallel_Common::bcast_double(scf_thr);
    Parallel_Common::bcast_int(scf_thr_type);
    Parallel_Common::bcast_int(scf_nmax);
    Parallel_Common::bcast_int(this->relax_nmax);
    Parallel_Common::bcast_bool(out_stru); // mohan add 2012-03-23

    // Parallel_Common::bcast_string( occupations );
    Parallel_Common::bcast_string(smearing_method);
    Parallel_Common::bcast_double(smearing_sigma);

    Parallel_Common::bcast_string(mixing_mode);
    Parallel_Common::bcast_double(mixing_beta);
    Parallel_Common::bcast_int(mixing_ndim);
    Parallel_Common::bcast_double(mixing_gg0); // mohan add 2014-09-27
    Parallel_Common::bcast_double(mixing_beta_mag);
    Parallel_Common::bcast_double(mixing_gg0_mag);
    Parallel_Common::bcast_double(mixing_gg0_min);
    Parallel_Common::bcast_double(mixing_angle);
    Parallel_Common::bcast_bool(mixing_tau);
    Parallel_Common::bcast_bool(mixing_dftu);

    Parallel_Common::bcast_string(read_file_dir);
    Parallel_Common::bcast_string(init_wfc);
    Parallel_Common::bcast_bool(psi_initializer);

    Parallel_Common::bcast_int(mem_saver);
    Parallel_Common::bcast_int(printe);
    Parallel_Common::bcast_string(init_chg);
    Parallel_Common::bcast_string(chg_extrap); // xiaohui modify 2015-02-01
    Parallel_Common::bcast_int(out_freq_elec);
    Parallel_Common::bcast_int(out_freq_ion);
    Parallel_Common::bcast_bool(out_chg);
    Parallel_Common::bcast_bool(out_dm);
    Parallel_Common::bcast_bool(out_dm1);
    Parallel_Common::bcast_bool(out_bandgap); // for bandgap printing

    Parallel_Common::bcast_bool(deepks_out_labels); // caoyu added 2020-11-24, mohan modified 2021-01-03
    Parallel_Common::bcast_bool(deepks_scf);
    Parallel_Common::bcast_bool(deepks_bandgap);
    Parallel_Common::bcast_bool(deepks_out_unittest);
    Parallel_Common::bcast_string(deepks_model);

    Parallel_Common::bcast_int(out_pot);
    Parallel_Common::bcast_int(out_wfc_pw);
    Parallel_Common::bcast_bool(out_wfc_r);
    Parallel_Common::bcast_int(out_dos);
    if(GlobalV::MY_RANK != 0) out_band.resize(2); /* If this line is absent, will cause segmentation fault in io_input_test_para */
    Parallel_Common::bcast_int(out_band.data(), 2);
    Parallel_Common::bcast_bool(out_proj_band);
    if(GlobalV::MY_RANK != 0) out_mat_hs.resize(2); /* If this line is absent, will cause segmentation fault in io_input_test_para */
    Parallel_Common::bcast_int(out_mat_hs.data(), 2);
    Parallel_Common::bcast_bool(out_mat_hs2); // LiuXh add 2019-07-15
    Parallel_Common::bcast_bool(out_mat_t);
    Parallel_Common::bcast_bool(out_mat_dh);
    Parallel_Common::bcast_bool(out_mat_xc);
    Parallel_Common::bcast_bool(out_mat_r); // jingan add 2019-8-14
    Parallel_Common::bcast_int(out_wfc_lcao);
    Parallel_Common::bcast_bool(out_alllog);
    Parallel_Common::bcast_bool(out_element_info);
    Parallel_Common::bcast_bool(out_app_flag);
    Parallel_Common::bcast_int(out_ndigits);
    Parallel_Common::bcast_int(out_interval);

    Parallel_Common::bcast_double(dos_emin_ev);
    Parallel_Common::bcast_double(dos_emax_ev);
    Parallel_Common::bcast_double(dos_edelta_ev);
    Parallel_Common::bcast_double(dos_scale);
    Parallel_Common::bcast_bool(dos_setemin);
    Parallel_Common::bcast_bool(dos_setemax);
    Parallel_Common::bcast_int(dos_nche);
    Parallel_Common::bcast_double(dos_sigma);

    // mohan add 2009-11-11
    Parallel_Common::bcast_double(lcao_ecut);
    Parallel_Common::bcast_double(lcao_dk);
    Parallel_Common::bcast_double(lcao_dr);
    Parallel_Common::bcast_double(lcao_rmax);
    // zheng daye add 2014/5/5
    Parallel_Common::bcast_string(mdp.md_type);
    Parallel_Common::bcast_string(mdp.md_thermostat);
    Parallel_Common::bcast_int(mdp.md_nstep);
    Parallel_Common::bcast_double(mdp.md_dt);
    Parallel_Common::bcast_int(mdp.md_tchain);
    Parallel_Common::bcast_double(mdp.msst_qmass);
    Parallel_Common::bcast_double(mdp.md_tfirst);
    Parallel_Common::bcast_double(mdp.md_tlast);
    Parallel_Common::bcast_int(mdp.md_dumpfreq);
    Parallel_Common::bcast_int(mdp.md_restartfreq);
    Parallel_Common::bcast_int(mdp.md_seed);
    Parallel_Common::bcast_int(mdp.md_prec_level);
    Parallel_Common::bcast_bool(mdp.md_restart);
    Parallel_Common::bcast_double(mdp.lj_rcut);
    Parallel_Common::bcast_double(mdp.lj_epsilon);
    Parallel_Common::bcast_double(mdp.lj_sigma);
    Parallel_Common::bcast_int(mdp.msst_direction);
    Parallel_Common::bcast_double(mdp.msst_vel);
    Parallel_Common::bcast_double(mdp.msst_vis);
    Parallel_Common::bcast_double(mdp.msst_tscale);
    Parallel_Common::bcast_double(mdp.md_tfreq);
    Parallel_Common::bcast_double(mdp.md_damp);
    Parallel_Common::bcast_string(mdp.pot_file);
    Parallel_Common::bcast_int(mdp.md_nraise);
    Parallel_Common::bcast_bool(cal_syns);
    Parallel_Common::bcast_double(dmax);
    Parallel_Common::bcast_double(mdp.md_tolerance);
    Parallel_Common::bcast_string(mdp.md_pmode);
    Parallel_Common::bcast_string(mdp.md_pcouple);
    Parallel_Common::bcast_int(mdp.md_pchain);
    Parallel_Common::bcast_double(mdp.md_pfirst);
    Parallel_Common::bcast_double(mdp.md_plast);
    Parallel_Common::bcast_double(mdp.md_pfreq);
    Parallel_Common::bcast_bool(mdp.dump_force);
    Parallel_Common::bcast_bool(mdp.dump_vel);
    Parallel_Common::bcast_bool(mdp.dump_virial);
    // Yu Liu add 2022-05-18
    Parallel_Common::bcast_bool(efield_flag);
    Parallel_Common::bcast_bool(dip_cor_flag);
    Parallel_Common::bcast_int(efield_dir);
    Parallel_Common::bcast_double(efield_pos_max);
    Parallel_Common::bcast_double(efield_pos_dec);
    Parallel_Common::bcast_double(efield_amp);
    // Yu Liu add 2022-09-13
    Parallel_Common::bcast_bool(gate_flag);
    Parallel_Common::bcast_double(zgate);
    Parallel_Common::bcast_bool(relax);
    Parallel_Common::bcast_bool(block);
    Parallel_Common::bcast_double(block_down);
    Parallel_Common::bcast_double(block_up);
    Parallel_Common::bcast_double(block_height);
    /* 	// Peize Lin add 2014-04-07
        Parallel_Common::bcast_bool( vdwD2 );
        Parallel_Common::bcast_double( vdwD2_scaling );
        Parallel_Common::bcast_double( vdwD2_d );
        Parallel_Common::bcast_string( vdwD2_C6_file );
        Parallel_Common::bcast_string( vdwD2_C6_unit );
        Parallel_Common::bcast_string( vdwD2_R0_file );
        Parallel_Common::bcast_string( vdwD2_R0_unit );
        Parallel_Common::bcast_string( vdwD2_model );
        Parallel_Common::bcast_int( vdwD2_period.x );
        Parallel_Common::bcast_int( vdwD2_period.y );
        Parallel_Common::bcast_int( vdwD2_period.z );
        Parallel_Common::bcast_double( vdwD2_radius );
        Parallel_Common::bcast_string( vdwD2_radius_unit ); */
    // jiyy add 2019-08-04
    Parallel_Common::bcast_string(vdw_method);
    Parallel_Common::bcast_string(vdw_s6);
    Parallel_Common::bcast_string(vdw_s8);
    Parallel_Common::bcast_string(vdw_a1);
    Parallel_Common::bcast_string(vdw_a2);
    Parallel_Common::bcast_double(vdw_d);
    Parallel_Common::bcast_bool(vdw_abc);
    Parallel_Common::bcast_string(vdw_cutoff_radius);
    Parallel_Common::bcast_string(vdw_radius_unit);
    Parallel_Common::bcast_double(vdw_cn_thr);
    Parallel_Common::bcast_string(vdw_cn_thr_unit);
    Parallel_Common::bcast_string(vdw_C6_file);
    Parallel_Common::bcast_string(vdw_C6_unit);
    Parallel_Common::bcast_string(vdw_R0_file);
    Parallel_Common::bcast_string(vdw_R0_unit);
    Parallel_Common::bcast_string(vdw_cutoff_type);
    Parallel_Common::bcast_int(vdw_cutoff_period.x);
    Parallel_Common::bcast_int(vdw_cutoff_period.y);
    Parallel_Common::bcast_int(vdw_cutoff_period.z);
    // Fuxiang He add 2016-10-26
    Parallel_Common::bcast_double(td_force_dt);
    Parallel_Common::bcast_bool(td_vext);
    Parallel_Common::bcast_string(td_vext_dire);
    Parallel_Common::bcast_int(propagator);
    Parallel_Common::bcast_int(td_stype);
    Parallel_Common::bcast_string(td_ttype);
    Parallel_Common::bcast_int(td_tstart);
    Parallel_Common::bcast_int(td_tend);
    Parallel_Common::bcast_double(td_lcut1);
    Parallel_Common::bcast_double(td_lcut2);
    Parallel_Common::bcast_string(td_gauss_freq);
    Parallel_Common::bcast_string(td_gauss_phase);
    Parallel_Common::bcast_string(td_gauss_sigma);
    Parallel_Common::bcast_string(td_gauss_t0);
    Parallel_Common::bcast_string(td_gauss_amp);
    Parallel_Common::bcast_string(td_trape_freq);
    Parallel_Common::bcast_string(td_trape_phase);
    Parallel_Common::bcast_string(td_trape_t1);
    Parallel_Common::bcast_string(td_trape_t2);
    Parallel_Common::bcast_string(td_trape_t3);
    Parallel_Common::bcast_string(td_trape_amp);
    Parallel_Common::bcast_string(td_trigo_freq1);
    Parallel_Common::bcast_string(td_trigo_freq2);
    Parallel_Common::bcast_string(td_trigo_phase1);
    Parallel_Common::bcast_string(td_trigo_phase2);
    Parallel_Common::bcast_string(td_trigo_amp);
    Parallel_Common::bcast_string(td_heavi_t0);
    Parallel_Common::bcast_string(td_heavi_amp);
    // Parallel_Common::bcast_string(td_hhg_freq1);
    // Parallel_Common::bcast_string(td_hhg_freq2);
    // Parallel_Common::bcast_string(td_hhg_amp1);
    // Parallel_Common::bcast_string(td_hhg_amp2);
    // Parallel_Common::bcast_string(td_hhg_phase1);
    // Parallel_Common::bcast_string(td_hhg_phase2);
    // Parallel_Common::bcast_string(td_hhg_freq1);
    // Parallel_Common::bcast_string(td_hhg_freq2);
    // Parallel_Common::bcast_string(td_hhg_t0);
    // Parallel_Common::bcast_string(td_hhg_sigma);
    Parallel_Common::bcast_bool(out_dipole);
    Parallel_Common::bcast_bool(out_efield);
    Parallel_Common::bcast_double(td_print_eij);
    Parallel_Common::bcast_int(td_edm);
    Parallel_Common::bcast_bool(test_skip_ewald);
    Parallel_Common::bcast_bool(ocp);
    Parallel_Common::bcast_string(ocp_set);
    Parallel_Common::bcast_bool(out_mul); // qifeng add 2019/9/10

    // Peize Lin add 2018-06-20
    Parallel_Common::bcast_string(exx_hybrid_alpha);
    Parallel_Common::bcast_double(exx_hse_omega);
    Parallel_Common::bcast_bool(exx_separate_loop);
    Parallel_Common::bcast_int(exx_hybrid_step);
    Parallel_Common::bcast_double(exx_lambda);
    Parallel_Common::bcast_double(exx_mixing_beta);
    Parallel_Common::bcast_string(exx_real_number);
    Parallel_Common::bcast_double(exx_pca_threshold);
    Parallel_Common::bcast_double(exx_c_threshold);
    Parallel_Common::bcast_double(exx_v_threshold);
    Parallel_Common::bcast_double(exx_dm_threshold);
    Parallel_Common::bcast_double(exx_schwarz_threshold);
    Parallel_Common::bcast_double(exx_cauchy_threshold);
    Parallel_Common::bcast_double(exx_c_grad_threshold);
    Parallel_Common::bcast_double(exx_v_grad_threshold);
    Parallel_Common::bcast_double(exx_cauchy_force_threshold);
    Parallel_Common::bcast_double(exx_cauchy_stress_threshold);
    Parallel_Common::bcast_double(exx_ccp_threshold);
    Parallel_Common::bcast_string(exx_ccp_rmesh_times);
    Parallel_Common::bcast_string(exx_distribute_type);
    Parallel_Common::bcast_int(exx_opt_orb_lmax);
    Parallel_Common::bcast_double(exx_opt_orb_ecut);
    Parallel_Common::bcast_double(exx_opt_orb_tolerence);

    Parallel_Common::bcast_bool(noncolin);
    Parallel_Common::bcast_bool(lspinorb);
    Parallel_Common::bcast_double(soc_lambda);

    // Parallel_Common::bcast_int( epsilon0_choice );
    Parallel_Common::bcast_double(cell_factor); // LiuXh add 20180619
    Parallel_Common::bcast_bool(restart_save);  // Peize Lin add 2020.04.04
    Parallel_Common::bcast_bool(restart_load);  // Peize Lin add 2020.04.04

    Parallel_Common::bcast_int(band_print_num);
    if(GlobalV::MY_RANK != 0)
    {
        bands_to_print.resize(band_print_num);
    }

    for(int i = 0; i < band_print_num; i++)
    {
        Parallel_Common::bcast_int(bands_to_print[i]);
    }

    //-----------------------------------------------------------------------------------
    // DFT+U (added by Quxin 2020-10-29)
    //-----------------------------------------------------------------------------------
    Parallel_Common::bcast_bool(dft_plus_u);
    Parallel_Common::bcast_bool(yukawa_potential);
    Parallel_Common::bcast_int(omc);
    Parallel_Common::bcast_double(yukawa_lambda);
    if (GlobalV::MY_RANK != 0)
    {
        hubbard_u = new double[this->ntype];
        orbital_corr = new int[this->ntype];
    }

    for (int i = 0; i < this->ntype; i++)
    {
        Parallel_Common::bcast_double(hubbard_u[i]);
        Parallel_Common::bcast_int(orbital_corr[i]);
    }

    //-----------------------------------------------------------------------------------
    // DFT+DMFT (added by Quxin 2020-08)
    //-----------------------------------------------------------------------------------
    Parallel_Common::bcast_bool(dft_plus_dmft);

    //-----------------------------------------------------------------------------------
    // RPA
    //-----------------------------------------------------------------------------------
    Parallel_Common::bcast_bool(rpa);
    Parallel_Common::bcast_bool(GlobalV::rpa_setorb);

    //----------------------------------------------------------------------------------
    //    implicit solvation model        (sunml added on 2022-04-04)
    //----------------------------------------------------------------------------------
    Parallel_Common::bcast_bool(imp_sol);
    Parallel_Common::bcast_double(eb_k);
    Parallel_Common::bcast_double(tau);
    Parallel_Common::bcast_double(sigma_k);
    Parallel_Common::bcast_double(nc_k);

    //----------------------------------------------------------------------------------
    //    OFDFT sunliang added on 2022-05-05
    //----------------------------------------------------------------------------------
    Parallel_Common::bcast_string(of_kinetic);
    Parallel_Common::bcast_string(of_method);
    Parallel_Common::bcast_string(of_conv);
    Parallel_Common::bcast_double(of_tole);
    Parallel_Common::bcast_double(of_tolp);
    Parallel_Common::bcast_double(of_tf_weight);
    Parallel_Common::bcast_double(of_vw_weight);
    Parallel_Common::bcast_double(of_wt_alpha);
    Parallel_Common::bcast_double(of_wt_beta);
    Parallel_Common::bcast_double(of_wt_rho0);
    Parallel_Common::bcast_bool(of_hold_rho0);
    Parallel_Common::bcast_double(of_lkt_a);
    Parallel_Common::bcast_bool(of_full_pw);
    Parallel_Common::bcast_int(of_full_pw_dim);
    Parallel_Common::bcast_bool(of_read_kernel);
    Parallel_Common::bcast_string(of_kernel_file);

    //==========================================================
    // spherical bessel  Peize Lin added on 2022-12-15
    //==========================================================
    Parallel_Common::bcast_bool(bessel_nao_smooth);
    Parallel_Common::bcast_double(bessel_nao_sigma);
    Parallel_Common::bcast_string(bessel_nao_ecut);
    /* newly support vector/list input of bessel_nao_rcut */
    int nrcut = bessel_nao_rcuts.size();
    Parallel_Common::bcast_int(nrcut);
    if (nrcut != 0) /* as long as its value is really given, bcast, otherwise not */
    {
        bessel_nao_rcuts.resize(nrcut);
        Parallel_Common::bcast_double(bessel_nao_rcuts.data(), nrcut);
    }
    /* end */
    Parallel_Common::bcast_double(bessel_nao_rcut);
    Parallel_Common::bcast_double(bessel_nao_tolerence);
    Parallel_Common::bcast_int(bessel_descriptor_lmax);
    Parallel_Common::bcast_bool(bessel_descriptor_smooth);
    Parallel_Common::bcast_double(bessel_descriptor_sigma);
    Parallel_Common::bcast_string(bessel_descriptor_ecut);
    Parallel_Common::bcast_double(bessel_descriptor_rcut);
    Parallel_Common::bcast_double(bessel_descriptor_tolerence);
    //----------------------------------------------------------------------------------
    //    device control denghui added on 2022-11-05
    //----------------------------------------------------------------------------------
    Parallel_Common::bcast_string(device);
    /**
     *  Deltaspin variables
     */
    Parallel_Common::bcast_bool(sc_mag_switch);
    Parallel_Common::bcast_bool(decay_grad_switch);
    Parallel_Common::bcast_double(sc_thr);
    Parallel_Common::bcast_int(nsc);
    Parallel_Common::bcast_int(nsc_min);
    Parallel_Common::bcast_int(sc_scf_nmin);
    Parallel_Common::bcast_string(sc_file);
    Parallel_Common::bcast_double(alpha_trial);
    Parallel_Common::bcast_double(sccut);

    Parallel_Common::bcast_bool(qo_switch);
    Parallel_Common::bcast_string(qo_basis);
    Parallel_Common::bcast_double(qo_thr);
    /* broadcasting std::vector is sometime a annorying task... */
    if (ntype != 0) /* ntype has been broadcasted before */
    {
        qo_strategy.resize(ntype); 
        Parallel_Common::bcast_string(qo_strategy.data(), ntype);
        qo_screening_coeff.resize(ntype);
        Parallel_Common::bcast_double(qo_screening_coeff.data(), ntype);
    }
    return;
}
#endif

void Input::Check(void)
{
    ModuleBase::TITLE("Input", "Check");

    if (ecutrho / ecutwfc < 4 - 1e-8)
    {
        ModuleBase::WARNING_QUIT("Input", "ecutrho/ecutwfc must >= 4");
    }

    if (ndx < nx || ndy < ny || ndz < nz)
    {
        ModuleBase::WARNING_QUIT("Input", "smooth grids is denser than dense grids");
    }

    if (nbands < 0)
        ModuleBase::WARNING_QUIT("Input", "NBANDS must >= 0");
    //	if(nbands_istate < 0) ModuleBase::WARNING_QUIT("Input","NBANDS_ISTATE must > 0");
    if (nb2d < 0)
        ModuleBase::WARNING_QUIT("Input", "nb2d must > 0");
    if (ntype <= 0)
        ModuleBase::WARNING_QUIT("Input", "ntype must > 0");

        // std::cout << "diago_proc=" << diago_proc << std::endl;
        // std::cout << " NPROC=" << GlobalV::NPROC << std::endl;
#ifndef USE_PAW
    if (use_paw)
    {
        ModuleBase::WARNING_QUIT("Input", "to use PAW, compile with USE_PAW");
        if (basis_type != "pw")
        {
            ModuleBase::WARNING_QUIT("Input", "PAW is for pw basis only");
        }
    }
#endif

    if (diago_proc > 1 && basis_type == "lcao" && diago_proc != GlobalV::NPROC)
    {
        ModuleBase::WARNING_QUIT("Input", "please don't set diago_proc with lcao base");
    }
    int kspacing_zero_num = 0;
    for (int i = 0; i < 3; i++)
    {
        if (kspacing[i] < 0.0)
        {
            ModuleBase::WARNING_QUIT("Input", "kspacing must > 0");
        }
        else if (kspacing[i] == 0.0)
        {
            kspacing_zero_num++;
        }
    }
    if (kspacing_zero_num > 0 && kspacing_zero_num < 3)
    {
        std::cout << "kspacing: " << kspacing[0] << " " << kspacing[1] << " " << kspacing[2] << std::endl;
        ModuleBase::WARNING_QUIT("Input", "kspacing must > 0");
    }

    if (nelec < 0.0)
    {
        ModuleBase::WARNING_QUIT("Input", "nelec < 0 is not allowed !");
    }

    if (dip_cor_flag && !efield_flag)
    {
        ModuleBase::WARNING_QUIT("Input", "dipole correction is not active if efield_flag=false !");
    }

    if (gate_flag && efield_flag && !dip_cor_flag)
    {
        ModuleBase::WARNING_QUIT("Input", "gate field cannot be used with efield if dip_cor_flag=false !");
    }

    if (ref_cell_factor < 1.0)
    {
        ModuleBase::WARNING_QUIT("Input", "ref_cell_factor must not be less than 1.0");
    }

    //----------------------------------------------------------
    // main parameters / electrons / spin ( 1/16 )
    //----------------------------------------------------------
    if (calculation == "nscf" || calculation == "get_S")
    {
        if (out_dos == 3 && symmetry == "1")
        {
            ModuleBase::WARNING_QUIT("Input::Check",
                                     "symmetry can't be used for out_dos==3(Fermi Surface Plotting) by now.");
        }
    }
    else if (calculation == "get_pchg")
    {
        if (basis_type == "pw") // xiaohui add 2013-09-01
        {
            ModuleBase::WARNING_QUIT("Input::Check", "calculate = get_pchg is only availble for LCAO.");
        }
    }
    else if (calculation == "get_wf")
    {
        if (basis_type == "pw") // xiaohui add 2013-09-01
        {
            ModuleBase::WARNING_QUIT("Input::Check", "calculate = get_wf is only availble for LCAO.");
        }
    }
    else if (calculation == "md") // mohan add 2011-11-04
    {
        // deal with input parameters , 2019-04-30
        if (mdp.md_dt < 0)
            ModuleBase::WARNING_QUIT("Input::Check", "time interval of MD calculation should be set!");
        if (mdp.md_type == "npt" && mdp.md_pfirst < 0)
            ModuleBase::WARNING_QUIT("Input::Check", "pressure of MD calculation should be set!");
        if (mdp.md_type == "msst")
        {
            if (mdp.msst_qmass <= 0)
            {
                ModuleBase::WARNING_QUIT("Input::Check", "msst_qmass must be greater than 0!");
            }
        }
        if (esolver_type == "dp")
        {
            if (access(mdp.pot_file.c_str(), 0) == -1)
            {
                ModuleBase::WARNING_QUIT("Input::Check", "Can not find DP model !");
            }
        }
    }
    else if (calculation == "gen_bessel")
    {
        if (basis_type != "pw")
        {
            ModuleBase::WARNING_QUIT("Input", "to generate descriptors, please use pw basis");
        }
    }
    // else if (calculation == "ofdft") // sunliang added on 2022-05-05
    // {
    //     if (pseudo_type != "blps")
    //     {
    //         ModuleBase::WARNING_QUIT("Input::Check", "pseudo_type in ofdft should be set as blps");
    //     }
    // }
    else if (calculation != "scf" && calculation != "relax" && calculation != "cell-relax"
             && calculation != "test_memory" && calculation != "test_neighbour")
    {
        ModuleBase::WARNING_QUIT("Input", "check 'calculation' !");
    }
    if (init_chg != "atomic" && init_chg != "file")
    {
        ModuleBase::WARNING_QUIT("Input", "wrong 'init_chg',not 'atomic', 'file',please check");
    }
    if (gamma_only_local == 0)
    {
        if (out_dm == 1)
        {
            ModuleBase::WARNING_QUIT("Input", "out_dm with k-point algorithm is not implemented yet.");
        }
    }
    else
    {
        if (out_dm1 == 1)
        {
            ModuleBase::WARNING_QUIT("Input", "out_dm1 is only for multi-k");
        }
    }

    if (chg_extrap == "dm" && basis_type == "pw") // xiaohui add 2013-09-01, xiaohui modify 2015-02-01
    {
        ModuleBase::WARNING_QUIT(
            "Input",
            "wrong 'chg_extrap=dm' is only available for local orbitals."); // xiaohui modify 2015-02-01
    }

    if ((init_wfc != "atomic") && (init_wfc != "random") && (init_wfc != "atomic+random") && (init_wfc != "nao")
        && (init_wfc != "nao+random") && (init_wfc != "file"))
    {
        if (psi_initializer)
        {
            ModuleBase::WARNING_QUIT(
                "Input",
                "wrong init_wfc, please use 'random', 'atomic(+random)', 'nao(+random)' or 'file' ");
        }
        else
        {
            ModuleBase::WARNING_QUIT("Input", "wrong init_wfc, please use 'atomic' or 'random' or 'file' ");
        }
    }

    if (nbands > 100000)
    {
        ModuleBase::WARNING_QUIT("Input", "nbnd >100000, out of range");
    }
    if (nelec > 0 && nbands > 0 && nelec > 2 * nbands)
    {
        ModuleBase::WARNING_QUIT("Input", "nelec > 2*nbnd , bands not enough!");
    }
    if (nspin != 1 && nspin != 2 && nspin != 4)
    {
        ModuleBase::WARNING_QUIT("Input", "nspin does not equal to 1, 2, or 4!");
    }
    if (basis_type == "pw") // xiaohui add 2013-09-01
    {
        if (ks_solver == "genelpa") // yshen add 2016-07-20
        {
            ModuleBase::WARNING_QUIT("Input", "genelpa can not be used with plane wave basis.");
        }
        else if (ks_solver == "scalapack_gvx") // Peize Lin add 2020.11.14
        {
            ModuleBase::WARNING_QUIT("Input", "scalapack_gvx can not be used with plane wave basis.");
        }
        else if (ks_solver == "lapack")
        {
            ModuleBase::WARNING_QUIT("Input", "lapack can not be used with plane wave basis.");
        }
        else if (ks_solver != "default" && ks_solver != "cg" && ks_solver != "dav" && ks_solver != "bpcg")
        {
            ModuleBase::WARNING_QUIT("Input", "please check the ks_solver parameter!");
        }

        if (gamma_only)
        {
            ModuleBase::WARNING_QUIT("Input", "gamma_only not implemented for plane wave now.");
        }

        if (out_proj_band == 1)
        {
            ModuleBase::WARNING_QUIT("Input", "out_proj_band not implemented for plane wave now.");
        }

        if (out_dos == 3)
        {
            ModuleBase::WARNING_QUIT("Input", "Fermi Surface Plotting not implemented for plane wave now.");
        }

        if (sc_mag_switch)
        {
            ModuleBase::WARNING_QUIT("Input", "Non-colliner Spin-constrained DFT not implemented for plane wave now.");
        }
    }
    else if (basis_type == "lcao")
    {
        if (ks_solver == "cg")
        {
            ModuleBase::WARNING_QUIT("Input", "not ready for cg method in lcao ."); // xiaohui add 2013-09-04
        }
        else if (ks_solver == "cg_in_lcao")
        {
            GlobalV::ofs_warning << "cg_in_lcao is under testing" << std::endl;
        }
        else if (ks_solver == "genelpa")
        {
#ifndef __MPI
            ModuleBase::WARNING_QUIT("Input", "genelpa can not be used for series version.");
#endif
#ifndef __ELPA
            ModuleBase::WARNING_QUIT(
                "Input",
                "Can not use genelpa if abacus is not compiled with ELPA. Please change ks_solver to scalapack_gvx.");
#endif
        }
        else if (ks_solver == "scalapack_gvx")
        {
#ifdef __MPI
            GlobalV::ofs_warning << "scalapack_gvx is under testing" << std::endl;
#else
            ModuleBase::WARNING_QUIT("Input", "scalapack_gvx can not be used for series version.");
#endif
        }
        else if (ks_solver == "lapack")
        {
#ifdef __MPI
            ModuleBase::WARNING_QUIT("Input",
                                     "ks_solver=lapack is not an option for parallel version of ABACUS (try genelpa).");
#else
            GlobalV::ofs_warning << " It's ok to use lapack." << std::endl;
#endif
        }
        else if (ks_solver == "cusolver")
        {
#ifndef __MPI
            ModuleBase::WARNING_QUIT("Input", "Cusolver can not be used for series version.");
#endif
        }
        else if (ks_solver != "default")
        {
            ModuleBase::WARNING_QUIT("Input", "please check the ks_solver parameter!");
        }

        if (kpar > 1)
        {
            ModuleBase::WARNING_QUIT("Input", "kpar > 1 has not been supported for lcao calculation.");
        }

        if (out_wfc_lcao != 0 && out_wfc_lcao != 1 && out_wfc_lcao != 2)
        {
            ModuleBase::WARNING_QUIT("Input", "out_wfc_lcao must be 0, 1, or 2");
        }
    }
    else if (basis_type == "lcao_in_pw")
    {
        if (ks_solver != "lapack")
        {
            ModuleBase::WARNING_QUIT("Input", "LCAO in plane wave can only done with lapack.");
        }
    }
    else
    {
        ModuleBase::WARNING_QUIT("Input", "please check the basis_type parameter!");
    }
    /*
    if (basis_type == "lcao" && !gamma_only_local) // xiaohui add 2013-09-01. Attention! Maybe there is some problem.
    {
        ModuleBase::WARNING("Input", "gamma_only_local algorithm is not used.");
    }
    */

    /* comment out because code cannot reach here anyway
    if (GlobalV::NPROC > 1 && ks_solver == "lapack") // xiaohui add 2013-09-01
    {
        if (basis_type != "lcao_in_pw") // xiaohui add 2013-09-01
        {
            ModuleBase::WARNING_QUIT("Input", "lapack can not be used when nproc > 1");
        }
    }
    */
    // pengfei add 13-8-10 a new method cg to bfgs
    if (relax_method != "sd" && relax_method != "cg" && relax_method != "bfgs" && relax_method != "cg_bfgs")
    {
        ModuleBase::WARNING_QUIT("Input", "relax_method can only be sd, cg, bfgs or cg_bfgs.");
    }

    if (bx > 10 || by > 10 || bz > 10)
    {
        ModuleBase::WARNING_QUIT("Input", "bx, or by, or bz is larger than 10!");
    }
    // jiyy add 2019-08-04
    if (vdw_method == "d2" || vdw_method == "d3_0" || vdw_method == "d3_bj")
    {
        if ((vdw_C6_unit != "Jnm6/mol") && (vdw_C6_unit != "eVA6"))
        {
            ModuleBase::WARNING_QUIT("Input", "vdw_C6_unit must be Jnm6/mol or eVA6");
        }
        if ((vdw_R0_unit != "A") && (vdw_R0_unit != "Bohr"))
        {
            ModuleBase::WARNING_QUIT("Input", "vdw_R0_unit must be A or Bohr");
        }
        if ((vdw_cutoff_type != "radius") && (vdw_cutoff_type != "period"))
        {
            ModuleBase::WARNING_QUIT("Input", "vdw_cutoff_type must be radius or period");
        }
        if ((vdw_cutoff_period.x <= 0) || (vdw_cutoff_period.y <= 0) || (vdw_cutoff_period.z <= 0))
        {
            ModuleBase::WARNING_QUIT("Input", "vdw_cutoff_period <= 0 is not allowd");
        }
        if (std::stod(vdw_cutoff_radius) <= 0)
        {
            ModuleBase::WARNING_QUIT("Input", "vdw_cutoff_radius <= 0 is not allowd");
        }
        if ((vdw_radius_unit != "A") && (vdw_radius_unit != "Bohr"))
        {
            ModuleBase::WARNING_QUIT("Input", "vdw_radius_unit must be A or Bohr");
        }
        if (vdw_cn_thr <= 0)
        {
            ModuleBase::WARNING_QUIT("Input", "vdw_cn_thr <= 0 is not allowd");
        }
        if ((vdw_cn_thr_unit != "A") && (vdw_cn_thr_unit != "Bohr"))
        {
            ModuleBase::WARNING_QUIT("Input", "vdw_cn_thr_unit must be A or Bohr");
        }
    }

    std::string dft_functional_lower = dft_functional;
    std::transform(dft_functional.begin(), dft_functional.end(), dft_functional_lower.begin(), tolower);
    if (dft_functional_lower == "hf" || dft_functional_lower == "pbe0" || dft_functional_lower == "hse"
        || dft_functional_lower == "scan0")
    {
        const double exx_hybrid_alpha_value = std::stod(exx_hybrid_alpha);
        if (exx_hybrid_alpha_value < 0 || exx_hybrid_alpha_value > 1)
        {
            ModuleBase::WARNING_QUIT("INPUT", "must 0 <= exx_hybrid_alpha <= 1");
        }
        if (exx_hybrid_step <= 0)
        {
            ModuleBase::WARNING_QUIT("INPUT", "must exx_hybrid_step > 0");
        }
        const double exx_ccp_rmesh_times_value = std::stod(exx_ccp_rmesh_times);
        if (exx_ccp_rmesh_times_value < 1)
        {
            ModuleBase::WARNING_QUIT("INPUT", "must exx_ccp_rmesh_times >= 1");
        }
        if (exx_distribute_type != "htime" && exx_distribute_type != "kmeans2" && exx_distribute_type != "kmeans1"
            && exx_distribute_type != "order")
        {
            ModuleBase::WARNING_QUIT("INPUT", "exx_distribute_type must be htime or kmeans2 or kmeans1");
        }
    }
    if (dft_functional_lower == "opt_orb")
    {
        if (exx_opt_orb_lmax < 0)
        {
            ModuleBase::WARNING_QUIT("INPUT", "exx_opt_orb_lmax must >=0");
        }
        if (exx_opt_orb_ecut < 0)
        {
            ModuleBase::WARNING_QUIT("INPUT", "exx_opt_orb_ecut must >=0");
        }
        if (exx_opt_orb_tolerence < 0)
        {
            ModuleBase::WARNING_QUIT("INPUT", "exx_opt_orb_tolerence must >=0");
        }
    }

    if (berry_phase)
    {
        if (basis_type != "pw" && basis_type != "lcao")
        {
            ModuleBase::WARNING_QUIT("Input", "calculate berry phase, please set basis_type = pw or lcao");
        }
        if (calculation != "nscf")
        {
            ModuleBase::WARNING_QUIT("Input", "calculate berry phase, please set calculation = nscf");
        }
        if (!(gdir == 1 || gdir == 2 || gdir == 3))
        {
            ModuleBase::WARNING_QUIT("Input", "calculate berry phase, please set gdir = 1 or 2 or 3");
        }
    }

    if (towannier90)
    {
        if (basis_type == "lcao_in_pw")
        {
            /*
                Developer's notes: on the repair of lcao_in_pw

                lcao_in_pw is a special basis_type, for scf calculation, it follows workflow of pw,
                but for nscf the toWannier90 calculation, the interface is in ESolver_KS_LCAO_elec,
                therefore lcao_in_pw for towannier90 calculation follows lcao.

                In the future lcao_in_pw will have its own ESolver.

                2023/12/22 use new psi_initializer to expand numerical atomic orbitals, ykhuang
            */
            basis_type = "lcao";
            wannier_method = 1; // it is the way to call toWannier90_lcao_in_pw
#ifdef __ELPA
            ks_solver = "genelpa";
#else
            ks_solver = "scalapack_gvx";
#endif
        }
        else
        {
            if ((basis_type != "pw") && (basis_type != "lcao"))
            {
                ModuleBase::WARNING_QUIT("Input", "to use towannier90, please set basis_type = pw, lcao or lcao_in_pw");
            }
        }
        if (calculation != "nscf")
        {
            ModuleBase::WARNING_QUIT("Input", "to use towannier90, please set calculation = nscf");
        }
        if (nspin == 2)
        {
            if (!(wannier_spin == "up" || wannier_spin == "down"))
            {
                ModuleBase::WARNING_QUIT("Input", "to use towannier90, please set wannier_spin = up or down");
            }
        }
    }

    if (read_file_dir != "auto")
    {
        const std::string ss = "test -d " + read_file_dir;
        if (system(ss.c_str()))
        {
            ModuleBase::WARNING_QUIT("Input", "please set right files directory for reading in.");
        }
    }

    if (true) // Numerical_Basis::output_overlap()
    {
        if (std::stod(bessel_nao_ecut) < 0)
        {
            ModuleBase::WARNING_QUIT("INPUT", "bessel_nao_ecut must >=0");
        }
        if (bessel_nao_rcut < 0)
        {
            ModuleBase::WARNING_QUIT("INPUT", "bessel_nao_rcut must >=0");
        }
    }
    if (true) // Numerical_Descriptor::output_descriptor()
    {
        if (std::stod(bessel_descriptor_ecut) < 0)
        {
            ModuleBase::WARNING_QUIT("INPUT", "bessel_descriptor_ecut must >=0");
        }
        if (bessel_descriptor_rcut < 0)
        {
            ModuleBase::WARNING_QUIT("INPUT", "bessel_descriptor_rcut must >=0");
        }
    }

    // Deltaspin variables checking
    if (sc_mag_switch)
    {
        if (sc_file == "none")
        {
            ModuleBase::WARNING_QUIT("INPUT", "sc_file (json format) must be set when sc_mag_switch > 0");
        }
        else
        {
            const std::string ss = "test -f " + sc_file;
            if (system(ss.c_str()))
            {
                ModuleBase::WARNING_QUIT("INPUT", "sc_file does not exist");
            }
        }
        if (nspin != 4 && nspin != 2)
        {
            ModuleBase::WARNING_QUIT("INPUT", "nspin must be 2 or 4 when sc_mag_switch > 0");
        }
        if (calculation != "scf")
        {
            ModuleBase::WARNING_QUIT("INPUT", "calculation must be scf when sc_mag_switch > 0");
        }
        if (sc_thr <= 0)
        {
            ModuleBase::WARNING_QUIT("INPUT", "sc_thr must > 0");
        }
        if (nsc <= 0)
        {
            ModuleBase::WARNING_QUIT("INPUT", "nsc must > 0");
        }
        if (nsc_min <= 0)
        {
            ModuleBase::WARNING_QUIT("INPUT", "nsc_min must > 0");
        }
        if (sc_scf_nmin < 2)
        {
            ModuleBase::WARNING_QUIT("INPUT", "sc_scf_nmin must >= 2");
        }
        if (alpha_trial <= 0)
        {
            ModuleBase::WARNING_QUIT("INPUT", "alpha_trial must > 0");
        }
        if (sccut <= 0)
        {
            ModuleBase::WARNING_QUIT("INPUT", "sccut must > 0");
        }
    }
    if(qo_switch)
    {
        /* first about rationality of parameters */
        if(qo_basis == "pswfc")
        {
            for(auto screen_coeff: qo_screening_coeff)
            {
                if(screen_coeff < 0)
                {
                    ModuleBase::WARNING_QUIT("INPUT", "screening coefficient must >= 0 to tune the pswfc decay");
                }
                if(std::fabs(screen_coeff) < 1e-6)
                {
                    ModuleBase::WARNING("INPUT", "every low screening coefficient might yield very high computational cost");
                }
            }
        }
        else if(qo_basis == "hydrogen")
        {
            if(qo_thr > 1e-6)
            {
                ModuleBase::WARNING("INPUT", "too high the convergence threshold might yield unacceptable result");
            }
        }
        /* then size of std::vector<> parameters */
        if(qo_screening_coeff.size() != ntype) ModuleBase::WARNING_QUIT("INPUT", "qo_screening_coeff.size() != ntype");
        if(qo_strategy.size() != ntype) ModuleBase::WARNING_QUIT("INPUT", "qo_strategy.size() != ntype");
    }

    return;
}

void Input::close_log(void) const
{

    ModuleBase::Global_File::close_all_log(GlobalV::MY_RANK, this->out_alllog);
}

void Input::read_bool(std::ifstream& ifs, bool& var)
{
    std::string str;
    ifs >> str;
    for (auto& i: str)
    {
        i = tolower(i);
    }
    if (str == "true")
    {
        var = true;
    }
    else if (str == "false")
    {
        var = false;
    }
    else if (str == "1")
    {
        var = true;
    }
    else if (str == "0")
    {
        var = false;
    }
    else if (str == "t")
    {
        var = true;
    }
    else if (str == "f")
    {
        var = false;
    }
    else
    {
        std::string warningstr = "Bad boolean parameter ";
        warningstr.append(str);
        warningstr.append(", please check the input parameters in file INPUT");
        ModuleBase::WARNING_QUIT("Input", warningstr);
    }
    ifs.ignore(150, '\n');
    return;
}

void Input::strtolower(char* sa, char* sb)
{
    char c;
    int len = strlen(sa);
    for (int i = 0; i < len; i++)
    {
        c = sa[i];
        sb[i] = tolower(c);
    }
    sb[len] = '\0';
}

template <typename T>
void Input::read_value2stdvector(std::ifstream& ifs, std::vector<T>& var)
{
    // reset var
    var.clear(); var.shrink_to_fit();
    std::string line;
    std::getline(ifs, line); // read the whole rest of line
    line = (line.find('#') == std::string::npos) ? line : line.substr(0, line.find('#')); // remove comments
    std::vector<std::string> tmp;
    std::string::size_type start = 0, end = 0;
    while ((start = line.find_first_not_of(" \t\n", end)) != std::string::npos) // find the first not of delimiters but not reaches the end
    {
        end = line.find_first_of(" \t\n", start); // find the first of delimiters starting from start pos
        tmp.push_back(line.substr(start, end - start)); // push back the substring
    }
    var.resize(tmp.size());
    // capture "this"'s member function cast_string and iterate from tmp.begin() to tmp.end(), transform to var.begin()
    std::transform(tmp.begin(), tmp.end(), var.begin(), [this](const std::string& s) { return cast_string<T>(s); });
}
template void Input::read_value2stdvector(std::ifstream& ifs, std::vector<int>& var);
template void Input::read_value2stdvector(std::ifstream& ifs, std::vector<double>& var);
template void Input::read_value2stdvector(std::ifstream& ifs, std::vector<std::string>& var);

// Conut how many types of atoms are listed in STRU
int Input::count_ntype(const std::string& fn)
{
    // Only RANK0 core can reach here, because this function is called during Input::Read.
    assert(GlobalV::MY_RANK == 0);

    std::ifstream ifa(fn.c_str(), std::ios::in);
    if (!ifa)
    {
        GlobalV::ofs_warning << fn;
        ModuleBase::WARNING_QUIT("Input::count_ntype", "Can not find the file containing atom positions.!");
    }

    int ntype_stru = 0;
    std::string temp;
    if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifa, "ATOMIC_SPECIES"))
    {
        while (true)
        {
            ModuleBase::GlobalFunc::READ_VALUE(ifa, temp);
            if (temp == "LATTICE_CONSTANT" || temp == "NUMERICAL_ORBITAL" || temp == "NUMERICAL_DESCRIPTOR"
                || temp == "PAW_FILES" || ifa.eof())
            {
                break;
            }
            else if (isalpha(temp[0]))
            {
                ntype_stru += 1;
            }
        }
    }
    return ntype_stru;
}
