//==========================================================
// Author: Lixin He,mohan
// DATE : 2008-11-6
//==========================================================
//#include "global.h"
#include "input.h"

#include "module_base/global_file.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/timer.h"
#include "src_parallel/parallel_common.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include <vector>
Input INPUT;

void Input::Init(const std::string &fn)
{
    ModuleBase::timer::tick("Input", "Init");
    this->Default();

    bool success = this->Read(fn);
    this->Default_2();

// xiaohui add 2015-09-16
#ifdef __MPI
    Parallel_Common::bcast_bool(input_error);
#endif
    if (input_error == 1)
    {
        ModuleBase::WARNING_QUIT("Input", "Bad parameter, please check the input parameters in file INPUT");
    }

#ifdef __MPI
    Parallel_Common::bcast_bool(success);
#endif
    if (!success)
    {
        ModuleBase::WARNING_QUIT("Input::Init", "Error during readin parameters.");
    }
#ifdef __MPI
    Bcast();
#endif

    // mohan move forward 2011-02-26
    //----------------------------------------------------------
    // OTHRE CLASS MEMBER FUNCTION :
    // NAME : Run::make_dir( dir name : OUT.suffix)
    //----------------------------------------------------------
    ModuleBase::Global_File::make_dir_out(this->suffix,
                                          this->calculation,
                                          this->out_mat_hs2,
                                          GlobalV::MY_RANK,
                                          this->mdp.md_restart,
                                          this->out_alllog); // xiaohui add 2013-09-01
    Check();

    time_t time_now = time(NULL);
    GlobalV::ofs_running << "                                                                                     "
                         << std::endl;
    GlobalV::ofs_running << "                             WELCOME TO ABACUS                                       "
                         << std::endl;
    GlobalV::ofs_running << "                                                                                     "
                         << std::endl;
    GlobalV::ofs_running << "               'Atomic-orbital Based Ab-initio Computation at UStc'                  "
                         << std::endl;
    GlobalV::ofs_running << "                                                                                     "
                         << std::endl;
    GlobalV::ofs_running << "                     Website: http://abacus.ustc.edu.cn/                             "
                         << std::endl;
    GlobalV::ofs_running << "                                                                                     "
                         << std::endl;

    GlobalV::ofs_running << std::setiosflags(ios::right);

#ifdef __MPI
    // GlobalV::ofs_running << "    Version: Parallel, under ALPHA test" << std::endl;
    GlobalV::ofs_running << "    Version: Parallel, in development" << std::endl;
    GlobalV::ofs_running << "    Processor Number is " << GlobalV::NPROC << std::endl;
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

    GlobalV::ofs_running << std::setiosflags(ios::left);
    std::cout << std::setiosflags(ios::left);

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
    stru_file = ""; // xiaohui modify 2015-02-01
    kpoint_file = ""; // xiaohui modify 2015-02-01
    pseudo_dir = "";
    orbital_dir = ""; // liuyu add 2021-08-14
    read_file_dir = "auto";
    // pseudo_type = "auto"; // mohan add 2013-05-20 (xiaohui add 2013-06-23)
    wannier_card = "";
    latname = "test";
    // xiaohui modify 2015-09-15, relax -> scf
    // calculation = "relax";
    calculation = "scf";
    pseudo_rcut = 15.0; // qianrui add this parameter 2021-5
    pseudo_mesh = false; // qianrui add this pararmeter
    ntype = 0;
    nbands = 0;
    nbands_sto = 256;
    nbands_istate = 5;
    pw_seed = 1;
    nche_sto = 100;
    seed_sto = 0;
    bndpar = 1;
    kpar = 1;
    initsto_freq = 0;
    method_sto = 2;
    npart_sto = 1;
    cal_cond = false;
    dos_nche = 100;
    cond_nche = 20;
    cond_dw = 0.1;
    cond_wcut = 10;
    cond_wenlarge = 10;
    cond_fwhm = 0.3;
    berry_phase = false;
    gdir = 3;
    towannier90 = false;
    NNKP = "seedname.nnkp";
    wannier_spin = "up";
    kspacing = 0.0;
    min_dist_coef = 0.2;
    //----------------------------------------------------------
    // electrons / spin
    //----------------------------------------------------------
    dft_functional = "default";
    xc_temperature = 0.0;
    nspin = 1;
    nelec = 0.0;
    lmaxmax = 2;
    tot_magnetization = 0.0;
    //----------------------------------------------------------
    // new function
    //----------------------------------------------------------
    // local_basis=0; xiaohui modify 2013-09-01
    // linear_scaling=false; xiaohui modify 2013-09-01
    basis_type = "pw"; // xiaohui add 2013-09-01
    ks_solver = "default"; // xiaohui add 2013-09-01
    search_radius = -1.0; // unit: a.u. -1.0 has no meaning.
    search_pbc = true;
    symmetry = 0;
    init_vel = false;
    symmetry_prec = 1.0e-5; // LiuXh add 2021-08-12, accuracy for symmetry
    cal_force = 0;
    out_force = false;
    force_thr = 1.0e-3;
    force_thr_ev2 = 0;
    stress_thr = 1.0e-2; // LiuXh add 20180515
    press1 = 0.0;
    press2 = 0.0;
    press3 = 0.0;
    cal_stress = false;
    fixed_axes = "None"; // pengfei 2018-11-9
    relax_method = "cg"; // pengfei  2014-10-13
    relax_cg_thr = 0.5; // pengfei add 2013-08-15
    out_level = "ie";
    out_md_control = false;
    relax_bfgs_w1 = 0.01; // mohan add 2011-03-13
    relax_bfgs_w2 = 0.5;
    relax_bfgs_rmax = 0.8; // bohr
    relax_bfgs_rmin = 1e-5;
    relax_bfgs_init = 0.5; // bohr
    nbspline = -1;
    //----------------------------------------------------------
    // ecutwfc
    //----------------------------------------------------------
    // gamma_only = false;
    gamma_only = false;
    gamma_only_local = false;
    ecutwfc = 0.0;
    ecutrho = 0.0;
    ncx = 0;
    ncy = 0;
    ncz = 0;
    nx = 0;
    ny = 0;
    nz = 0;
    bx = 2;
    by = 2;
    bz = 2;
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
    scf_thr = 1.0e-9;
    scf_nmax = 40;
    relax_nmax = 0;
    out_stru = 0;
    //----------------------------------------------------------
    // occupation
    //----------------------------------------------------------
    occupations = "smearing"; // pengfei 2014-10-13
    smearing_method = "fixed";
    smearing_sigma = 0.01;
    //----------------------------------------------------------
    //  charge mixing
    //----------------------------------------------------------
    mixing_mode = "pulay";
    mixing_beta = 0.7;
    mixing_ndim = 8;
    mixing_gg0 = 0.00; // used in kerker method. mohan add 2014-09-27
    //----------------------------------------------------------
    // potential / charge / wavefunction / energy
    //----------------------------------------------------------
    init_wfc = "atomic";
    mem_saver = 0;
    printe = 100; // must > 0
    init_chg = "atomic";
    chg_extrap = "atomic"; // xiaohui modify 2015-02-01
    out_freq_elec = 0;
    out_freq_ion = 0;
    out_chg = 0;
    out_dm = 0;

    deepks_out_labels = 0; // caoyu added 2020-11-24, mohan added 2021-01-03
    deepks_scf = 0;
    deepks_bandgap = 0;
    deepks_out_unittest = 0;
    bessel_lmax = 2; // mohan added 2021-01-03
    bessel_rcut = 6.0;
    bessel_tol = 1.0e-12;

    out_pot = 0;
    out_wfc_pw = 0;
    out_wfc_r = 0;
    out_dos = 0;
    out_band = 0;
    out_proj_band = 0;
    out_mat_hs = 0;
    out_mat_hs2 = 0; // LiuXh add 2019-07-15
    out_mat_r = 0; // jingan add 2019-8-14
    out_wfc_lcao = false;
    out_alllog = false;
    dos_emin_ev = -15; //(ev)
    dos_emax_ev = 15; //(ev)
    dos_edelta_ev = 0.01; //(ev)
    dos_scale = 0.01;
    b_coef = 0.07;
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
    efield_pos_max = 0.5;
    efield_pos_dec = 0.1;
    efield_amp  = 0.0;
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

    exx_hybrid_alpha = 0.25;
    exx_hse_omega = 0.11;

    exx_separate_loop = true;
    exx_hybrid_step = 100;

    exx_lambda = 0.3;

    exx_pca_threshold = 0;
    exx_c_threshold = 0;
    exx_v_threshold = 0;
    exx_dm_threshold = 0;
    exx_schwarz_threshold = 0;
    exx_cauchy_threshold = 0;
    exx_ccp_threshold = 1E-8;
    exx_ccp_rmesh_times = 10;

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
    tddft = 0;
    td_scf_thr = 1e-9;
    td_dt = 0.02;
    td_force_dt = 0.02;
    td_val_elec_01 = 1;
    td_val_elec_02 = 1;
    td_val_elec_03 = 1;
    td_vext = 0;
    td_vext_dire = 1;

    td_timescale = 0.5;
    td_vexttype = 1;
    td_vextout = 0;
    td_dipoleout = 0;

    //----------------------------------------------------------			//Fuxiang He add 2016-10-26
    // constrained DFT
    //----------------------------------------------------------
    GlobalV::ocp = 0;
    // ocp_n = 0;
    GlobalV::ocp_set = "none";
    // for(int i=0; i<10000; i++)
    // {
    // GlobalV::ocp_kb[i] = 0.0;
    // }

    cell_factor = 1.2; // LiuXh add 20180619

    out_mul = 0; // qi feng add 2019/9/10

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
    dft_plus_u = false; // 1:DFT+U correction; 0：standard DFT calcullation
    yukawa_potential = false;
    yukawa_lambda = -1.0;
    double_counting = 1;
    omc = false;
    dftu_type = 2;

    //==========================================================
    //    DFT+DMFT     Xin Qu added on 2020-08
    //==========================================================
    dft_plus_dmft = false;

    //==========================================================
    //    implicit solvation model       sunml added on 2022-04-04
    //==========================================================
    imp_sol = 0;
    eb_k = 80.0;
    tau = 1.0798 * 1e-5;
    sigma_k = 0.6;
    nc_k = 0.00037;

    return;
}

bool Input::Read(const std::string &fn)
{
    ModuleBase::TITLE("Input", "Read");

    if (GlobalV::MY_RANK != 0)
        return false;

    std::ifstream ifs(fn.c_str(), ios::in); // "in_datas/input_parameters"

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
        std::cout << " Error parameter list." << std::endl;
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
            read_value(ifs, pseudo_mesh);
        }
        else if (strcmp("calculation", word) == 0) // which type calculation
        {
            read_value(ifs, calculation);
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
            read_value(ifs, nbands_sto);
        }
        else if (strcmp("kspacing", word) == 0)
        {
            read_value(ifs, kspacing);
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
            read_value(ifs, cal_cond);
        }
        else if (strcmp("cond_nche", word) == 0)
        {
            read_value(ifs, cond_nche);
        }
        else if (strcmp("cond_dw", word) == 0)
        {
            read_value(ifs, cond_dw);
        }
        else if (strcmp("cond_wcut", word) == 0)
        {
            read_value(ifs, cond_wcut);
        }
        else if (strcmp("cond_wenlarge", word) == 0)
        {
            read_value(ifs, cond_wenlarge);
        }
        else if (strcmp("cond_fwhm", word) == 0)
        {
            read_value(ifs, cond_fwhm);
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
            read_value(ifs, berry_phase);
        }
        else if (strcmp("gdir", word) == 0) // berry phase calculation
        {
            read_value(ifs, gdir);
        }
        else if (strcmp("towannier90", word) == 0) // add by jingan for wannier90
        {
            read_value(ifs, towannier90);
        }
        else if (strcmp("nnkpfile", word) == 0) // add by jingan for wannier90
        {
            read_value(ifs, NNKP);
        }
        else if (strcmp("wannier_spin", word) == 0) // add by jingan for wannier90
        {
            read_value(ifs, wannier_spin);
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
        else if (strcmp("lmaxmax", word) == 0)
        {
            read_value(ifs, lmaxmax);
        }

        else if (strcmp("tot_magnetization", word) == 0)
        {
            read_value(ifs, tot_magnetization);
        }
        //----------------------------------------------------------
        // new function
        //----------------------------------------------------------
        // else if (strcmp("local_basis", word) == 0)
        //{
        //    read_value(ifs, local_basis);
        //} xiaohui modify 2013-09-01
        else if (strcmp("basis_type", word) == 0)
        {
            read_value(ifs, basis_type);
        } // xiaohui add 2013-09-01
        // else if (strcmp("linear_scaling", word) == 0)
        //{
        //     read_value(ifs, linear_scaling);
        // } xiaohui modify 2013-09-01
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
            read_value(ifs, search_pbc);
        }
        else if (strcmp("symmetry", word) == 0)
        {
            read_value(ifs, symmetry);
        }
        else if (strcmp("init_vel", word) == 0)
        {
            read_value(ifs, init_vel);
        }
        else if (strcmp("symmetry_prec", word) == 0) // LiuXh add 2021-08-12, accuracy for symmetry
        {
            read_value(ifs, symmetry_prec);
        }
        else if (strcmp("cal_force", word) == 0)
        {
            read_value(ifs, cal_force);
        }
        else if (strcmp("out_force", word) == 0)
        {
            read_value(ifs, out_force);
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
            read_value(ifs, cal_stress);
        }
        else if (strcmp("fixed_axes", word) == 0)
        {
            read_value(ifs, fixed_axes);
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
        //		else if (strcmp("gauss_pao_flag", word) == 0)
        //		else if (strcmp("gauss_pao_flag", word) == 0)
        //		else if (strcmp("gauss_pao_flag", word) == 0)
        //		{
        //			read_value(ifs, gauss_PAO_flag);
        //		}
        //----------------------------------------------------------
        // plane waves
        //----------------------------------------------------------
        else if (strcmp("gamma_only", word) == 0)
        {
            read_value(ifs, gamma_only);
        }
        else if (strcmp("ecutwfc", word) == 0)
        {
            read_value(ifs, ecutwfc);
            ecutrho = 4.0 * ecutwfc;
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
            read_value(ifs, colour);
        }
        else if (strcmp("nbspline", word) == 0)
        {
            read_value(ifs, nbspline);
        }
        else if (strcmp("t_in_h", word) == 0)
        {
            read_value(ifs, t_in_h);
        }
        else if (strcmp("vl_in_h", word) == 0)
        {
            read_value(ifs, vl_in_h);
        }
        else if (strcmp("vnl_in_h", word) == 0)
        {
            read_value(ifs, vnl_in_h);
        }
        else if (strcmp("vh_in_h", word) == 0)
        {
            read_value(ifs, vh_in_h);
        }
        else if (strcmp("vion_in_h", word) == 0)
        {
            read_value(ifs, vion_in_h);
        }
        else if (strcmp("test_force", word) == 0)
        {
            read_value(ifs, test_force);
        }
        else if (strcmp("test_stress", word) == 0)
        {
            read_value(ifs, test_stress);
        }
        //----------------------------------------------------------
        // iteration
        //----------------------------------------------------------
        else if (strcmp("scf_thr", word) == 0)
        {
            read_value(ifs, scf_thr);
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
            read_value(ifs, out_stru);
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
            // 2015-06-15, xiaohui
            if (mixing_mode == "pulay-kerker")
            {
                mixing_gg0 = 1.5;
            }
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
            read_value(ifs, out_chg);
        }
        else if (strcmp("out_dm", word) == 0)
        {
            read_value(ifs, out_dm);
        }
        else if (strcmp("deepks_out_labels", word) == 0) // caoyu added 2020-11-24, mohan modified 2021-01-03
        {
            read_value(ifs, deepks_out_labels);
        }
        else if (strcmp("deepks_scf", word) == 0) // caoyu added 2020-11-24, mohan modified 2021-01-03
        {
            read_value(ifs, deepks_scf);
        }
        else if (strcmp("deepks_bandgap", word) == 0) // caoyu added 2020-11-24, mohan modified 2021-01-03
        {
            read_value(ifs, deepks_bandgap);
        }
        else if (strcmp("deepks_out_unittest", word) == 0) // mohan added 2021-01-03
        {
            read_value(ifs, deepks_out_unittest);
        }
        else if (strcmp("deepks_model", word) == 0) // caoyu added 2021-06-03
        {
            read_value(ifs, deepks_model);
        }
        else if (strcmp("bessel_lmax", word) == 0) // QO added 2021-12-15
        {
            read_value(ifs, bessel_lmax);
        }
        else if (strcmp("bessel_rcut", word) == 0) // QO added 2021-12-15
        {
            read_value(ifs, bessel_rcut);
        }
        else if (strcmp("bessel_tol", word) == 0) // QO added 2021-12-15
        {
            read_value(ifs, bessel_tol);
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
            read_value(ifs, out_wfc_r);
        }
        // mohan add 20090909
        else if (strcmp("out_dos", word) == 0)
        {
            read_value(ifs, out_dos);
        }
        else if (strcmp("out_band", word) == 0)
        {
            read_value(ifs, out_band);
        }
        else if (strcmp("out_proj_band", word) == 0)
        {
            read_value(ifs, out_proj_band);
        }

        else if (strcmp("out_mat_hs", word) == 0)
        {
            read_value(ifs, out_mat_hs);
        }
        // LiuXh add 2019-07-15
        else if (strcmp("out_mat_hs2", word) == 0)
        {
            read_value(ifs, out_mat_hs2);
        }
        else if (strcmp("out_mat_r", word) == 0)
        {
            read_value(ifs, out_mat_r);
        }
        else if (strcmp("out_wfc_lcao", word) == 0)
        {
            read_value(ifs, out_wfc_lcao);
        }
        else if (strcmp("out_alllog", word) == 0)
        {
            read_value(ifs, out_alllog);
        }
        else if (strcmp("out_element_info", word) == 0)
        {
            read_value(ifs, out_element_info);
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
            read_value(ifs, b_coef);
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
        // about molecular dynamics
        // added begin by zheng daye
        else if (strcmp("md_type", word) == 0)
        {
            read_value(ifs, mdp.md_type);
        }
        else if (strcmp("md_nstep", word) == 0)
        {
            read_value(ifs, mdp.md_nstep);
        }
        else if (strcmp("md_dt", word) == 0)
        {
            read_value(ifs, mdp.md_dt);
        }
        else if (strcmp("md_mnhc", word) == 0)
        {
            read_value(ifs, mdp.md_mnhc);
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
        else if (strcmp("md_restart", word) == 0)
        {
            read_value(ifs, mdp.md_restart);
        }
        // added by zheng daye
        //----------------------------------------------------------
        // Classic MD
        // Yu Liu add 2021-07-30
        //----------------------------------------------------------
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
        else if (strcmp("md_ensolver", word) == 0)
        {
            read_value(ifs, mdp.md_ensolver);
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
        //----------------------------------------------------------
        // efield and dipole correction
        // Yu Liu add 2022-05-18
        //----------------------------------------------------------
        else if (strcmp("efield_flag", word) == 0)
        {
            read_value(ifs, efield_flag);
        }
        else if (strcmp("dip_cor_flag", word) == 0)
        {
            read_value(ifs, dip_cor_flag);
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
            read_value(ifs, efield_amp );
        }
        //----------------------------------------------------------
        // gatefield (compensating charge)
        // Yu Liu add 2022-09-13
        //----------------------------------------------------------
        else if (strcmp("gate_flag", word) == 0)
        {
            read_value(ifs, gate_flag);
        }
        else if (strcmp("zgate", word) == 0)
        {
            read_value(ifs, zgate);
        }
        else if (strcmp("relax", word) == 0)
        {
            read_value(ifs, relax);
        }
        else if (strcmp("block", word) == 0)
        {
            read_value(ifs, block);
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
        else if (strcmp("tddft", word) == 0)
        {
            read_value(ifs, tddft);
        }
        else if (strcmp("td_scf_thr", word) == 0)
        {
            read_value(ifs, td_scf_thr);
        }
        else if (strcmp("td_dt", word) == 0)
        {
            read_value(ifs, td_dt);
        }
        else if (strcmp("td_force_dt", word) == 0)
        {
            read_value(ifs, td_force_dt);
        }
        else if (strcmp("td_val_elec_01", word) == 0)
        {
            read_value(ifs, td_val_elec_01);
        }
        else if (strcmp("td_val_elec_02", word) == 0)
        {
            read_value(ifs, td_val_elec_02);
        }
        else if (strcmp("td_val_elec_03", word) == 0)
        {
            read_value(ifs, td_val_elec_03);
        }
        else if (strcmp("td_vext", word) == 0)
        {
            read_value(ifs, td_vext);
        }
        else if (strcmp("td_vext_dire", word) == 0)
        {
            read_value(ifs, td_vext_dire);
        }
        else if (strcmp("td_timescale", word) == 0)
        {
            read_value(ifs, td_timescale);
        }
        else if (strcmp("td_vexttype", word) == 0)
        {
            read_value(ifs, td_vexttype);
        }
        else if (strcmp("td_vextout", word) == 0)
        {
            read_value(ifs, td_vextout);
        }
        else if (strcmp("td_dipoleout", word) == 0)
        {
            read_value(ifs, td_dipoleout);
        }
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
            read_value(ifs, vdw_abc);
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
            read_value(ifs, restart_save);
        }
        else if (strcmp("restart_load", word) == 0)
        {
            read_value(ifs, restart_load);
        }
        else if (strcmp("ocp", word) == 0)
        {
            read_value(ifs, GlobalV::ocp);
        }
        else if (strcmp("ocp_set", word) == 0)
        {
            getline(ifs, GlobalV::ocp_set);
            //			ifs.ignore(150, '\n');
        }
        else if (strcmp("out_mul", word) == 0)
        {
            read_value(ifs, out_mul);
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
            read_value(ifs, exx_separate_loop);
        }
        else if (strcmp("exx_hybrid_step", word) == 0)
        {
            read_value(ifs, exx_hybrid_step);
        }
        else if (strcmp("exx_lambda", word) == 0)
        {
            read_value(ifs, exx_lambda);
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
            read_value(ifs, noncolin);
        }
        else if (strcmp("lspinorb", word) == 0)
        {
            read_value(ifs, lspinorb);
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
            read_value(ifs, test_skip_ewald);
        }
        //--------------
        //----------------------------------------------------------------------------------
        //         Xin Qu added on 2020-10-29 for DFT+U
        //----------------------------------------------------------------------------------
        else if (strcmp("dft_plus_u", word) == 0)
        {
            ifs >> dft_plus_u;
        }
        else if (strcmp("dftu_type", word) == 0)
            ifs.ignore(150, '\n');
        else if (strcmp("yukawa_potential", word) == 0)
            ifs.ignore(150, '\n');
        else if (strcmp("hubbard_u", word) == 0)
            ifs.ignore(150, '\n');
        else if (strcmp("hund_j", word) == 0)
            ifs.ignore(150, '\n');
        else if (strcmp("double_counting", word) == 0)
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
            ifs >> dft_plus_dmft;
        }
        //----------------------------------------------------------------------------------
        //    implicit solvation model       sunml added on 2022-04-04
        //----------------------------------------------------------------------------------
        else if (strcmp("imp_sol", word) == 0)
        {
            read_value(ifs, imp_sol);
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

    //----------------------------------------------------------
    //       DFT+U    Xin Qu  added on 2020-10-29
    //----------------------------------------------------------
    hubbard_u = new double[ntype];
    for (int i = 0; i < ntype; i++)
    {
        hubbard_u[i] = 0.0;
    }

    hund_j = new double[ntype];
    for (int i = 0; i < ntype; i++)
    {
        hund_j[i] = 0.0;
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
            if(ifs.eof() != 0) break;
            strtolower(word1, word); // convert uppercase std::string to lower case; word1 --> word

            if (strcmp("dftu_type", word) == 0)
            {
                ifs >> dftu_type;
            }
            else if (strcmp("yukawa_potential", word) == 0)
            {
                ifs >> yukawa_potential;
            }
            else if (strcmp("yukawa_lambda", word) == 0)
            {
                ifs >> yukawa_lambda;
            }
            else if (strcmp("double_counting", word) == 0)
            {
                ifs >> double_counting;
            }
            else if (strcmp("hubbard_u", word) == 0)
            {
                for (int i = 0; i < ntype; i++)
                {
                    ifs >> hubbard_u[i];
                    hubbard_u[i] /= ModuleBase::Ry_to_eV;
                }
            }
            else if (strcmp("hund_j", word) == 0)
            {
                for (int i = 0; i < ntype; i++)
                {
                    ifs >> hund_j[i];
                    hund_j[i] /= ModuleBase::Ry_to_eV;
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
                ifs >> omc;
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

            if (hund_j[i] < -1.0e-3)
            {
                std::cout << " WRONG ARGUMENTS OF hund_j " << std::endl;
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

        if (strcmp("genelpa", ks_solver.c_str()) != 0 && strcmp(ks_solver.c_str(), "scalapack_gvx") != 0)
        {
            std::cout
                << " WRONG ARGUMENTS OF ks_solver in DFT+U routine, only genelpa and scalapack_gvx are supportted "
                << std::endl;
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
            else if (strcmp("hund_j", word) == 0)
            {
                for (int i = 0; i < ntype; i++)
                {
                    ifs >> hund_j[i];
                    hund_j[i] /= ModuleBase::Ry_to_eV;
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

            if (hund_j[i] < -1.0e-3)
            {
                std::cout << " WRONG ARGUMENTS OF hund_j " << std::endl;
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

        if (strcmp("genelpa", ks_solver.c_str()) != 0)
        {
            std::cout << " WRONG ARGUMENTS OF ks_solver in DFT+DMFT routine, only genelpa is support " << std::endl;
            exit(0);
        }
    }

    if (basis_type == "pw" && gamma_only !=0) // pengfei Li add 2015-1-31
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
    }
    if((out_mat_r || out_mat_hs2) && gamma_only_local)
    {
        ModuleBase::WARNING_QUIT("Input", "printing of H(R)/S(R)/r(R) is not available for gamma only calculations");
    }

    return true;
} // end read_parameters

void Input::Default_2(void) // jiyy add 2019-08-04
{
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
    if(calculation.substr(0,3) != "sto")    bndpar = 1;
    if(bndpar > GlobalV::NPROC) bndpar = GlobalV::NPROC;
    if(method_sto != 1 && method_sto != 2) 
    {
        method_sto = 2;
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
    Parallel_Common::bcast_double(pseudo_rcut);
    Parallel_Common::bcast_bool(pseudo_mesh);
    Parallel_Common::bcast_int(ntype);
    Parallel_Common::bcast_int(nbands);
    Parallel_Common::bcast_int(nbands_sto);
    Parallel_Common::bcast_int(nbands_istate);
    Parallel_Common::bcast_double(kspacing);
    Parallel_Common::bcast_double(min_dist_coef);
    Parallel_Common::bcast_int(nche_sto);
    Parallel_Common::bcast_int(seed_sto);
    Parallel_Common::bcast_int(pw_seed);
    Parallel_Common::bcast_double(emax_sto);
    Parallel_Common::bcast_double(emin_sto);
    Parallel_Common::bcast_int(initsto_freq);
    Parallel_Common::bcast_int(method_sto);
    Parallel_Common::bcast_int(npart_sto);
    Parallel_Common::bcast_bool(cal_cond);
    Parallel_Common::bcast_int(cond_nche);
    Parallel_Common::bcast_double(cond_dw);
    Parallel_Common::bcast_double(cond_wcut);
    Parallel_Common::bcast_int(cond_wenlarge);
    Parallel_Common::bcast_double(cond_fwhm);
    Parallel_Common::bcast_int(bndpar);
    Parallel_Common::bcast_int(kpar);
    Parallel_Common::bcast_bool(berry_phase);
    Parallel_Common::bcast_int(gdir);
    Parallel_Common::bcast_bool(towannier90);
    Parallel_Common::bcast_string(NNKP);
    Parallel_Common::bcast_string(wannier_spin);

    Parallel_Common::bcast_string(dft_functional);
    Parallel_Common::bcast_double(xc_temperature);
    Parallel_Common::bcast_int(nspin);
    Parallel_Common::bcast_double(nelec);
    Parallel_Common::bcast_int(lmaxmax);

    Parallel_Common::bcast_double(tot_magnetization);

    Parallel_Common::bcast_string(basis_type); // xiaohui add 2013-09-01
    Parallel_Common::bcast_string(ks_solver); // xiaohui add 2013-09-01
    Parallel_Common::bcast_double(search_radius);
    Parallel_Common::bcast_bool(search_pbc);
    Parallel_Common::bcast_double(search_radius);
    Parallel_Common::bcast_int(symmetry);
    Parallel_Common::bcast_bool(init_vel); // liuyu 2021-07-14
    Parallel_Common::bcast_double(symmetry_prec); // LiuXh add 2021-08-12, accuracy for symmetry
    Parallel_Common::bcast_int(cal_force);
    Parallel_Common::bcast_bool(out_force);
    Parallel_Common::bcast_double(force_thr);
    Parallel_Common::bcast_double(force_thr_ev2);
    Parallel_Common::bcast_double(stress_thr); // LiuXh add 20180515
    Parallel_Common::bcast_double(press1);
    Parallel_Common::bcast_double(press2);
    Parallel_Common::bcast_double(press3);
    Parallel_Common::bcast_bool(cal_stress);
    Parallel_Common::bcast_string(fixed_axes);
    Parallel_Common::bcast_string(relax_method);
    Parallel_Common::bcast_double(relax_cg_thr); // pengfei add 2013-08-15
    Parallel_Common::bcast_string(out_level);
    Parallel_Common::bcast_bool(out_md_control);
    Parallel_Common::bcast_double(relax_bfgs_w1);
    Parallel_Common::bcast_double(relax_bfgs_w2);
    Parallel_Common::bcast_double(relax_bfgs_rmax);
    Parallel_Common::bcast_double(relax_bfgs_rmin);
    Parallel_Common::bcast_double(relax_bfgs_init);

    Parallel_Common::bcast_bool(gamma_only);
    Parallel_Common::bcast_bool(gamma_only_local);
    Parallel_Common::bcast_double(ecutwfc);
    Parallel_Common::bcast_double(ecutrho);
    Parallel_Common::bcast_int(ncx);
    Parallel_Common::bcast_int(ncy);
    Parallel_Common::bcast_int(ncz);
    Parallel_Common::bcast_int(nx);
    Parallel_Common::bcast_int(ny);
    Parallel_Common::bcast_int(nz);
    Parallel_Common::bcast_int(bx);
    Parallel_Common::bcast_int(by);
    Parallel_Common::bcast_int(bz);

    Parallel_Common::bcast_int(diago_proc); // mohan add 2012-01-03
    Parallel_Common::bcast_int(pw_diag_nmax);
    Parallel_Common::bcast_int(diago_cg_prec);
    Parallel_Common::bcast_int(pw_diag_ndim);
    Parallel_Common::bcast_double(pw_diag_thr);
    Parallel_Common::bcast_int(nb2d);
    Parallel_Common::bcast_int(nurse);
    Parallel_Common::bcast_bool(colour);
    Parallel_Common::bcast_int(nbspline);
    Parallel_Common::bcast_int(t_in_h);
    Parallel_Common::bcast_int(vl_in_h);
    Parallel_Common::bcast_int(vnl_in_h);
    Parallel_Common::bcast_int(vh_in_h);
    Parallel_Common::bcast_int(vion_in_h);

    Parallel_Common::bcast_int(test_force);
    Parallel_Common::bcast_int(test_stress);

    Parallel_Common::bcast_double(scf_thr);
    Parallel_Common::bcast_int(scf_nmax);
    Parallel_Common::bcast_int(this->relax_nmax);
    Parallel_Common::bcast_int(out_stru); // mohan add 2012-03-23

    // Parallel_Common::bcast_string( occupations );
    Parallel_Common::bcast_string(smearing_method);
    Parallel_Common::bcast_double(smearing_sigma);

    Parallel_Common::bcast_string(mixing_mode);
    Parallel_Common::bcast_double(mixing_beta);
    Parallel_Common::bcast_int(mixing_ndim);
    Parallel_Common::bcast_double(mixing_gg0); // mohan add 2014-09-27

    Parallel_Common::bcast_string(read_file_dir);
    Parallel_Common::bcast_string(init_wfc);
    Parallel_Common::bcast_int(mem_saver);
    Parallel_Common::bcast_int(printe);
    Parallel_Common::bcast_string(init_chg);
    Parallel_Common::bcast_string(chg_extrap); // xiaohui modify 2015-02-01
    Parallel_Common::bcast_int(out_freq_elec);
    Parallel_Common::bcast_int(out_freq_ion);
    Parallel_Common::bcast_int(out_chg);
    Parallel_Common::bcast_int(out_dm);

    Parallel_Common::bcast_bool(deepks_out_labels); // caoyu added 2020-11-24, mohan modified 2021-01-03
    Parallel_Common::bcast_bool(deepks_scf);
    Parallel_Common::bcast_bool(deepks_bandgap);
    Parallel_Common::bcast_bool(deepks_out_unittest);
    Parallel_Common::bcast_string(deepks_model);
    Parallel_Common::bcast_int(bessel_lmax);
    Parallel_Common::bcast_double(bessel_rcut);
    Parallel_Common::bcast_double(bessel_tol);
    
    Parallel_Common::bcast_int(out_pot);
    Parallel_Common::bcast_int(out_wfc_pw);
    Parallel_Common::bcast_int(out_wfc_r);
    Parallel_Common::bcast_int(out_dos);
    Parallel_Common::bcast_int(out_band);
    Parallel_Common::bcast_int(out_proj_band);
    Parallel_Common::bcast_int(out_mat_hs);
    Parallel_Common::bcast_int(out_mat_hs2); // LiuXh add 2019-07-15
    Parallel_Common::bcast_int(out_mat_r); // jingan add 2019-8-14
    Parallel_Common::bcast_bool(out_wfc_lcao);
    Parallel_Common::bcast_bool(out_alllog);
    Parallel_Common::bcast_bool(out_element_info);

    Parallel_Common::bcast_double(dos_emin_ev);
    Parallel_Common::bcast_double(dos_emax_ev);
    Parallel_Common::bcast_double(dos_edelta_ev);
    Parallel_Common::bcast_double(dos_scale);
    Parallel_Common::bcast_bool(dos_setemin);
    Parallel_Common::bcast_bool(dos_setemax);
    Parallel_Common::bcast_int(dos_nche);
    Parallel_Common::bcast_double(b_coef);

    // mohan add 2009-11-11
    Parallel_Common::bcast_double(lcao_ecut);
    Parallel_Common::bcast_double(lcao_dk);
    Parallel_Common::bcast_double(lcao_dr);
    Parallel_Common::bcast_double(lcao_rmax);
    /*
        // mohan add 2011-11-07
        Parallel_Common::bcast_double( mdp.dt );
        Parallel_Common::bcast_int( md_restart );
        Parallel_Common::bcast_double( md_tolv );
        Parallel_Common::bcast_string( md_thermostat );
        Parallel_Common::bcast_double( md_temp0 );
        Parallel_Common::bcast_int( md_tstep );
        Parallel_Common::bcast_double( md_delt );
    */
    // zheng daye add 2014/5/5
    Parallel_Common::bcast_int(mdp.md_type);
    Parallel_Common::bcast_int(mdp.md_nstep);
    Parallel_Common::bcast_double(mdp.md_dt);
    Parallel_Common::bcast_int(mdp.md_mnhc);
    Parallel_Common::bcast_double(mdp.msst_qmass);
    Parallel_Common::bcast_double(mdp.md_tfirst);
    Parallel_Common::bcast_double(mdp.md_tlast);
    Parallel_Common::bcast_int(mdp.md_dumpfreq);
    Parallel_Common::bcast_int(mdp.md_restartfreq);
    Parallel_Common::bcast_int(mdp.md_seed);
    Parallel_Common::bcast_bool(mdp.md_restart);
    Parallel_Common::bcast_double(mdp.lj_rcut);
    Parallel_Common::bcast_double(mdp.lj_epsilon);
    Parallel_Common::bcast_double(mdp.lj_sigma);
    Parallel_Common::bcast_string(mdp.md_ensolver);
    Parallel_Common::bcast_int(mdp.msst_direction);
    Parallel_Common::bcast_double(mdp.msst_vel);
    Parallel_Common::bcast_double(mdp.msst_vis);
    Parallel_Common::bcast_double(mdp.msst_tscale);
    Parallel_Common::bcast_double(mdp.md_tfreq);
    Parallel_Common::bcast_double(mdp.md_damp);
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
    Parallel_Common::bcast_int(tddft);
    Parallel_Common::bcast_int(td_val_elec_01);
    Parallel_Common::bcast_int(td_val_elec_02);
    Parallel_Common::bcast_int(td_val_elec_03);
    Parallel_Common::bcast_double(td_scf_thr);
    Parallel_Common::bcast_double(td_dt);
    Parallel_Common::bcast_double(td_force_dt);
    Parallel_Common::bcast_int(td_vext);
    Parallel_Common::bcast_int(td_vext_dire);
    Parallel_Common::bcast_double(td_timescale);
    Parallel_Common::bcast_int(td_vexttype);
    Parallel_Common::bcast_int(td_vextout);
    Parallel_Common::bcast_int(td_dipoleout);
    Parallel_Common::bcast_bool(test_skip_ewald);
    Parallel_Common::bcast_int(GlobalV::ocp);
    Parallel_Common::bcast_string(GlobalV::ocp_set);
    Parallel_Common::bcast_int(out_mul); // qifeng add 2019/9/10

    // Peize Lin add 2018-06-20
    Parallel_Common::bcast_double(exx_hybrid_alpha);
    Parallel_Common::bcast_double(exx_hse_omega);
    Parallel_Common::bcast_bool(exx_separate_loop);
    Parallel_Common::bcast_int(exx_hybrid_step);
    Parallel_Common::bcast_double(exx_lambda);
    Parallel_Common::bcast_double(exx_pca_threshold);
    Parallel_Common::bcast_double(exx_c_threshold);
    Parallel_Common::bcast_double(exx_v_threshold);
    Parallel_Common::bcast_double(exx_dm_threshold);
    Parallel_Common::bcast_double(exx_schwarz_threshold);
    Parallel_Common::bcast_double(exx_cauchy_threshold);
    Parallel_Common::bcast_double(exx_ccp_threshold);
    Parallel_Common::bcast_double(exx_ccp_rmesh_times);
    Parallel_Common::bcast_string(exx_distribute_type);
    Parallel_Common::bcast_int(exx_opt_orb_lmax);
    Parallel_Common::bcast_double(exx_opt_orb_ecut);
    Parallel_Common::bcast_double(exx_opt_orb_tolerence);

    Parallel_Common::bcast_bool(noncolin);
    Parallel_Common::bcast_bool(lspinorb);
    Parallel_Common::bcast_double(soc_lambda);

    // Parallel_Common::bcast_int( epsilon0_choice );
    Parallel_Common::bcast_double(cell_factor); // LiuXh add 20180619
    Parallel_Common::bcast_bool(restart_save); // Peize Lin add 2020.04.04
    Parallel_Common::bcast_bool(restart_load); // Peize Lin add 2020.04.04

    //-----------------------------------------------------------------------------------
    // DFT+U (added by Quxin 2020-10-29)
    //-----------------------------------------------------------------------------------
    Parallel_Common::bcast_bool(dft_plus_u);
    Parallel_Common::bcast_bool(yukawa_potential);
    Parallel_Common::bcast_bool(omc);
    Parallel_Common::bcast_int(dftu_type);
    Parallel_Common::bcast_int(double_counting);
    Parallel_Common::bcast_double(yukawa_lambda);
    if (GlobalV::MY_RANK != 0)
    {
        hubbard_u = new double[this->ntype];
        hund_j = new double[this->ntype];
        orbital_corr = new int[this->ntype];
    }

    for (int i = 0; i < this->ntype; i++)
    {
        Parallel_Common::bcast_double(hubbard_u[i]);
        Parallel_Common::bcast_double(hund_j[i]);
        Parallel_Common::bcast_int(orbital_corr[i]);
    }

    //-----------------------------------------------------------------------------------
    // DFT+DMFT (added by Quxin 2020-08)
    //-----------------------------------------------------------------------------------
    Parallel_Common::bcast_bool(dft_plus_dmft);

    //----------------------------------------------------------------------------------
    //    implicit solvation model        (sunml added on 2022-04-04)
    //----------------------------------------------------------------------------------
    Parallel_Common::bcast_bool(imp_sol);
    Parallel_Common::bcast_double(eb_k);
    Parallel_Common::bcast_double(tau);
    Parallel_Common::bcast_double(sigma_k);
    Parallel_Common::bcast_double(nc_k);

    return;
}
#endif

void Input::Check(void)
{
    ModuleBase::TITLE("Input", "Check");

    if (nbands < 0)
        ModuleBase::WARNING_QUIT("Input", "NBANDS must >= 0");
    //	if(nbands_istate < 0) ModuleBase::WARNING_QUIT("Input","NBANDS_ISTATE must > 0");
    if (nb2d < 0)
        ModuleBase::WARNING_QUIT("Input", "nb2d must > 0");
    if (ntype <= 0)
        ModuleBase::WARNING_QUIT("Input", "ntype must > 0");

    // std::cout << "diago_proc=" << diago_proc << std::endl;
    // std::cout << " NPROC=" << GlobalV::NPROC << std::endl;
    if (diago_proc > 1 && basis_type == "lcao")
    {
        ModuleBase::WARNING_QUIT("Input", "please don't set diago_proc with lcao base");
    }
    if (diago_proc <= 0)
    {
        diago_proc = GlobalV::NPROC;
    }
    else if (diago_proc > GlobalV::NPROC)
    {
        diago_proc = GlobalV::NPROC;
    }

    if (kspacing < 0.0)
    {
        ModuleBase::WARNING_QUIT("Input", "kspacing must > 0");
    }

    if (nelec < 0.0)
    {
        ModuleBase::WARNING_QUIT("Input", "nelec < 0 is not allowed !");
    }

    if(dip_cor_flag && !efield_flag)
    {
        ModuleBase::WARNING_QUIT("Input", "dipole correction is not active if efield_flag=false !");
    }

    if(gate_flag && efield_flag && !dip_cor_flag)
    {
        ModuleBase::WARNING_QUIT("Input", "gate field cannot be used with efield if dip_cor_flag=false !");
    }

    //----------------------------------------------------------
    // main parameters / electrons / spin ( 1/16 )
    //----------------------------------------------------------
    if (calculation == "scf")
    {
        if (mem_saver == 1)
        {
            mem_saver = 0;
            ModuleBase::GlobalFunc::AUTO_SET("mem_savre", "0");
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
    else if (calculation == "sto-scf") // qianrui 2021-2-20
    {
        if (mem_saver == 1)
        {
            mem_saver = 0;
            ModuleBase::GlobalFunc::AUTO_SET("mem_savre", "0");
        }
        this->relax_nmax = 1;
    }
    else if (calculation == "relax") // pengfei 2014-10-13
    {
        if (mem_saver == 1)
        {
            mem_saver = 0;
            ModuleBase::GlobalFunc::AUTO_SET("mem_savre", "0");
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

        // if (local_basis == 0 && linear_scaling == 0) xiaohui modify 2013-09-01
        if (basis_type == "pw" && calculation == "get_S") // xiaohui add 2013-09-01. Attention! maybe there is some problem
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
        if (out_dos == 3 && symmetry)
        {
            ModuleBase::WARNING_QUIT("Input::Check",
                                     "symmetry can't be used for out_dos==3(Fermi Surface Plotting) by now.");
        }
    }
    else if (calculation == "istate")
    {
        GlobalV::CALCULATION = "istate";
        this->relax_nmax = 1;
        out_stru = 0;
        out_dos = 0;
        out_band = 0;
        out_proj_band = 0;
        cal_force = 0;
        init_wfc = "file";
        init_chg = "atomic"; // useless,
        chg_extrap = "atomic"; // xiaohui modify 2015-02-01
        out_chg = 1; // this leads to the calculation of state charge.
        out_dm = 0;
        out_pot = 0;

        // if(!local_basis || !linear_scaling) xiaohui modify 2013-09-01
        if (basis_type == "pw") // xiaohui add 2013-09-01
        {
            ModuleBase::WARNING_QUIT("Input::Check", "calculate = istate is only availble for LCAO.");
        }
    }
    else if (calculation == "ienvelope")
    {
        GlobalV::CALCULATION = "ienvelope"; // mohan fix 2011-11-04
        this->relax_nmax = 1;
        out_stru = 0;
        out_dos = 0;
        out_band = 0;
        out_proj_band = 0;
        cal_force = 0;
        init_wfc = "file";
        init_chg = "atomic";
        chg_extrap = "atomic"; // xiaohui modify 2015-02-01
        out_chg = 1;
        out_dm = 0;
        out_pot = 0;
        // if(!local_basis || !linear_scaling) xiaohui modify 2013-09-01
        if (basis_type == "pw") // xiaohui add 2013-09-01
        {
            ModuleBase::WARNING_QUIT("Input::Check", "calculate = istate is only availble for LCAO.");
        }
    }
    else if (calculation == "md" || calculation == "sto-md") // mohan add 2011-11-04
    {
        GlobalV::CALCULATION = calculation;
        symmetry = 0;
        cal_force = 1;
        if (mdp.md_nstep == 0)
        {
            GlobalV::ofs_running << "md_nstep should be set. Autoset md_nstep to 50!" << endl;
            mdp.md_nstep = 50;
        }
        if (!out_md_control)
            out_level = "m"; // zhengdy add 2019-04-07

        // deal with input parameters , 2019-04-30
        // if(basis_type == "pw" ) ModuleBase::WARNING_QUIT("Input::Check","calculate = MD is only availble for LCAO.");
        if (mdp.md_dt < 0)
            ModuleBase::WARNING_QUIT("Input::Check", "time interval of MD calculation should be set!");
        if (mdp.md_tfirst < 0 && tddft==0)
            ModuleBase::WARNING_QUIT("Input::Check", "temperature of MD calculation should be set!");
        if (mdp.md_tlast < 0.0)
            mdp.md_tlast = mdp.md_tfirst;

        if(mdp.md_tfreq == 0)
        {
            mdp.md_tfreq = 1.0/40.0/mdp.md_dt;
        }
        if(mdp.md_restart) 
        {
            init_vel = 1;
        }
        if(mdp.md_ensolver == "LJ" || mdp.md_ensolver == "DP" || mdp.md_type == 4)
        {
            cal_stress = 1;
        }
        if(mdp.md_type == 4)
        {
            if(mdp.msst_qmass <= 0)
            {
                ModuleBase::WARNING_QUIT("Input::Check", "msst_qmass must be greater than 0!");
            }
        }
        // if(mdp.md_tfirst!=mdp.md_tlast)
        // {
        //     std::ifstream file1;
        //     file1.open("ChangeTemp.dat");
        //     if(!file1)                      // Peize Lin fix bug 2016-08-06
        //    {
        //         std::ofstream file;
        //         file.open("ChangeTemp.dat");
        //         for(int ii=0;ii<30;ii++)
        //         {
        //             file<<mdp.md_tfirst+(mdp.md_tlast-mdp.md_tfirst)/double(30)*double(ii+1)<<" ";
        //         }
        //         file.close();
        //     }
        //     else
        //         file1.close();
        // }
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
    else if(calculation == "test_neighbour")
    {
        this->relax_nmax = 1;
    }
    else if(calculation == "gen_bessel")
    {
        this->relax_nmax = 1;
        if(basis_type != "pw")
        {
            ModuleBase::WARNING_QUIT("Input","to generate descriptors, please use pw basis");
        }
    }
    else
    {
        ModuleBase::WARNING_QUIT("Input", "check 'calculation' !");
    }
    if (init_chg != "atomic" && init_chg != "file")
    {
        ModuleBase::WARNING_QUIT("Input", "wrong 'init_chg',not 'atomic', 'file',please check");
    }
    // xiaohui modify 2014-05-10, chg_extrap value changes to 0~7
    // if (chg_extrap <0 ||chg_extrap > 7)
    //{
    //	ModuleBase::WARNING_QUIT("Input","wrong 'chg_extrap',neither 0~7.");
    // }xiaohui modify 2015-02-01
    if (gamma_only_local == 0)
    {
        if (out_dm == 1)
        {
            ModuleBase::WARNING_QUIT("Input", "out_dm with k-point algorithm is not implemented yet.");
        }
    }

    // if(chg_extrap==4 && local_basis==0) xiaohui modify 2013-09-01
    if (chg_extrap == "dm" && basis_type == "pw") // xiaohui add 2013-09-01, xiaohui modify 2015-02-01
    {
        ModuleBase::WARNING_QUIT(
            "Input",
            "wrong 'chg_extrap=dm' is only available for local orbitals."); // xiaohui modify 2015-02-01
    }

    if (chg_extrap == "dm" || cal_force > 1)
    {
        // if(out_dm==0) out_dm = 10000;//at least must output the density matrix at the last electron iteration step.
    }
    // if(chg_extrap != "dm")//xiaohui add 2015-02-01
    //{
    //	if(calculation=="relax")//xiaohui add 2015-02-01
    //	{
    //		chg_extrap = "first-order";
    //	}
    //	if(calculation=="md")//xiaohui add 2015-02-01
    //	{
    //		chg_extrap = "second-order";
    //	}
    // }

    if (GlobalV::CALCULATION == "nscf" && init_chg != "file")
    {
        init_chg = "file";
        ModuleBase::GlobalFunc::AUTO_SET("init_chg", init_chg);
    }

    if (init_wfc != "atomic" && init_wfc != "atomic+random" && init_wfc != "random" && init_wfc != "file")
    {
        ModuleBase::WARNING_QUIT("Input", "wrong init_wfc, please use 'atomic' or 'random' or 'file' ");
    }

    if (nbands > 100000)
    {
        ModuleBase::WARNING_QUIT("Input", "nbnd >100000, out of range");
    }
    if (nelec > 0 && nbands > 0 && nelec > 2 * nbands)
    {
        ModuleBase::WARNING_QUIT("Input", "nelec > 2*nbnd , bands not enough!");
    }
    if (nspin < 1 || nspin > 4)
    {
        ModuleBase::WARNING_QUIT("Input", "nspin out of range!");
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
        }
        else if (ks_solver == "dav")
        {
            GlobalV::ofs_warning << " It's ok to use dav." << std::endl;
        }
        else if (ks_solver == "genelpa") // yshen add 2016-07-20
        {
            ModuleBase::WARNING_QUIT("Input", "genelpa can not be used with plane wave basis.");
        }
        else if (ks_solver == "scalapack_gvx") // Peize Lin add 2020.11.14
        {
            ModuleBase::WARNING_QUIT("Input", "scalapack_gvx can not be used with plane wave basis.");
        }
        else if (ks_solver == "hpseps")
        {
            ModuleBase::WARNING_QUIT("Input", "hpseps can not be used with plane wave basis."); // xiaohui add
                                                                                                // 2013-09-04
        }
        else if (ks_solver == "selinv")
        {
            ModuleBase::WARNING_QUIT("Input", "selinv can not be used with plane wave basis."); // xiaohui add
                                                                                                // 2013-09-04
        }
        else if (ks_solver == "lapack")
        {
            ModuleBase::WARNING_QUIT("Input", "lapack can not be used with plane wave basis.");
        }
        else
        {
            ModuleBase::WARNING_QUIT("Input", "please check the ks_solver parameter!");
        }
    }
    else if (basis_type == "lcao")
    {
        if (ks_solver == "default")
        {
            ks_solver = "genelpa";
            ModuleBase::GlobalFunc::AUTO_SET("ks_solver", "genelpa");
        }
        else if (ks_solver == "cg")
        {
            ModuleBase::WARNING_QUIT("Input", "not ready for cg method in lcao ."); // xiaohui add 2013-09-04
        }
        else if (ks_solver == "genelpa")
        {
#ifdef __MPI
//				GlobalV::ofs_warning << "genelpa is under testing" << std::endl;
#else
            ModuleBase::WARNING_QUIT("Input", "genelpa can not be used for series version.");
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
        else if (ks_solver == "hpseps")
        {
#ifdef __MPI
            GlobalV::ofs_warning << "It's not a good choice to use hpseps!" << std::endl;
            if (gamma_only)
                ModuleBase::WARNING_QUIT("Input", "hpseps can not be used for gamma_only.");
#else
            ModuleBase::WARNING_QUIT("Input", "hpseps can not be used for series version.");
#endif
        }
        else if (ks_solver == "lapack")
        {
#ifdef __MPI
            ModuleBase::WARNING_QUIT("Input",
                                     "ks_solver=lapack is not an option for parallel version of ABACUS (try hpseps).");
#else
            GlobalV::ofs_warning << " It's ok to use lapack." << std::endl;
#endif
        }
        else if (ks_solver == "selinv")
        {
            ModuleBase::WARNING_QUIT("Input", "not ready for selinv method in lcao .");
        }
        else if (ks_solver == "linear_scaling")
        {
            ModuleBase::WARNING_QUIT("Input", "not ready for linear_scaling method in lcao .");
        }
        else if (ks_solver == "cusolver")
        {
#ifndef __MPI
            ModuleBase::WARNING_QUIT("Input","Cusolver can not be used for series version.");
#endif
        }
        else
        {
            ModuleBase::WARNING_QUIT("Input", "please check the ks_solver parameter!");
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

    if (basis_type == "pw" && gamma_only)
    {
        ModuleBase::WARNING_QUIT("Input", "gamma_only not implemented for plane wave now.");
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

    if (basis_type == "lcao" && !gamma_only_local) // xiaohui add 2013-09-01. Attention! Maybe there is some problem.
    {
        ModuleBase::WARNING("Input", "gamma_only_local algorithm is not used.");
    }

    if (basis_type == "lcao" && kpar > 1)
    {
        ModuleBase::WARNING_QUIT("Input", "kpar > 1 has not been supported for lcao calculation.");
    }
    // new rule, mohan add 2012-02-11
    // otherwise, there need wave functions transfers
    // if(diago_type=="cg") xiaohui modify 2013-09-01
    if (ks_solver == "cg") // xiaohui add 2013-09-01
    {
        if (diago_proc != GlobalV::NPROC)
        {
            ModuleBase::WARNING("Input", "when CG is used for diago, diago_proc==GlobalV::NPROC");
            diago_proc = GlobalV::NPROC;
        }
    }

    if (GlobalV::NPROC > 1 && ks_solver == "lapack") // xiaohui add 2013-09-01
    {
        // if(local_basis ==4 && linear_scaling==0) xiaohui modify 2013-09-01
        if (basis_type == "lcao_in_pw") // xiaohui add 2013-09-01
        {
        }
        else
        {
            ModuleBase::WARNING_QUIT("Input", "lapack can not be used when nproc > 1");
        }
    }

    // pengfei add 13-8-10 a new method cg to bfgs
    if (relax_method != "sd" && relax_method != "cg" && relax_method != "bfgs" && relax_method != "cg_bfgs")
    {
        ModuleBase::WARNING_QUIT("Input", "relax_method can only be sd, cg, bfgs or cg_bfgs.");
    }

    if (basis_type == "pw")
    {
        bx = 1;
        by = 1;
        bz = 1;
    }
    else if (bx > 10)
    {
        ModuleBase::WARNING_QUIT("Input", "bx is too large!");
    }
    else if (by > 10)
    {
        ModuleBase::WARNING_QUIT("Input", "by is too large!");
    }
    else if (bz > 10)
    {
        ModuleBase::WARNING_QUIT("Input", "bz is too large!");
    }

    if (basis_type == "lcao")
    {
        if (lcao_ecut == 0)
        {
            lcao_ecut = ecutwfc;
            ModuleBase::GlobalFunc::AUTO_SET("lcao_ecut", ecutwfc);
        }
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

    if (dft_functional == "hf" || dft_functional == "pbe0" || dft_functional == "hse" || dft_functional == "scan0")
    {
        if (exx_hybrid_alpha < 0 || exx_hybrid_alpha > 1)
        {
            ModuleBase::WARNING_QUIT("INPUT", "must 0 < exx_hybrid_alpha < 1");
        }
        if (exx_hybrid_step <= 0)
        {
            ModuleBase::WARNING_QUIT("INPUT", "must exx_hybrid_step > 0");
        }
        if (exx_ccp_rmesh_times < 1)
        {
            ModuleBase::WARNING_QUIT("INPUT", "must exx_ccp_rmesh_times >= 1");
        }
        if (exx_distribute_type != "htime" && exx_distribute_type != "kmeans2" && exx_distribute_type != "kmeans1"
            && exx_distribute_type != "order")
        {
            ModuleBase::WARNING_QUIT("INPUT", "exx_distribute_type must be htime or kmeans2 or kmeans1");
        }
    }
    if (dft_functional == "opt_orb")
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

    // 2015-06-15, xiaohui
    if (mixing_mode == "pulay" && mixing_gg0 > 0.0)
    {
        ModuleBase::WARNING("Input", "To use pulay-kerker mixing method, please set mixing_type=pulay-kerker");
    }

    if (berry_phase)
    {
        if (basis_type == "pw")
        {
            if (!(calculation == "nscf"))
                ModuleBase::WARNING_QUIT("Input", "calculate berry phase, please set calculation = nscf");
        }
        else if (basis_type == "lcao" && ks_solver == "genelpa" || ks_solver == "scalapack_gvx" || ks_solver == "cusolver")
        {
            if (!(calculation == "nscf"))
                ModuleBase::WARNING_QUIT("Input", "calculate berry phase, please set calculation = nscf");
        }
        else
        {
            ModuleBase::WARNING_QUIT("Input", "calculate berry phase, please set basis_type = pw or lcao");
        }

        if (!(gdir == 1 || gdir == 2 || gdir == 3))
        {
            ModuleBase::WARNING_QUIT("Input", "calculate berry phase, please set gdir = 1 or 2 or 3");
        }
    }

    if (towannier90)
    {
        if (basis_type == "pw" || basis_type == "lcao")
        {
            if (!(calculation == "nscf"))
                ModuleBase::WARNING_QUIT("Input", "to use towannier90, please set calculation = nscf");
        }
        else
        {
            ModuleBase::WARNING_QUIT("Input", "to use towannier90, please set basis_type = pw or lcao");
        }

        if (nspin == 2)
        {
            if (!(wannier_spin == "up" || wannier_spin == "down"))
            {
                ModuleBase::WARNING_QUIT("Input", "to use towannier90, please set wannier_spin = up or down");
            }
        }
    }

    const std::string ss = "test -d " + read_file_dir;
    if (read_file_dir == "auto")
    {
        GlobalV::global_readin_dir = GlobalV::global_out_dir;
    }
    else if (system(ss.c_str()))
    {
        ModuleBase::WARNING_QUIT("Input", "please set right files directory for reading in.");
    }
    else
    {
        GlobalV::global_readin_dir = read_file_dir + '/';
    }

    return;
}

void Input::close_log(void) const
{

    ModuleBase::Global_File::close_all_log(GlobalV::MY_RANK, this->out_alllog);
}

void Input::readbool(std::ifstream &ifs, bool &var)
{
    std::string str;
    ifs >> str;
    if (str == "true")
    {
        var = true;
    }
    else
    {
        var = false;
    }
    ifs.ignore(100, '\n');
    return;
}

void Input::strtolower(char *sa, char *sb)
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
