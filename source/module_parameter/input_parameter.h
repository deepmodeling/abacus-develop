#ifndef INPUT_PARAMETER_H
#define INPUT_PARAMETER_H
#include "md_parameter.h"
#include "module_base/vector3.h"

#include <string>
#include <vector>

// It stores all input parameters both defined in INPUT file and not defined in
// INPUT file
struct Input_para
{
    // ---------------------------------------------------------------
    // --------------       INPUT  Parameters         ----------------
    // ---------------------------------------------------------------
    // ==============   #Parameters (1.System)  =====================
    std::string suffix = "ABACUS";      ///< suffix of out put dir
    int ntype = 0;                      ///< number of atom types
    std::string calculation = "scf";    ///< "scf" : self consistent calculation.
                                        ///< "nscf" : non-self consistent calculation.
                                        ///< "relax" : cell relaxations
    std::string esolver_type = "ksdft"; ///< the energy solver: ksdft, sdft, ofdft, tddft, lj, dp
    /* symmetry level:
      -1, no symmetry at all;
      0, only basic time reversal would be considered;
      1, point group symmetry would be considered*/
    std::string symmetry = "default";
    double symmetry_prec = 1.0e-6;  ///< LiuXh add 2021-08-12, accuracy for symmetry
    bool symmetry_autoclose = true; ///< whether to close symmetry automatically
                                    ///< when error occurs in symmetry analysis
    bool cal_force = false;         ///< calculate the force
    bool cal_stress = false;        ///< calculate the stress
    int kpar = 1;                   ///< ecch pool is for one k point
    int bndpar = 1;                 ///< parallel for stochastic/deterministic bands
    std::string latname = "none";   ///< lattice name
    double ecutwfc = 0;            ///< energy cutoff for wavefunctions
    double ecutrho = 0;             ///< energy cutoff for charge/potential

    int nx = 0, ny = 0, nz = 0;    ///< three dimension of FFT wavefunc
    int ndx = 0, ndy = 0, ndz = 0; ///< three dimension of FFT smooth charge density

    double cell_factor = 1.2;           ///< LiuXh add 20180619
    double erf_ecut = 0;                ///< the value of the constant energy cutoff
    double erf_height = 0;              ///< the height of the energy step for reciprocal vectors
    double erf_sigma = 0.1;             ///< the width of the energy step for reciprocal vectors
    int fft_mode = 0;                   ///< fftw mode 0: estimate, 1: measure, 2: patient, 3: exhaustive
    std::string init_wfc = "atomic";    ///< "file","atomic","random"
    bool psi_initializer = false;       ///< whether use psi_initializer to initialize wavefunctions
    int pw_seed = 0;                    ///< random seed for initializing wave functions
    std::string init_chg = "atomic";    ///< "file","atomic"
    bool dm_to_rho = false;             ///< read density matrix from npz format and calculate charge density
    std::string chg_extrap = "default"; ///< xiaohui modify 2015-02-01
    bool init_vel = false;              ///< read velocity from STRU or not  liuyu 2021-07-14
    
    std::string input_file = "INPUT";   ///< input file name
    std::string stru_file = "STRU";     ///< file contains atomic positions --
                                        ///< xiaohui modify 2015-02-01
    std::string kpoint_file = "KPT";    ///< file contains k-points -- xiaohui modify 2015-02-01
    std::string pseudo_dir = "";        ///< directory of pseudopotential
    std::string orbital_dir = "";       ///< directory of orbital file
    std::string read_file_dir = "auto"; ///< directory of files for reading
    bool restart_load = false;
    std::string wannier_card = "none";              ///< input card for wannier functions.
    int mem_saver = 0;                              ///< 1: save psi when nscf calculation.
    int diago_proc = 0;                             ///< the number of procs used to diag. mohan add 2012-01-13
    int nbspline = -1;                              ///< the order of B-spline basis(>=0) if it is -1 (default)
    std::vector<double> kspacing = {0.0, 0.0, 0.0}; ///< kspacing for k-point generation
    double min_dist_coef = 0.2;                     ///< allowed minimum distance between two atoms

    std::string device = "auto";
    std::string precision = "double";

    // ==============   #Parameters (2.Electronic structure) ===========================
    std::string ks_solver = "default"; ///< xiaohui add 2013-09-01
    std::string basis_type = "pw";     ///< xiaohui add 2013-09-01, for structural adjustment
    bool use_paw = false;              ///< whether to use PAW in pw calculation
    int nbands = 0;                    ///< number of bands
    double nelec = 0.0;                ///< total number of electrons
    double nelec_delta = 0.0;          ///< change in the number of total electrons
    double nupdown = 0.0;
    std::string dft_functional = "default"; ///< input DFT functional.
    double xc_temperature = 0.0;            ///< only relevant if finite temperature functional is used
    double pseudo_rcut = 15.0;              ///< cut-off radius for calculating msh
    bool pseudo_mesh = false;               ///< 0: use msh to normalize radial wave functions; 1:
                                            ///< use mesh, which is used in QE.
    int nspin = 1;                          ///< LDA ; LSDA ; non-linear spin
    int pw_diag_nmax = 50;
    double pw_diag_thr = 0.01; ///< used in cg method
    int pw_diag_ndim = 4;      ///< dimension of workspace for Davidson diagonalization
    int diago_cg_prec = 1;     ///< mohan add 2012-03-31

    std::string smearing_method = "gauss"; ///< "gauss",
                                           ///< "mp","methfessel-paxton"
                                           ///< "mv","marzari-vanderbilt","cold"
                                           ///< "fd","fermi-dirac"
    double smearing_sigma = 0.015;         ///<

    std::string mixing_mode = "broyden"; ///< "plain","broyden",...
    double mixing_beta = -1.0;           ///< 0 : no_mixing
    int mixing_ndim = 8;                 ///< used in Broyden method
    double mixing_restart = 0.0;         ///< mixing will restart once if drho is
                                         ///< smaller than mixing_restart
    double mixing_gg0 = 1.0;             ///< used in kerker method
    double mixing_beta_mag = -10.0;
    double mixing_gg0_mag = 0.0;
    double mixing_gg0_min = 0.1;
    double mixing_angle = -10.0;
    bool mixing_tau = false;  ///< whether to mix tau in mgga
    bool mixing_dftu = false; ///< whether to mix locale in DFT+U
    bool mixing_dmr = false;  ///< whether to mix real space density matrix

    bool gamma_only = false;   ///< for plane wave.
    int scf_nmax = 100;        ///< number of max elec iter
    double scf_thr = -1.0;     ///< \sum |rhog_out - rhog_in |^2
    double scf_ene_thr = -1.0; ///< energy threshold for scf convergence, in eV
    int scf_thr_type = -1;     ///< type of the criterion of scf_thr, 1: reci drho, 2: real drho
    bool final_scf = false;    ///< whether to do final scf
    bool scf_os_stop = false;  ///< whether to stop scf when oscillation is detected
    double scf_os_thr = -0.01;  ///< drho threshold for oscillation
    int scf_os_ndim = 0;       ///< number of old iterations used for oscillation detection

    bool lspinorb = false;   ///< consider the spin-orbit interaction
    bool noncolin = false;   ///< using non-collinear-spin
    double soc_lambda = 1.0; ///< The fraction of averaged SOC pseudopotential
                             ///< is given by (1-soc_lambda)

    // ==============   #Parameters (3.LCAO) ===========================
    int nb2d = 0;                              ///< matrix 2d division.
    int lmaxmax = 2;                           ///< maximum of l channels used
    double lcao_ecut = 0.0;                    ///< ecut of two center integral
    double lcao_dk = 0.01;                     ///< delta k used in two center integral
    double lcao_dr = 0.01;                     ///< dr used in two center integral
    double lcao_rmax = 30.0;                   ///< rmax(a.u.) to make table.
    double search_radius = -1.0;               ///< 11.1
    bool search_pbc = true;                    ///< 11.2
    int bx = 0, by = 0, bz = 0;                ///< big mesh ball. 0: auto set bx/by/bz
    int elpa_num_thread = -1;                  ///< Number of threads need to use in elpa
    int nstream = 4;                           ///< Number of streams in CUDA as per input data
    std::string bessel_nao_ecut = "default";   ///< energy cutoff for spherical bessel functions(Ry)
    double bessel_nao_tolerence = 1e-12;       ///< tolerance for spherical bessel root
    std::vector<double> bessel_nao_rcuts = {}; ///< No specific values provided for bessel_nao_rcuts
    bool bessel_nao_smooth = true;             ///< spherical bessel smooth or not
    double bessel_nao_sigma = 0.1;             ///< spherical bessel smearing_sigma
    // ==========================================================
    //  spherical bessel  Peize Lin added on 2022-12-15
    // ==========================================================
    //  the following are used when generating orb_matrix.dat
    //  int		bessel_nao_lmax;		///< lmax used in descriptor

    // ==============   #Parameters (4.Relaxation) ===========================
    std::string relax_method = "cg"; ///< methods to move_ion: sd, bfgs, cg...
    bool relax_new = true;
    bool relax = false; ///< allow relaxation along the specific direction
    double relax_scale_force = 0.5;
    int relax_nmax = 0;        ///< number of max ionic iter
    double relax_cg_thr = 0.5; ///< threshold when cg to bfgs, pengfei add 2011-08-15
    double force_thr = -1;     ///< threshold of force in unit (Ry/Bohr)
    double force_thr_ev = -1;  ///< threshold of force in unit (eV/Angstrom)
    double force_thr_ev2 = 0;  ///< invalid force threshold, mohan add 2011-04-17
    double stress_thr = 0.5;   ///< Pengfei Li 2017-11-01 ///<LiuXh update 20180515
    double press1 = 0;
    double press2 = 0;
    double press3 = 0;
    double relax_bfgs_w1 = 0.01;     ///< wolfe condition 1
    double relax_bfgs_w2 = 0.5;      ///< wolfe condition 2
    double relax_bfgs_rmax = 0.8;    ///< trust radius max
    double relax_bfgs_rmin = 1e-05;  ///< trust radius min
    double relax_bfgs_init = 0.5;    ///< initial move
    std::string fixed_axes = "None"; ///< which axes are fixed
    bool fixed_ibrav = false;        ///< whether to keep type of lattice; must be used
                                     ///< along with latname
    bool fixed_atoms = false;        ///< whether to fix atoms during vc-relax

    // ==============   #Parameters (5.Molecular dynamics) ===========================
    MD_para mdp;
    double ref_cell_factor = 1; ///< construct a reference cell bigger than the
                                ///< initial cell liuyu 2023-03-21
    bool cal_syns = false;      ///< calculate asynchronous S matrix to output
    double dmax = 0.01;         ///< maximum displacement of all atoms in one step (bohr)

    // ==============   #Parameters (6.OFDFT) ===========================
    // OFDFT  sunliang added on 2022-05-05
    std::string of_kinetic = "wt";               ///< Kinetic energy functional, such as TF, VW, WT, TF+
    std::string of_method = "tn";                ///< optimization method, include cg1, cg2, tn (default), bfgs
    std::string of_conv = "energy";              ///< select the convergence criterion,
                                                 ///< potential, energy (default), or both
    double of_tole = 1e-06;                      ///< tolerance of the energy change (in Ry) for
                                                 ///< determining the convergence, default=2e-6 Ry
    double of_tolp = 1e-05;                      ///< tolerance of potential for determining the
                                                 ///< convergence, default=1e-5 in a.u.
    double of_tf_weight = 1.0;                   ///< weight of TF KEDF
    double of_vw_weight = 1.0;                   ///< weight of vW KEDF
    double of_wt_alpha = 5. / 6.;                ///< parameter alpha of WT KEDF
    double of_wt_beta = 5. / 6.;                 ///< parameter beta of WT KEDF
    double of_wt_rho0 = 0.0;                     ///< set the average density of system, in Bohr^-3
    bool of_hold_rho0 = false;                   ///< If set to 1, the rho0 will be fixed even if the volume of
                                                 ///< system has changed, it will be set to 1 automatically if
                                                 ///< of_wt_rho0 is not zero.
    double of_lkt_a = 1.3;                       ///< parameter a of LKT KEDF
    bool of_full_pw = true;                      ///< If set to 1, ecut will be ignored while collecting
                                                 ///< planewaves, so that all planewaves will be used.
    int of_full_pw_dim = 0;                      ///< If of_full_pw = 1, the dimension of FFT will be restricted to
                                                 ///< be (0) either odd or even; (1) odd only; (2) even only.
    bool of_read_kernel = false;                 ///< If set to 1, the kernel of WT KEDF will be
                                                 ///< filled from file of_kernel_file, not from
                                                 ///< formula. Only usable for WT KEDF.
    std::string of_kernel_file = "WTkernel.txt"; ///< The name of WT kernel file.

    // ==============   #Parameters (7.stochastic DFT) ===========================
    int method_sto = 2;        ///< different methods for sdft, 1: slow, less memory 2:
                               ///< fast, more memory
    int npart_sto = 1;         ///< for method_sto = 2, reduce memory
    int nbands_sto = 256;      ///< number of stochastic bands //qianrui 2021-2-5
    int nche_sto = 100;        ///< number of orders for Chebyshev expansion in
                               ///< stochastic DFT ///<qinarui 2021-2-5
    double emin_sto = 0;       ///< Emin to normalize H
    double emax_sto = 0;       ///< Emax to normalize H
    int seed_sto = 0;          ///< random seed for sDFT
    double initsto_ecut = 0.0; ///< maximum ecut to init stochastic bands
    int initsto_freq = 0;      ///< frequency to init stochastic orbitals when running md

    // ==============   #Parameters (8.DeepKS) ===========================
    //==========================================================
    // DeepKS -- added by caoyu and mohan
    //==========================================================
    bool deepks_out_labels = false;   ///< (need libnpy) prints energy and force labels and
                                      ///< descriptors for training, wenfei 2022-1-12
    bool deepks_scf = false;          ///< (need libnpy and libtorch) if set to true, a trained model
                                      ///< would be needed to calculate V_delta and F_delta
    bool deepks_bandgap = false;      ///< for bandgap label. QO added 2021-12-15
    int deepks_v_delta = 0;           ///< for v_delta label. xuan added
    bool deepks_equiv = false;        ///< whether to use equivariant version of DeePKS
    bool deepks_out_unittest = false; ///< if set to true, prints intermediate quantities that shall
                                      ///< be used for making unit test
    std::string deepks_model = "None";              ///< needed when deepks_scf=1
    
    int bessel_descriptor_lmax = 2;                 ///< lmax used in descriptor
    std::string bessel_descriptor_ecut = "default"; ///< energy cutoff for spherical bessel functions(Ry)
    double bessel_descriptor_tolerence = 1e-12;     ///< tolerance for spherical bessel root
    double bessel_descriptor_rcut = 6.0;            ///< radial cutoff for spherical bessel functions(a.u.)
    bool bessel_descriptor_smooth = true;           ///< spherical bessel smooth or not
    double bessel_descriptor_sigma = 0.1;           ///< spherical bessel smearing_sigma

    // ==============   #Parameters (9.rt-tddft) ===========================
    double td_force_dt = 0.02;      ///<"fs"
    bool td_vext = false;           ///< add extern potential or not
    std::string td_vext_dire = "1"; ///< vext direction

    bool init_vecpot_file = false; ///< initialize the vector potential, though file or integral
    double td_print_eij = -1.0;    ///< threshold to output Eij elements
    int td_edm = 0;                ///< 0: new edm method   1: old edm method
    int propagator = 0;            ///< method of propagator
    int td_stype = 0;              ///< type of space domain  0 : length gauge  1: velocity gauge
    std::string td_ttype = "0";    ///< type of time domain
    ///<  0  Gauss type function.
    ///<  1  trapezoid type function.
    ///<  2  Trigonometric functions, sin^2.
    ///<  3  heaviside function.
    ///<  4  HHG function.
    int td_tstart = 1;
    int td_tend = 1000;

    ///< space domain parameters
    ///< length gauge
    double td_lcut1 = 0.05;
    double td_lcut2 = 0.95;

    ///< time domain parameters
    ///< Gauss
    std::string td_gauss_freq = "22.13"; ///< time(fs)^-1
    std::string td_gauss_phase = "0.0";
    std::string td_gauss_sigma = "30.0"; ///< time(fs)
    std::string td_gauss_t0 = "100.0";
    std::string td_gauss_amp = "0.25"; ///< V/A

    ///< trapezoid
    std::string td_trape_freq = "1.60"; ///< time(fs)^-1
    // Trapezoidal
    std::string td_trape_phase = "0.0";
    std::string td_trape_t1 = "1875.0";
    std::string td_trape_t2 = "5625.0";
    std::string td_trape_t3 = "7500.0";
    std::string td_trape_amp = "2.74"; // V/A

    // Trigonometric
    std::string td_trigo_freq1 = "1.164656"; // time(fs)^-1
    std::string td_trigo_freq2 = "0.029116"; // time(fs)^-1
    std::string td_trigo_phase1 = "0.0";
    std::string td_trigo_phase2 = "0.0";
    std::string td_trigo_amp = "2.74"; // V/A

    // Heaviside
    std::string td_heavi_t0 = "100.0";
    std::string td_heavi_amp = "1.0"; // V/A

    bool ocp = false;
    // std::string ocp_set = "";
    std::vector<double> ocp_kb = {}; ///< OCP kb values

    // ==============   #Parameters (10.lr-tddft) ===========================
    int lr_nstates = 1; ///< the number of 2-particle states to be solved
    std::vector<std::string> lr_init_xc_kernel = {};    ///< The method to initalize the xc kernel
    int nocc = -1;      ///< the number of occupied orbitals to form the 2-particle basis
    int nvirt = 1;      ///< the number of virtual orbitals to form the 2-particle basis (nocc + nvirt <= nbands)
    std::string xc_kernel = "LDA"; ///< exchange correlation (XC) kernel for LR-TDDFT
    std::string lr_solver = "dav"; ///< the eigensolver for LR-TDDFT
    double lr_thr = 1e-2;          ///< convergence threshold of the LR-TDDFT eigensolver
    bool out_wfc_lr = false; ///< whether to output the eigenvectors (excitation amplitudes) in the particle-hole basis
    bool lr_unrestricted = false; ///< whether to use the unrestricted construction for LR-TDDFT
    std::vector<double> abs_wavelen_range = {}; ///< the range of wavelength(nm) to output the absorption spectrum
    double abs_broadening = 0.01;                     ///< the broadening (eta) for LR-TDDFT absorption spectrum
    std::string ri_hartree_benchmark = "none"; ///< whether to use the RI approximation for the Hartree potential in LR-TDDFT for benchmark (with FHI-aims/ABACUS read-in style)
    std::vector<int> aims_nbasis = {};  ///< the number of basis functions for each atom type used in FHI-aims (for benchmark)
    // ==============   #Parameters (11.Output) ===========================
    bool out_stru = false;                ///< outut stru file each ion step
    int out_freq_elec = 0;                ///< the frequency ( >= 0) of electronic iter to output charge
                                          ///< 0: output only when converged
    int out_freq_ion = 0;                 ///< the frequency ( >= 0 ) of ionic step to output charge density;
                                          ///< 0: output only when ion steps are finished
    std::vector<int> out_chg = {0, 3};    ///< output charge density. 0: no; 1: yes
    int out_pot = 0;                      ///< yes or no
    int out_wfc_pw = 0;                   ///< 0: no; 1: txt; 2: dat
    bool out_wfc_r = false;               ///< 0: no; 1: yes
    int printe = 0;                       ///< Print out energy for each band for every printe step, default is scf_nmax
    std::vector<int> out_band = {0, 8};   ///< band calculation pengfei 2014-10-13
    int out_dos = 0;                      ///< dos calculation. mohan add 20090909
    bool out_mul = false;                 ///< qifeng add 2019-9-10
    bool out_proj_band = false;           ///< projected band structure calculation jiyy add 2022-05-11
    std::string out_level = "ie";         ///< control the output information.
    bool out_dm = false;                  ///< output density matrix. (gamma point)
    bool out_dm1 = false;                 ///< output density matrix. (multi-k points)
    bool out_bandgap = false;             ///< QO added for bandgap printing
    std::vector<int> out_mat_hs = {0, 8}; ///< output H matrix and S matrix in local basis.
    std::vector<int> out_mat_tk = {0, 8}; ///< output T(k) matrix in local basis.
    bool out_mat_hs2 = false;             ///< LiuXh add 2019-07-16, output H(R) matrix and
                                          ///< S(R) matrix in local basis.
    bool out_mat_dh = false;
    bool out_mat_xc = false;      ///< output exchange-correlation matrix in
                                  ///< KS-orbital representation.
    bool out_eband_terms = false; ///< output the band energy terms separately
    bool out_hr_npz = false;      ///< output exchange-correlation matrix in
                                  ///< KS-orbital representation.
    bool out_dm_npz = false;

    int out_interval = 1;
    bool out_app_flag = true; ///< whether output r(R), H(R), S(R), T(R), and dH(R) matrices
                              ///< in an append manner during MD liuyu 2023-03-20
    int out_ndigits = 8;      ///< Assuming 8 digits precision is needed for matrices output
    bool out_mat_t = false;
    bool out_element_info = false;        ///< output information of all elements
    bool out_mat_r = false;               ///< jingan add 2019-8-14, output r(R) matrix.
    int out_wfc_lcao = 0;                 ///< output the wave functions in local basis.
    bool out_dipole = false;              ///< output the dipole or not
    bool out_efield = false;              ///< output the efield or not
    bool out_current = false;             ///< output the current or not
    bool out_current_k = false;           ///< output tddft current for all k points
    bool out_vecpot = false;              ///< output the vector potential or not
    bool restart_save = false;            ///< restart //Peize Lin add 2020-04-04
    bool rpa = false;                     ///< rpa calculation
    int nbands_istate = 5;                ///< number of bands around fermi level for get_pchg calculation.
    std::vector<int> bands_to_print = {}; ///< specify the bands to be calculated for partial charge
    std::vector<int> out_pchg = {};       ///< specify the bands to be calculated for partial charge
    std::vector<int> out_wfc_norm = {};   ///< specify the bands to be calculated for norm of wfc
    std::vector<int> out_wfc_re_im = {};  ///< specify the bands to be calculated for real and imaginary parts of wfc
    bool if_separate_k = false; ///< whether to write partial charge for all k-points to individual files or merge them
    std::vector<int> out_elf = {0, 3};    ///< output the electron localization function (ELF). 0: no; 1: yes

    // ==============   #Parameters (12.Postprocess) ===========================
    double dos_emin_ev = -15.0;
    double dos_emax_ev = 15.0;
    double dos_edelta_ev = 0.01;
    double dos_scale = 0.01;
    double dos_sigma = 0.07; ///< pengfei 2014-10-13
    int dos_nche = 100;      ///< orders of Chebyshev expansions for dos

    bool cal_cond = false;      ///< calculate electronic conductivities
    double cond_che_thr = 1e-8; ///< control the error of Chebyshev expansions
                                ///< for conductivities
    double cond_dw = 0.1;       ///< d\omega for conductivities
    double cond_wcut = 10;      ///< cutoff \omega for conductivities
    double cond_dt = 0.02;      ///< dt to integrate conductivities
    int cond_dtbatch = 0;       ///< exp(iH*dt*cond_dtbatch) is expanded with Chebyshev expansion.
    int cond_smear = 1;         ///< smearing method for conductivities 1: Gaussian 2: Lorentzian
    double cond_fwhm = 0.4;     ///< FWHM for conductivities
    bool cond_nonlocal = true;  ///< if calculate nonlocal effects

    bool berry_phase = false; ///< berry phase calculation: calculate berry phase or not
    int gdir = 3;             ///< berry phase calculation: calculate the polarization in
                              ///< the direction of the lattice vector

    ///<==========================================================
    ///< Wannier functions
    ///<==========================================================
    bool towannier90 = false;               ///< add by jingan for wannier90: use wannier90
                                            ///< code interface or not
    std::string nnkpfile = "seedname.nnkp"; ///< add by jingan for wannier90: the wannier90 code
                                            ///< nnkp file name
    std::string wannier_spin = "up";        ///< add by jingan for wannier90: calculate
                                            ///< spin in wannier90 code interface
    int wannier_method = 1;                 ///< different implementation methods under Lcao basis set
    bool out_wannier_mmn = true;            ///< add by renxi for wannier90: output .mmn file or not
    bool out_wannier_amn = true;            ///< output .amn file or not
    bool out_wannier_unk = false;           ///< output UNK. file or not
    bool out_wannier_eig = true;            ///< output .eig file or not
    bool out_wannier_wvfn_formatted = true; ///< output UNK. file in text format or in binary format

    // ==============   #Parameters (13.Model) ===========================
    // ==========================================================
    //  efield and dipole correction
    //  Yu Liu add 2022-05-18
    // ==========================================================
    bool efield_flag = false;   ///< add electric field
    bool dip_cor_flag = false;  ///< dipole correction
    int efield_dir = 2;         ///< the direction of the electric field or dipole correction
    double efield_pos_max = -1; ///< position of the maximum of the saw-like
                                ///< potential along crystal axis efield_dir
    double efield_pos_dec = -1; ///< zone in the unit cell where the saw-like potential decreases
    double efield_amp = 0;      ///< amplitude of the electric field

    // ==========================================================
    //  gatefield (compensating charge)
    //  Yu Liu add 2022-09-13
    // ==========================================================
    bool gate_flag = false;    ///< compensating charge or not
    double zgate = 0.5;        ///< position of charged plate
    bool block = false;        ///< add a block potential or not
    double block_down = 0.45;  ///< low bound of the block
    double block_up = 0.55;    ///< high bound of the block
    double block_height = 0.1; ///< height of the block

    //    implicit solvation model       Menglin Sun added on 2022-04-04
    bool imp_sol = false;    ///< true: implicit solvation correction; false:
                             ///< vacuum calculation(default)
    double eb_k = 80;        ///< the relative permittivity of the bulk solvent
    double tau = 1.0798e-05; ///< the effective surface tension parameter
    double sigma_k = 0.6;    ///< the width of the diffuse cavity
    double nc_k = 0.00037;   ///< the cut-off charge density

    // ==============  #Parameters (14.vdW Correction) ===========================
    // ==========================================================
    //  vdw
    //  Peize Lin add 2014-03-31, jiyy update 2019-08-01
    // ==========================================================
    std::string vdw_method = "none";                        ///< the method of calculating vdw (none; d2; d3_0; d3_bj)
    std::string vdw_s6 = "default";                         ///< scale parameter of d2/d3_0/d3_bj
    std::string vdw_s8 = "default";                         ///< scale parameter of d3_0/d3_bj
    std::string vdw_a1 = "default";                         ///< damping parameter of d3_0/d3_bj
    std::string vdw_a2 = "default";                         ///< damping parameter of d3_bj
    double vdw_d = 20.0;                                    ///< damping parameter of d2
    bool vdw_abc = false;                                   ///< third-order term?
    std::string vdw_C6_file = "default";                    ///< filename of C6
    std::string vdw_C6_unit = "Jnm6/mol";                   ///< unit of C6, Jnm6/mol or eVA6
    std::string vdw_R0_file = "default";                    ///< filename of R0
    std::string vdw_R0_unit = "A";                          ///< unit of R0, A or Bohr
    std::string vdw_cutoff_type = "radius";                 ///< expression model of periodic
                                                            ///< structure, radius or period
    std::string vdw_cutoff_radius = "default";              ///< radius cutoff for periodic structure
    std::string vdw_radius_unit = "Bohr";                   ///< unit of radius cutoff for periodic structure
    double vdw_cn_thr = 40.0;                               ///< radius cutoff for cn
    std::string vdw_cn_thr_unit = "Bohr";                   ///< unit of cn_thr, Bohr or Angstrom
    ModuleBase::Vector3<int> vdw_cutoff_period = {3, 3, 3}; ///< periods of periodic structure

    // ==============   #Parameters (15.exx) ====================
    // ==========================================================
    //  exx
    //  Peize Lin add 2018-06-20
    // ==========================================================
    std::string exx_hybrid_alpha = "default";   ///< fraction of Fock exchange in hybrid functionals
    double exx_hse_omega = 0.11;                ///< range-separation parameter in HSE functional
    bool exx_separate_loop = true;              ///< if 1, a two-step method is employed, else it will start
                                                ///< with a GGA-Loop, and then Hybrid-Loop
    int exx_hybrid_step = 100;                  ///< the maximal electronic iteration number in
                                                ///< the evaluation of Fock exchange
    double exx_mixing_beta = 1.0;               ///< mixing_beta for outer-loop when exx_separate_loop=1
    double exx_lambda = 0.3;                    ///< used to compensate for divergence points at G=0 in the
                                                ///< evaluation of Fock exchange using lcao_in_pw method
    std::string exx_real_number = "default";          ///< exx calculated in real or complex
    double exx_pca_threshold = 0.0001;          ///< threshold to screen on-site ABFs in exx
    double exx_c_threshold = 0.0001;            ///< threshold to screen C matrix in exx
    double exx_v_threshold = 0.1;               ///< threshold to screen C matrix in exx
    double exx_dm_threshold = 0.0001;           ///< threshold to screen density matrix in exx
    double exx_schwarz_threshold = 0;           ///< threshold to screen exx using Cauchy-Schwartz inequality
    double exx_cauchy_threshold = 1e-07;        ///< threshold to screen exx using Cauchy-Schwartz inequality
    double exx_c_grad_threshold = 0.0001;       ///< threshold to screen nabla C matrix in exx
    double exx_v_grad_threshold = 0.1;          ///< threshold to screen nabla V matrix in exx
    double exx_c_grad_r_threshold = 0.0001;     ///< threshold to screen nabla C matrix in exx
    double exx_v_grad_r_threshold = 0.1;        ///< threshold to screen nabla V matrix in exx
    double exx_cauchy_force_threshold = 1e-07;  ///< threshold to screen exx force using Cauchy-Schwartz
                                                ///< inequality
    double exx_cauchy_stress_threshold = 1e-07; ///< threshold to screen exx stress using Cauchy-Schwartz
                                                ///< inequality
    std::string exx_ccp_rmesh_times = "default";      ///< how many times larger the radial mesh required for
                                                ///< calculating Columb potential is to that of atomic orbitals
    std::string exx_distribute_type = "htime";  ///< distribute type (assuming default as no specific value
                                                ///< provided)
    int exx_opt_orb_lmax = 0;                   ///< the maximum l of the spherical Bessel functions for opt ABFs
    double exx_opt_orb_ecut = 0.0;              ///< the cut-off of plane wave expansion for opt ABFs
    double exx_opt_orb_tolerence = 0.0;         ///< the threshold when solving for the zeros of spherical Bessel
                                                ///< functions for opt ABFs
    bool exx_symmetry_realspace = true; ///< whether to reduce the real-space sector in when using symmetry=1 in EXX calculation
    double rpa_ccp_rmesh_times = 10.0;          ///< how many times larger the radial mesh required for
                                                ///< calculating Columb potential is to that of atomic orbitals
    bool out_ri_cv = false;   ///<Whether to output the coefficient tensor C and ABFs-representation Coulomb matrix V
    // ==============   #Parameters (16.dft+u) ======================
    //    DFT+U       Xin Qu added on 2020-10-29
    int dft_plus_u = 0;                    ///< 0: standard DFT calculation (default)
    bool dft_plus_dmft = false;            ///< true:DFT+DMFT; false: standard DFT calcullation(default)
    bool yukawa_potential = false;         ///< default: false
    double yukawa_lambda = -1.0;           ///< default: -1.0, which means we calculate lambda
    double uramping_eV = -1.0;             ///< U-Ramping method (eV)
    int omc = 0;                           ///< the mode of occupation matrix control
    double onsite_radius = 0.0;            ///< radius of the sphere for onsite projection (Bohr)
    std::vector<double> hubbard_u_eV = {}; ///< Hubbard Coulomb interaction parameter U(ev)
    std::vector<int> orbital_corr = {};    ///< which correlated orbitals need corrected ; d:2 ,f:3, do not
                                           ///< need correction:-1

    // ==============   #Parameters (17.non-collinear spin-constrained DFT) =========
    /**
     * 0: none spin-constrained DFT;
     * 1: constrain atomic spin;
     */
    bool sc_mag_switch = false;     ///< the switch to open the DeltaSpin function, 0: no
                                    ///< spin-constrained DFT; 1: constrain atomic magnetization
    bool decay_grad_switch = false; ///< the switch to use the local approximation of gradient
                                    ///< decay, 0: no local approximation; 1: apply the method
    double sc_thr = 1e-06;          ///< threshold for spin-constrained DFT in uB
    int nsc = 100;                  ///< maximum number of inner lambda loop
    int nsc_min = 2;                ///< minimum number of inner lambda loop
    int sc_scf_nmin = 2;            ///< minimum number of outer scf loop before initial lambda loop
    double alpha_trial = 0.01;      ///< initial trial step size for lambda in eV/uB^2
    double sccut = 3.0;             ///< restriction of step size in eV/uB
    double sc_scf_thr = 1e-3;       ///< minimum number of outer scf loop before initial lambda loop
    double sc_drop_thr = 1e-3;      ///< threshold for lambda-loop threshold cutoff in spin-constrained DFT

    // ==============   #Parameters (18.Quasiatomic Orbital analysis) =========
    ///<==========================================================
    ///< variables for Quasiatomic Orbital analysis
    ///<==========================================================
    bool qo_switch = false;       ///< 0: no QO analysis; 1: QO analysis
    std::string qo_basis = "szv"; ///< type of QO basis function: hydrogen: hydrogen-like basis,
                                  ///< pswfc: read basis from pseudopotential
    double qo_thr = 1e-06;
    std::vector<std::string> qo_strategy = {};   ///< No specific values provided for qo_strategy
    std::vector<double> qo_screening_coeff = {}; ///< No specific values provided for qo_screening_coeff

    // ==============   #Parameters (19.PEXSI) =====================
    int pexsi_npole = 40;
    bool pexsi_inertia = true;
    int pexsi_nmax = 80;
    bool pexsi_comm = true;
    bool pexsi_storage = true;
    int pexsi_ordering = 0;
    int pexsi_row_ordering = 1;
    int pexsi_nproc = 1;
    bool pexsi_symm = true;
    bool pexsi_trans = false;
    int pexsi_method = 1;
    int pexsi_nproc_pole = 1;
    double pexsi_temp = 0.015;
    double pexsi_gap = 0;
    double pexsi_delta_e = 20.0;
    double pexsi_mu_lower = -10;
    double pexsi_mu_upper = 10;
    double pexsi_mu = 0.0;
    double pexsi_mu_thr = 0.05;
    double pexsi_mu_expand = 0.3;
    double pexsi_mu_guard = 0.2;
    double pexsi_elec_thr = 0.001;
    double pexsi_zero_thr = 1e-10;

    // ==============   #Parameters (20.Test) ====================
    bool out_alllog = false;         ///< output all logs.
    int nurse = 0;                   ///< used for debug.
    bool t_in_h = true;              ///< calculate the T or not.
    bool vl_in_h = true;             ///< calculate the vloc or not.
    bool vnl_in_h = true;            ///< calculate the vnl or not.
    bool vh_in_h = true;             ///< calculate the hartree potential or not
    bool vion_in_h = true;           ///< calculate the local ionic potential or not
                                     ///< //only relevant when vl_in_h = 1
    bool test_force = false;         ///< test the force.
    bool test_stress = false;        ///< test the stress.
    bool test_skip_ewald = false;    ///< variables for test only
    int  test_atom_input = false;    ///< variables for test_atom_input only
    int  test_symmetry = false;      ///< variables for test_lattice only
    int  test_wf = 0;                ///< variables for test_wf only
    int  test_grid = false;          ///< variables for test_grid only
    int  test_charge = false;        ///< variables for test_vloc only
    int  test_energy = false;        ///< variables for test_energy only
    int  test_gridt = false;         ///< variables for test_gridt only
    int  test_pseudo_cell = false;   ///< variables for test_pseudo_cell only
    int  test_pp = 0;                ///< variables for test_pp only
    int  test_relax_method = false;  ///< variables for test_relax_method only
    int  test_deconstructor = false; ///< variables for test_deconstructor only
};
#endif
