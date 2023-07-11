#ifndef INPUT_H
#define INPUT_H

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "module_base/vector3.h"
#include "module_md/md_para.h"

class Input
{
  public:
    ~Input()
    {
        delete[] hubbard_u;
        delete[] orbital_corr;
    }
    void Init(const std::string &fn);

    void Print(const std::string &fn) const;

    void close_log(void) const;

    //==========================================================
    // directories of files
    //==========================================================

    std::string suffix; // suffix of out put dir
    std::string stru_file; // file contains atomic positions -- xiaohui modify 2015-02-01
    std::string pseudo_dir; // directory of pseudopotential
    std::string orbital_dir; // directory of orbital file
    std::string read_file_dir; // directory of files for reading
    // std::string pseudo_type; // the type of pseudopotential, mohan add 2013-05-20, ABACUS supports
    //                          // UPF format (default) and vwr format. (xiaohui add 2013-06-23)
    std::string kpoint_file; // file contains k-points -- xiaohui modify 2015-02-01
    std::string wannier_card; // input card for wannier functions.
    std::string latname; // lattice name

    std::string calculation; // "scf" : self consistent calculation.
                             // "nscf" : non-self consistent calculation.
                             // "relax" : cell relaxations
    std::string esolver_type;    // the energy solver: ksdft, sdft, ofdft, tddft, lj, dp
    double pseudo_rcut; // cut-off radius for calculating msh
    bool pseudo_mesh; // 0: use msh to normalize radial wave functions;  1: use mesh, which is used in QE.
    int ntype; // number of atom types
    int nbands; // number of bands
    int nbands_istate; // number of bands around fermi level for get_pchg calculation.
    int pw_seed; // random seed for initializing wave functions qianrui 2021-8-12

    bool init_vel;             // read velocity from STRU or not  liuyu 2021-07-14
    double ref_cell_factor;    // construct a reference cell bigger than the initial cell  liuyu 2023-03-21

    /* symmetry level: 
      -1, no symmetry at all; 
      0, only basic time reversal would be considered; 
      1, point group symmetry would be considered*/
    std::string symmetry; 
    double symmetry_prec; // LiuXh add 2021-08-12, accuracy for symmetry
    int kpar; // ecch pool is for one k point

    bool berry_phase; // berry phase calculation
    int gdir; // berry phase calculation
    double kspacing[3];
    double min_dist_coef;
    //==========================================================
    // Wannier functions
    //==========================================================
    bool towannier90; // add by jingan for wannier90
    std::string nnkpfile; // add by jingan for wannier90
    std::string wannier_spin; // add by jingan for wannier90

    //==========================================================
    // Stochastic DFT
    //==========================================================
    int nche_sto; // number of orders for Chebyshev expansion in stochastic DFT //qinarui 2021-2-5
    int nbands_sto;			// number of stochastic bands //qianrui 2021-2-5
    std::string nbndsto_str; // string parameter for stochastic bands
    int seed_sto; // random seed for sDFT
    double emax_sto; // Emax & Emin to normalize H
    double emin_sto;
    int bndpar; //parallel for stochastic/deterministic bands
    int initsto_freq; //frequency to init stochastic orbitals when running md
    int method_sto; //different methods for sdft, 1: slow, less memory  2: fast, more memory
    int npart_sto; //for method_sto = 2, reduce memory
    bool cal_cond; //calculate electronic conductivities
    int cond_nche; //orders of Chebyshev expansions for conductivities
    double cond_dw; //d\omega for conductivities
    double cond_wcut; //cutoff \omega for conductivities
    double cond_dt;  //dt to integrate conductivities
    int cond_dtbatch; //exp(iH*dt*cond_dtbatch) is expanded with Chebyshev expansion.
    double cond_fwhm; //FWHM for conductivities 
    bool cond_nonlocal; //if calculate nonlocal effects

    //==========================================================
    // electrons / spin
    //==========================================================
    std::string dft_functional; // input DFT functional.
    double xc_temperature; // only relevant if finite temperature functional is used
    int nspin; // LDA ; LSDA ; non-linear spin
    double nupdown = 0.0;
    double nelec; // total number of electrons
    int lmaxmax;

    //==========================================================
    // LCAO parameters
    //==========================================================
    std::string basis_type; // xiaohui add 2013-09-01, for structural adjustment
    std::string ks_solver; // xiaohui add 2013-09-01

    //==========================================================
    // Forces
    //==========================================================
    bool cal_force;
    double force_thr; // threshold of force in unit (Ry/Bohr)
    double force_thr_ev2; // invalid force threshold, mohan add 2011-04-17

    //==========================================================
    // Stress
    //==========================================================
    double stress_thr; // Pengfei Li 2017-11-01 //LiuXh update 20180515
    double press1;
    double press2;
    double press3;
    bool cal_stress; // calculate the stress

    std::string fixed_axes; // which axes are fixed
    bool fixed_ibrav; //whether to keep type of lattice; must be used along with latname
    bool fixed_atoms; //whether to fix atoms during vc-relax
    std::string relax_method; // methods to move_ion: sd, bfgs, cg...

    //For now, this is only relevant if we choose to use
    //CG relaxation method. If set to true, then the new
    //implementation will be used; if set to false, then
    //the original implementation will be used
    //Default is true
    bool relax_new;

    double relax_cg_thr; // threshold when cg to bfgs, pengfei add 2011-08-15

    double relax_bfgs_w1; // wolfe condition 1
    double relax_bfgs_w2; // wolfe condition 2

    double relax_bfgs_rmax; // trust radius max
    double relax_bfgs_rmin; // trust radius min
    double relax_bfgs_init; // initial move

    double relax_scale_force;

    //==========================================================
    // Planewave
    //==========================================================
    bool gamma_only; // for plane wave.
    bool gamma_only_local; // for local orbitals.

    double ecutwfc; // energy cutoff for wavefunctions
    double ecutrho; // energy cutoff for charge/potential

    int ncx, ncy, ncz; // three dimension of FFT charge/grid
    int nx, ny, nz; // three dimension of FFT wavefunc
    int bx, by, bz; // big mesh ball. mohan add 2011-04-21

    //==========================================================
    // technique
    //==========================================================
    int diago_proc; // the number of procs used to diag. mohan add 2012-01-13
    int pw_diag_nmax;
    int diago_cg_prec; // mohan add 2012-03-31
    int pw_diag_ndim;
    double pw_diag_thr; // used in cg method

    int nb2d; // matrix 2d division.

    int nurse; // used for debug.
    int nbspline; // the order of B-spline basis(>=0) if it is -1 (default), B-spline for Sturcture Factor isnot used.

    bool colour; // used for fun.

    bool t_in_h; // calculate the T or not.
    bool vl_in_h; // calculate the vloc or not.
    bool vnl_in_h; // calculate the vnl or not.

    bool vh_in_h; // calculate the hartree potential or not
    bool vion_in_h; // calculate the local ionic potential or not
    // only relevant when vl_in_h = 1

    bool test_force; // test the force.
    bool test_stress; // test the stress.

    //==========================================================
    // iteration
    //==========================================================
    double scf_thr; // \sum |rhog_out - rhog_in |^2
    int scf_thr_type; // type of the criterion of scf_thr, 1: reci drho, 2: real drho
    int scf_nmax; // number of max elec iter
    int relax_nmax; // number of max ionic iter
    bool out_stru; // outut stru file each ion step
    std::string out_level; // control the output information.
    bool out_md_control; // internal parameter , added by zhengdy 2019-04-07

    //==========================================================
    // occupation
    //==========================================================
    std::string occupations; // "fixed","smearing","tetrahedra","from_input"

    std::string smearing_method; // "gaussian",
                                 // "mp","methfessel-paxton"
                                 // "mv","marzari-vanderbilt","cold"
                                 // "fd","fermi-dirac"
    double smearing_sigma; //

    //==========================================================
    // charge mixing
    //==========================================================
    std::string mixing_mode; // "plain","broyden",...
    double mixing_beta; // 0 : no_mixing
    int mixing_ndim; // used in Broyden method
    double mixing_gg0; // used in kerker method. mohan add 2014-09-27
    bool mixing_tau; // whether to mix tau in mgga
    bool mixing_dftu; //whether to mix locale in DFT+U

    //==========================================================
    // potential / charge / wavefunction / energy
    //==========================================================

    std::string init_wfc; // "file","atomic","random"
    std::string init_chg; // "file","atomic"

    std::string chg_extrap; // xiaohui modify 2015-02-01

    int mem_saver; // 1: save psi when nscf calculation.

    int printe; // mohan add 2011-03-16
    int out_freq_elec;  // the frequency ( >= 0) of electronic iter to output charge density and wavefunction. 0: output only when converged
    int out_freq_ion;  // the frequency ( >= 0 ) of ionic step to output charge density and wavefunction. 0: output only when ion steps are finished
    bool out_chg; // output charge density. 0: no; 1: yes
    bool out_dm; // output density matrix.
    bool out_dm1;
    int out_pot; // yes or no
    int out_wfc_pw; // 0: no; 1: txt; 2: dat
    bool out_wfc_r; // 0: no; 1: yes
    int out_dos; // dos calculation. mohan add 20090909
    bool out_band; // band calculation pengfei 2014-10-13
    bool out_proj_band; // projected band structure calculation jiyy add 2022-05-11
    bool out_mat_hs; // output H matrix and S matrix in local basis.
    bool cal_syns; // calculate asynchronous S matrix to output
    double dmax; // maximum displacement of all atoms in one step (bohr)
    bool out_mat_hs2; // LiuXh add 2019-07-16, output H(R) matrix and S(R) matrix in local basis.
    bool out_mat_dh;
    int out_interval;
    bool out_app_flag;    // whether output r(R), H(R), S(R), T(R), and dH(R) matrices in an append manner during MD  liuyu 2023-03-20
    bool out_mat_t;
    bool out_mat_r; // jingan add 2019-8-14, output r(R) matrix.
    bool out_wfc_lcao; // output the wave functions in local basis.
    bool out_alllog; // output all logs.
    bool out_element_info; // output infomation of all element

    bool out_bandgap; // QO added for bandgap printing
    
    double dos_emin_ev;
    double dos_emax_ev;
    double dos_edelta_ev;
    double dos_scale;
    int dos_nche; //orders of Chebyshev expansions for dos
    bool dos_setemin = false; //true: emin is set
    bool dos_setemax = false; //true: emax is set

    double dos_sigma; //  pengfei 2014-10-13

    //==========================================================
    // two center integrals in LCAO
    // mohan add 2009-11-11
    //==========================================================
    double lcao_ecut; // ecut of two center integral
    double lcao_dk; // delta k used in two center integral
    double lcao_dr; // dr used in two center integral
    double lcao_rmax; // rmax(a.u.) to make table.
    double search_radius; // 11.1
    bool search_pbc; // 11.2

    //==========================================================
    // molecular dynamics
    // added by Daye Zheng
    //==========================================================
    MD_para mdp;

    //==========================================================
    // efield and dipole correction
    // Yu Liu add 2022-05-18
    //==========================================================
    bool efield_flag;        // add electric field
    bool dip_cor_flag;        // dipole correction
    int efield_dir;           // the direction of the electric field or dipole correction
    double efield_pos_max;     // position of the maximum of the saw-like potential along crystal axis efield_dir
    double efield_pos_dec;      // zone in the unit cell where the saw-like potential decreases
    double efield_amp ;        // amplitude of the electric field

    //==========================================================
    // gatefield (compensating charge)
    // Yu Liu add 2022-09-13
    //==========================================================
    bool gate_flag;                 // compensating charge or not
    double zgate;                   // position of charged plate
    bool relax;                     // allow relaxation along the specific direction
    bool block;                     // add a block potential or not
    double block_down;              // low bound of the block
    double block_up;                // high bound of the block
    double block_height;            // height of the block

    //==========================================================
    // vdw
    // Peize Lin add 2014-03-31, jiyy update 2019-08-01
    //==========================================================
    std::string vdw_method; // the method of vdw calculation
    std::string vdw_s6; // scale parameter
    std::string vdw_s8; // scale parameter
    std::string vdw_a1; // damping function parameter
    std::string vdw_a2; // damping function parameter
    double vdw_d; // damping function parameter d
    bool vdw_abc; // third-order term?
    std::string vdw_cutoff_radius; // cutoff radius for std::pair interactions
    std::string vdw_radius_unit; //"Bohr" or "Angstrom"
    double vdw_cn_thr; // cutoff radius for calculating the coordination number
    std::string vdw_cn_thr_unit; //"Bohr" or "Angstrom"
    std::string vdw_C6_file;
    std::string vdw_C6_unit; //"Bohr" or "Angstrom"
    std::string vdw_R0_file;
    std::string vdw_R0_unit; //"Bohr" or "Angstrom"
    std::string vdw_cutoff_type; //"period" or "radius"
    ModuleBase::Vector3<int> vdw_cutoff_period;

    bool ocp;
    std::string ocp_set;
    bool out_mul; // qifeng add 2019-9-10
    // added by zhengdy-soc
    bool noncolin;
    bool lspinorb;
    double soc_lambda;

    //==========================================================
    // exx
    // Peize Lin add 2018-06-20
    //==========================================================
    std::string exx_hybrid_alpha;
    double exx_hse_omega;

    bool exx_separate_loop; // 0 or 1
    int exx_hybrid_step;
    double exx_mixing_beta; // only for exx_separate_loop=1

    double exx_lambda;

    std::string exx_real_number;
    double exx_pca_threshold;
    double exx_c_threshold;
    double exx_v_threshold;
    double exx_dm_threshold;
    double exx_schwarz_threshold;
    double exx_cauchy_threshold;
    double exx_c_grad_threshold;
    double exx_v_grad_threshold;
    double exx_cauchy_force_threshold;
    double exx_cauchy_stress_threshold;
    double exx_ccp_threshold;
    std::string exx_ccp_rmesh_times;

    std::string exx_distribute_type;

    int exx_opt_orb_lmax;
    double exx_opt_orb_ecut;
    double exx_opt_orb_tolerence;

    //==========================================================
    // tddft
    // Fuxiang He add 2016-10-26
    //==========================================================
    double td_force_dt; //"fs"
    bool td_vext; // add extern potential or not
    std::string td_vext_dire; // vext direction
    bool out_dipole; // output the dipole or not
    bool out_efield; // output the efield or not

    double td_print_eij; // threshold to output Eij elements
    int td_edm; //0: new edm method   1: old edm method

    int propagator; // method of propagator

    int td_stype ; //type of space domain  0 : length gauge  1: velocity gauge

    std::string td_ttype ; //type of time domain
    //  0  Gauss type function.
    //  1  trapezoid type function.
    //  2  Trigonometric functions, sin^2.
    //  3  heaviside function.
    //  4  HHG function.

    int td_tstart ;
    int td_tend ;

    // space domain parameters

    //length gauge
    double td_lcut1;
    double td_lcut2;

    // time domain parameters

    // Gauss
    std::string td_gauss_freq ; // time(fs)^-1 
    std::string td_gauss_phase ;
    std::string td_gauss_sigma ; // time(fs)
    std::string td_gauss_t0 ;
    std::string td_gauss_amp ;  // V/A

    // trapezoid
    std::string td_trape_freq ; // time(fs)^-1 
    std::string td_trape_phase ;
    std::string td_trape_t1 ;
    std::string td_trape_t2 ;
    std::string td_trape_t3 ;
    std::string td_trape_amp ; // V/A

    // Trigonometric
    std::string td_trigo_freq1 ; // time(fs)^-1 
    std::string td_trigo_freq2 ; // time(fs)^-1 
    std::string td_trigo_phase1 ;
    std::string td_trigo_phase2 ;
    std::string td_trigo_amp ; // V/A

    // Heaviside
    std::string td_heavi_t0;
    std::string td_heavi_amp ; // V/A

    // HHG
    // std::string td_hhg_amp1; // V/A
    // std::string td_hhg_amp2; // V/A
    // std::string td_hhg_freq1; // time(fs)^-1 
    // std::string td_hhg_freq2; // time(fs)^-1 
    // std::string td_hhg_phase1;
    // std::string td_hhg_phase2;
    // std::string td_hhg_t0;
    // std::string td_hhg_sigma; // time(fs)

    //==========================================================
    // restart
    // Peize Lin add 2020-04-04
    //==========================================================
    bool restart_save;
    bool restart_load;
    // xiaohui add 2015-09-16
    bool input_error;
    double cell_factor; // LiuXh add 20180619

    //==========================================================
    //    DFT+U       Xin Qu added on 2020-10-29
    //==========================================================
    bool dft_plus_u;             ///< true:DFT+U correction; false: standard DFT calculation(default)
    int* orbital_corr = nullptr; ///< which correlated orbitals need corrected ; d:2 ,f:3, do not need correction:-1
    double* hubbard_u = nullptr; ///< Hubbard Coulomb interaction parameter U(ev)
    int omc;                     ///< whether turn on occupation matrix control method or not
    bool yukawa_potential;       ///< default:false
    double yukawa_lambda;        ///< default:-1.0, which means we calculate lambda

    //==========================================================
    //    DFT+DMFT       Xin Qu added on 2021-08
    //==========================================================
    bool dft_plus_dmft; // true:DFT+DMFT; false: standard DFT calcullation(default)

    //==========================================================
    //    RPA           Rong Shi added on 2022-04
    //==========================================================
    bool rpa;
    std::string coulomb_type;

    //==========================================================
    // DeepKS -- added by caoyu and mohan
    //==========================================================
    bool deepks_out_labels; // (need libnpy) prints energy and force labels and descriptors for training, wenfei
                            // 2022-1-12
    bool deepks_scf; //(need libnpy and libtorch) if set 1, a trained model would be needed to cal V_delta and F_delta
    bool deepks_bandgap; // for bandgap label. QO added 2021-12-15

    bool deepks_out_unittest; // if set 1, prints intermediate quantities that shall be used for making unit test

    std::string deepks_model; // needed when deepks_scf=1

    //==========================================================
    //    implicit solvation model       Menglin Sun added on 2022-04-04
    //==========================================================
    bool imp_sol; // true:implicit solvation correction; false: vacuum calculation(default)
    double eb_k;
    double tau;
    double sigma_k;
    double nc_k;

    //==========================================================
    // OFDFT  sunliang added on 2022-05-05
    //==========================================================
    std::string of_kinetic; // Kinetic energy functional, such as TF, VW, WT, TF+
    std::string of_method;  // optimization method, include cg1, cg2, tn (default), bfgs
    std::string of_conv;    // select the convergence criterion, potential, energy (default), or both
    double of_tole;    // tolerance of the energy change (in Ry) for determining the convergence, default=2e-6 Ry
    double of_tolp;    // tolerance of potential for determining the convergence, default=1e-5 in a.u.
    double of_tf_weight;  // weight of TF KEDF
    double of_vw_weight;  // weight of vW KEDF
    double of_wt_alpha;   // parameter alpha of WT KEDF
    double of_wt_beta;    // parameter beta of WT KEDF
    double of_wt_rho0;    // set the average density of system, in Bohr^-3
    bool of_hold_rho0;  // If set to 1, the rho0 will be fixed even if the volume of system has changed, it will be set to 1 automaticly if of_wt_rho0 is not zero.
    double of_lkt_a;    // parameter a of LKT KEDF
    bool of_full_pw;    // If set to 1, ecut will be ignored while collecting planewaves, so that all planewaves will be used.
    int of_full_pw_dim; // If of_full_pw = 1, the dimention of FFT will be testricted to be (0) either odd or even; (1) odd only; (2) even only.
    bool of_read_kernel; // If set to 1, the kernel of WT KEDF will be filled from file of_kernel_file, not from formula. Only usable for WT KEDF.
    std::string of_kernel_file; // The name of WT kernel file.

    //==========================================================
    // spherical bessel  Peize Lin added on 2022-12-15
    //==========================================================
    // the following are used when generating orb_matrix.dat
		//int		bessel_nao_lmax;		// lmax used in descriptor
	bool	bessel_nao_smooth;		// spherical bessel smooth or not
	double	bessel_nao_sigma;		// spherical bessel smearing_sigma
	std::string	bessel_nao_ecut;		// energy cutoff for spherical bessel functions(Ry)
	double	bessel_nao_rcut;		// radial cutoff for spherical bessel functions(a.u.)
	double	bessel_nao_tolerence;	// tolerence for spherical bessel root
    // the following are used when generating jle.orb
	int		bessel_descriptor_lmax;			// lmax used in descriptor
	bool	bessel_descriptor_smooth;		// spherical bessel smooth or not
	double	bessel_descriptor_sigma;		// spherical bessel smearing_sigma
	std::string	bessel_descriptor_ecut;			// energy cutoff for spherical bessel functions(Ry)
	double	bessel_descriptor_rcut;			// radial cutoff for spherical bessel functions(a.u.)
	double	bessel_descriptor_tolerence;	// tolerence for spherical bessel root

    //==========================================================
    //    device control denghui added on 2022-11-15
    //==========================================================
    std::string device;
    //==========================================================
    //    precision control denghui added on 2023-01-01
    //==========================================================
    std::string precision;

    //==========================================================
    // variables for test only
    //==========================================================
    bool test_skip_ewald = false;

  private:
    //==========================================================
    // MEMBER FUNCTIONS :
    // NAME : Read()
    // NAME : Default(set the default value for variables)
    // NAME : Check(check the values)
    // NAME : Bcast(only used in paralle case,only read in data
    //        from first cpu, and distribute the data to the
    //        other processors)
    //==========================================================

    bool Read(const std::string &fn);

    void Default(void);

    void Default_2(void); // jiyy add 2019-08-04

    void Check(void);

#ifdef __MPI
    void Bcast(void);
#endif

    int count_ntype(const std::string &fn); // sunliang add 2022-12-06

  public:
    template <class T> static void read_value(std::ifstream &ifs, T &var)
    {
        ifs >> var;
        std::string line;
        getline(ifs, line);
        return;
    }
    void read_kspacing(std::ifstream &ifs)
    {
        std::string s;
        std::getline(ifs, s);
        std::stringstream ss(s);
        // read 3 values
        int count = 0;
        while ((ss >> kspacing[count]) && count < 3)
        {
            count++;
        }
        // if not read even one value, or read two values, the input is invalid.
        if (count == 0 || count == 2)
        {
            std::cout << "kspacing can only accept one or three double values." << std::endl;
            ifs.setstate(std::ios::failbit);
        }
        // if only read one value, set all to kspacing[0]
        if (count == 1)
        {
            kspacing[1] = kspacing[0];
            kspacing[2] = kspacing[0];
        }
        // std::cout << "count: " << count << " kspacing: " << kspacing[0] << " " << kspacing[1] << " " << kspacing[2]
        // << std::endl;
    };

    void strtolower(char *sa, char *sb);
    void read_bool(std::ifstream &ifs, bool &var);
};

extern Input INPUT;
#endif // INPUT
