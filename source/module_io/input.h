#ifndef INPUT_H
#define INPUT_H

#include <fstream>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include "input_conv.h"
#include "module_base/vector3.h"
#include "module_parameter/md_parameter.h"

class Input
{
  public:
    ~Input()
    {
        delete[] hubbard_u;
        delete[] orbital_corr;
    }

    //==========================================================
    // directories of files
    //==========================================================

    std::string stru_file;     // file contains atomic positions -- xiaohui modify
                               // 2015-02-01
    std::string pseudo_dir;    // directory of pseudopotential
    std::string orbital_dir;   // directory of orbital file
    std::string read_file_dir; // directory of files for reading
    // std::string pseudo_type; // the type of pseudopotential, mohan add
    // 2013-05-20, ABACUS supports
    //                          // UPF format (default) and vwr format. (xiaohui
    //                          add 2013-06-23)
    std::string kpoint_file;  // file contains k-points -- xiaohui modify 2015-02-01
    std::string wannier_card; // input card for wannier functions.
    std::string latname;      // lattice name

    std::string calculation;  // "scf" : self consistent calculation.
                              // "nscf" : non-self consistent calculation.
                              // "relax" : cell relaxations
    std::string esolver_type; // the energy solver: ksdft, sdft, ofdft, tddft, lj, dp
    double pseudo_rcut;       // cut-off radius for calculating msh
    bool pseudo_mesh;         // 0: use msh to normalize radial wave functions;  1: use
                              // mesh, which is used in QE.
    int ntype;                // number of atom types
    int nbands;               // number of bands
    int nbands_istate;        // number of bands around fermi level for get_pchg
                              // calculation.

    int pw_seed; // random seed for initializing wave functions qianrui
                 // 2021-8-12

    bool init_vel;          // read velocity from STRU or not  liuyu 2021-07-14
    double ref_cell_factor; // construct a reference cell bigger than the
                            // initial cell  liuyu 2023-03-21

    int kpar; // ecch pool is for one k point

    bool berry_phase; // berry phase calculation
    int gdir;         // berry phase calculation
    double kspacing[3];
    double min_dist_coef;
    //==========================================================
    // Wannier functions
    //==========================================================
    bool towannier90;         // add by jingan for wannier90
    std::string nnkpfile;     // add by jingan for wannier90
    std::string wannier_spin; // add by jingan for wannier90
    int wannier_method;       // different implementation methods under Lcao basis set
    bool out_wannier_mmn;     // add by renxi for wannier90
    bool out_wannier_amn;
    bool out_wannier_unk;
    bool out_wannier_eig;
    bool out_wannier_wvfn_formatted;

    //==========================================================
    // Stochastic DFT
    //==========================================================
    int nche_sto;              // number of orders for Chebyshev expansion in stochastic DFT
                               // //qinarui 2021-2-5
    int nbands_sto;            // number of stochastic bands //qianrui 2021-2-5
    std::string nbndsto_str;   // string parameter for stochastic bands
    int seed_sto;              // random seed for sDFT
    double initsto_ecut = 0.0; // maximum ecut to init stochastic bands
    double emax_sto;           // Emax & Emin to normalize H
    double emin_sto;
    int bndpar;          // parallel for stochastic/deterministic bands
    int initsto_freq;    // frequency to init stochastic orbitals when running md
    int method_sto;      // different methods for sdft, 1: slow, less memory  2:
                         // fast, more memory
    int npart_sto;       // for method_sto = 2, reduce memory
    bool cal_cond;       // calculate electronic conductivities
    double cond_che_thr; // control the error of Chebyshev expansions for
                         // conductivities
    int cond_smear;      // smearing method for conductivities 1: Gaussian 2:
                         // Lorentzian
    double cond_dw;      // d\omega for conductivities
    double cond_wcut;    // cutoff \omega for conductivities
    double cond_dt;      // dt to integrate conductivities
    int cond_dtbatch;    // exp(iH*dt*cond_dtbatch) is expanded with Chebyshev
                         // expansion.
    double cond_fwhm;    // FWHM for conductivities
    bool cond_nonlocal;  // if calculate nonlocal effects

    //==========================================================
    // electrons / spin
    //==========================================================
    std::string dft_functional; // input DFT functional.
    double xc_temperature;      // only relevant if finite temperature functional is
                                // used
    int nspin;                  // LDA ; LSDA ; non-linear spin
    int lmaxmax;

    //==========================================================
    // LCAO parameters
    //==========================================================
    std::string basis_type; // xiaohui add 2013-09-01, for structural adjustment
    std::string ks_solver;  // xiaohui add 2013-09-01

    //==========================================================
    // Forces
    //==========================================================
    bool cal_force;
    double force_thr; // threshold of force in unit (Ry/Bohr)

    //==========================================================
    // Stress
    //==========================================================
    double stress_thr; // Pengfei Li 2017-11-01 //LiuXh update 20180515
    double press1;
    double press2;
    double press3;
    bool cal_stress; // calculate the stress
    int nstream;
    std::string fixed_axes;   // which axes are fixed
    bool fixed_ibrav;         // whether to keep type of lattice; must be used along
                              // with latname
    bool fixed_atoms;         // whether to fix atoms during vc-relax
    std::string relax_method; // methods to move_ion: sd, bfgs, cg...

    // For now, this is only relevant if we choose to use
    // CG relaxation method. If set to true, then the new
    // implementation will be used; if set to false, then
    // the original implementation will be used
    // Default is true
    bool relax_new;

    double relax_cg_thr; // threshold when cg to bfgs, pengfei add 2011-08-15

    double relax_scale_force;

    //==========================================================
    // Planewave
    //==========================================================
    bool gamma_only;       // for plane wave.
    bool gamma_only_local; // for local orbitals.
    int fft_mode = 0;      // fftw mode 0: estimate, 1: measure, 2: patient, 3: exhaustive

    double ecutwfc; // energy cutoff for wavefunctions
    double ecutrho; // energy cutoff for charge/potential

    double erf_ecut;   // the value of the constant energy cutoff
    double erf_height; // the height of the energy step for reciprocal vectors
    double erf_sigma;  // the width of the energy step for reciprocal vectors

    int ncx, ncy, ncz; // three dimension of FFT charge/grid
    int nx, ny, nz;    // three dimension of FFT wavefunc
    int bx, by, bz;    // big mesh ball. mohan add 2011-04-21
    int ndx, ndy, ndz; // three dimension of FFT smooth charge density

    //==========================================================
    // technique
    //==========================================================
    int diago_proc; // the number of procs used to diag. mohan add 2012-01-13
    int pw_diag_nmax;
    int diago_cg_prec; // mohan add 2012-03-31
    int pw_diag_ndim;
    double pw_diag_thr; // used in cg method

    int nb2d; // matrix 2d division.

    int nurse;    // used for debug.
    int nbspline; // the order of B-spline basis(>=0) if it is -1 (default),
                  // B-spline for Sturcture Factor isnot used.

    bool colour; // used for fun.

    bool t_in_h;   // calculate the T or not.
    bool vl_in_h;  // calculate the vloc or not.
    bool vnl_in_h; // calculate the vnl or not.

    bool vh_in_h;   // calculate the hartree potential or not
    bool vion_in_h; // calculate the local ionic potential or not
    // only relevant when vl_in_h = 1

    bool test_force;  // test the force.
    bool test_stress; // test the stress.

    //==========================================================
    // iteration
    //==========================================================
    double scf_thr;        // \sum |rhog_out - rhog_in |^2
    int scf_thr_type;      // type of the criterion of scf_thr, 1: reci drho, 2: real
                           // drho
    int scf_nmax;          // number of max elec iter
    int relax_nmax;        // number of max ionic iter
    std::string out_level; // control the output information.
    bool out_md_control;   // internal parameter , added by zhengdy 2019-04-07

    //==========================================================
    // occupation
    //==========================================================

    std::string smearing_method; // "gaussian",
                                 // "mp","methfessel-paxton"
                                 // "mv","marzari-vanderbilt","cold"
                                 // "fd","fermi-dirac"
    double smearing_sigma;       //

    //==========================================================
    // potential / charge / wavefunction / energy
    //==========================================================

    std::string init_wfc; // "file","atomic","random"
    std::string init_chg; // "file","atomic"
    bool psi_initializer; // whether use psi_initializer to initialize
                          // wavefunctions

    std::string chg_extrap; // xiaohui modify 2015-02-01

    int mem_saver; // 1: save psi when nscf calculation.

    int printe;                // mohan add 2011-03-16
    int out_freq_elec;         // the frequency ( >= 0) of electronic iter to output
                               // charge density and wavefunction. 0: output only when
                               // converged
    int out_freq_ion;          // the frequency ( >= 0 ) of ionic step to output charge
                               // density and wavefunction. 0: output only when ion steps
                               // are finished
    int out_chg;               // output charge density. 0: no; 1: yes
    int out_pot;               // yes or no
    int out_wfc_pw;            // 0: no; 1: txt; 2: dat
    bool out_wfc_r;            // 0: no; 1: yes
    int out_dos;               // dos calculation. mohan add 20090909
    std::vector<int> out_band; // band calculation pengfei 2014-10-13
    bool out_proj_band;        // projected band structure calculation jiyy add
                               // 2022-05-11
    bool cal_syns;             // calculate asynchronous S matrix to output
    double dmax;               // maximum displacement of all atoms in one step (bohr)
    int out_interval;
    bool out_app_flag; // whether output r(R), H(R), S(R), T(R), and dH(R)
                       // matrices in an append manner during MD  liuyu
                       // 2023-03-20
    int out_ndigits;
    bool out_mat_r;        // jingan add 2019-8-14, output r(R) matrix.
    bool out_alllog;       // output all logs.
    bool out_element_info; // output infomation of all element

    bool out_bandgap; // QO added for bandgap printing

    double dos_emin_ev;
    double dos_emax_ev;
    double dos_edelta_ev;
    double dos_scale;
    int dos_nche;             // orders of Chebyshev expansions for dos
    bool dos_setemin = false; // true: emin is set
    bool dos_setemax = false; // true: emax is set

    double dos_sigma; //  pengfei 2014-10-13

    //==========================================================
    // two center integrals in LCAO
    // mohan add 2009-11-11
    //==========================================================
    double lcao_ecut;     // ecut of two center integral
    double lcao_dk;       // delta k used in two center integral
    double lcao_dr;       // dr used in two center integral
    double lcao_rmax;     // rmax(a.u.) to make table.
    double search_radius; // 11.1
    bool search_pbc;      // 11.2
    double onsite_radius; // the radius of on-site orbitals

    //==========================================================
    // molecular dynamics
    // added by Daye Zheng
    //==========================================================
    MD_para mdp;

    //==========================================================
    // vdw
    // Peize Lin add 2014-03-31, jiyy update 2019-08-01
    //==========================================================
    std::string vdw_method;        // the method of vdw calculation
    std::string vdw_s6;            // scale parameter
    std::string vdw_s8;            // scale parameter
    std::string vdw_a1;            // damping function parameter
    std::string vdw_a2;            // damping function parameter
    double vdw_d;                  // damping function parameter d
    bool vdw_abc;                  // third-order term?
    std::string vdw_cutoff_radius; // cutoff radius for std::pair interactions
    std::string vdw_radius_unit;   //"Bohr" or "Angstrom"
    double vdw_cn_thr;             // cutoff radius for calculating the coordination number
    std::string vdw_cn_thr_unit;   //"Bohr" or "Angstrom"
    std::string vdw_C6_file;
    std::string vdw_C6_unit; //"Bohr" or "Angstrom"
    std::string vdw_R0_file;
    std::string vdw_R0_unit;     //"Bohr" or "Angstrom"
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
    // tddft
    // Fuxiang He add 2016-10-26
    //==========================================================

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
    int dft_plus_u;              ///< 1:DFT+U correction; 2:old DFT+U method; 0:standard DFT
                                 ///< calculation(default)
    int* orbital_corr = nullptr; ///< which correlated orbitals need corrected ;
                                 ///< d:2 ,f:3, do not need correction:-1
    double* hubbard_u = nullptr; ///< Hubbard Coulomb interaction parameter U(ev)
    int omc;                     ///< whether turn on occupation matrix control method or not
    bool yukawa_potential;       ///< default:false
    double yukawa_lambda;        ///< default:-1.0, which means we calculate lambda
    double uramping;             ///< default:-1.0, which means we do not use U-Ramping method

    //==========================================================
    //    DFT+DMFT       Xin Qu added on 2021-08
    //==========================================================
    bool dft_plus_dmft; // true:DFT+DMFT; false: standard DFT
                        // calcullation(default)

    //==========================================================
    //    RPA           Rong Shi added on 2022-04
    //==========================================================
    bool rpa;
    std::string coulomb_type;
    //==========================================================
    // OFDFT  sunliang added on 2022-05-05
    //==========================================================
    std::string of_kinetic;     // Kinetic energy functional, such as TF, VW, WT, TF+
    std::string of_method;      // optimization method, include cg1, cg2, tn (default), bfgs
    std::string of_conv;        // select the convergence criterion, potential, energy
                                // (default), or both
    double of_tole;             // tolerance of the energy change (in Ry) for determining
                                // the convergence, default=2e-6 Ry
    double of_tolp;             // tolerance of potential for determining the convergence,
                                // default=1e-5 in a.u.
    double of_tf_weight;        // weight of TF KEDF
    double of_vw_weight;        // weight of vW KEDF
    double of_wt_alpha;         // parameter alpha of WT KEDF
    double of_wt_beta;          // parameter beta of WT KEDF
    double of_wt_rho0;          // set the average density of system, in Bohr^-3
    bool of_hold_rho0;          // If set to 1, the rho0 will be fixed even if the volume
                                // of system has changed, it will be set to 1 automaticly
                                // if of_wt_rho0 is not zero.
    double of_lkt_a;            // parameter a of LKT KEDF
    bool of_full_pw;            // If set to 1, ecut will be ignored while collecting
                                // planewaves, so that all planewaves will be used.
    int of_full_pw_dim;         // If of_full_pw = 1, the dimention of FFT will be
                                // testricted to be (0) either odd or even; (1) odd
                                // only; (2) even only.
    bool of_read_kernel;        // If set to 1, the kernel of WT KEDF will be filled
                                // from file of_kernel_file, not from formula. Only
                                // usable for WT KEDF.
    std::string of_kernel_file; // The name of WT kernel file.

    //==========================================================
    // spherical bessel  Peize Lin added on 2022-12-15
    //==========================================================
    // the following are used when generating orb_matrix.dat
    // int		bessel_nao_lmax;		// lmax used in descriptor
    bool bessel_nao_smooth;      // spherical bessel smooth or not
    double bessel_nao_sigma;     // spherical bessel smearing_sigma
    std::string bessel_nao_ecut; // energy cutoff for spherical bessel functions(Ry)
    double bessel_nao_rcut;      // radial cutoff for spherical bessel functions(a.u.)
    std::vector<double> bessel_nao_rcuts;
    double bessel_nao_tolerence;        // tolerence for spherical bessel root
                                        // the following are used when generating jle.orb
    int bessel_descriptor_lmax;         // lmax used in descriptor
    bool bessel_descriptor_smooth;      // spherical bessel smooth or not
    double bessel_descriptor_sigma;     // spherical bessel smearing_sigma
    std::string bessel_descriptor_ecut; // energy cutoff for spherical bessel
                                        // functions(Ry)
    double bessel_descriptor_rcut;      // radial cutoff for spherical bessel
                                        // functions(a.u.)
    double bessel_descriptor_tolerence; // tolerence for spherical bessel root

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

    //==========================================================
    // variables for PAW
    //==========================================================
    bool use_paw = false;

    //==========================================================
    // variables for elpa
    //==========================================================

    bool check_input = false;

    std::time_t start_time;
    std::time_t get_start_time(void) const
    {
        return start_time;
    }
    
    //==========================================================
    //Beyond DFT
    //==========================================================
    int lr_nstates;   // the number of 2-particle states to be solved
    int nocc;         // the number of occupied orbitals to form the 2-particle basis
    int nvirt;        // the number of virtual orbitals to form the 2-particle basis (nocc + nvirt <= nbands)
    std::string xc_kernel; // xc kernel for LR-TDDFT
    std::string lr_solver; // the solver for LR-TDDFT
    double lr_thr;  // convergence threshold of the LR-TDDFT eigensolver
    std::vector<double> abs_wavelen_range;  // the range of wavelength(nm) to output the absorption spectrum 
    bool out_wfc_lr;
    double abs_broadening;   // the broadening (eta) for LR-TDDFT absorption spectrum

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

    // start time

    bool Read(const std::string& fn);

    void Default(void);

    void Default_2(void); // jiyy add 2019-08-04

    void Check(void);

#ifdef __MPI
    void Bcast(void);
#endif

    int count_ntype(const std::string& fn); // sunliang add 2022-12-06

    std::string bands_to_print_; // specify the bands to be calculated in the get_pchg
                                 // calculation, formalism similar to ocp_set.

  public:
    template <class T>
    static void read_value(std::ifstream& ifs, T& var)
    {
        ifs >> var;
        std::string line;
        getline(ifs, line); // read the rest of the line, directly discard it.
        return;
    }
    void read_kspacing(std::ifstream& ifs)
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
        // std::cout << "count: " << count << " kspacing: " << kspacing[0] << "
        // " << kspacing[1] << " " << kspacing[2]
        // << std::endl;
    };

    /* I hope this function would be more and more useful if want to support
    vector/list of input */
    template <typename T>
    void read_value2stdvector(std::ifstream& ifs, std::vector<T>& var);
    template <typename T>
    typename std::enable_if<std::is_same<T, double>::value, T>::type cast_string(const std::string& str)
    {
        return std::stod(str);
    }
    template <typename T>
    typename std::enable_if<std::is_same<T, int>::value, T>::type cast_string(const std::string& str)
    {
        if (str == "true" || str == "1")
            return 1;
        else if (str == "false" || str == "0")
            return 0;
        else
            return std::stoi(str);
    }
    template <typename T>
    typename std::enable_if<std::is_same<T, bool>::value, T>::type cast_string(const std::string& str)
    {
        return (str == "true" || str == "1");
    }
    template <typename T>
    typename std::enable_if<std::is_same<T, std::string>::value, T>::type cast_string(const std::string& str)
    {
        return str;
    }
    void strtolower(char* sa, char* sb);
    void read_bool(std::ifstream& ifs, bool& var);

    // Return the const string pointer of private member bands_to_print_
    // Not recommended to use this function directly, use get_out_band_kb()
    // instead
    const std::string* get_bands_to_print() const
    {
        return &bands_to_print_;
    }
    // Return parsed bands_to_print_ as a vector of integers
    std::vector<int> get_out_band_kb() const
    {
        std::vector<int> out_band_kb;
        Input_Conv::parse_expression(bands_to_print_, out_band_kb);
        return out_band_kb;
    }
};

extern Input INPUT;
#endif // INPUT