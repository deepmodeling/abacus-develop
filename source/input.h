#ifndef INPUT_H
#define INPUT_H

#include "module_base/vector3.h"
#include <fstream>
#include <string>
#include "module_md/MD_parameters.h"

using namespace std;

class Input
{
	public:
    
	Input();
    ~Input();

    void Init(const string &fn);

    void Print(const string &fn)const;

	void close_log(void)const;

//==========================================================
// directories of files 
//==========================================================

    string suffix;			// suffix of out put dir
    string atom_file;		// file contains atomic positions -- xiaohui modify 2015-02-01
    string pseudo_dir;      // directory of pseudopotential
	string read_file_dir;   // directory of files for reading
    string pseudo_type;     // the type of pseudopotential, mohan add 2013-05-20, ABACUS supports
			    			// UPF format (default) and vwr format. (xiaohui add 2013-06-23)
    string kpoint_file;		// file contains k-points -- xiaohui modify 2015-02-01
	string wannier_card;	// input card for wannier functions.
    string latname;			// lattice name


    string calculation;		// "scf" : self consistent calculation.
						    // "nscf" : non-self consistent calculation.
							// "relax" : cell relaxations
    double pseudo_rcut;     // cut-off radius for calculating msh
	bool renormwithmesh;     // 0: use msh to normalize radial wave functions;  1: use mesh, which is used in QE.
    int ntype;				// number of atom types
    int nbands;				// number of bands
	int nbands_sto;			// number of stochastic bands //qianrui 2021-2-5
	int nbands_istate;		// number of bands around fermi level for istate calculation.
	int nche_sto;			// number of orders for Chebyshev expansion in stochastic DFT //qinarui 2021-2-5
	int seed_sto;  // random seed for sDFT
	double emax_sto;		//Emax & Emin to normalize H
	double emin_sto;
	string stotype;

	bool set_vel;           // read velocity from STRU or not  liuyu 2021-07-14

    bool symmetry;			// turn on symmetry or not
	int npool; 				// ecch pool is for one k point

    bool berry_phase;		// berry phase calculation
	int gdir;               // berry phase calculation

//==========================================================
// Wannier functions 
//==========================================================
	bool towannier90;       // add by jingan for wannier90
	string NNKP;            // add by jingan for wannier90
	string wannier_spin;    // add by jingan for wannier90

	bool mlwf_flag; 		// add by mohan

//==========================================================
// E field 
//==========================================================
    int efield;				// add electrical field
	int edir;
	double emaxpos;
	double eopreg;
	double eamp;

//==========================================================
// Optical properties 
//==========================================================
	bool opt_epsilon2;		// true : calculate the dielectric functions
	int  opt_nbands;		// number of bands for optical transition matrix
    bool lda_plus_u;		// true : lda plus u calculation

//==========================================================
// electrons / spin
//==========================================================
    string dft_functional;	// input DFT functional.
	int nspin;				// LDA ; LSDA ; non-linear spin
    double nelec;			// total number of electrons
    int lmaxmax;
    double tot_magnetization;

//==========================================================
// LCAO parameters 
//==========================================================
	string basis_type; 			//xiaohui add 2013-09-01, for structural adjustment
	string ks_solver; 			//xiaohui add 2013-09-01

//==========================================================
// Forces 
//==========================================================
    int force;
    bool force_set;
    double force_thr;		// threshold of force in unit (Ry/Bohr)
	double force_thr_ev2;	// invalid force threshold, mohan add 2011-04-17

//==========================================================
// Stress 
//==========================================================
    double stress_thr;      // Pengfei Li 2017-11-01 //LiuXh update 20180515
    double press1;
    double press2;
    double press3;
	bool stress;			// calculate the stress

	string fixed_axes;      //which axes are fixed
	string ion_dynamics;	// methods to move_ion: sd, bfgs, cg...

    double cg_threshold;    // threshold when cg to bfgs, pengfei add 2011-08-15

	double bfgs_w1;			// wolfe condition 1	
	double bfgs_w2;			// wolfe condition 2

	double trust_radius_max;	// trust radius max
	double trust_radius_min;	// trust radius min
	double trust_radius_ini;	// initial move

//==========================================================
// Planewave
//==========================================================
    bool gamma_only;		// for plane wave.
	bool gamma_only_local;	// for local orbitals.

    double ecutwfc;			// energy cutoff for wavefunctions
    double ecutrho;			// energy cutoff for charge/potential

    int ncx,ncy,ncz;		// three dimension of FFT charge/grid
    int nx,ny,nz;			// three dimension of FFT wavefunc
	int bx,by,bz;			// big mesh ball. mohan add 2011-04-21 

//==========================================================
// technique 
//==========================================================
	int diago_proc;			// the number of procs used to diag. mohan add 2012-01-13
    int diago_cg_maxiter;
	int diago_cg_prec;		// mohan add 2012-03-31
    int diago_david_ndim;
    double ethr;			// used in cg method

	int nb2d;				// matrix 2d division.

	int nurse;				// used for debug.

	bool colour;			// used for fun.

	int t_in_h;				// calculate the T or not.
	int vl_in_h;			// calculate the vloc or not.
	int vnl_in_h;			// calculate the vnl or not.

	int vh_in_h;			// calculate the hartree potential or not
	int vxc_in_h;			// calculate the xc potential or not
	int vion_in_h;			// calculate the local ionic potential or not
	//only relevant when vl_in_h = 1

	int test_force;			// test the force.
	int test_stress;		// test the stress.

//==========================================================
// iteration
//==========================================================
    double dr2;				// \sum |rhog_out - rhog_in |^2
    int niter;				// number of max elec iter
    int nstep;				// number of max ionic iter
	int out_stru;			// outut stru file each ion step
	string out_level;		// control the output information.
    bool out_md_control;    // internal parameter , added by zhengdy 2019-04-07

//==========================================================
// occupation
//==========================================================
    string occupations;		// "fixed","smearing","tetrahedra","from_input"

    string smearing;		// "gaussian",
						    // "mp","methfessel-paxton"
						    // "mv","marzari-vanderbilt","cold"
						    // "fd","fermi-dirac"
    double degauss;			//

//==========================================================
// charge mixing
//==========================================================
    string mixing_mode;		// "plain","broden",...
    double mixing_beta;		// 0 : no_mixing
    int mixing_ndim;		// used in Broden method
	double mixing_gg0;      // used in kerker method. mohan add 2014-09-27

//==========================================================
// potential / charge / wavefunction / energy
//==========================================================
    string restart_mode;	//

    string start_wfc;		// "file","atomic","random"
    string start_pot;		// "file","atomic"

	string charge_extrap;	// xiaohui modify 2015-02-01

	int mem_saver;			// 1: save psi when nscf calculation.

	int printe;				// mohan add 2011-03-16
    int out_charge;		// output charge density.
	int out_dm; // output density matrix.
	int out_potential;		// yes or no
    bool out_wf;			// yes or no
	int out_dos;			// dos calculation. mohan add 20090909
    int out_band;                   // band calculation pengfei 2014-10-13
	int out_hs;			// output H matrix and S matrix in local basis.
	int out_hs2;			//LiuXh add 2019-07-16, output H(R) matrix and S(R) matrix in local basis.
	int out_r_matrix;   // jingan add 2019-8-14, output r(R) matrix.
	bool out_lowf;			// output the wave functions in local basis.
	bool out_alllog; 		// output all logs.

	double dos_emin_ev;
	double dos_emax_ev;
	double dos_edelta_ev;
	double dos_scale;

    double b_coef;           //  pengfei 2014-10-13

//==========================================================
// two center integrals in LCAO
// mohan add 2009-11-11 
//==========================================================
	double lcao_ecut;		// ecut of two center integral
	double lcao_dk;			// delta k used in two center integral
	double lcao_dr;			// dr used in two center integral
	double lcao_rmax;		// rmax(a.u.) to make table.
    double search_radius;	//11.1
	bool search_pbc; // 11.2


//==========================================================
// selected inversion method 
//==========================================================
	// selinv method parameter (cooperate with LinLin)
	int selinv_npole;
	double selinv_temp;
	double selinv_gap;
	double selinv_deltae;
	double selinv_mu;
	double selinv_threshold;
	int selinv_niter;

//==========================================================
// molecular dynamics 
// added by Daye Zheng
//==========================================================
/*    int md_mdtype;                   //choose ensemble
	double md_tauthermo;
	double md_taubaro;
	double md_dt;                    //time step
	int md_nresn;                     //parameter during integrater
	int md_nyosh;                      //parameter during integrater
	double md_qmass;                   //mass of thermostat
	double md_tfirst;                    //temperature begin
	double md_tlast;                    //temperature end
	int md_dumpmdfred;                  //The period to dump MD information for monitoring and restarting MD
	string md_mdoutpath;                //output path for md
	bool md_domsd;                   //whether compute <r(t)-r(0)>
	bool md_domsdatom;                //whether compute msd for each atom
	int md_rstmd;                    //whether restart;
	int md_outputstressperiod;      //period to output stress
	int md_fixtemperature;          //period to change temperature
	double md_ediff;             //parameter for constraining total energy change
	double md_ediffg;             //parameter for constraining max force change
	int md_msdstartTime;            //choose which step that msd be calculated */
	MD_parameters mdp;

//==========================================================
// vdw
// Peize Lin add 2014-03-31, jiyy update 2019-08-01
//==========================================================
    string vdw_method;          //the method of vdw calculation                 
    string vdw_s6;              //scale parameter 
	string vdw_s8;              //scale parameter 
	string vdw_a1;             //damping function parameter 
	string vdw_a2;             //damping function parameter 
	double vdw_d;               //damping function parameter d
	bool vdw_abc;               //third-order term?
    string vdw_radius;          //cutoff radius for pair interactions
	string vdw_radius_unit;	    //"Bohr" or "Angstrom"
	double vdw_cn_thr;          //cutoff radius for calculating the coordination number
	string vdw_cn_thr_unit;     //"Bohr" or "Angstrom"
	string vdw_C6_file;
	string vdw_C6_unit;		    //"Bohr" or "Angstrom"
	string vdw_R0_file;
	string vdw_R0_unit;		    //"Bohr" or "Angstrom"
	string vdw_model;			//"period" or "radius"
	Vector3<int> vdw_period;

//==========================================================
// Spectrum
// pengfei Li add 2016-11-23
//==========================================================    
	//bool     epsilon;               // calculate epsilon or not
	string   spectral_type;          // the type of the calculated spectrum
	int      spectral_method;        // 0: tddft(linear response)
	int      eels_method;            // 0: hilbert_transform method; 1: standard method
	int      absorption_method;      // 0: vasp's method  1: pwscf's method
	//int		 epsilon_choice;         // 0: hilbert_transform method; 1: standard method
	string   kernel_type;           // the kernel type: rpa, tdlda ...
	string system_type;                 // bulk or surface
	double  eta;                   // unit(Ry)
	double  domega;                // unit(Ry)
	int     nomega;
	int     ecut_chi;                   // the dimension of G 
	//int     oband;                 // the number of "occupied" bands  
	double  q_start[3];            // the position of the first q point in direct coordinate
	double  q_direct[3];             // the q direction
	//int     start_q;               // the serial number of the start qpoint
	//int     interval_q;            // the interval of the qpoints
	int     nq;                    // the total number of qpoints for calculation
	bool     out_epsilon;           // output epsilon or not
	bool     out_chi;               // output chi or not
	bool     out_chi0;              // output chi0 or not
	double  fermi_level;            // the change the fermi level(Ry)
	bool     coulomb_cutoff;         // turn on or off the Coulomb_cutoff 0/1
	//bool     epsilon0;              // calculate the macroscopic dielectric constant or not
	//double   intersmear;            // eta
	double   intrasmear;            // Eta
	double   shift;
	bool     metalcalc;             // metal or not
	double   eps_degauss;            // degauss
	//int	epsilon0_choice;             // 0: vasp's method  1: pwscf's method
	bool     kmesh_interpolation;          // calculting <i,0|j,R>
	double  qcar[100][3];          // the Cartesian position of q points(unit: 2*PI/lat0) 
	int ocp;
	//int ocp_n;
	string ocp_set;
	//double  ocp_kb[10000];
	int     lcao_box[3];           // the scale for searching the existence of the overlap <i,0|j,R>
	int    mulliken;//qifeng add 2019-9-10
	//added by zhengdy-soc
	bool noncolin;
	bool lspinorb;
	double soc_lambda;

	//	bool starting_spin_angle;
	double *angle1;
	double *angle2;
		

//==========================================================
// exx
// Peize Lin add 2018-06-20
//==========================================================		
	string exx_hybrid_type;		// "no", "hf", "pbe0", "hse"

	double exx_hybrid_alpha;
	double exx_hse_omega;
	
	bool exx_separate_loop;		// 0 or 1
	int exx_hybrid_step;
	
	double exx_lambda;
	
	double exx_pca_threshold;
	double exx_c_threshold;
	double exx_v_threshold;
	double exx_dm_threshold;
	double exx_schwarz_threshold;
	double exx_cauchy_threshold;
	double exx_ccp_threshold;
	double exx_ccp_rmesh_times;
	
	string exx_distribute_type;
	
	int exx_opt_orb_lmax;
	double exx_opt_orb_ecut;
	double exx_opt_orb_tolerence;
	
//==========================================================
// tddft
// Fuxiang He add 2016-10-26
//==========================================================
	int tddft;			//calculate tddft or not
	double td_dr2;			//threshold for electronic iteration of tddft
	double td_dt;			//"fs"
	double td_force_dt;			//"fs"
	int td_val_elec_01;			//valence electron 01
	int td_val_elec_02;			//valence electron 02
	int td_val_elec_03;			//valence electron 03
	int td_vext;			//add extern potential or not
	int td_vext_dire;			//vext direction
	double td_timescale;			//"fs"
	int td_vexttype;
	int td_vextout; 			// output the electronic potential or not
	int td_dipoleout;			// output the dipole or not


	
//==========================================================
// restart
// Peize Lin add 2020-04-04
//==========================================================
	bool restart_save;
	bool restart_load;
	//xiaohui add 2015-09-16
	bool input_error;
    double cell_factor; //LiuXh add 20180619
	
//==========================================================
// new DM algorithm test
// add by Shen Yu @ 2019/5/9
// newDM values:
//  0: not use new DM algorithm; 
//  1: use new DM algorithm and show no debug information
//  2: use new DM algorithm and only show key debug information
//  3: use new DM algorithm and show all detail debug information
//==========================================================
    int new_dm;

//==========================================================
//    DFT+U       Xin Qu added on 2020-10-29
//==========================================================
    bool dft_plus_u;                //true:DFT+U correction; false：standard DFT calcullation(default)	
	int *orbital_corr;              // which correlated orbitals need corrected ; d:2 ,f:3, do not need correction:-1
    double *hubbard_u;              //Hubbard Coulomb interaction parameter U(ev)
	double *hund_j;                 //Hund exchange parameter J(ev) 	
	bool omc;                       //whether turn on occupation matrix control method or not
	bool yukawa_potential;          //default:false
	double yukawa_lambda;            //default:0.0

	//The two parameters below are not usable currently

	int dftu_type;                  //1:rotationally invarient formalism; 2:simplified form(default)
	int double_counting;            // 1:FLL(fully localized limit)(default); 2:AMF(around mean field)

//==========================================================
// DeepKS -- added by caoyu and mohan 
//==========================================================
    int out_descriptor; // output descritpor for deepks. caoyu added 2020-11-24, mohan modified 2021-01-03
	int lmax_descriptor; // lmax used in descriptor, mohan added 2021-01-03
	int deepks_scf;	//if set 1, a trained model would be needed to cal V_delta and F_delta
	string model_file;		//needed when deepks_scf=1

//==========================================================
// variables for test only
//==========================================================
	bool test_just_neighbor;
	
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

    bool Read(const string &fn);

    void Default(void);

	void Default_2(void);    //jiyy add 2019-08-04											  

    void Check(void);

#ifdef __MPI
    void Bcast(void);
#endif

	public:

    template <class T>
    static void read_value(ifstream &ifs, T &var)
    {
        ifs >> var;
        ifs.ignore(150, '\n');
        return;
    }

	void strtolower(char *sa, char *sb);
	void readbool(ifstream &ifs, bool &var);
};

extern Input INPUT;
#endif //INPUT
