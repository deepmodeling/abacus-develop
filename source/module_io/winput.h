//==========================================================
// AUTHOR : mohan
// DATE : 2008-11-08
// Last Update: 2010-06-07
//==========================================================
#ifndef WINPUT_H
#define WINPUT_H

#include <string>
#include <fstream>

//==========================================================
// CLASS :
// NAME : winput
// ( Readin wannier parameters.
//   Check wannier parameters.
//   Print wannier parameters )
//==========================================================
class winput 
{
public:

	winput();
	~winput();

//==========================================================
// MEMBER FUNCTIONS : 
// NAME : target( no use now)
// NAME : before_iter ( call wannier::runnning )
// NAME : after_iter ( call wannier::running )
// NAME : begin_stop_flag ( only use readin information ,
// 		  stop at the very beginning ).
// NAME : end_flag ( output data, etc. )
//==========================================================
	static std::string target;
	static bool before_iter;
	static bool after_iter;
	static bool begin_stop_flag;
	static bool end_flag;

//==========================================================
// MEMBER FUNCTIONS : 
// NAME : rcut( PAO cutoff )
// NAME : trunc_ao ( PAO cutoff )
// NAME : trunc_wlmr ( 1D Radial wave function cutoff )
// NAME : trunc_wan ( 3D wannier function cutoff )
//==========================================================
	static double rcut;
	static double trunc_ao;
	static double trunc_wlmr;
	static double trunc_wan;

//==========================================================
// MEMBER FUNCTIONS : 
// NAME : wlmr_dir( if reconstruction , this is wlmr adress)
// NAME : wf_type ( type of wannier functions)
// NAME : build_wf
// NAME : imp_pao
// NAME : b_out_wf
// NAME : b_fftwan
// NAME : b_plot_build
// NAME : b_plot_atomic
//==========================================================
	static std::string wlmr_dir;
	static std::string wf_type;
	static bool build_wf;
	static int imp_pao;
	static bool b_out_wf;
	static bool b_fftwan;//add 2008-07-20
	static bool b_plot_build;//add 2008-06-04
	static bool b_plot_atomic;//add 2008-06-04

//==========================================================
// MEMBER FUNCTIONS : 
// NAME : trial ( trial wave functions , "atomic" or "gauss")
// NAME : bs(parameters for gauss orbit)
// NAME : bp
// NAME : px
// NAME : g1
// NAME : g2
//==========================================================
	static std::string trial;//"atomic" or "gauss" 
	static double bs;//parameters for gauss orbit
	static double bp;
	static double px;
	static double g1;
	static double g2;

	static int bloch_begin;
	static int bloch_end;

	static double fermi_t;
	static double clm2_lowest;

	static int sph_proj;//"1" spherical project,"2": first minus atomic orbitals
	static bool sph_type;//0:Rewrite 1:Skip

	static bool b_recon;//"1" reconstruction of wannier function
	static bool b_mix_wf;// add 2008-06-15
	static double mix_wf;//add 2008-06-13
	static bool recon_wanq;

	static bool speed_mode;

	static bool b_near_atom;
	static double range0;
	static double range1;

	static int L_start;
	static int L_end;
	static int atom_start;
	static int atom_end;

	static bool plot_wanq;//add 2008-01-26
	static std::string plot_option;//(110),[110] etc.
	static int n_unitcell;//number of unitcell to plot
	static bool out_all;
	static bool out_chg;
	static std::string charge_type;
	static bool cal_bands; //for wan  	   wan basis + wan charge
	static bool cal_bands2;//for semi-wan ;pw basis + wan charge add 2008-4-11
	static bool cal_dos;
	static int mesh;
	static double dr;

	static bool no_center;
	static int sum_lm;
	static bool compare_atomic;

	static int out_spillage; // output spillage file.
	static std::string spillage_outdir;
	
	static void Init(const std::string &fn);
	static void Print(const std::string &fn);
private:

	template <class T>
	static void read_value(std::ifstream &ifs, T &var)
	{
		ifs >> var;
		ifs.ignore(75, '\n');
		return;
	}

	static void Read(const std::string &fn);
	static void Default();
	static void Check(void);
#ifdef __MPI
	static void Bcast();
#endif
};

#endif
