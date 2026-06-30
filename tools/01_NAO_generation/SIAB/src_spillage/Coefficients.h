#ifndef COEFFICIENTS_H
#define COEFFICIENTS_H

#include "common.h"
class Coefficients
{
	public:
	//============================
	// C4(type, L, N, ie)
	//============================
	realArray C4;
	realArray C4_old;
	realArray C4_accumulate;
	//IntArray acc_count;
	IntArray accept_number;
	realArray accept_rate;

	// constructor, available for SZ, DZ, DZP...
	Coefficients();
	~Coefficients();
	
	void init( const int &ntype, const int &lmax, const int &namx, 
			const double &ecut, const double &ecut_jlq, 
			const double &rcut, const double &tolerence);
	void trial_c4( const int &t, const int &l, const int &n, const int &ie );

	void update_c4( const int &t, const int &l, const int &n, const int &ie );
	void go_back_c4( const int &t, const int &l, const int &n, const int &ie );
	void copy_c4_to_old(void); // for normalize
	void kinetic_energy( const int &t, const int &l, const int &n, const int &ie);// mohan add 2009-07-23
	void accumulate_num_zero(void);
	void accumulating_C4( const int &t, const int &l, const int &n, const int &ie );

	int ntype;
	int lmax;
	int nmax;
	int enumber;
	double rcut;
	double tolerence;
	int fix;

	string *output_c4_name;
	static bool write_c4_flag;

	private:
	void allocate_C4();
	double ecut;
	double ecut_jlq; // mohan add 2009-07-23
	int enumber_jlq; // mohan add 2009-07-23
	int test;
	int accumulate_num;
};

#endif
