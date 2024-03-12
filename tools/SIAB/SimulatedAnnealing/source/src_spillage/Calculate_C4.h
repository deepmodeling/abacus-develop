#ifndef CALCULATE_C4_H
#define CALCULATE_C4_H

#include "common.h"

class Calculate_C4
{
	public:

	Calculate_C4();
	~Calculate_C4();

	void init( ifstream &ifs, const double &ecut, const double &rcut, const double &tolerence, const int &lmax );
	void Test( ifstream &ifs, const double &ecut, const double &tolerence);
	void plot( ifstream &ifs );

	// Peize Lin add accept_rate 2015-12-24
	static void norm_ic( realArray &C4, realArray &accept_rate, const int &ntype, const int &lmax, 
		const int &nmax, const int &enumber, const double &tolerence,
		const double &rcut, const double &dr);
	
	static void Find_Eigenvalues(const int &Enumber, matrix &en,
        const int &lmax, const double &rcut, const double &tolerence);

	private:

	void Test_Orthogonal(const double &rcut, const double &dr, 
			const int &L, const matrix &en, const int &Enumber, const int &test_ei);

	void key_solver(const string &name, const int &L, 
			const matrix &en, const int &Enumber, double *C4);
	
	static int test;

	
	
};

#endif
