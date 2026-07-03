#ifndef PSI_SECOND_H
#define PSI_SECOND_H

#include "common.h"

class Psi_Second
{
	public:
	Psi_Second();
	~Psi_Second();

	// be called in MultiZeta.cpp.
	void init(ifstream &ifs, 
		const int &ne_in, 
		const int &lmax_in,
		const double &rcut_in, 
		const int &total_nchi_in, 
		const double &tolerence);

	int calculate_oscillation(const int &L, const int &N, const double *c4, const int &ic);	
	
	int oscillation;
	
	ofstream ofso;
	ofstream ofsk;
	
	bool* below_max_os;
	
	int *os_number;

	void norm_c4( double *c4, const int &L);
	// new method, mohan 2010-04-14
	double get_ecut( double *c4, const int &L);

	matrix eigenvalue;
	matrix jjnorm; //mohan add 2010-04-18
	double*** jle;
	
	private:

	double dr; // about mesh
	int mesh; // about mesh
	double* r; // about mesh
	double* rab; // about mesh
	double rcut; // about mesh
	
	double* eigen1;

	int total_nchi;
	int lmax;
	int ne;
	int test;
	int start;
	int max_os;

	double* psi; // psi
	double* psi_first;//first derivative
	double* psi_second; //second derivative
	double* os_position;
	int count;
};

#endif
