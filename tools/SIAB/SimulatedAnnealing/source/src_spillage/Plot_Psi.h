#ifndef PLOT_PSI_H
#define PLOT_PSI_H

#include "common.h"

class Plot_Psi
{
	public:
	Plot_Psi();
	~Plot_Psi();

	void radial_wave_function(void);

	private:

	// (1) about uniform mesh.
	double dr;
	int meshU; // U stands for uniform.
	int meshL;// L stands for log.
	double* rU;// uniform mesh.
	double* rL;// log mesh.
	int meshK;
	double dk;//classic bug, change dk to 'int' type. then see .... 

	// (2) about logmesh.
	double xmin;
	double zed;
	double dx;
	double xmax;

	// (3) about the wave functions.
	matrix psiU; // sum c*j(l,x)
	matrix psiL; //
	matrix psiUK;

	int total_nchi;
	double rcut;
	int lmax;
	int ne;
	int test;

	void allocate(void);

	// mohan add smooth, sigam and rcut 2009-08-28
	void get_psi( const int &mesh, const double *r,const int &l, 
		double *psi1, double *c4, const double *eigen1, bool logmesh)const;
	void orthogonalization( const int &mesh, double *r, const matrix &psi );
	
	void print_orbitals_plot(void)const;
	void print_orbitalsU_used(void)const;
	void print_orbitalsL_used(void)const;

	// mohan add 2010-04-13
	void establish_ecut( ofstream &ofs, const int &ik, double *psi, const int &L);
};

#endif
