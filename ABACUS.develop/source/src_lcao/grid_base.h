//=========================================================
//AUTHOR : mohan
//DATE : 2008-09-16
//=========================================================
#ifndef GRID_BASE_H
#define GRID_BASE_H

#include "../src_pw/tools.h"
#include "ORB_atomic_lm.h"

class Grid_Base
{
public:
	void init( 
		const Matrix3 &latvec_in, 
		const double &lat0_in, 
		const int &nx_in,
		const int &ny_in,
		const int &nz_in,
		const int &nxyz_in);

protected:

	Grid_Base();
	~Grid_Base();
	
//==========================================================
// EXPLAIN : Integral On 3D Real Space For Local Potential
// MEMBER FUNCTION :
//===========================================================
	// define pointer
	const double* vlocal;

	void get_rcut_max(void);
	void edge_grid_points(
		const Vector3<double> &R10,
		const Vector3<double> &R20,
		const Vector3<double> &max_direct_coodinate,
		const Vector3<double> &min_direct_coodinate);
	
	void get_small_box( 
		const Vector3<double> &tau, 
		const int &T,
		Vector3<double> &tau_max_direct,
		Vector3<double> &tau_min_direct);

	int* ijk_index;

	Matrix3 latvec,latvec0;
	Vector3<double> a1, a2, a3;
	double a1_len, a2_len, a3_len;
	Vector3<double> da_d;
	double da1, da2, da3;
	int nx, ny, nz, nxyz;
	Vector3<double> *cartesian;
	double lat0;

	int test;
	double *Rcut_max;
	Vector3<double> *Rcut_max_direct;
	int grid_number, grid_number_last;
	
	double *norm1, *norm2;
	
	double Rcut1, Rcut2;
	
	Vector3<double> *dR1, *dR2;

	const Numerical_Orbital_Lm* pointer1;
	const Numerical_Orbital_Lm* pointer2;

	int iw1_all,iw2_all;
	int index1, index2;
	Vector3<int> edge_min, edge_max;

	enum cal_type{ cal_charge, cal_local } job;
	
	double** yy1;
	double** yy2;
	int n1,n2; // (lmax+1)^2
	int n1_last, n2_last;
	int lmax1, lmax2;

};

#endif
