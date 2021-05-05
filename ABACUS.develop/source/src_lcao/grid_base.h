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

	Matrix3 latvec;
	Matrix3 latvec0;

	Vector3<double> a1;
	Vector3<double> a2;
	Vector3<double> a3;

	double a1_len;
	double a2_len;
	double a3_len;

	Vector3<double> da_d;

	double da1;
	double da2;
	double da3;

	int nx;
	int ny;
	int nz;
	int nxyz;

	Vector3<double> *cartesian;

	double lat0;

	int test;

	double *Rcut_max;

	Vector3<double> *Rcut_max_direct;

	int grid_number;
	int grid_number_last;
	
	double *norm1;
	double *norm2;
	
	double Rcut1;
	double Rcut2;
	
	Vector3<double> *dR1;
	Vector3<double> *dR2;

	const Numerical_Orbital_Lm* pointer1;
	const Numerical_Orbital_Lm* pointer2;

	int iw1_all;
	int iw2_all;
	int index1;
	int index2;

	Vector3<int> edge_min;
	Vector3<int> edge_max;

	enum cal_type{ cal_charge, cal_local } job;
	
	double** yy1;
	double** yy2;

	int n1; // (lmax+1)^2
	int n2;
	int n1_last;
	int n2_last;
	int lmax1;
	int lmax2;

};

#endif
