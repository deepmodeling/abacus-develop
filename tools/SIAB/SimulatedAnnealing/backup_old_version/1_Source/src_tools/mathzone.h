#ifndef MATHZONE_H
#define MATHZONE_H

#include <cstdlib>
#include <cmath>

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cassert>

using namespace std;

class Mathzone 
{
public:

	Mathzone();
	~Mathzone();

	//========================================================
	// Spherical Bessel
	//========================================================
	static void Spherical_Bessel
	(
		const int &msh,	//number of grid points
		const double *r,//radial grid
		const double &q,	//
		const int &l,	//angular momentum
		double *jl		//jl(1:msh) = j_l(q*r(i)),spherical bessel function
	);
	static void Spherical_Bessel_Roots
	(
		const int &num,
		const int &l,
		const double &epsilon,
		double* eigenvalue,
		const double &rcut
	);
	static void Simpson_Integral(
		const int mesh,
		const double *func,
		const double *rab,
		double &asum
	);
	static double PI;

private:
	static double Spherical_Bessel_7(const int n, const double &x);
	static double BESSJY(double x, double xnu, double *rj, double *ry, double *rjp, double *ryp);
	static void BESCHB(double x, double *gam1, double *gam2, double *gampl, double *gammi);
	static double CHEBEV(double a, double b, double c[], int m, double x);
	static int IMAX(int a, int b);

	
};

#endif
