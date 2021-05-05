#ifndef SPH_BESSEL_H
#define SPH_BESSEL_H

#include "realarray.h"
using namespace std;

class Sph_Bessel
{
public:

    Sph_Bessel();
    ~Sph_Bessel();
	
	void jlx
    (
        const int &msh,	//number of grid points
        const double *r,//radial grid
        const double &q,	//
        const int &l,	//angular momentum
        double *jl	//jl(1:msh) = j_l(q*r(i)),spherical bessel function
    );

private:
	
	double jlx7(const int &l, const double &x1);
	void BESSJY(double x, double xnu, double *rj, double *ry, double *rjp, double *ryp);
    void BESCHB(double x, double *gam1, double *gam2, double *gampl, double *gammi);
    double CHEBEV(double a, double b, double c[], int m, double x);
	int IMAX(int a, int b);

	double eps;
	double fpmin;
	double maxit;
	double xmin;
	double pi;

};

#endif
