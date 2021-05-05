#ifndef MATHZONE_H
#define MATHZONE_H

#include "realarray.h"
#include <vector>
#include <map>
#include <cassert>
using namespace std;

class Mathzone
{
public:

    Mathzone();
    ~Mathzone();

    template<class T>
    static T Max3(const T &a,const T &b,const T &c)
    {
        if (a>=b && a>=c) return a;
        else if (b>=a && b>=c) return b;
        else if (c>=a && c>=b) return c;
		else throw runtime_error(TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
    }

    // be careful, this can only be used for plane wave
    // during parallel calculation
    static void norm_pw(complex<double> *u, const int n);

    //========================================================
    // Polynomial_Interpolation
    //========================================================
    static void Polynomial_Interpolation
    (
        const realArray &table,
        const int &dim1,
        const int &dim2,
        realArray &y,
        const int &dim_y,
        const int &table_length,
        const double &table_interval,
        const double &x
    );
    static double Polynomial_Interpolation
    (
        const realArray &table,
        const int &dim1,
        const int &dim2,
        const int &table_length,
        const double &table_interval,
        const double &x				// input value
    );
	static double Polynomial_Interpolation             // pengfei Li 2018-3-23
    (
        const realArray &table,
        const int &dim1,
        const int &dim2,
		const int &dim3,
        const int &table_length,
        const double &table_interval,
        const double &x				// input value
    );
    static double Polynomial_Interpolation
    (
        const double *table,
        const int &table_length,
        const double &table_interval,
        const double &x				// input value
    );
    static double Polynomial_Interpolation_xy
    (
        const double *xpoint,
        const double *ypoint,
        const int table_length,
        const double &x             // input value
    );

    //========================================================
    // Spherical Bessel
    //========================================================
    static void Spherical_Bessel
    (
        const int &msh,	//number of grid points
        const double *r,//radial grid
        const double &q,	//
        const int &l,	//angular momentum
        double *jl	//jl(1:msh) = j_l(q*r(i)),spherical bessel function
    );

	static void Spherical_Bessel
	(           
	    const int &msh, //number of grid points
		const double *r,//radial grid
		const double &q,    //
		const int &l,   //angular momentum
		double *sj,     //jl(1:msh) = j_l(q*r(i)),spherical bessel function
		double *sjp
	);

	
    static void Spherical_Bessel_Roots
    (
        const int &num,
        const int &l,
        const double &epsilon,
        double* eigenvalue,
        const double &rcut
    );
private:

    static double Spherical_Bessel_7(const int n, const double &x);
    static void BESSJY(double x, double xnu, double *rj, double *ry, double *rjp, double *ryp);	// Peize Lin change double to void 2019-05-01
    static void BESCHB(double x, double *gam1, double *gam2, double *gampl, double *gammi);
    static double CHEBEV(double a, double b, double c[], int m, double x);
    static int IMAX(int a, int b);

public:

    static void Ylm_Real
    (
        const int lmax2, 			// lmax2 = (lmax+1)^2
        const int ng,				//
        const Vector3<double> *g, 	// g_cartesian_vec(x,y,z)
        matrix &ylm 				// output
    );
	
	static void Ylm_Real2
	(
    	const int lmax2, 			// lmax2 = (lmax+1)^2
    	const int ng,				//
    	const Vector3<double> *g, 	// g_cartesian_vec(x,y,z)
    	matrix &ylm 				// output
	);

	static void rlylm
	(
    	const int lmax, 	
    	const double& x,				
    	const double& y,
		const double& z, // g_cartesian_vec(x,y,z)
    	double* rly 	 // output
	);

    static long double Fact(const int n);
    static int Semi_Fact(const int n);

	// Peize Lin accelerate 2017-10-02
    static void Simpson_Integral
    (
        const int mesh,
        const double *func,
        const double *rab,
        double &asum
    );
	// Peize Lin accelerate 2017-10-02
	static void Simpson_Integral
	(
		const int mesh,
		const double *func,
		const double dr,
		double &asum
	);

    // Peize Lin add 2016-02-14
    static void Simpson_Integral_0toall
    (
        const int mesh,
        const double *func,
        const double *rab,
        double *asum
    );

    // Peize Lin add 2016-02-14
    static void Simpson_Integral_alltoinf
    (
        const int mesh,
        const double *func,
        const double *rab,
        double *asum
    );     
	
	// Peize Lin add 2016-08-03
	template< typename Type >
	static vector<Type> Pointwise_Product( const vector<Type> &f1, const vector<Type> &f2 )
	{
		assert(f1.size()==f2.size());
		vector<Type> f(f1.size());
		for( int ir=0; ir!=f.size(); ++ir )
			f[ir] = f1[ir] * f2[ir];
		return f;
	}

//==========================================================
// MEMBER FUNCTION :
// NAME : Direct_to_Cartesian
// use lattice vector matrix R
// change the direct vector (dx,dy,dz) to cartesuab vectir
// (cx,cy,cz)
// (dx,dy,dz) = (cx,cy,cz) * R
//
// NAME : Cartesian_to_Direct
// the same as above
// (cx,cy,cz) = (dx,dy,dz) * R^(-1)
//==========================================================
    static inline void Direct_to_Cartesian
    (
        const double &dx,const double &dy,const double &dz,
        const double &R11,const double &R12,const double &R13,
        const double &R21,const double &R22,const double &R23,
        const double &R31,const double &R32,const double &R33,
        double &cx,double &cy,double &cz)
    {
        static Matrix3 lattice_vector;
        static Vector3<double> direct_vec, cartesian_vec;
        lattice_vector.e11 = R11;
        lattice_vector.e12 = R12;
        lattice_vector.e13 = R13;
        lattice_vector.e21 = R21;
        lattice_vector.e22 = R22;
        lattice_vector.e23 = R23;
        lattice_vector.e31 = R31;
        lattice_vector.e32 = R32;
        lattice_vector.e33 = R33;

        direct_vec.x = dx;
        direct_vec.y = dy;
        direct_vec.z = dz;

        cartesian_vec = direct_vec * lattice_vector;
        cx = cartesian_vec.x;
        cy = cartesian_vec.y;
        cz = cartesian_vec.z;
        return;
    }

    static inline void Cartesian_to_Direct
    (
        const double &cx,const double &cy,const double &cz,
        const double &R11,const double &R12,const double &R13,
        const double &R21,const double &R22,const double &R23,
        const double &R31,const double &R32,const double &R33,
        double &dx,double &dy,double &dz)
    {
        static Matrix3 lattice_vector, inv_lat;
        lattice_vector.e11 = R11;
        lattice_vector.e12 = R12;
        lattice_vector.e13 = R13;
        lattice_vector.e21 = R21;
        lattice_vector.e22 = R22;
        lattice_vector.e23 = R23;
        lattice_vector.e31 = R31;
        lattice_vector.e32 = R32;
        lattice_vector.e33 = R33;

        inv_lat = lattice_vector.Inverse();

        static Vector3<double> direct_vec, cartesian_vec;
        cartesian_vec.x = cx;
        cartesian_vec.y = cy;
        cartesian_vec.z = cz;
        direct_vec = cartesian_vec * inv_lat;
        dx = direct_vec.x;
        dy = direct_vec.y;
        dz = direct_vec.z;
        return;
    }


    static void To_Polar_Coordinate
    (
        const double &x_cartesian,
        const double &y_cartesian,
        const double &z_cartesian,
        double &r,
        double &theta,
        double &phi
    );
	
	
	// coeff1 * x1 + (1-coeff1) * x2
	// Peize Lin add 2017-08-09
	template< typename T, typename T_coeff >
	static T Linear_Mixing( const T & x1, const T & x2, const T_coeff & coeff1 )
	{
		return coeff1 * x1 + (1-coeff1) * x2;
	}
	template< typename T, typename T_coeff >
	static vector<T> Linear_Mixing( const vector<T> & x1, const vector<T> & x2, const T_coeff & coeff1 )
	{
		assert(x1.size()==x2.size());
		vector<T> x;
		for( size_t i=0; i!=x1.size(); ++i )
			x.push_back( Linear_Mixing( x1[i], x2[i], coeff1 ) );
		return x;
	}
	template< typename T1, typename T2, typename T_coeff >
	static map<T1,T2> Linear_Mixing( const map<T1,T2> & x1, const map<T1,T2> & x2, const T_coeff & coeff1 )
	{
		map<T1,T2> x;
		for( const auto & x1i : x1 )
			x.insert( make_pair( x1i.first, Linear_Mixing( x1i.second, x2.at(x1i.first), coeff1 ) ) );
		return x;
	}	

};

#endif
