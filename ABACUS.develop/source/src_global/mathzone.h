#ifndef MATHZONE_H
#define MATHZONE_H

#include "realarray.h"
#include "vector3.h"
#include "../module_base/matrix3.h"
#include "global_function.h"
#include <vector>
#include <map>
#include <cassert>
#include <complex>
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

public:



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
