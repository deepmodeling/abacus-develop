#ifndef MATH_YLMREAL_H
#define MATH_YLMREAL_H

#include "vector3.h"
#include "../module_base/matrix.h"

using namespace std;

class YlmReal 
{
	public:

    YlmReal();
    ~YlmReal();

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

	private:

    static long double Fact(const int n);
    static int Semi_Fact(const int n);

};

#endif
