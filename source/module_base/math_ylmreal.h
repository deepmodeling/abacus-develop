#ifndef MATH_YLMREAL_H
#define MATH_YLMREAL_H

#include "vector3.h"
#include "matrix.h"

namespace ModuleBase
{

class YlmReal 
{
	public:

    YlmReal();
    ~YlmReal();

	/**
	 * @brief spherical harmonic function (real form) an array of vectors
	 * 
	 * @param lmax2 [in] lmax2 = (lmax + 1)^2 ; lmax = angular quantum number
	 * @param ng [in] the number of vectors
	 * @param g [in] an array of vectors
	 * @param ylm [out] Ylm; column index represent vector, row index represent Y00, Y10, Y11, Y1-1, Y20,Y21,Y2-1,Y22.Y2-2,...;
	 */
    static void Ylm_Real
    (
        const int lmax2, 		
        const int ng,				
        const ModuleBase::Vector3<double> *g, 
        matrix &ylm 	
    );

    /**
	 * @brief spherical harmonic function (real form) an array
	 *
	 * @param lmax2 [in] lmax2 = (lmax + 1)^2 ; lmax = angular quantum number
	 * @param ng [in] the number of vectors
	 * @param g [in] an array of vectors
	 * @param ylm [out] Ylm; column index represent vector, row index represent Y00, Y10, Y11, Y1-1, Y20,Y21,Y2-1,Y22.Y2-2,...;
	 */
	template <typename FPTYPE, typename Device>
    static void Ylm_Real(Device * ctx, const int lmax2, const int ng, const FPTYPE *g, FPTYPE * ylm);

	/**
	 * @brief gradient of spherical harmonic function (real form) an array of vectors
	 * 
	 * @param lmax2 [in] lmax2 = (lmax + 1)^2 ; lmax = angular quantum number
	 * @param ng [in] the number of vectors
	 * @param g [in] an array of vectors
	 * @param ylm [out] Ylm; column index represent vector, row index represent Y00, Y10, Y11, Y1-1, Y20,Y21,Y2-1,Y22.Y2-2,...;
	 * @param dylmx/dylmy/dylmz [out] \nabla Ylm; column index represent vector, row index represent dY00/dxyz, dY10/dxyz,...;
	 */
    static void grad_Ylm_Real
    (
        const int lmax2, 		
        const int ng,				
        const ModuleBase::Vector3<double> *g, 
		matrix &ylm,
        matrix &dylmx,
		matrix &dylmy,
		matrix &dylmz	
    );
	
	/**
	 * @brief spherical harmonic function (Herglotz generating form) of an array of vectors
	 * 
	 * @param lmax2 [in] lmax2 = (lmax + 1)^2 ; lmax = angular quantum number
	 * @param ng [in] the number of vectors
	 * @param g [in] an array of vectors
	 * @param ylm [out] Ylm; column index represent vector, row index represent Y00, Y10, Y11, Y1-1, Y20,Y21,Y2-1,Y22.Y2-2,...;
	 */
	static void Ylm_Real2
	(
    	const int lmax2, 		
    	const int ng,				
    	const ModuleBase::Vector3<double> *g, 
    	matrix &ylm 		
	);

	/**
	 * @brief spherical harmonic function (Herglotz generating form) of a vector
	 * 
	 * @param lmax [in] maximum angular quantum number
	 * @param x [in] x part of the vector
	 * @param y [in] y part of the vector
	 * @param z [in] z part of the vector
	 * @param rly [in] Ylm, Y00, Y10, Y11, Y1-1, Y20,Y21,Y2-1,Y22.Y2-2,...
	 */
	static void rlylm
	(
    	const int lmax, 	
    	const double& x,				
    	const double& y,
		const double& z, 
    	double* rly 
	);

	private:

    static long double Fact(const int n);
    static int Semi_Fact(const int n);

};

}

#endif
