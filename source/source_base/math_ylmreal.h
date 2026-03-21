#ifndef MATH_YLMREAL_H
#define MATH_YLMREAL_H

#include "vector3.h"
#include "matrix.h"
#include "realarray.h"
#include <vector>
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
	// ------------------------------------------------------------------
    // 以下为 Ylm_Real 函数重构拆分出的辅助函数
    // ------------------------------------------------------------------
    
    // 根据 lmax2 计算 lmax，如果输入不合法返回 -1
    static int get_lmax(const int lmax2);

    // 计算极角 (cost = cos(theta)) 和方位角 (phi)
    static void compute_polar_angles(
        int ng, 
        const ModuleBase::Vector3<double>* g, 
        std::vector<double>& cost, 
        std::vector<double>& phi
    );

    // 递推计算连带勒让德多项式
    static void compute_Legendre_Polynomials(
        int lmax, 
        int ng, 
        const std::vector<double>& cost, 
        ModuleBase::realArray& p
    );

    // 组装最终的实球谐函数矩阵
    static void assemble_ylm(
        int lmax, 
        int ng, 
        const std::vector<double>& phi, 
        const ModuleBase::realArray& p, 
        matrix& ylm
		);

};

}

#endif
