#ifndef YLM_H
#define YLM_H

#include <vector>

#include "vector3.h"

namespace ModuleBase
{

class Ylm
{
	public:
	Ylm(){};
	~Ylm(){};

	static int nlm;


	/**
	 * @brief Get the ylm real object 
	 * 
	 * @param Lmax [in] maximum angular quantum number + 1
	 * @param vec [in] the vector to be calculated	
	 * @param ylmr [out] calculated Ylm, Y00, Y10, Y11, Y1-1, Y20, Y21, Y2-1, Y22, Y2-2...
	 */
	static void get_ylm_real( 
			const int &Lmax , 
			const ModuleBase::Vector3<double> &vec, 
			double ylmr[]);

	/**
	 * @brief Get the ylm real object and the gradient 
	 * 
	 * @param Lmax [in] maximum angular quantum number + l
	 * @param vec [in] the vector to be calculated
	 * @param ylmr [out] calculated Ylm, Y00, Y10, Y11, Y1-1, Y20, Y21, Y2-1, Y22, Y2-2...
	 * @param dylmdr [out] gradient of Ylm, [dY00/dx, dY00/dy, dY00/dz], [dY10/dx, dY10/dy, dY10/dz], [dY11/dx, dY11/dy, dY11/dz],...
	 */
	static void get_ylm_real( 
			const int &Lmax , 
			const ModuleBase::Vector3<double> &vec, 
			double ylmr[],
			double dylmdr[][3]);

	/**
	 * @brief Get the ylm real (solid) object (not used anymore)
	 * 
	 * @param Lmax [in] maximum angular quantum number + l
	 * @param x [in] x
	 * @param y [in] y 
	 * @param z [in] z
	 * @param rly [in] calculated Ylm, Y00, Y10, Y11, Y1-1, Y20, Y21, Y2-1, Y22, Y2-2...
	 */
	static void rlylm(
			const int& Lmax, 
			const double& x,
			const double& y,
			const double& z,
			double rly[]);

	/**
	 * @brief Get the ylm real (solid) object and the gradient (not used anymore)
	 * 
	 * @param Lmax [in] maximum angular quantum number + 1
	 * @param x [in] x
	 * @param y [in] y
	 * @param z [in] z
	 * @param rly [in] calculated Ylm, Y00, Y10, Y11, Y1-1, Y20, Y21, Y2-1, Y22, Y2-2...
	 * @param grly [out] gradient of Ylm, [dY00/dx, dY00/dy, dY00/dz], [dY10/dx, dY10/dy, dY10/dz], [dY11/dx, dY11/dy, dY11/dz],...
	 */
	static void rlylm(
			const int& Lmax, 
			const double& x,
			const double& y,
			const double& z,
			double rly[],
			double grly[][3]);

	/**
	 * @brief Get the ylm real object (used in grid integration)
	 * 
	 * @param Lmax [in] maximum angular quantum number
	 * @param xdr [in] x/r
	 * @param ydr [in] y/r
	 * @param zdr [in] z/r
	 * @param rly [in] calculated Ylm, Y00, Y10, Y11, Y1-1, Y20, Y21, Y2-1, Y22, Y2-2...
	 */
	static void sph_harm(
			const int& Lmax,
			const double& xdr,
			const double& ydr,
			const double& zdr,
			std::vector<double> &rly);
	
	/**
	 * @brief Get the ylm real object (used in getting overlap) 
	 * 
	 * @param Lmax [in] maximum angular quantum number
	 * @param x [in] x/r
	 * @param y [in] y/r
	 * @param z [in] z/r
	 * @param rly [in] calculated Ylm, Y00, Y10, Y11, Y1-1, Y20, Y21, Y2-1, Y22, Y2-2...
	 * @author Peize Lin
	 * @date 2016-08-26
	 */
	static void rl_sph_harm(
			const int& Lmax,
			const double& x,
			const double& y,
			const double& z,
			std::vector<double>& rly);
	
	/**
	 * @brief Get the ylm real object and the gradient (used in getting derivative of overlap)
	 * 
	 * @param Lmax [in] maximum angular quantum number
	 * @param x [in] x/r
	 * @param y [in] y/r
	 * @param z [in] z/r
	 * @param rly [in] calculated Ylm, Y00, Y10, Y11, Y1-1, Y20, Y21, Y2-1, Y22, Y2-2...
	 * @param grly [out] gradient of Ylm, [dY00/dx, dY00/dy, dY00/dz], [dY10/dx, dY10/dy, dY10/dz], [dY11/dx, dY11/dy, dY11/dz],...
	 */
	static void grad_rl_sph_harm(
			const int& Lmax,
			const double& x,
			const double& y,
			const double& z,
			std::vector<double>& rly,
			std::vector<std::vector<double>>& grly);

	/**
	 * @brief Get the hessian of r^l Ylm (used in getting derivative of overlap)
	 * 
	 * @param Lmax [in] maximum angular quantum number
	 * @param x [in] x
	 * @param y [in] y
	 * @param z [in] z
	 * @param hrly [out] hessian of Ylm, [dY00/dx2, dY00/dxy, dY00/dxz, dY00/dyy, dY00/dyz, dY00/dzz] , ...
	 */
	static void hes_rl_sph_harm(
			const int& Lmax,
			const double& x,
			const double& y,
			const double& z,
			std::vector<std::vector<double>>& hrly);

	//calculate the coefficient of Ylm, ylmcoef.		
	static void set_coefficients ();
	static std::vector<double> ylmcoef;
	
	//static void test();
	//static void test1();
	//static void test2();
	
	//set the first n elements of u to be 0.0
	static void ZEROS(double u[], const int& n);

	private:
	static long double Fact(const int n);
	static int Semi_Fact(const int n); 
	static double sgn(const double x);
};

}

#endif
