#ifndef YLM_H
#define YLM_H

#include "../module_base/vector3.h"
#include <vector>

class Ylm
{
	public:
	Ylm(){};
	~Ylm(){};

	static int nlm;

	// (1) for check
	static void get_ylm_real( 
			const int &Lmax , 
			const Vector3<double> &vec, 
			double ylmr[]);

	// (2) for check
	static void get_ylm_real( 
			const int &Lmax , 
			const Vector3<double> &vec, 
			double ylmr[],
			double dylmdr[][3]);

	// (3) not used anymore.
	static void rlylm(
			const int& Lmax, // max momentum of l +1
			const double& x,
			const double& y,
			const double& z,
			double rly[]);

	// (4) not used anymore.
	static void rlylm(
			const int& Lmax, // max momentum of l +1
			const double& x,
			const double& y,
			const double& z,
			double rly[],
			double grly[][3]);
	
	// (5) used in grid integration.
	static void sph_harm(
			const int& Lmax,
			const double& xdr,
			const double& ydr,
			const double& zdr,
			std::vector<double> &rly);
	
	// (6) used in getting overlap.
		// Peize Lin change rly 2016-08-26
	static void rl_sph_harm(
			const int& Lmax,
			const double& x,
			const double& y,
			const double& z,
			std::vector<double>& rly);
	
	// (6) used in getting derivative of overlap.
		// Peize Lin change rly, grly 2016-08-26
	static void grad_rl_sph_harm(
			const int& Lmax,
			const double& x,
			const double& y,
			const double& z,
			std::vector<double>& rly,
			std::vector<std::vector<double>>& grly);
			
	static void set_coefficients ();
	static std::vector<double> ylmcoef;
	
	static void test();
	static void test1();
	static void test2();

	static void ZEROS(double u[], const int& n);

	private:
	static long double Fact(const int n);
	static int Semi_Fact(const int n); 
	static double sgn(const double x);
};

#endif
