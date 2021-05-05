#ifndef INTEGRAL_H
#define INTEGRAL_H

class Integral
{
	public:
	Integral();
	~Integral();

	static double Gauss_Legendre
	(
		const double& a,
		const double& b,
	 	const double* F,
	 	const double* Rad,
	 	const int& Msh
	);

	private:

	static void gauleg();

	static int n_root ;
	static double* gauleg_w;
	static double* gauleg_x;
	static bool calc_wx;

};
#endif
