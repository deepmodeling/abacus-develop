#ifndef INTEGRAL_G_H
#define INTEGRAL_G_H

class Integral_G
{
	public:

	Integral_G();
	~Integral_G();

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
