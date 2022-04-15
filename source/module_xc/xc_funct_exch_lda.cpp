// This file contains realization of LDA exchange functionals
// Spin unpolarized ones:
//  1. slater: ordinary Slater exchange with alpha=2/3
//  2. slater1: Slater exchange with alpha=1
//  3. slater_rxc : Slater exchange with alpha=2/3 and Relativistic exchange
// And their spin polarized counterparts:
//  1. slater_spin
//  2. slater1_spin
//  3. slater_rxc_spin

#include "xc_functional.h"

//Slater exchange with alpha=2/3
void XC_Functional::slater(const double &rs, double &ex, double &vx)
{
	// f = -9/8*(3/2pi)^(2/3)
	const double f = -0.687247939924714e0;
	const double alpha = 2.00 / 3.00;
	ex = f * alpha / rs;
	vx = 4.0 / 3.0 * f * alpha / rs;
	return;
}

//Slater exchange with alpha=1, corresponding to -1.374/r_s Ry
//used to recover old results
void XC_Functional::slater1(const double &rs, double &ex, double &vx)
{
	const double f = -0.687247939924714e0;
	const double alpha = 1.0;
	ex = f * alpha / rs;
	vx = 4.0 / 3.0 * f * alpha / rs;
	return;
}

// Slater exchange with alpha=2/3 and Relativistic exchange
void XC_Functional::slater_rxc(const double &rs, double &ex, double &vx)
{
	static const double pi = 3.14159265358979;
    const double trd = 1.0 / 3.0;
    //const double ftrd = 4.0 / 3.0;
    //const double tftm = pow(2.0, ftrd) - 2.0;
    const double a0 = pow((4.0 / (9.0 * pi)), trd);
    // X-alpha parameter:
    const double alp = 2 * trd;

    double vxp = -3 * alp / (2 * pi * a0 * rs);
    double exp = 3 * vxp / 4;
    const double beta = 0.014 / rs;
    const double sb = sqrt(1 + beta * beta);
    const double alb = log(beta + sb);
    vxp = vxp * (-0.5 + 1.5 * alb / (beta * sb));
    double x = (beta * sb - alb) / (beta * beta);
    exp = exp * (1.0 - 1.5 * x * x);
    vx = vxp;
    ex = exp;
	return;	
}

// Slater exchange with alpha=2/3, spin-polarized case
void XC_Functional::slater_spin( const double &rho, const double &zeta, 
		double &ex, double &vxup, double &vxdw)
{
    const double f = - 1.107838149573033610;
    const double alpha = 2.00 / 3.00;
    // f = -9/8*(3/pi)^(1/3)
    const double third = 1.0 / 3.0;
    const double p43 = 4.0 / 3.0;

    double rho13 = pow(((1.0 + zeta) * rho) , third);
    double exup = f * alpha * rho13;
    vxup = p43 * f * alpha * rho13;
    rho13 = pow(((1.0 - zeta) * rho) , third);
    double exdw = f * alpha * rho13;
    vxdw = p43 * f * alpha * rho13;
    ex = 0.50 * ((1.0 + zeta) * exup + (1.0 - zeta) * exdw);

    return;
}

// Slater exchange with alpha=2/3, spin-polarized case
void XC_Functional::slater1_spin( const double &rho, const double &zeta, double &ex, double &vxup, double &vxdw)
{
    const double f = - 1.107838149573033610;
	const double alpha = 1.00;
	const double third = 1.0 / 3.0;
	const double p43 = 4.0 / 3.0;
    // f = -9/8*(3/pi)^(1/3)

    double rho13 = pow(((1.0 + zeta) * rho) , third);
    double exup = f * alpha * rho13;
    vxup = p43 * f * alpha * rho13;
    rho13 = pow(((1.0 - zeta) * rho) , third);
    double exdw = f * alpha * rho13;
    vxdw = p43 * f * alpha * rho13;
    ex = 0.50 * ((1.0 + zeta) * exup + (1.0 - zeta) * exdw);

    return;
} // end subroutine slater1_spin

// Slater exchange with alpha=2/3, relativistic exchange case
void XC_Functional::slater_rxc_spin( const double &rho, const double &z, 
		double &ex, double &vxup, double &vxdw)
{
    if (rho <= 0.0)
    {
        ex = vxup = vxdw = 0.0;
        return;
    }

	const double pi = 3.141592653589790;

	const double trd = 1.0 / 3.0;
	const double ftrd = 4.0 / 3.0;
	double tftm = pow(2.0, ftrd) - 2;
	double a0 = pow((4 / (9 * pi)), trd);

	double alp = 2 * trd;

	double fz = (pow((1 + z), ftrd) + pow((1 - z), ftrd) - 2) / tftm;
	double fzp = ftrd * (pow((1 + z), trd) - pow((1 - z), trd)) / tftm;

    double rs = pow((3 / (4 * pi * rho)), trd);
    double vxp = -3 * alp / (2 * pi * a0 * rs);
    double exp = 3 * vxp / 4;

    double beta = 0.014 / rs;
    double sb = sqrt(1 + beta * beta);
    double alb = log(beta + sb);
    vxp = vxp * (-0.5 + 1.5 * alb / (beta * sb));
    exp = exp * (1.0- 1.5*( (beta*sb-alb) / (beta*beta) )
			        * 1.5*( (beta*sb-alb) / (beta*beta) ));

    double x = (beta * sb - alb) / (beta * beta);

    exp = exp * (1.0 - 1.5 * x * x);

    double vxf = pow(2.0, trd) * vxp;

    double exf = pow(2.0, trd) * exp;

    vxup  = vxp + fz * (vxf - vxp) + (1 - z) * fzp * (exf - exp);

    vxdw  = vxp + fz * (vxf - vxp) - (1 + z) * fzp * (exf - exp);

    ex    = exp + fz * (exf - exp);

    return;
}
