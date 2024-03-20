// This file contains realization of LDA correlation functionals
// Spin unpolarized ones:
//  1. pw : Perdew-Wang LDA correlation
//  2. pz : Perdew-Zunger LDA correlation
//  3. lyp : Lee-Yang-Parr correlation
//  4. vwn : Vosko-Wilk-Nusair LDA correlation
//  5. wigner : Wigner
//  6. hl : Hedin-Lunqvist
//  7. gl : Gunnarson-Lunqvist
// And some of their spin polarized counterparts:
//  1. pw_spin
//  2. pz_spin, which calls pz_polarized

#include "xc_functional.h"

void XC_Functional::pw(const double &rs, const int &iflag, double &ec, double &vc)
{
    const double a = 0.0310910;
    const double b1 = 7.59570;
    const double b2 = 3.58760;
    const double c0 = a;
    const double c1 = 0.0466440;
    const double c2 = 0.006640;
    const double c3 = 0.010430;
    const double d0 = 0.43350;
    const double d1 = 1.44080;
    double lnrs, rs12, rs32, rs2, om, dom, olog;
    double a1[2] = { 0.213700, 0.0264810 };
	double b3[2] = { 1.63820, -0.466470 };
	double b4[2] = { 0.492940, 0.133540 };

	// high- and low-density formulae implemented but not used in PW case
	// (reason: inconsistencies in PBE/PW91 functionals)
	
   	if (rs < 1 && iflag == 1)
    {
        // high density formula
        lnrs = log(rs);
        ec = c0 * lnrs - c1 + c2 * rs * lnrs - c3 * rs;
        vc = c0 * lnrs - (c1 + c0 / 3.0) + 2.0 / 3.0 * c2 * rs *
             lnrs - (2.0 * c3 + c2) / 3.0 * rs;
    }
    else if (rs > 100.0 && iflag == 1)
    {
        // low density formula
        ec = - d0 / rs + d1 / pow(rs, 1.50);
        vc = - 4.0 / 3.0 * d0 / rs + 1.50 * d1 / pow(rs, 1.50);
    }
    else
    {
        // interpolation formula
        rs12 = sqrt(rs);
        rs32 = rs * rs12;
        rs2 = rs * rs;
        om = 2.0 * a * (b1 * rs12 + b2 * rs + b3 [iflag] * rs32 + b4 [
                            iflag] * rs2);
        dom = 2.0 * a * (0.50 * b1 * rs12 + b2 * rs + 1.50 * b3 [
                             iflag] * rs32 + 2.0 * b4 [iflag] * rs2);
        olog = log(1.0 + 1.0 / om);
        ec = - 2.0 * a * (1.0 + a1 [iflag] * rs) * olog;
        vc = - 2.0 * a * (1.0 + 2.0 / 3.0 * a1 [iflag] * rs)
             * olog - 2.0 / 3.0 * a * (1.0 + a1 [iflag] * rs) * dom /
             (om * (om + 1.0));
    }
	return;
}

//LDA parameterization form Monte Carlo data
//iflag=0: J.P. Perdew and A. Zunger, PRB 23, 5048 (1981)
//iflag=1: G. Ortiz and P. Ballone, PRB 50, 1391 (1994)
void XC_Functional::pz(const double &rs, const int &iflag, double &ec, double &vc)
{
//	ModuleBase::TITLE("XC_Functional","pz");
	const double a[2] = {0.0311, 0.031091};
    const double b[2] = { -0.048, -0.046644 };
    const double c[2] = { 0.0020, 0.00419 };
    const double d[2] = { -0.0116, -0.00983 };
    const double gc[2] = { -0.1423, -0.103756};
    const double b1[2] = { 1.0529, 0.56371 };
    const double b2[2] = { 0.3334, 0.27358 };

    if (rs < 1.0)
    {
        // high density formula
        const double lnrs = log(rs);
        ec = a [iflag] * lnrs + b [iflag] + c [iflag] * rs * lnrs + d [iflag] * rs;
        vc = a [iflag] * lnrs + (b [iflag] - a [iflag] / 3.0) + 2.0 /
             3.0 * c [iflag] * rs * lnrs + (2.0 * d [iflag] - c [iflag])
             / 3.0 * rs;
    }
    else
    {
        // interpolation formula
        const double rs12 = sqrt(rs);
        const double ox = 1.0 + b1 [iflag] * rs12 + b2 [iflag] * rs;
        const double dox = 1.0 + 7.0 / 6.0 * b1 [iflag] * rs12 + 4.0 / 3.0 *
              b2 [iflag] * rs;
        ec = gc [iflag] / ox;
        vc = ec * dox / ox;
    }
    return;
}

// C. Lee, W. Yang, and R.G. Parr, PRB 37, 785 (1988)
// LDA part only
void XC_Functional::lyp(const double &rs, double &ec, double &vc)
{
    const double a = 0.04918e0;
    const double b = 0.1320 * 2.87123400018819108e0;
    // pi43 = (4pi/3)^(1/3)
    const double pi43 = 1.61199195401647e0;
    const double c = 0.2533e0 * pi43;
    const double d = 0.349e0 * pi43;
    double ecrs, ox;

    ecrs = b * exp(- c * rs);
    ox = 1.0 / (1.0 + d * rs);
    ec = - a * ox * (1.0 + ecrs);
    vc = ec - rs / 3.0 * a * ox * (d * ox + ecrs * (d * ox + c));

    return;
}

// S.H. Vosko, L. Wilk, and M. Nusair, Can. J. Phys. 58, 1200 (1980)
void XC_Functional::vwn(const double &rs, double &ec, double &vc)
{
    const double a = 0.0310907;
    const double b = 3.72744;
    const double c = 12.9352;
    const double x0 = -0.10498;
    double q, f1, f2, f3, rs12, fx, qx, tx, tt;

    q = sqrt(4.0 * c - b * b);
    f1 = 2.0 * b / q;
    f2 = b * x0 / (x0 * x0 + b * x0 + c);
    f3 = 2.0 * (2.0 * x0 + b) / q;
    rs12 = sqrt(rs);
    fx = rs + b * rs12 + c;
    qx = atan(q / (2.0 * rs12 + b));

    double x;
    x = (rs12 - x0);
    ec = a * (log(rs / fx) + f1 * qx - f2 * (log((x * x) /
              fx) + f3 * qx));
    tx = 2.0 * rs12 + b;
    tt = tx * tx + q * q;

    vc = ec - rs12 * a / 6.0 * (2.0 / rs12 - tx / fx - 4.0 * b /
	tt - f2 * (2.0 / (rs12 - x0) - tx / fx - 4.0 * (2.0 * x0 + b)/ tt));

    return;
}

void XC_Functional::wigner( const double &rs, double &ec, double &vc)
{
    const double pi34 = 0.6203504908994e0;
    // pi34=(3/4pi)^(1/3), rho13=rho^(1/3)

	double rho13 = pi34 / rs;
    double x = 1.0 + 12.570 * rho13;
    vc = - rho13 * ((0.9436560 + 8.89630 * rho13) / (x * x));
    ec = - 0.7380 * rho13 * (0.9590 / (1.0 + 12.570 * rho13));
    return;
}

// L. Hedin and  B.I. Lundqvist,  J. Phys. C 4, 2064 (1971)
void XC_Functional::hl( const double &rs, double &ec, double &vc)
{
    double a, x;
    a = log(1.00 + 21.0 / rs);
    x = rs / 21.00;
    ec = a + (pow(x, 3) * a - x * x) + x / 2.0 - 1.00 / 3.00;
    ec = - 0.02250 * ec;
    vc = - 0.02250 * a;
    return;
}

// O. Gunnarsson and B. I. Lundqvist, PRB 13, 4274 (1976)
void XC_Functional::gl( const double &rs, double &ec, double &vc)
{
    const double c = 0.0333;
    const double r = 11.4;
    // c=0.0203, r=15.9 for the paramagnetic case
    double x = rs / r;
    vc = - c * log(1.0 + 1.0 / x);
    ec = - c * ((1.0 + pow(x, 3)) * log(1.0 + 1.0 / x) - 1.00 /
                3.00 + x * (0.50 - x));
	return;	
}

// J.P. Perdew and Y. Wang, PRB 45, 13244 (1992)
void XC_Functional::pw_spin( const double &rs, const double &zeta, 
		double &ec, double &vcup, double &vcdw)
{
	const double a = 0.0310910;
    const double a1 = 0.213700;
    const double b1 = 7.59570;
    const double b2 = 3.58760;
    const double b3 = 1.63820;
    const double b4 = 0.492940;
    //const double c0 = a;
    //const double c1 = 0.0466440;
    //const double c2 = 0.006640;
    //const double c3 = 0.010430;
    //const double d0 = 0.43350;
    //const double d1 = 1.44080;
	// polarized parameters
    const double ap = 0.0155450;
    const double a1p = 0.205480;
    const double b1p = 14.11890;
    const double b2p = 6.19770;
    const double b3p = 3.36620;
    const double b4p = 0.625170;
    //const double c0p = ap;
    //const double c1p = 0.0255990;
    //const double c2p = 0.003190;
    //const double c3p = 0.003840;
    //const double d0p = 0.32870;
    //const double d1p = 1.76970;
	// antiferro
    const double aa = 0.0168870;
    const double a1a = 0.111250;
    const double b1a = 10.3570;
    const double b2a = 3.62310;
    const double b3a = 0.880260;
    const double b4a = 0.496710;
    //const double c0a = aa;
    //const double c1a = 0.0354750;
    //const double c2a = 0.001880;
    //const double c3a = 0.005210;
    //const double d0a = 0.22400;
    //const double d1a = 0.39690;

	const double fz0 = 1.7099210;
	double rs12, rs32, rs2, zeta2, zeta3, zeta4, fz, dfz;
	double om, dom, olog, epwc, vpwc;
    double omp, domp, ologp, epwcp, vpwcp;
    double oma, doma, ologa, alpha, vpwca;

	//     if(rs.lt.0.5d0) then
    // high density formula (not implemented)
    //
    //     else if(rs.gt.100.d0) then
    // low density formula  (not implemented)
    //
    //     else
    // interpolation formula
    zeta2 = zeta * zeta;
    zeta3 = zeta2 * zeta;
    zeta4 = zeta3 * zeta;
    rs12 = sqrt(rs);
    rs32 = rs * rs12;
    rs2 = rs * rs;
    // unpolarised
    om = 2.0 * a * (b1 * rs12 + b2 * rs + b3 * rs32 + b4 * rs2);
    dom = 2.0 * a * (0.50 * b1 * rs12 + b2 * rs + 1.50 * b3 * rs32
                     + 2.0 * b4 * rs2);
    olog = log(1.0 + 1.00 / om);
    epwc = - 2.0 * a * (1.0 + a1 * rs) * olog;
    vpwc = - 2.0 * a * (1.0 + 2.0 / 3.0 * a1 * rs) * olog - 2.0 /
           3.0 * a * (1.0 + a1 * rs) * dom / (om * (om + 1.0));
    // polarized
    omp = 2.0 * ap * (b1p * rs12 + b2p * rs + b3p * rs32 + b4p * rs2);
    domp = 2.0 * ap * (0.50 * b1p * rs12 + b2p * rs + 1.50 * b3p *
                       rs32 + 2.0 * b4p * rs2);
    ologp = log(1.0 + 1.00 / omp);
    epwcp = - 2.0 * ap * (1.0 + a1p * rs) * ologp;
    vpwcp = - 2.0 * ap * (1.0 + 2.0 / 3.0 * a1p * rs) * ologp -
            2.0 / 3.0 * ap * (1.0 + a1p * rs) * domp / (omp * (omp + 1.0));
    // antiferro
    oma = 2.0 * aa * (b1a * rs12 + b2a * rs + b3a * rs32 + b4a * rs2);
    doma = 2.0 * aa * (0.50 * b1a * rs12 + b2a * rs + 1.50 * b3a *
                       rs32 + 2.0 * b4a * rs2);
    ologa = log(1.0 + 1.00 / oma);
    alpha = 2.0 * aa * (1.0 + a1a * rs) * ologa;
    vpwca = + 2.0 * aa * (1.0 + 2.0 / 3.0 * a1a * rs) * ologa +
            2.0 / 3.0 * aa * (1.0 + a1a * rs) * doma / (oma * (oma + 1.0));

	fz = (pow((1.0 + zeta) , (4.0 / 3.0)) + pow((1.0 - zeta) , (4.0 /
            3.0)) - 2.0) / (pow(2.0, (4.0 / 3.0)) - 2.0);
    dfz = (pow((1.0 + zeta) , (1.0 / 3.0)) - pow((1.0 - zeta) , (1.0 /
            3.0))) * 4.0 / (3.0 * (pow(2.0 , (4.0 / 3.0)) - 2.0));

    ec = epwc + alpha * fz * (1.0 - zeta4) / fz0 + (epwcp - epwc)
         * fz * zeta4;

    vcup = vpwc + vpwca * fz * (1.0 - zeta4) / fz0 
	+ (vpwcp - vpwc) * fz * zeta4 
	+ (alpha / fz0 * (dfz * (1.0 - zeta4) - 4.0 * fz * zeta3) 
	+ (epwcp - epwc) * (dfz * zeta4 + 4.0 * fz * zeta3)) * (1.0 - zeta);

    vcdw = vpwc + vpwca * fz * (1.0 - zeta4) / fz0 
	+ (vpwcp - vpwc) * fz * zeta4 
	- (alpha / fz0 * (dfz * (1.0 - zeta4) - 4.0 * fz * zeta3) 
	+ (epwcp - epwc) * (dfz * zeta4 + 4.0 * fz * zeta3)) * (1.0 + zeta);
	return;
}

// J.P. Perdew and Y. Wang, PRB 45, 13244 (1992)
void XC_Functional::pz_spin( const double &rs, const double &zeta, 
	double &ec, double &vcup, double &vcdw)
{
//	ModuleBase::TITLE("XC_Functional","pz_spin");
	double ecu, vcu, ecp, vcp ;
    static const double p43 = 4.00 / 3.0;
    static const double third = 1.0 / 3.0;

    // unpolarized part (Perdew-Zunger formula)
    XC_Functional::pz(rs, 0, ecu, vcu);

    // polarization contribution
	XC_Functional::pz_polarized(rs, ecp, vcp);

    double fz = (pow((1.00 + zeta) , p43) + pow((1.0 - zeta) , p43) - 2.0) /
         (pow(2.0 , p43) - 2.0);
    double dfz = p43 * (pow((1.00 + zeta) , third) - pow((1.0 - zeta) , third))
          / (pow(2.0, p43) - 2.0);

	//correlation energy
    ec = ecu + fz * (ecp - ecu);
	
	// spin up correlation potential
    vcup = vcu + fz * (vcp - vcu) + (ecp - ecu) * dfz * (1.0 - zeta);
	// spin down correlation potential
    vcdw = vcu + fz * (vcp - vcu) + (ecp - ecu) * dfz * (- 1.0 - zeta);

    return;
}

// J.P. Perdew and A. Zunger, PRB 23, 5048 (1981)
void XC_Functional::pz_polarized( const double &rs, double &ec, double &vc)
{
//	ModuleBase::TITLE("XC_Functional","pz_polarized");
	static const double a = 0.015550;
	static const double b = -0.02690;
	static const double c = 0.00070;
	static const double d = -0.00480;
	static const double gc = -0.08430;
	static const double b1 = 1.39810;
	static const double b2 = 0.26110;
	
	if (rs < 1.00)
	{
		double lnrs = log(rs);
		ec = a * lnrs + b + c * rs * lnrs + d * rs;
		vc = a * lnrs + (b - a / 3.0) + 2.0 / 3.0 * c * rs * lnrs +
			(2.0 * d - c) / 3.0 * rs;
	}
	else
	{
		// interpolation formula
		double rs12 = sqrt(rs);
		double ox = 1.0 + b1 * rs12 + b2 * rs;
		double dox = 1.0 + 7.0 / 6.0 * b1 * rs12 + 4.0 / 3.0 * b2 * rs;
		ec = gc / ox;
		vc = ec * dox / ox;
	}
	return;
}