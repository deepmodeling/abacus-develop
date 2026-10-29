#ifndef ABACUS_XC_MATH_H
#define ABACUS_XC_MATH_H

#include <cmath>

#if defined(__CUDACC__)
#define ABACUS_XC_HOST_DEVICE __host__ __device__
#else
#define ABACUS_XC_HOST_DEVICE
#endif

// Scalar exchange-correlation formulas shared by host and device evaluation.
namespace xc_math {

// Slater exchange with alpha = 2/3 for spin-polarized density.
ABACUS_XC_HOST_DEVICE inline void slater_spin( const double &rho, const double &zeta, 
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


ABACUS_XC_HOST_DEVICE inline void pw_spin( const double &rs, const double &zeta, 
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

ABACUS_XC_HOST_DEVICE inline void pbex(const double &rho, const double &grho, const int &iflag, 
double &sx, double &v1x, double &v2x)
{
    // PBE exchange (without Slater exchange):
    // iflag=0  J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)
    // iflag=1  "revised' PBE: Y. Zhang et al., PRL 80, 890 (1998)

    // input: charge and squared gradient
    // output: energy
    // output: potential
    // (3*pi2*|rho|)^(1/3)
    // |grho|
    // |grho|/(2*kf*|rho|)
    // s^2
    // n*ds/dn
    // n*ds/d(gn)
    // exchange energy LDA part
    // exchange energy gradient part
    
	// numerical coefficients (NB: c2=(3 pi^2)^(1/3) )
    const double third = 1.0 / 3.0;
    const double pi = 3.14159265358979323846;
    const double c1 = 0.750 / pi;
    const double c2 = 3.0936677262801360;
    const double c5 = 4.0 * third;
    // parameters of the functional
    double k[3] = { 0.8040, 1.24500, 0.8040 };
    const double mu[3] = {0.2195149727645171, 0.2195149727645171, 0.12345679012345679} ;//modified by zhengdy, to ensure the same parameters with another dft code.

    const double agrho = sqrt(grho);
    const double kf = c2 * pow(rho, third);
    const double dsg = 0.50 / kf;
    const double s1 = agrho * dsg / rho;
    const double s2 = s1 * s1;
    const double ds = - c5 * s1;

    // Energy
    const double f1 = s2 * mu[iflag] / k [iflag];
    const double f2 = 1.0 + f1;
    const double f3 = k [iflag] / f2;
    const double fx = k [iflag] - f3;
    const double exunif = - c1 * kf;
    sx = exunif * fx;

    // Potential
    const double dxunif = exunif * third;
    const double dfx1 = f2 * f2;
    const double dfx = 2.0 * mu[iflag] * s1 / dfx1;
    
	v1x = sx + dxunif * fx + exunif * dfx * ds;
    v2x = exunif * dfx * dsg / agrho;
    sx = sx * rho;

	return;
}

ABACUS_XC_HOST_DEVICE inline void pbec_spin(double rho, double zeta, double grho, const int &iflag, double &sc,
               double &v1cup, double &v1cdw, double &v2c)
{
    //-----------------------------------------------------------

    // PBE correlation (without LDA part) - spin-polarized
    // J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996).

    // USE kinds
    // implicit none
    // real(kind=DP) :: rho, zeta, grho, sc, v1cup, v1cdw, v2c
    double ga, be[3];//mohan add
    // parameter :
    ga = 0.0310910;
    be[1] = 0.06672455060314922;//zhengdy add 2019-09-12, ensure same parameter with another dft code.
	be[2] = 0.0460000;//mohan add 2012-05-28
    double third, pi34, xkf, xks;
    // parameter :
    third = 1.0 / 3.0;
    pi34 = 0.62035049089940;
    // parameter :
    xkf = 1.9191582926775130;
    xks = 1.1283791670955130;
    // pi34=(3/4pi)^(1/3), xkf=(9 pi/4)^(1/3), xks= sqrt(4/pi)
    double kf, ks, rs, ec, vcup, vcdw, t, expe, af, y, xy, qy,
    s1, h0, ddh0;
    double fz, fz2, fz3, dfz, bfup, bfdw, dh0up, dh0dw, dh0zup, dh0zdw; // fz4

    rs = pi34 / pow(rho, third);
	pw_spin(rs, zeta, ec, vcup, vcdw); //mohan fix bug 2012-05-28
    kf = xkf / rs;
    ks = xks * sqrt(kf);
    fz = 0.50 * (pow((1.0 + zeta) , (2.0 / 3.0)) + pow((1.0 - zeta) , (
                     2.0 / 3.0)));
    fz2 = fz * fz;
    fz3 = fz2 * fz;
    //fz4 = fz3 * fz;
    dfz = (pow((1.0 + zeta) , (- 1.0 / 3.0)) - pow((1.0 - zeta) , (-
            1.0 / 3.0))) / 3.0;
    t = sqrt(grho) / (2.0 * fz * ks * rho);
    expe = exp(- ec / (fz3 * ga));
    af = be[iflag] / ga * (1.0 / (expe - 1.0));
    bfup = expe * (vcup - ec) / fz3;
    bfdw = expe * (vcdw - ec) / fz3;
    y = af * t * t;
    xy = (1.0 + y) / (1.0 + y + y * y);
    qy = y * y * (2.0 + y) / pow( (1.0 + y + y * y),2 );
//    qy *= qy; //bug
    s1 = 1.0 + be[iflag] / ga * t * t * xy;
    h0 = fz3 * ga * log(s1);
    dh0up = be[iflag] * t * t * fz3 / s1 * (- 7.0 / 3.0 * xy - qy *
                                     (af * bfup / be[iflag] - 7.0 / 3.0));
    dh0dw = be[iflag] * t * t * fz3 / s1 * (- 7.0 / 3.0 * xy - qy *
                                     (af * bfdw / be[iflag] - 7.0 / 3.0));
    dh0zup = (3.0 * h0 / fz - be[iflag] * t * t * fz2 / s1 * (2.0 * xy -
              qy * (3.0 * af * expe * ec / fz3 / be[iflag] + 2.0))) * dfz * (1.0 -
                      zeta);
    dh0zdw = - (3.0 * h0 / fz - be[iflag] * t * t * fz2 / s1 * (2.0 * xy -
                qy * (3.0 * af * expe * ec / fz3 / be[iflag] + 2.0))) * dfz * (1.0 +
                        zeta);
    ddh0 = be[iflag] * fz / (2.0 * ks * ks * rho) * (xy - qy) / s1;
    sc = rho * h0;
    v1cup = h0 + dh0up + dh0zup;
    v1cdw = h0 + dh0dw + dh0zdw;
    v2c = ddh0;
    return;
}

} // namespace xc_math

#undef ABACUS_XC_HOST_DEVICE
#endif
