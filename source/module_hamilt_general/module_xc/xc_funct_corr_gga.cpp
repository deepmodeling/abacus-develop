// This file contains realizations of gradient correction to correlation part
// Spin unpolarized ones:
//  1. perdew86 : P86
//  2. ggac : PW91
//  3. pbec
//  4. glyp
// And some of their spin polarized counterparts:
//  1. perdew86_spin
//  2. ggac_spin
//  3. pbec_spin

#include "xc_functional.h"

void XC_Functional::perdew86(const double rho, const double grho, double &sc, double &v1c, double &v2c)
{
    // Perdew gradient correction on correlation: PRB 33, 8822 (1986)

    // USE kinds
    // implicit none
    // real(kind=DP) :: rho, grho, sc, v1c, v2c
    double p1, p2, p3, p4, pc1, pc2, pci;
    // parameter :
    p1 = 0.0232660;
    p2 = 7.389e-6;
    p3 = 8.7230;
    p4 = 0.4720;
    // parameter :
    pc1 = 0.0016670;
    pc2 = 0.0025680;
    pci = pc1 + pc2;
    double third, pi34;
    // parameter :
    third = 1.0 / 3.0;
    pi34 = 0.62035049089940;  // pi34=(3/4pi)^(1/3)
    double rho13, rho43, rs, rs2, rs3, cna, cnb, cn, drs;
    double dcna, dcnb, dcn, phi, ephi;

    rho13 = pow(rho, third);
    rho43 = pow(rho13, 4);
    rs = pi34 / rho13;
    rs2 = rs * rs;
    rs3 = rs * rs2;
    cna = pc2 + p1 * rs + p2 * rs2;
    cnb = 1.0 + p3 * rs + p4 * rs2 + 1.e4 * p2 * rs3;
    cn = pc1 + cna / cnb;
    drs = - third * pi34 / rho43;
    dcna = (p1 + 2.0 * p2 * rs) * drs;
    dcnb = (p3 + 2.0 * p4 * rs + 3.e4 * p2 * rs2) * drs;
    dcn = dcna / cnb - cna / (cnb * cnb) * dcnb;
    phi = 0.1920 * pci / cn * sqrt(grho) * pow(rho, (- 7.0 / 6.0));
    // SdG: in the original paper 1.745*0.11=0.19195 is used
    ephi = exp(- phi);
    sc = grho / rho43 * cn * ephi;
    v1c = sc * ((1.0 + phi) * dcn / cn - ((4.0 / 3.0) - (7.0 /
                                          6.0) * phi) / rho);
    v2c = cn * ephi / rho43 * (2.0 - phi);

    return;
} //end subroutine perdew86

void XC_Functional::ggac(const double &rho,const double &grho, double &sc, double &v1c, double &v2c)
{
    // Perdew-Wang GGA (PW91) correlation part
    double al, pa, pb, pc, pd, cx, cxc0, cc0;
    // parameter :
    al = 0.090;
    pa = 0.0232660;
    pb = 7.389e-6;
    pc = 8.7230;
    pd = 0.4720;
    // parameter :
    cx = -0.0016670;
    cxc0 = 0.0025680;
    cc0 = - cx + cxc0;
    double third, pi34, nu, be, xkf, xks;
    // parameter :
    third = 1.0 / 3.0;
    pi34 = 0.62035049089940;
    // parameter :
    nu = 15.7559203494831440;
    be = nu * cc0;
    // parameter :
    xkf = 1.9191582926775130;
    xks = 1.1283791670955130;
    // pi34=(3/4pi)^(1/3),  nu=(16/pi)*(3 pi^2)^(1/3)
    // xkf=(9 pi/4)^(1/3), xks= sqrt(4/pi)
    double kf, ks, rs, rs2, rs3, ec, vc, t, expe, af, bf, y, xy,
    qy, s1;
    double h0, dh0, ddh0, ee, cn, dcn, cna, dcna, cnb, dcnb, h1,
    dh1, ddh1;

	rs = pi34 / pow(rho, third);
    rs2 = rs * rs;
    rs3 = rs * rs2;

	// call function: pw
	XC_Functional::pw(rs, 0, ec, vc);
    kf = xkf / rs;
    ks = xks * sqrt(kf);
    t = sqrt(grho) / (2.0 * ks * rho);
    expe = exp(- 2.0 * al * ec / (be * be));
    af = 2.0 * al / be * (1.0 / (expe - 1.0));
    bf = expe * (vc - ec);
    y = af * t * t;
    xy = (1.0 + y) / (1.0 + y + y * y);

    double x;
    x = 1.0 + y + y * y;
    qy = y * y * (2.0 + y) / (x * x);
    s1 = 1.0 + 2.0 * al / be * t * t * xy;
    h0 = be * be / (2.0 * al) * log(s1);
    dh0 = be * t * t / s1 * (- 7.0 / 3.0 * xy - qy * (af * bf /
                             be - 7.0 / 3.0));
    ddh0 = be / (2.0 * ks * ks * rho) * (xy - qy) / s1;

    x = ks / kf * t;
    ee = - 100.0 * (x * x);
    cna = cxc0 + pa * rs + pb * rs2;
    dcna = pa * rs + 2.0 * pb * rs2;
    cnb = 1.0 + pc * rs + pd * rs2 + 1.e4 * pb * rs3;
    dcnb = pc * rs + 2.0 * pd * rs2 + 3.e4 * pb * rs3;
    cn = cna / cnb - cx;
    dcn = dcna / cnb - cna * dcnb / (cnb * cnb);
    h1 = nu * (cn - cc0 - 3.0 / 7.0 * cx) * t * t * exp(ee);
    dh1 = - third * (h1 * (7.0 + 8.0 * ee) + nu * t * t * exp(ee)
                     * dcn);
    ddh1 = 2.0 * h1 * (1.0 + ee) * rho / grho;
    sc = rho * (h0 + h1);
    v1c = h0 + h1 + dh0 + dh1;
    v2c = ddh0 + ddh1;

    return;
} //end subroutine ggac

void XC_Functional::pbec(const double &rho, const double &grho, const int &iflag, double &sc, double &v1c, double &v2c)
{
	// PBE correlation (without LDA part)
	// iflag=0: J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996).
	// iflag=1: J.P.Perdew et al., PRL 100, 136406 (2008).
	const double ga = 0.0310906908696548950;
	const double be[2] = {0.06672455060314922, 0.046};

	const double third = 1.0 / 3.0;
	const double pi34 = 0.62035049089940;
	const double xkf = 1.9191582926775130;
	const double xks = 1.1283791670955130;

	// pi34=(3/4pi)^(1/3), xkf=(9 pi/4)^(1/3), xks= sqrt(4/pi)
	double ec, vc;

	const double rs = pi34 / pow(rho, third);

	XC_Functional::pw(rs, 0, ec, vc);

	const double kf = xkf / rs;
	const double ks = xks * sqrt(kf);
	const double t = sqrt(grho) / (2.0 * ks * rho);
	const double expe = exp(- ec / ga);
	const double af = be[iflag] / ga * (1.0 / (expe - 1.0));
	const double bf = expe * (vc - ec);
	const double y = af * t * t;
	const double xy = (1.0 + y) / (1.0 + y + y * y);

	const double x = 1.0 + y + y * y;
	const double qy = y * y * (2.0 + y) / (x * x);
	const double s1 = 1.0 + be[iflag] / ga * t * t * xy;
	const double h0 = ga * log(s1);
	const double dh0 = be[iflag] * t * t / s1 * (- 7.0 / 3.0 * xy - qy * (af * bf /be[iflag] - 7.0 / 3.0));
	const double ddh0 = be[iflag] / (2.0 * ks * ks * rho) * (xy - qy) / s1;

	sc = rho * h0;
	v1c = h0 + dh0;
    v2c = ddh0;

	return;
}

void XC_Functional::glyp(const double &rho, const double &grho, double &sc, double &v1c, double &v2c)
{
    //-----------------------------------------------------------------------
    // Lee Yang Parr: gradient correction part

    double a, b, c, d;
    a = 0.049180;
    b = 0.1320;
    c = 0.25330;
    d = 0.3490;
    double rhom13, rhom43, rhom53, om, xl, ff, dom, dxl;

    rhom13 = pow(rho, (- 1.0 / 3.0));
    om = exp(- c * rhom13) / (1.0 + d * rhom13);
    xl = 1.0 + (7.0 / 3.0) * (c * rhom13 + d * rhom13 / (1.0 + d *
                              rhom13));
    ff = a * b * grho / 24.0;
    rhom53 = pow(rhom13, 5);
    sc = ff * rhom53 * om * xl;
    dom = - om * (c + d + c * d * rhom13) / (1.0 + d * rhom13);

    double x;
    x = 1.0 + d * rhom13;
    dxl = (7.0 / 3.0) * (c + d + 2.0 * c * d * rhom13 + c * d * d *
                         rhom13 * rhom13) / (x * x);
    rhom43 = pow(rhom13, 4);
    v1c = - ff * rhom43 / 3.0 * (5.0 * rhom43 * om * xl + rhom53 *
                                 dom * xl + rhom53 * om * dxl);
    v2c = 2.0 * sc / grho;

    return;
} // end subroutine glyp

//-----------------------------------------------------------------------
void XC_Functional::perdew86_spin(double rho, double zeta, double grho, double &sc,
                   double &v1cup, double &v1cdw, double &v2c)
{
    //-------------------------------------------------------------------
    // Perdew gradient correction on correlation: PRB 33, 8822 (1986)
    // spin-polarized case

    // USE kinds
    // implicit none
    // real(kind=DP) :: rho, zeta, grho, sc, v1cup, v1cdw, v2c
    double p1, p2, p3, p4, pc1, pc2, pci;
    // parameter :
    p1 = 0.0232660;
    p2 = 7.389e-6;
    p3 = 8.7230;
    p4 = 0.4720;
    // parameter :
    pc1 = 0.0016670;
    pc2 = 0.0025680;
    pci = pc1 + pc2;
    double third, pi34;
    // parameter :
    third = 1.0 / 3.0;
    pi34 = 0.62035049089940;
    // pi34=(3/4pi)^(1/3)

    double rho13, rho43, rs, rs2, rs3, cna, cnb, cn, drs;
    double dcna, dcnb, dcn, phi, ephi, dd, ddd;

    rho13 = pow(rho, third);
    rho43 = pow(rho13, 4);
    rs = pi34 / rho13;
    rs2 = rs * rs;
    rs3 = rs * rs2;
    cna = pc2 + p1 * rs + p2 * rs2;
    cnb = 1.0 + p3 * rs + p4 * rs2 + 1.e4 * p2 * rs3;
    cn = pc1 + cna / cnb;
    drs = - third * pi34 / rho43;
    dcna = (p1 + 2.0 * p2 * rs) * drs;
    dcnb = (p3 + 2.0 * p4 * rs + 3.e4 * p2 * rs2) * drs;
    dcn = dcna / cnb - cna / (cnb * cnb) * dcnb;
    phi = 0.1920 * pci / cn * sqrt(grho) * pow(rho, (- 7.0 / 6.0));
    //SdG: in the original paper 1.745*0.11=0.19195 is used
    dd = pow((2.0) , third) * sqrt(pow(((1.0 + zeta) * 0.50) , (5.0 /
                                       3.0)) + pow(((1.0 - zeta) * 0.50) , (5.0 / 3.0)));
    ddd = pow((2.0) , (- 4.0 / 3.0)) * 5.0 * (pow(((1.0 + zeta)
            * 0.50) , (2.0 / 3.0)) - pow(((1.0 - zeta) * 0.50) , (2.0 /
                                         3.0))) / (3.0 * dd);
    ephi = exp(- phi);
    sc = grho / rho43 * cn * ephi / dd;
    v1cup = sc * ((1.0 + phi) * dcn / cn - ((4.0 / 3.0) -
                                            (7.0 / 6.0) * phi) / rho) - sc * ddd / dd * (1.0 - zeta)
            / rho;
    v1cdw = sc * ((1.0 + phi) * dcn / cn - ((4.0 / 3.0) -
                                            (7.0 / 6.0) * phi) / rho) + sc * ddd / dd * (1.0 + zeta)
            / rho;
    v2c = cn * ephi / rho43 * (2.0 - phi) / dd;

    return;
} //end subroutine perdew86_spin

// There seems to be something wrong with it
// lots of terms evaluates to inf / nan in unit test
/*
void XC_Functional::ggac_spin(double rho, double zeta, double grho, double &sc,
               double &v1cup, double &v1cdw, double &v2c)
{
    //-------------------------------------------------------------------
    // Perdew-Wang GGA (PW91) correlation part - spin-polarized

    // USE kinds
    // implicit none
    // real(kind=DP) :: rho, zeta, grho, sc, v1cup, v1cdw, v2c
    double al, pa, pb, pc, pd, cx, cxc0, cc0;
    // parameter :
    al = 0.090;
    pa = 0.0232660;
    pb = 7.389e-6;
    pc = 8.7230;
    pd = 0.4720;
    // parameter :
    cx = - 0.0016670;
    cxc0 = 0.0025680;
    cc0 = - cx + cxc0;

    double third, pi34, nu, be, xkf, xks;
    // parameter :
    third = 1.0 / 3.0;
    pi34 = 0.6203504908994e0l;
    // parameter :
    nu = 15.7559203494831440;
    be = nu * cc0;
    // parameter :
    xkf = 1.9191582926775130;
    xks = 1.1283791670955130;
    // pi34=(3/4pi)^(1/3),  nu=(16/pi)*(3 pi^2)^(1/3)
    // xkf=(9 pi/4)^(1/3), xks= sqrt(4/pi)
    double kf, ks, rs, rs2, rs3, ec, vcup, vcdw, t, expe, af, y,
    xy, qy, s1, h0, ddh0, ee, cn, dcn, cna, dcna, cnb, dcnb, h1, dh1,
    ddh1, fz, fz2, fz3, fz4, dfz, bfup, dh0up, dh0zup, //bfdw, dh0dw,
    dh0zdw, dh1zup, dh1zdw;

    rs = pi34 / pow(rho, third);
    rs2 = rs * rs;
    rs3 = rs * rs2;
	XC_Functional::pw_spin(rs, zeta, ec, vcup, vcdw);
    kf = xkf / rs;
    ks = xks * sqrt(kf);
    fz = 0.50 * (pow((1.0 + zeta) , (2.0 / 3.0)) + pow((1.0 - zeta) , (
                     2.0 / 3.0)));
    fz2 = fz * fz;
    fz3 = fz2 * fz;
    fz4 = fz3 * fz;
    dfz = (pow((1.0 + zeta) , (- 1.0 / 3.0)) - pow((1.0 - zeta) , (-
            1.0 / 3.0))) / 3.0;
    t = sqrt(grho) / (2.0 * fz * ks * rho);
    expe = exp(- 2.0 * al * ec / (fz3 * be * be));
    af = 2.0 * al / be * (1.0 / (expe - 1.0));
    bfup = expe * (vcup - ec) / fz3;
    //bfdw = expe * (vcdw - ec) / fz3;
    y = af * t * t;
    xy = (1.0 + y) / (1.0 + y + y * y);
    qy = y * y * (2.0 + y) / (1.0 + y + y * y) ;	// **2;
    qy *= qy;
    s1 = 1.0 + 2.0 * al / be * t * t * xy;
    h0 = fz3 * be * be / (2.0 * al) * log(s1);
    dh0up = be * t * t * fz3 / s1 * (- 7.0 / 3.0 * xy - qy *
                                     (af * bfup / be - 7.0 / 3.0));
    //dh0dw = be * t * t * fz3 / s1 * (- 7.0 / 3.0 * xy - qy *
    //                                 (af * bfdw / be - 7.0 / 3.0));
    dh0zup = (3.0 * h0 / fz - be * t * t * fz2 / s1 * (2.0 * xy -
              qy * (3.0 * af * expe * ec / fz3 / be + 2.0))) * dfz * (1.0 -
                      zeta);
    dh0zdw = - (3.0 * h0 / fz - be * t * t * fz3 / s1 * (2.0 * xy -
                qy * (3.0 * af * expe * ec / fz3 / be + 2.0))) * dfz * (1.0 +
                        zeta);
    ddh0 = be * fz / (2.0 * ks * ks * rho) * (xy - qy) / s1;
    ee = - 100.0 * fz4 * (ks / kf * t);	// **2;
    ee *= ee;
    cna = cxc0 + pa * rs + pb * rs2;
    dcna = pa * rs + 2.0 * pb * rs2;
    cnb = 1.0 + pc * rs + pd * rs2 + 1.e4 * pb * rs3;
    dcnb = pc * rs + 2.0 * pd * rs2 + 3.e4 * pb * rs3;
    cn = cna / cnb - cx;
    dcn = dcna / cnb - cna * dcnb / (cnb * cnb);
    h1 = nu * (cn - cc0 - 3.0 / 7.0 * cx) * fz3 * t * t * exp(ee);
    dh1 = - third * (h1 * (7.0 + 8.0 * ee) + fz3 * nu * t * t * exp
                     (ee) * dcn);
    ddh1 = 2.0 * h1 * (1.0 + ee) * rho / grho;
    dh1zup = (1.0 - zeta) * dfz * h1 * (1.0 + 2.0 * ee / fz);
    dh1zdw = - (1.0 + zeta) * dfz * h1 * (1.0 + 2.0 * ee / fz);
    sc = rho * (h0 + h1);
    v1cup = h0 + h1 + dh0up + dh1 + dh0zup + dh1zup;
    v1cdw = h0 + h1 + dh0up + dh1 + dh0zdw + dh1zdw;
    v2c = ddh0 + ddh1;
    return;
} // end subroutine ggac_spin
*/

//---------------------------------------------------------------
void XC_Functional::pbec_spin(double rho, double zeta, double grho, const int &iflag, double &sc,
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
	XC_Functional::pw_spin(rs, zeta, ec, vcup, vcdw); //mohan fix bug 2012-05-28
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
} // end subroutine pbec_spin