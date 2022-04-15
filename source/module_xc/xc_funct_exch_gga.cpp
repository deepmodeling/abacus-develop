// This file contains realizations of gradient correction to exchange part
// Spin unpolarized ones:
//  1. becke88 : Becke88 exchange
//  2. ggax : PW91 exchange
//  3. pbex : PBE exchange (and revPBE)
//  4. optx : OPTX, Handy et al.
//  5. wcx : Wu-Cohen exchange
// And some of their spin polarized counterparts:
//  1. becke88_spin

#include "xc_functional.h"

void XC_Functional::becke88(const double &rho, const double &grho, double &sx, double &v1x, double &v2x)
{
    //-----------------------------------------------------------------------
    // Becke exchange: A.D. Becke, PRA 38, 3098 (1988)
    // only gradient-corrected part, no Slater term included
    //
    double beta, third, two13;
    beta = 0.00420;
    third = 1.0 / 3.0;
    two13 = 1.2599210498948730;
    double rho13, rho43, xs, xs2, sa2b8, shm1, dd, dd2, ee;

    rho13 = pow(rho, third);
    rho43 = pow(rho13, 4);
    xs = two13 * sqrt(grho) / rho43;
    xs2 = xs * xs;
    sa2b8 = sqrt(1.00 + xs2);
    shm1 = log(xs + sa2b8);
    dd = 1.00 + 6.00 * beta * xs * shm1;
    dd2 = dd * dd;
    ee = 6.00 * beta * xs2 / sa2b8 - 1.0;
    sx = two13 * grho / rho43 * (- beta / dd);
    v1x = - (4.0 / 3.0) / two13 * xs2 * beta * rho13 * ee / dd2;
    v2x = two13 * beta * (ee - dd) / (rho43 * dd2);

    return;
} // end subroutine becke88

void XC_Functional::ggax(const double &rho, const double &grho, double &sx, double &v1x, double &v2x)
{
    //-----------------------------------------------------------------------
    // Perdew-Wang GGA (PW91), exchange part:
    // J.P. Perdew et al.,PRB 46, 6671 (1992)
    const double f1 = 0.196450;
    const double f2 = 7.79560;
    const double f3 = 0.27430;
    const double f4 = 0.15080;
    const double f5 = 0.0040;
    const double fp1 = -0.0192920212964260;
    const double fp2 = 0.1616204596739950;
    // fp1 = -3/(16 pi)*(3 pi^2)^(-1/3)
    double rhom43, s, s2, s3, s4, exps, as, sa2b8, shm1, bs, das,
    dbs, dls;

    rhom43 = pow(rho, (- 4.0 / 3.0));
    s = fp2 * sqrt(grho) * rhom43;
    s2 = s * s;
    s3 = s2 * s;
    s4 = s2 * s2;
    exps = f4 * exp(- 100.0 * s2);
    as = f3 - exps - f5 * s2;
    sa2b8 = sqrt(1.00 + f2 * f2 * s2);
    shm1 = log(f2 * s + sa2b8);
    bs = 1.0 + f1 * s * shm1 + f5 * s4;
    das = (200.0 * exps - 2.0 * f5) * s;
    dbs = f1 * (shm1 + f2 * s / sa2b8) + 4.0 * f5 * s3;
    dls = (das / as - dbs / bs);
    
	sx = fp1 * grho * rhom43 * as / bs;
	v1x = - 4.0 / 3.0 * sx / rho * (1.0 + s * dls);
    v2x = fp1 * rhom43 * as / bs * (2.0 + s * dls);

    return;
} //end subroutine ggax

void XC_Functional::pbex(const double &rho, const double &grho, const int &iflag, 
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

void XC_Functional::optx(const double rho, const double grho, double &sx, double &v1x, double &v2x)
{
    //     OPTX, Handy et al. JCP 116, p. 5411 (2002) and refs. therein
    //     Present release: Mauro Boero, Tsukuba, 10/9/2002
    //--------------------------------------------------------------------------
    //     rhoa = rhob = 0.5 * rho in LDA implementation
    //     grho is the SQUARE of the gradient of rho// --> gr=sqrt(grho)
    //     sx  : total exchange correlation energy at point r
    //     v1x : d(sx)/drho
    //     v2x : 1/gr*d(sx)/d(gr)
    //--------------------------------------------------------------------------
    // use kinds, only: DP
    // implicit none
    // real(kind=DP) :: rho, grho, sx, v1x, v2x

    // parameter :
    double small = 1.e-30;
    double smal2 = 1.e-10;
    //.......coefficients and exponents....................
    // parameter :
    double o43 = 4.00 / 3.00,
                 two13 = 1.2599210498948730,
                         two53 = 3.1748021039363990,
                                 gam = 0.0060,
                                       a1cx = 0.97845711702844210,
                                              a2 = 1.431690;
    double gr, rho43, xa, gamx2, uden, uu;
    //.......OPTX in compact form..........................

    if (rho <= small)
    {
        sx = 0.00;
        v1x = 0.00;
        v2x = 0.00;
    }
    else
    {
        gr = (grho > smal2) ? grho : smal2;	//max()
        rho43 = pow(rho, o43);
        xa = two13 * sqrt(gr) / rho43;
        gamx2 = gam * xa * xa;
        uden = 1.e+00 / (1.e+00 + gamx2);
        uu = a2 * gamx2 * gamx2 * uden * uden;
        uden = rho43 * uu * uden;
        sx = -rho43 * (a1cx + uu) / two13;
        v1x = o43 * (sx + two53 * uden) / rho;
        v2x = -two53 * uden / gr;
    } //endif

    return;
} // end subroutine optx

void XC_Functional::wcx(const double &rho,const double &grho, double &sx, double &v1x, double &v2x)
{
  double kf, agrho, s1, s2, es2, ds, dsg, exunif, fx;
  // (3*pi2*|rho|)^(1/3)
  // |grho|
  // |grho|/(2*kf*|rho|)
  // s^2
  // n*ds/dn
  // n*ds/d(gn)
  // exchange energy LDA part
  // exchange energy gradient part
  double dxunif, dfx, f1, f2, f3, dfx1, x1, x2, x3, dxds1, dxds2, dxds3;
  // numerical coefficients (NB: c2=(3 pi^2)^(1/3) )
  double third, c1, c2, c5, teneightyone;	// c6
  const double pi = 3.14159265358979323846;

  third = 1.0/3.0;
  c1 = 0.75/pi;
  c2 = 3.093667726280136;
  c5 = 4.0 * third;
  teneightyone = 0.123456790123;
  // parameters of the functional
  double k, mu, cwc;
  k = 0.804;
  mu = 0.2195149727645171;
  cwc = 0.00793746933516;
  //
  agrho = sqrt (grho);
  kf = c2 * pow(rho,third);
  dsg = 0.5 / kf;
  s1 = agrho * dsg / rho;
  s2 = s1 * s1;
  es2 = exp(-s2);
  ds = - c5 * s1;
  //
  //   Energy
  //
  // x = 10/81 s^2 + (mu - 10/81) s^2 e^-s^2 + ln (1 + c s^4)
  x1 = teneightyone * s2;
  x2 = (mu - teneightyone) * s2 * es2;
  x3 = log(1.0 + cwc * s2 * s2);
  f1 = (x1 + x2 + x3) / k;
  f2 = 1.0 + f1;
  f3 = k / f2;
  fx = k - f3;
  exunif = - c1 * kf;
  sx = exunif * fx;
  //
  //   Potential
  //
  dxunif = exunif * third;
  dfx1 = f2 * f2;
  dxds1 = teneightyone;
  dxds2 = (mu - teneightyone) * es2 * (1.0 - s2);
  dxds3 = 2.0 * cwc * s2 / (1.0 + cwc * s2 *s2);
  dfx = 2.0 * s1 * (dxds1 + dxds2 + dxds3) / dfx1;
  v1x = sx + dxunif * fx + exunif * dfx * ds;
  v2x = exunif * dfx * dsg / agrho;

  sx = sx * rho;
  return;
}

void XC_Functional::becke88_spin(double rho, double grho, double &sx, double &v1x, double &v2x)
{
    //-------------------------------------------------------------------
    // Becke exchange: A.D. Becke, PRA 38, 3098 (1988) - Spin polarized case

    // USE kinds
    // implicit none
    // real(kind=DP) :: rho, grho, sx, v1x, v2x
    // input: charge
    // input: gradient
    // output: the up and down energies
    // output: first part of the potential
    // output: the second part of the potential

    double beta, third;
    // parameter :
    beta = 0.00420;
    third = 1.0 / 3.0;
    double rho13, rho43, xs, xs2, sa2b8, shm1, dd, dd2, ee;

    rho13 = pow(rho, third);
    rho43 = pow(rho13, 4);
    xs = sqrt(grho) / rho43;
    xs2 = xs * xs;
    sa2b8 = sqrt(1.00 + xs2);
    shm1 = log(xs + sa2b8);
    dd = 1.00 + 6.00 * beta * xs * shm1;
    dd2 = dd * dd;
    ee = 6.00 * beta * xs2 / sa2b8 - 1.0;
    sx = grho / rho43 * (- beta / dd);
    v1x = - (4.0 / 3.0) * xs2 * beta * rho13 * ee / dd2;
    v2x = beta * (ee - dd) / (rho43 * dd2);

    return;
} //end subroutine becke88_spin
