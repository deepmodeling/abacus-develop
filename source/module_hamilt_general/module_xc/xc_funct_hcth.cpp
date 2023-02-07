// This file contains realizations of the HCTH GGA functional

#include "xc_functional.h"

void XC_Functional::hcth(const double rho, const double grho, double &sx, double &v1x, double &v2x)
{
    //     ===============================================================
    //     HCTH/120, JCP 109, p. 6264 (1998)
    //     Parameters set-up after N.L. Doltsisnis & M. Sprik (1999)
    //     Present release: Mauro Boero, Tsukuba, 11/05/2004
    //--------------------------------------------------------------------------
    //     rhoa = rhob = 0.5 * rho
    //     grho is the SQUARE of the gradient of rho// --> gr=sqrt(grho)
    //     sx  : total exchange correlation energy at point r
    //     v1x : d(sx)/drho  (eq. dfdra = dfdrb in original)
    //     v2x : 1/gr*d(sx)/d(gr) (eq. 0.5 * dfdza = 0.5 * dfdzb in original)
    //--------------------------------------------------------------------------
    // USE kinds
    // implicit none
    // real(kind=DP) :: rho, grho, sx, v1x, v2x

    // parameter :
    const double pi = 3.14159265358979323846;
    double o3 = 1.00 / 3.00;
    double o34 = 4.00 / 3.00;
    double fr83 = 8.0 / 3.0;
    //double cg0[6], cg1[6], caa[6], cab[6], cx[6];
    double cg0[7], cg1[7], caa[7], cab[7], cx[7];	//mohan modify 2007-10-13
    double r3q2, r3pi, gr, rho_o3, rho_o34, xa, xa2, ra, rab,
    dra_drho, drab_drho, g, dg, era1, dera1_dra, erab0, derab0_drab,
    ex, dex_drho, uaa, uab, ux, ffaa, ffab,  dffaa_drho, dffab_drho,
    denaa, denab, denx, f83rho, bygr, gaa, gab, gx, taa, tab, txx,
    dgaa_drho, dgab_drho, dgx_drho, dgaa_dgr, dgab_dgr, dgx_dgr;

    r3q2 = std::pow(2.0, (-o3));
    r3pi = std::pow((3.0 / pi), o3);
    //.....coefficients for pwf correlation......................................
    cg0[1] = 0.0310910;
    cg0[2] = 0.2137000;
    cg0[3] = 7.5957000;
    cg0[4] = 3.5876000;
    cg0[5] = 1.6382000;
    cg0[6] = 0.4929400;

    cg1[1] = 0.0155450;
    cg1[2] = 0.2054800;
    cg1[3] = 14.1189000;
    cg1[4] = 6.1977000;
    cg1[5] = 3.3662000;
    cg1[6] = 0.6251700;
    //......hcth-19-4.....................................
    caa[1] =  0.489508e+00;
    caa[2] = -0.260699e+00;
    caa[3] =  0.432917e+00;
    caa[4] = -0.199247e+01;
    caa[5] =  0.248531e+01;
    caa[6] =  0.200000e+00;

    cab[1] =  0.514730e+00;
    cab[2] =  0.692982e+01;
    cab[3] = -0.247073e+02;
    cab[4] =  0.231098e+02;
    cab[5] = -0.113234e+02;
    cab[6] =  0.006000e+00;

    cx[1] =  0.109163e+01;
    cx[2] = -0.747215e+00;
    cx[3] =  0.507833e+01;
    cx[4] = -0.410746e+01;
    cx[5] =  0.117173e+01;
    cx[6] =   0.004000e+00;
    //...........................................................................
    gr = sqrt(grho);
    rho_o3 = pow(rho, (o3));
    rho_o34 = pow(rho, (o34));
    xa = 1.259921050 * gr / rho_o34;
    xa2 = xa * xa;
    ra = 0.7815926420 / rho_o3;
    rab = r3q2 * ra;
    dra_drho = -0.2605308810 / rho_o34;
    drab_drho = r3q2 * dra_drho;
    XC_Functional::pwcorr(ra, cg1, g, dg);
    era1 = g;
    dera1_dra = dg;
    XC_Functional::pwcorr(rab, cg0, g, dg);
    erab0 = g;
    derab0_drab = dg;
    ex = -0.750 * r3pi * rho_o34;
    dex_drho = -r3pi * rho_o3;
    uaa = caa[6] * xa2;
    uaa = uaa / (1.00 + uaa);
    uab = cab[6] * xa2;
    uab = uab / (1.00 + uab);
    ux = cx[6] * xa2;
    ux = ux / (1.00 + ux);
    ffaa = rho * era1;
    ffab = rho * erab0 - ffaa;
    dffaa_drho = era1 + rho * dera1_dra * dra_drho;
    dffab_drho = erab0 + rho * derab0_drab * drab_drho - dffaa_drho;
    // mb-> i-loop removed
    denaa = 1.0 / (1.00 + caa[6] * xa2);
    denab = 1.0 / (1.00 + cab[6] * xa2);
    denx = 1.0 / (1.00 + cx[6] * xa2);
    f83rho = fr83 / rho;
    bygr = 2.00 / gr;
    gaa = caa[1] + uaa * (caa[2] + uaa * (caa[3] + uaa * (caa[4] + uaa * caa[5])));
    gab = cab[1] + uab * (cab[2] + uab * (cab[3] + uab * (cab[4] + uab * cab[5])));
    gx = cx[1] + ux * (cx[2] + ux * (cx[3] + ux * (cx[4] + ux * cx[5])));
    taa = denaa * uaa * (caa[2] + uaa * (2.0 * caa[3] + uaa
                                         * (3.0 * caa[4] + uaa * 4.0 * caa[5])));
    tab = denab * uab * (cab[2] + uab * (2.0 * cab[3] + uab
                                         * (3.0 * cab[4] + uab * 4.0 * cab[5])));
    txx = denx * ux * (cx[2] + ux * (2.0 * cx[3] + ux
                                     * (3.0 * cx[4] + ux * 4.0 * cx[5])));
    dgaa_drho = -f83rho * taa;
    dgab_drho = -f83rho * tab;
    dgx_drho = -f83rho * txx;
    dgaa_dgr = bygr * taa;
    dgab_dgr = bygr * tab;
    dgx_dgr = bygr * txx;
    // mb
    sx = ex * gx + ffaa * gaa + ffab * gab;
    v1x = dex_drho * gx + ex * dgx_drho
          + dffaa_drho * gaa + ffaa * dgaa_drho
          + dffab_drho * gab + ffab * dgab_drho;
    v2x = (ex * dgx_dgr + ffaa * dgaa_dgr + ffab * dgab_dgr) / gr;
    return;
} //end subroutine hcth

void XC_Functional::pwcorr(const double r, const double c[], double &g, double &dg)
{
    // USE kinds
    // implicit none
    // real(kind=DP) :: r, g, dg, c(6)
    double r12, r32, r2, rb, drb, sb;

    r12 = sqrt(r);
    r32 = r * r12;
    r2 = r * r;
    rb = c[3] * r12 + c[4] * r + c[5] * r32 + c[6] * r2;	//c[i] i=0--5;
    sb = 1.00 + 1.00 / (2.00 * c[1] * rb);
    g = -2.00 * c[1] * (1.00 + c[2] * r) * log(sb);
    drb = c[3] / (2.00 * r12) + c[4] + 1.50 * c[5] * r12 + 2.00 * c[6] * r;
    dg = (1.00 + c[2] * r) * drb / (rb * rb * sb) - 2.00 * c[1] * c[2] * log(sb);

    return;
} //end subroutine pwcorr