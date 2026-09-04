#include "vnl_pw.h"

#include "source_io/module_parameter/parameter.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/math_sphbes.h"
#include "source_base/math_integral.h"
#include "source_base/math_polyint.h"
#include "source_base/timer.h"

#include <cmath>
#include <complex>
#include <vector>

#ifdef __LCAO

/**
 * @file vnl_pw_alpha.cpp
 * @brief LCAO-mode alpha-channel VNL initialization and Gaunt coefficient helpers.
 *
 * This file contains:
 * - Cal_C(): compute the C coefficient for alpha channel
 * - CG(): wrapper for Gaunt coefficient lookup via MGT
 * - init_vnl_alpha(): build the tab_alpha interpolation table including L index
 *
 * All functions are guarded by #ifdef __LCAO.
 */

/**
 * @brief Compute the C coefficient for alpha channel.
 *
 * Combines spherical harmonics with Gaunt coefficients to produce
 * the alpha-resolved Kleinman-Bylander projector coefficient.
 *
 * @param alpha  channel index (0, 1, 2)
 * @param lu     orbital angular momentum l
 * @param mu     magnetic quantum number m
 * @param L      total angular momentum L
 * @param M      total magnetic quantum number M
 * @return complex coefficient
 */
std::complex<double> pseudopot_cell_vnl::Cal_C(int alpha, int lu, int mu, int L, int M) // pengfei Li  2018-3-23
{
    std::complex<double> cf;
    if (alpha == 0)
    {
        cf = -sqrt(4 * ModuleBase::PI / 3) * CG(lu, mu, 1, 1, L, M);
    }
    else if (alpha == 1)
    {
        cf = -sqrt(4 * ModuleBase::PI / 3) * CG(lu, mu, 1, 2, L, M);
    }
    else if (alpha == 2)
    {
        cf = sqrt(4 * ModuleBase::PI / 3) * CG(lu, mu, 1, 0, L, M);
    }
    else
    {
        ModuleBase::WARNING_QUIT("pseudopot_cell_vnl_alpha", "alpha must be 0~2");
    }

    return cf;
}

/**
 * @brief Wrapper for Gaunt coefficient lookup via MGT.
 *
 * Converts (l, m) pairs to linear indices and delegates to
 * ORB_gaunt_table::Gaunt_Coefficients.
 *
 * @param l1 first angular momentum
 * @param m1 first magnetic quantum number
 * @param l2 second angular momentum
 * @param m2 second magnetic quantum number
 * @param L  total angular momentum
 * @param M  total magnetic quantum number
 * @return Gaunt coefficient value
 */
double pseudopot_cell_vnl::CG(int l1, int m1, int l2, int m2, int L, int M) // pengfei Li 2018-3-23
{
    int dim = L * L + M;
    int dim1 = l1 * l1 + m1;
    int dim2 = l2 * l2 + m2;

    return MGT.Gaunt_Coefficients(dim1, dim2, dim);
}

/**
 * @brief Build the alpha-channel interpolation table tab_alpha.
 *
 * For each atom type and each beta function, compute the radial
 * integral of betar * j_L(q*r) * r^2 for all L channels, using
 * Simpson integration on the pseudopotential radial grid.
 *
 * @param ucell unit cell (read-only)
 */
void pseudopot_cell_vnl::init_vnl_alpha(const UnitCell& ucell) // pengfei Li 2018-3-23
{
    if (PARAM.inp.test_pp) {
        ModuleBase::TITLE("pseudopot_cell_vnl", "init_vnl_alpha");
    }
    ModuleBase::timer::start("ppcell_vnl", "init_vnl_alpha");

    for (int it = 0; it < ucell.ntype; it++)
    {
        int BetaIndex = 0;
        for (int ib = 0; ib < ucell.atoms[it].ncpp.nbeta; ib++)
        {
            const int l = ucell.atoms[it].ncpp.lll[ib];
            for (int m = 0; m < 2 * l + 1; m++)
            {
                this->nhtol(it, BetaIndex) = l;
                this->nhtolm(it, BetaIndex) = l * l + m;
                this->indv(it, BetaIndex) = ib;
                ++BetaIndex;
            }
        }
    }

    // max number of beta functions
    const int nbrx = 10;

    const double pref = ModuleBase::FOUR_PI / sqrt(ucell.omega);
    this->tab_alpha.create(ucell.ntype, nbrx, lmaxkb + 2, PARAM.globalv.nqx);
    this->tab_alpha.zero_out();
    GlobalV::ofs_running << "\n Init Non-Local PseudoPotential table( including L index) : ";
    for (int it = 0; it < ucell.ntype; it++)
    {
        const int nbeta = ucell.atoms[it].ncpp.nbeta;
        int kkbeta = ucell.atoms[it].ncpp.kkbeta;

        // mohan modify 2008-3-31
        // mohan add kkbeta>0 2009-2-27
        if ((kkbeta % 2 == 0) && kkbeta > 0)
        {
            kkbeta--;
        }

        std::vector<double> jl(kkbeta);
        std::vector<double> aux(kkbeta);

        for (int ib = 0; ib < nbeta; ib++)
        {
            for (int L = 0; L <= lmaxkb + 1; L++)
            {
                for (int iq = 0; iq < PARAM.globalv.nqx; iq++)
                {
                    const double q = iq * PARAM.globalv.dq;
                    ModuleBase::Sphbes::Spherical_Bessel(kkbeta, ucell.atoms[it].ncpp.r.data(), q, L, jl.data());

                    for (int ir = 0; ir < kkbeta; ir++)
                    {
                        aux[ir] = ucell.atoms[it].ncpp.betar(ib, ir) * jl[ir]
                                  * ucell.atoms[it].ncpp.r[ir] * ucell.atoms[it].ncpp.r[ir];
                    }
                    double vqint = 0.0;
                    ModuleBase::Integral::Simpson_Integral(kkbeta, aux.data(), ucell.atoms[it].ncpp.rab.data(), vqint);
                    this->tab_alpha(it, ib, L, iq) = vqint * pref;
                }
            }
        }
    }
    ModuleBase::timer::end("ppcell_vnl", "init_vnl_alpha");
    GlobalV::ofs_running << "\n Init Non-Local-Pseudopotential done(including L)." << std::endl;
}

#endif // __LCAO
