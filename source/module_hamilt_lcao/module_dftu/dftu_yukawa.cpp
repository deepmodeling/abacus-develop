//==========================================================
// Author:Xin Qu
// DATE : 2019-12-10
//==========================================================
#include "module_base/constants.h"
#include "module_base/global_function.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"
#include "dftu.h"

#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string.h>

namespace ModuleDFTU
{

void DFTU::cal_yukawa_lambda(double** rho, const int& nrxx)
{
    ModuleBase::TITLE("DFTU", "cal_yukawa_lambda");

    if (INPUT.yukawa_lambda > 0)
    {
        this->lambda = INPUT.yukawa_lambda;
        return;
    }

    double sum_rho = 0.0;
    double sum_rho_lambda = 0.0;
    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        for (int ir = 0; ir < nrxx; ir++)
        {
            double rho_ir = rho[is][ir];
            sum_rho += rho_ir;

            double lambda_ir = 2 * pow(3 * rho_ir / ModuleBase::PI, (double)1.0 / 6.0);
            sum_rho_lambda += lambda_ir * rho_ir;
        }
    }

    double val1 = 0.0;
    double val2 = 0.0;

#ifdef __MPI
    MPI_Allreduce(&sum_rho, &val1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&sum_rho_lambda, &val2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

    this->lambda = val2 / val1;

    // rescaling
    this->lambda /= 1.6;

    return;
}

void DFTU::cal_slater_Fk(const int L, const int T)
{
    ModuleBase::TITLE("DFTU", "cal_slater_Fk");

    if (Yukawa)
    {
        for (int chi = 0; chi < GlobalC::ucell.atoms[T].l_nchi[L]; chi++)
        {
            //	if(chi!=0) continue;
            const int mesh = GlobalC::ORB.Phi[T].PhiLN(L, chi).getNr();

            for (int k = 0; k <= L; k++)
            {
                for (int ir0 = 1; ir0 < mesh; ir0++)
                {
                    double r0 = GlobalC::ORB.Phi[T].PhiLN(L, chi).getRadial(ir0);
                    const double rab0 = GlobalC::ORB.Phi[T].PhiLN(L, chi).getRab(ir0);
                    const double R_L0 = GlobalC::ORB.Phi[T].PhiLN(L, chi).getPsi(ir0);

                    for (int ir1 = 1; ir1 < mesh; ir1++)
                    {
                        double bslval, hnkval;
                        double r1 = GlobalC::ORB.Phi[T].PhiLN(L, chi).getRadial(ir1);
                        const double rab1 = GlobalC::ORB.Phi[T].PhiLN(L, chi).getRab(ir1);
                        const double R_L1 = GlobalC::ORB.Phi[T].PhiLN(L, chi).getPsi(ir1);

                        int l = 2 * k;
                        if (ir0 < ir1) // less than
                        {
                            bslval = this->spherical_Bessel(l, r0, lambda);
                            hnkval = this->spherical_Hankel(l, r1, lambda);
                        }
                        else // greater than
                        {
                            bslval = this->spherical_Bessel(l, r1, lambda);
                            hnkval = this->spherical_Hankel(l, r0, lambda);
                        }
                        this->Fk[T][L][chi][k] -= (4 * k + 1) * lambda * pow(R_L0, 2) * bslval * hnkval * pow(R_L1, 2)
                                                  * pow(r0, 2) * pow(r1, 2) * rab0 * rab1;
                    }
                }
            }
        }
    }

    return;
}

void DFTU::cal_slater_UJ(double** rho, const int& nrxx)
{
    ModuleBase::TITLE("DFTU", "cal_slater_UJ");
    if (!Yukawa)
        return;

    this->cal_yukawa_lambda(rho, nrxx);

    for (int it = 0; it < GlobalC::ucell.ntype; it++)
    {
        const int NL = GlobalC::ucell.atoms[it].nwl + 1;

        for (int l = 0; l < NL; l++)
        {
            int N = GlobalC::ucell.atoms[it].l_nchi[l];
            for (int n = 0; n < N; n++)
            {
                ModuleBase::GlobalFunc::ZEROS(ModuleBase::GlobalFunc::VECTOR_TO_PTR(this->Fk[it][l][n]), l + 1);
            }
        }
    }

    for (int T = 0; T < GlobalC::ucell.ntype; T++)
    {
        const int NL = GlobalC::ucell.atoms[T].nwl + 1;

        for (int L = 0; L < NL; L++)
        {
            const int N = GlobalC::ucell.atoms[T].l_nchi[L];

            if (L >= INPUT.orbital_corr[T] && INPUT.orbital_corr[T] != -1)
            {
                if (L != INPUT.orbital_corr[T])
                    continue;
                this->cal_slater_Fk(L, T);

                for (int n = 0; n < N; n++)
                {
                    if (n != 0)
                        continue;

                    switch (L)
                    {
                    case 1: // p electrons
                        this->U_Yukawa[T][L][n] = this->Fk[T][L][n][0];
                        this->J_Yukawa[T][L][n] = this->Fk[T][L][n][1] / 5.0;
                        break;

                    case 2: // d electrons
                        this->U_Yukawa[T][L][n] = this->Fk[T][L][n][0];
                        this->J_Yukawa[T][L][n] = (this->Fk[T][L][n][1] + this->Fk[T][L][n][2]) / 14.0;
                        break;

                    case 3: // f electrons
                        if (Yukawa)
                            this->U_Yukawa[T][L][n] = this->Fk[T][L][n][0];
                        this->J_Yukawa[T][L][n] = (286.0 * this->Fk[T][L][n][1] + 195.0 * this->Fk[T][L][n][2]
                                                   + 250.0 * this->Fk[T][L][n][3])
                                                  / 6435.0;
                        break;
                    }

                    // Hartree to Rydeberg
                    this->U_Yukawa[T][L][n] *= 2.0;
                    this->J_Yukawa[T][L][n] *= 2.0;
                } // end n
            } // end if
        } // end L
    } // end T

    return;
}

double DFTU::spherical_Bessel(const int k, const double r, const double lambda)
{
    ModuleBase::TITLE("DFTU", "spherical_Bessel");

    double val;
    double x = r * lambda;
    if (k == 0)
    {
        if (x < 1.0e-3)
            val = 1 + pow(x, 2) / 6.0;
        else
            val = sinh(x) / x;
    }
    else if (k == 2)
    {
        if (x < 1.0e-2)
            val = -pow(x, 2) / 15.0 - pow(x, 4) / 210.0 - pow(x, 6) / 7560.0;
        else
            val = 3 * cosh(x) / pow(x, 2) + (-3 - pow(x, 2)) * sinh(x) / pow(x, 3);
    }
    else if (k == 4)
    {
        if (x < 5.0e-1)
            val = pow(x, 4) / 945.0 + pow(x, 6) / 20790.0 + pow(x, 8) / 1081080.0 + pow(x, 10) / 97297200.0;
        else
            val = -5 * (21 + 2 * pow(x, 2)) * cosh(x) / pow(x, 4)
                  + (105 + 45 * pow(x, 2) + pow(x, 4)) * sinh(x) / pow(x, 5);
    }
    else if (k == 6)
    {
        if (x < 9.0e-1)
            val = -pow(x, 6) / 135135.0 - pow(x, 8) / 4054050.0 - pow(x, 10) / 275675400.0;
        else
            val = 21 * (495 + 60 * pow(x, 2) + pow(x, 4)) * cosh(x) / pow(x, 6)
                  + (-10395 - 4725 * pow(x, 2) - 210 * pow(x, 4) - pow(x, 6)) * sinh(x) / pow(x, 7);
    }
    return val;
}

double DFTU::spherical_Hankel(const int k, const double r, const double lambda)
{
    ModuleBase::TITLE("DFTU", "spherical_Bessel");

    double val;
    double x = r * lambda;
    if (k == 0)
    {
        if (x < 1.0e-3)
            val = -1 / x + 1 - x / 2.0 + pow(x, 2) / 6.0;
        else
            val = -exp(-x) / x;
    }
    else if (k == 2)
    {
        if (x < 1.0e-2)
            val = 3 / pow(x, 3) - 1 / (2 * x) + x / 8 - pow(x, 2) / 15.0 + pow(x, 3) / 48.0;
        else
            val = exp(-x) * (3 + 3 * x + pow(x, 2)) / pow(x, 3);
    }
    else if (k == 4)
    {
        if (x < 5.0e-1)
            val = -105 / pow(x, 5) + 15 / (2 * pow(x, 3)) - 3 / (8 * x) + x / 48 - pow(x, 3) / 384.0
                  + pow(x, 4) / 945.0;
        else
            val = -exp(-x) * (105 + 105 * x + 45 * pow(x, 2) + 10 * pow(x, 3) + pow(x, 4)) / pow(x, 5);
    }
    else if (k == 6)
    {
        if (x < 9.0e-1)
            val = 10395 / pow(x, 7) - 945 / (2 * pow(x, 5)) + 105 / (8 * pow(x, 3)) - 5 / (16 * x) + x / 128.0
                  - pow(x, 3) / 3840.0 + pow(x, 5) / 46080.0 - pow(x, 6) / 135135.0;
        else
            val = exp(-x)
                  * (10395 + 10395 * x + 4725 * pow(x, 2) + 1260 * pow(x, 3) + 210 * pow(x, 4) + 21 * pow(x, 5)
                     + pow(x, 6))
                  / pow(x, 7);
    }
    return val;
}

} // namespace ModuleDFTU