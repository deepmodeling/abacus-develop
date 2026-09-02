#ifdef __LCAO
#include "source_base/constants.h"
#include "source_base/global_function.h"
#include "dftu_lcao.h"
#include "dftu_yukawa.h"
#include "source_io/module_parameter/parameter.h"

#include <cmath>
#include <limits>
#include <sstream>


void DFTU_LCAO::cal_yukawa_lambda(Plus_U& dftu, double** rho, const int& nrxx)
{
    ModuleBase::TITLE("DFTU_LCAO", "cal_yukawa_lambda");

    // read from the global PARAM.inp.nspin instead of a Plus_U member;
    // the member indirection is being removed during the refactor
    const int nspin = PARAM.inp.nspin;

    if (dftu.get_yukawa_lambda() > 0)
    {
        dftu.set_lambda(dftu.get_yukawa_lambda());
        return;
    }

    // Track sum_rho, sum_rho_lambda, and the global min of rho. min_rho is
    // kept only for the diagnostic message: negative grid-point values are
    // expected (atomic-density superposition on a coarse FFT grid oscillates
    // near atom boundaries) and are clamped to 0 before being raised to the
    // 1/6 power below, so they no longer produce NaN here. The remaining
    // failure mode is sum_rho == 0 (rho not populated, or grid/normalization
    // mismatch), which still makes lambda = sum_rho_lambda/sum_rho = 0/0.
    double sum_rho = 0.0;
    double sum_rho_lambda = 0.0;
    double min_rho = std::numeric_limits<double>::max();
    for (int is = 0; is < nspin; is++)
    {
        if(nspin == 4 && is > 0)
        {
            continue;// for non-collinear spin case, first spin contains the charge density
        }
        for (int ir = 0; ir < nrxx; ir++)
        {
            const double rho_ir = rho[is][ir];
            if (rho_ir < min_rho)
            {
                min_rho = rho_ir;
            }
            sum_rho += rho_ir;

            // Clamp rho to 0 before pow(3*rho/PI, 1/6): pow(negative, 1/6)
            // returns NaN in C++ and would silently poison sum_rho_lambda,
            // lambda, U/J and the Hamiltonian. Negative rho here is a
            // numerical artifact of superposing atomic densities onto the
            // FFT grid; the physical charge density is non-negative, and
            // the contribution of these near-boundary points to the
            // Thomas-Fermi-like integral is negligible.
            const double rho_safe = std::max(rho_ir, 0.0);
            const double lambda_ir = 2 * pow(3 * rho_safe / ModuleBase::PI, (double)1.0 / 6.0);
            sum_rho_lambda += lambda_ir * rho_ir;
        }
    }

    double val1 = 0.0;
    double val2 = 0.0;
    double min_rho_global = min_rho;

#ifdef __MPI
    MPI_Allreduce(&sum_rho, &val1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&sum_rho_lambda, &val2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&min_rho, &min_rho_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#else
    val1 = sum_rho;
    val2 = sum_rho_lambda;
#endif

    // Fail loud instead of writing NaN into lambda -> U/J -> Hamiltonian.
    // After the clamp above, NaN can only slip in if sum_rho == 0 (0/0);
    // min_rho is still reported for diagnosis when this fires.
    if (!std::isfinite(val1) || !std::isfinite(val2) || val1 <= 0.0)
    {
        std::ostringstream oss;
        oss << "rho sanity check failed: sum_rho=" << val1
            << " sum_rho_lambda=" << val2
            << " min_rho=" << min_rho_global
            << " (need finite sum_rho > 0); nspin=" << nspin
            << " nrxx=" << nrxx;
        ModuleBase::WARNING_QUIT("DFTU_LCAO::cal_yukawa_lambda", oss.str());
    }

    dftu.set_lambda(val2 / val1);

    // rescaling
    dftu.set_lambda(dftu.get_lambda() / 1.6);

    return;
}

void DFTU_LCAO::cal_slater_Fk(Plus_U& dftu, const UnitCell& ucell, const int L, const int T)
{
    ModuleBase::TITLE("DFTU_LCAO", "cal_slater_Fk");

    if (dftu.use_yukawa())
    {
        const LCAO_Orbitals* orb = dftu.get_ptr_orb();
        const double lambda_val = dftu.get_lambda();
        auto& Fk = dftu.get_Fk_data();

        for (int chi = 0; chi < ucell.atoms[T].l_nchi[L]; chi++)
        {
            //    if(chi!=0) continue;
            const int mesh = orb->Phi[T].PhiLN(L, chi).getNr();

            for (int k = 0; k <= L; k++)
            {
                for (int ir0 = 1; ir0 < mesh; ir0++)
                {
                    double r0 = orb->Phi[T].PhiLN(L, chi).getRadial(ir0);
                    const double rab0 = orb->Phi[T].PhiLN(L, chi).getRab(ir0);
                    const double R_L0 = orb->Phi[T].PhiLN(L, chi).getPsi(ir0);

                    for (int ir1 = 1; ir1 < mesh; ir1++)
                    {
                        double bslval, hnkval;
                        double r1 = orb->Phi[T].PhiLN(L, chi).getRadial(ir1);
                        const double rab1 = orb->Phi[T].PhiLN(L, chi).getRab(ir1);
                        const double R_L1 = orb->Phi[T].PhiLN(L, chi).getPsi(ir1);

                        int l = 2 * k;
                        if (ir0 < ir1) // less than
                        {
                            bslval = DFTU_LCAO::spherical_Bessel(l, r0, lambda_val);
                            hnkval = DFTU_LCAO::spherical_Hankel(l, r1, lambda_val);
                        }
                        else // greater than
                        {
                            bslval = DFTU_LCAO::spherical_Bessel(l, r1, lambda_val);
                            hnkval = DFTU_LCAO::spherical_Hankel(l, r0, lambda_val);
                        }
                        Fk[T][L][chi][k] -= (4 * k + 1) * lambda_val * pow(R_L0, 2) * bslval * hnkval * pow(R_L1, 2)
                                                  * pow(r0, 2) * pow(r1, 2) * rab0 * rab1;
                    }
                }
            }
        }
    }

    return;
}

void DFTU_LCAO::cal_slater_UJ(Plus_U& dftu, const UnitCell& ucell, double** rho, const int& nrxx)
{
    ModuleBase::TITLE("DFTU_LCAO", "cal_slater_UJ");
    if (!dftu.use_yukawa())
    {
        return;
    }

    cal_yukawa_lambda(dftu, rho, nrxx);

    auto& Fk = dftu.get_Fk_data();

    for (int it = 0; it < ucell.ntype; it++)
    {
        const int NL = ucell.atoms[it].nwl + 1;

        for (int l = 0; l < NL; l++)
        {
            int N = ucell.atoms[it].l_nchi[l];
            for (int n = 0; n < N; n++)
            {
                ModuleBase::GlobalFunc::ZEROS(ModuleBase::GlobalFunc::VECTOR_TO_PTR(Fk[it][l][n]), l + 1);
            }
        }
    }

    for (int T = 0; T < ucell.ntype; T++)
    {
        const int NL = ucell.atoms[T].nwl + 1;

        for (int L = 0; L < NL; L++)
        {
            const int N = ucell.atoms[T].l_nchi[L];

            if (L >= dftu.get_orbital_corr(T) && dftu.get_orbital_corr(T) != -1)
            {
                if (L != dftu.get_orbital_corr(T))
                {
                    continue;
                }
                cal_slater_Fk(dftu, ucell, L, T);


                if( L == 1)
                {
                    dftu.set_U_Yukawa(T, L, 0, Fk[T][L][0][0]);
                    dftu.set_J_Yukawa(T, L, 0, Fk[T][L][0][1] / 5.0);
                }
                else if( L == 2)
                {
                    dftu.set_U_Yukawa(T, L, 0, Fk[T][L][0][0]);
                    dftu.set_J_Yukawa(T, L, 0, (Fk[T][L][0][1] + Fk[T][L][0][2]) / 14.0);
                }
                else if( L == 3)
                {
                    dftu.set_U_Yukawa(T, L, 0, Fk[T][L][0][0]);
                    dftu.set_J_Yukawa(T, L, 0, (286.0 * Fk[T][L][0][1] + 195.0 * Fk[T][L][0][2]
                                                + 250.0 * Fk[T][L][0][3])
                                                / 6435.0);
                }

                // Hartree to Rydeberg
                dftu.set_U_Yukawa(T, L, 0, dftu.get_U_Yukawa(T, L, 0) * 2.0);
                dftu.set_J_Yukawa(T, L, 0, dftu.get_J_Yukawa(T, L, 0) * 2.0);
                // update current U with calculated U-J from Slater integrals
                dftu.set_u_current(T, dftu.get_U_Yukawa(T, L, 0) - dftu.get_J_Yukawa(T, L, 0));
            } // end if
        } // end L
    } // end T

    return;
}

double DFTU_LCAO::spherical_Bessel(const int k, const double r, const double lambda)
{
    ModuleBase::TITLE("DFTU_LCAO", "spherical_Bessel");

    double val=0.0;
    double x = r * lambda;
    if (k == 0)
    {
        if (x < 1.0e-3)
        {
            val = 1 + pow(x, 2) / 6.0;
        }
        else
        {
            val = sinh(x) / x;
        }
    }
    else if (k == 2)
    {
        if (x < 1.0e-2)
        {
            val = -pow(x, 2) / 15.0 - pow(x, 4) / 210.0 - pow(x, 6) / 7560.0;
        }
        else
        {
            val = 3 * cosh(x) / pow(x, 2) + (-3 - pow(x, 2)) * sinh(x) / pow(x, 3);
        }
    }
    else if (k == 4)
    {
        if (x < 5.0e-1)
        {
            val = pow(x, 4) / 945.0 + pow(x, 6) / 20790.0 + pow(x, 8) / 1081080.0 + pow(x, 10) / 97297200.0;
        }
        else
        {
            val = -5 * (21 + 2 * pow(x, 2)) * cosh(x) / pow(x, 4)
                  + (105 + 45 * pow(x, 2) + pow(x, 4)) * sinh(x) / pow(x, 5);
        }
    }
    else if (k == 6)
    {
        if (x < 9.0e-1)
        {
            val = -pow(x, 6) / 135135.0 - pow(x, 8) / 4054050.0 - pow(x, 10) / 275675400.0;
        }
        else
        {
            val = 21 * (495 + 60 * pow(x, 2) + pow(x, 4)) * cosh(x) / pow(x, 6)
                  + (-10395 - 4725 * pow(x, 2) - 210 * pow(x, 4) - pow(x, 6)) * sinh(x) / pow(x, 7);
        }
    }
    return val;
}

double DFTU_LCAO::spherical_Hankel(const int k, const double r, const double lambda)
{
    ModuleBase::TITLE("DFTU_LCAO", "spherical_Hankel");

    double val=0.0;
    double x = r * lambda;
    if (k == 0)
    {
        if (x < 1.0e-3)
        {
            val = -1 / x + 1 - x / 2.0 + pow(x, 2) / 6.0;
        }
        else
        {
            val = -exp(-x) / x;
        }
    }
    else if (k == 2)
    {
        if (x < 1.0e-2)
        {
            val = 3 / pow(x, 3) - 1 / (2 * x) + x / 8 - pow(x, 2) / 15.0 + pow(x, 3) / 48.0;
        }
        else
        {
            val = exp(-x) * (3 + 3 * x + pow(x, 2)) / pow(x, 3);
        }
    }
    else if (k == 4)
    {
        if (x < 5.0e-1)
        {
            val = -105 / pow(x, 5) + 15 / (2 * pow(x, 3)) - 3 / (8 * x) + x / 48 - pow(x, 3) / 384.0
                  + pow(x, 4) / 945.0;
        }
        else
        {
            val = -exp(-x) * (105 + 105 * x + 45 * pow(x, 2) + 10 * pow(x, 3) + pow(x, 4)) / pow(x, 5);
        }
    }
    else if (k == 6)
    {
        if (x < 9.0e-1)
        {
            val = 10395 / pow(x, 7) - 945 / (2 * pow(x, 5)) + 105 / (8 * pow(x, 3)) - 5 / (16 * x) + x / 128.0
                  - pow(x, 3) / 3840.0 + pow(x, 5) / 46080.0 - pow(x, 6) / 135135.0;
        }
        else
        {
            val = exp(-x)
                  * (10395 + 10395 * x + 4725 * pow(x, 2) + 1260 * pow(x, 3) + 210 * pow(x, 4) + 21 * pow(x, 5)
                     + pow(x, 6))
                  / pow(x, 7);
        }
    }
    return val;
}

#endif
