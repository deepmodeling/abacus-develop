#include "yukawa_screening.h"

#include "source_base/constants.h"
#include "source_base/parallel_reduce.h"
#include "source_base/tool_quit.h"
#include "source_base/tool_title.h"
#include "source_cell/unitcell.h"
#ifdef __LCAO
#include "source_basis/module_ao/orb_read.h"
#endif

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

void YukawaScreening::init(const UnitCell& cell,
                           const std::vector<int>& orbital_corr,
                           double yukawa_lambda_cfg)
{
    this->yukawa_lambda_cfg_ = yukawa_lambda_cfg;
    this->lambda_ = 0.0;
    this->orbital_corr_ = orbital_corr;

    this->Fk_.resize(cell.ntype);
    this->U_Yukawa_.resize(cell.ntype);
    this->J_Yukawa_.resize(cell.ntype);

    for (int it = 0; it < cell.ntype; it++)
    {
        const int NL = cell.atoms[it].nwl + 1;

        this->Fk_[it].resize(NL);
        this->U_Yukawa_[it].resize(NL);
        this->J_Yukawa_[it].resize(NL);

        for (int l = 0; l < NL; l++)
        {
            const int N = cell.atoms[it].l_nchi[l];

            this->Fk_[it][l].resize(N);
            for (int n = 0; n < N; n++)
            {
                this->Fk_[it][l][n].resize(l + 1, 0.0);
            }

            this->U_Yukawa_[it][l].resize(N, 0.0);
            this->J_Yukawa_[it][l].resize(N, 0.0);
        }
    }
}

void YukawaScreening::cal_lambda(double** rho, int nrxx, int nspin)
{
    ModuleBase::TITLE("YukawaScreening", "cal_lambda");

    if (this->yukawa_lambda_cfg_ > 0)
    {
        this->lambda_ = this->yukawa_lambda_cfg_;
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
        if (nspin == 4 && is > 0)
        {
            continue; // for non-collinear spin case, first spin contains the charge density
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
        ModuleBase::WARNING_QUIT("YukawaScreening::cal_lambda", oss.str());
    }

    this->lambda_ = val2 / val1;

    // rescaling
    this->lambda_ /= 1.6;
}

void YukawaScreening::cal_slater_Fk(const UnitCell& ucell, int L, int T, const LCAO_Orbitals* orb)
{
    ModuleBase::TITLE("YukawaScreening", "cal_slater_Fk");

#ifdef __LCAO
    const double lambda_val = this->lambda_;

    for (int chi = 0; chi < ucell.atoms[T].l_nchi[L]; chi++)
    {
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
                        bslval = spherical_Bessel(l, r0, lambda_val);
                        hnkval = spherical_Hankel(l, r1, lambda_val);
                    }
                    else // greater than
                    {
                        bslval = spherical_Bessel(l, r1, lambda_val);
                        hnkval = spherical_Hankel(l, r0, lambda_val);
                    }
                    this->Fk_[T][L][chi][k] -= (4 * k + 1) * lambda_val * pow(R_L0, 2) * bslval * hnkval
                                               * pow(R_L1, 2) * pow(r0, 2) * pow(r1, 2) * rab0 * rab1;
                }
            }
        }
    }
#else
    (void)ucell;
    (void)L;
    (void)T;
    (void)orb;
    ModuleBase::WARNING_QUIT("YukawaScreening::cal_slater_Fk",
                             "Slater integrals require numerical orbitals; compile with __LCAO");
#endif
}

void YukawaScreening::cal_slater_UJ(const UnitCell& ucell,
                                    double** rho,
                                    int nrxx,
                                    int nspin,
                                    const LCAO_Orbitals* orb)
{
    ModuleBase::TITLE("YukawaScreening", "cal_slater_UJ");

    this->cal_lambda(rho, nrxx, nspin);

    for (int it = 0; it < ucell.ntype; it++)
    {
        const int NL = ucell.atoms[it].nwl + 1;

        for (int l = 0; l < NL; l++)
        {
            int N = ucell.atoms[it].l_nchi[l];
            for (int n = 0; n < N; n++)
            {
                std::fill(this->Fk_[it][l][n].begin(), this->Fk_[it][l][n].end(), 0.0);
            }
        }
    }

    for (int T = 0; T < ucell.ntype; T++)
    {
        const int NL = ucell.atoms[T].nwl + 1;

        for (int L = 0; L < NL; L++)
        {
            if (L >= this->orbital_corr_[T] && this->orbital_corr_[T] != -1)
            {
                if (L != this->orbital_corr_[T])
                {
                    continue;
                }
                this->cal_slater_Fk(ucell, L, T, orb);

                if (L == 1)
                {
                    this->U_Yukawa_[T][L][0] = this->Fk_[T][L][0][0];
                    this->J_Yukawa_[T][L][0] = this->Fk_[T][L][0][1] / 5.0;
                }
                else if (L == 2)
                {
                    this->U_Yukawa_[T][L][0] = this->Fk_[T][L][0][0];
                    this->J_Yukawa_[T][L][0] = (this->Fk_[T][L][0][1] + this->Fk_[T][L][0][2]) / 14.0;
                }
                else if (L == 3)
                {
                    this->U_Yukawa_[T][L][0] = this->Fk_[T][L][0][0];
                    this->J_Yukawa_[T][L][0] = (286.0 * this->Fk_[T][L][0][1]
                                                + 195.0 * this->Fk_[T][L][0][2]
                                                + 250.0 * this->Fk_[T][L][0][3]) / 6435.0;
                }

                // Hartree to Rydeberg
                this->U_Yukawa_[T][L][0] *= 2.0;
                this->J_Yukawa_[T][L][0] *= 2.0;
            } // end if
        } // end L
    } // end T
}

double YukawaScreening::spherical_Bessel(const int k, const double r, const double lambda)
{
    double val = 0.0;
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

double YukawaScreening::spherical_Hankel(const int k, const double r, const double lambda)
{
    double val = 0.0;
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
