#include "source_pw/module_pwdft/dftu_output.h"

#include "source_pw/module_pwdft/dftu_base.h"
#include "source_base/constants.h"
#include "source_base/global_function.h"
#include "source_base/timer.h"
#include "source_io/module_parameter/parameter.h"

#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

// local inline helpers for eigenvalue calculation
// migrated from dftu_base.cpp, mohan 2025-11-08
inline void JacobiRotate(std::vector<std::vector<double>>& A, int p, int q, int n)
{
    if (std::abs(A[p][q]) > 1e-10)
    {
        double r = (A[q][q] - A[p][p]) / (2.0 * A[p][q]);
        double t = 0.0;
        if (r >= 0)
        {
            t = 1.0 / (r + sqrt(1.0 + r * r));
        }
        else
        {
            t = -1.0 / (-r + sqrt(1.0 + r * r));
        }
        double c = 1.0 / sqrt(1.0 + t * t);
        double s = t * c;

        A[p][p] -= t * A[p][q];
        A[q][q] += t * A[p][q];
        A[p][q] = A[q][p] = 0.0;

        for (int k = 0; k < n; k++)
        {
            if (k != p && k != q)
            {
                double Akp = c * A[k][p] - s * A[k][q];
                double Akq = s * A[k][p] + c * A[k][q];
                A[k][p] = A[p][k] = Akp;
                A[k][q] = A[q][k] = Akq;
            }
        }
    }
}

inline std::vector<double> CalculateEigenvalues(std::vector<std::vector<double>>& A, int n)
{
    std::vector<double> eigenvalues(n);
    while (true)
    {
        int p = 0, q = 1;
        for (int i = 0; i < n; i++)
        {
            for (int j = i + 1; j < n; j++)
            {
                if (std::abs(A[i][j]) > std::abs(A[p][q]))
                {
                    p = i;
                    q = j;
                }
            }
        }

        if (std::abs(A[p][q]) < 1e-10)
        {
            for (int i = 0; i < n; i++)
            {
                eigenvalues[i] = A[i][i];
            }
            break;
        }

        JacobiRotate(A, p, q, n);
    }
    return eigenvalues;
}


namespace dftu_io
{

void output(const Plus_U_Base& dftu,
            const UnitCell& ucell,
            bool out_chg,
            const std::string& global_out_dir,
            int nspin,
            int npol)
{
    ModuleBase::TITLE("dftu_io", "output");

    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    GlobalV::ofs_running << " | #DFT+U INFORMATION# |" << std::endl;
    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>" << std::endl;

    for (int T = 0; T < ucell.ntype; T++)
    {
        const int NL = ucell.atoms[T].nwl + 1;

        for (int L = 0; L < NL; L++)
        {
            const int N = ucell.atoms[T].l_nchi[L];

            if (L >= dftu.get_orbital_corr(T) && dftu.has_correlated_orbital(T))
            {
                if (L != dftu.get_orbital_corr(T))
                {
                    continue;
                }

                if (!dftu.use_yukawa)
                {
                    GlobalV::ofs_running << " Type=" << T+1 << " L=" << L << " ORBITAL=" << 0
                                         << " U=" << dftu.get_u_current(T) * ModuleBase::Ry_to_eV << " eV" << std::endl;
                }
                else
                {
                    for (int n = 0; n < N; n++)
                    {
                        if (n != 0)
                        {
                            continue;
                        }
                        double Ueff = (dftu.get_U_Yukawa(T, L, n) - dftu.get_J_Yukawa(T, L, n)) * ModuleBase::Ry_to_eV;
                        GlobalV::ofs_running << " Type=" << T+1 << " L=" << L << "  ORBITAL=" << n
                                             << " U=" << dftu.get_U_Yukawa(T, L, n) * ModuleBase::Ry_to_eV << " eV"
                                             << " J=" << dftu.get_J_Yukawa(T, L, n) * ModuleBase::Ry_to_eV << " eV"
                                             << std::endl;
                    }
                }
            }
        }
    }

    GlobalV::ofs_running << " Local Occupation Matrices for each atom" << std::endl;
    dftu_io::write_occup_m(dftu, ucell, GlobalV::ofs_running, true, nspin, npol);

    // Write dm_onsite.txt
    if (out_chg && GlobalV::MY_RANK == 0)
    {
        std::ofstream ofdftu;
        ofdftu.open(global_out_dir + "dm_onsite.txt");
        if (!ofdftu)
        {
            ModuleBase::WARNING_QUIT("dftu_io::output", "Can't create file dm_onsite.txt");
        }
        dftu_io::write_occup_m(dftu, ucell, ofdftu, false, nspin, npol);
        ofdftu.close();
    }

    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    GlobalV::ofs_running << " | # END DFT+U INFO    |" << std::endl;
    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>" << std::endl << std::endl;

    return;
}


void write_occup_m(const Plus_U_Base& dftu,
                   const UnitCell& ucell,
                   std::ofstream& ofs,
                   bool diag,
                   int nspin,
                   int npol)
{
    ModuleBase::TITLE("dftu_io", "write_occup_m");

    if (GlobalV::MY_RANK != 0)
    {
        return;
    }

    for (int T = 0; T < ucell.ntype; T++)
    {
        if (!dftu.has_correlated_orbital(T))
        {
            continue;
        }
        const int NL = ucell.atoms[T].nwl + 1;
        const int LC = dftu.get_orbital_corr(T);

        for (int I = 0; I < ucell.atoms[T].na; I++)
        {
            const int iat = ucell.itia2iat(T, I);

            for (int l = 0; l < NL; l++)
            {
                if (l != dftu.get_orbital_corr(T))
                {
                    continue;
                }

                const int N = ucell.atoms[T].l_nchi[l];

                for (int n = 0; n < N; n++)
                {
                    if (n != 0)
                    {
                        continue;
                    }

                    ofs << "\n Atom= " << iat+1;
                    ofs << " L= " << l;
                    ofs << " ORBITAL= " << n << std::endl;

                    if (nspin == 1 || nspin == 2)
                    {
                        double sum0[2];
                        for (int is = 0; is < 2; is++)
                        {
                            if (diag)
                            {
                                std::vector<std::vector<double>> A(2 * l + 1, std::vector<double>(2 * l + 1));
                                for (int m0 = 0; m0 < 2 * l + 1; m0++)
                                {
                                    for (int m1 = 0; m1 < 2 * l + 1; m1++)
                                    {
                                        A[m0][m1] = dftu.get_occ_mat(iat, l, n, is, m0, m1);
                                    }
                                }
                                std::vector<double> eigenvalues = CalculateEigenvalues(A, 2 * l + 1);
                                sum0[is] = 0.0;
                                ofs << " Eigenvalues for spin=" << is+1 << std::endl;
                                ofs << std::setprecision(8) << std::fixed;
                                for (int i = 0; i < 2 * l + 1; i++)
                                {
                                    ofs << std::setw(12) << eigenvalues[i];
                                    sum0[is] += eigenvalues[i];
                                }
                                ofs << std::endl;
                                ofs << " sum is " << std::setw(12) << sum0[is] << std::endl;
                            }
                            ofs << " spin= " << is+1 << std::endl;
                            ofs << std::setprecision(8) << std::fixed;
                            for (int m0 = 0; m0 < 2 * l + 1; m0++)
                            {
                                for (int m1 = 0; m1 < 2 * l + 1; m1++)
                                {
                                    ofs << std::setw(12)
                                        << dftu.get_occ_mat(iat, l, n, is, m0, m1);
                                }
                                ofs << std::endl;
                            }
                        }
                        if (diag)
                        {
                            ofs << std::setw(12) << std::setprecision(8)
                                << std::fixed << " Magnetism for atom " << iat+1 << ": " << sum0[0] - sum0[1]
                                << std::endl;
                        }
                    }
                    else if (nspin == 4) // SOC
                    {
                        if (diag)
                        {
                            double sum0[4];
                            std::vector<std::vector<double>> A(2 * l + 1, std::vector<double>(2 * l + 1));
                            int index = 0;
                            for (int is = 0; is < 4; is++)
                            {
                                for (int m0 = 0; m0 < 2 * l + 1; m0++)
                                {
                                    for (int m1 = 0; m1 < 2 * l + 1; m1++)
                                    {
                                        A[m0][m1] = dftu.get_occ_mat(iat, l, n, 0, m0, m1);
                                        index++;
                                    }
                                }
                                std::vector<double> eigenvalues = CalculateEigenvalues(A, 2 * l + 1);
                                sum0[is] = 0.0;
                                ofs << " Eigenvalues for is=" << is << std::endl;
                                ofs << std::setprecision(8) << std::fixed;
                                for (int i = 0; i < 2 * l + 1; i++)
                                {
                                    ofs << std::setw(12) << eigenvalues[i];
                                    sum0[is] += eigenvalues[i];
                                }
                                ofs << std::endl;
                                ofs << " sum is " << std::setw(12) << sum0[is] << std::endl;
                            }
                            ofs << std::setw(12) << std::setprecision(8)
                                << std::fixed << " Magnetism for atom " << iat + 1 << ": "
                                << sum0[1] << " " << sum0[2] << " " << sum0[3] << std::endl;
                        }
                        else
                        {
                            for (int m0 = 0; m0 < 2 * l + 1; m0++)
                            {
                                for (int ipol0 = 0; ipol0 < npol; ipol0++)
                                {
                                    const int m0_all = m0 + (2 * l + 1) * ipol0;

                                    for (int m1 = 0; m1 < 2 * l + 1; m1++)
                                    {
                                        for (int ipol1 = 0; ipol1 < npol; ipol1++)
                                        {
                                            int m1_all = m1 + (2 * l + 1) * ipol1;
                                            ofs << std::setw(12) << std::setprecision(8) << std::fixed
                                                << dftu.get_occ_mat(iat, l, n, 0, m0_all, m1_all);
                                        }
                                    }
                                    ofs << std::endl;
                                }
                            }
                        }
                    }
                } // n
            }     // l
        }         // I
    }             // T

    return;
}


} // namespace dftu_io
