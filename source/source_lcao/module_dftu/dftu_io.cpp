#include "dftu.h"
#include "source_base/timer.h"
#include "source_io/module_parameter/parameter.h"
#include <cstring>
#include <iomanip>


void Plus_U::output(const UnitCell& ucell,
                    bool out_chg,
                    const std::string& global_out_dir,
                    int nspin,
                    int npol)
{
    ModuleBase::TITLE("Plus_U", "output");

    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    GlobalV::ofs_running << " | #DFT+U INFORMATION# |" << std::endl;
    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>" << std::endl;

    for (int T = 0; T < ucell.ntype; T++)
    {
        const int NL = ucell.atoms[T].nwl + 1;

        for (int L = 0; L < NL; L++)
        {
            const int N = ucell.atoms[T].l_nchi[L];

            if (L >= get_orbital_corr(T) && has_correlated_orbital(T))
            {
                if (L != get_orbital_corr(T))
                {
                    continue;
                }

                if (!Yukawa)
                {
                    GlobalV::ofs_running << " Type=" << T+1 << " L=" << L << " ORBITAL=" << 0
                                         << " U=" << this->U[T] * ModuleBase::Ry_to_eV << " eV" << std::endl;
                }
                else
                {
                    for (int n = 0; n < N; n++)
                    {
                        if (n != 0)
                        {
                            continue;
                        }
                        double Ueff = (this->U_Yukawa[T][L][n] - this->J_Yukawa[T][L][n]) * ModuleBase::Ry_to_eV;
                        GlobalV::ofs_running << " Type=" << T+1 << " L=" << L << "  ORBITAL=" << n
                                             << " U=" << this->U_Yukawa[T][L][n] * ModuleBase::Ry_to_eV << " eV"
                                             << " J=" << this->J_Yukawa[T][L][n] * ModuleBase::Ry_to_eV << " eV"
                                             << std::endl;
                    }
                }
            }
        }
    }

    GlobalV::ofs_running << " Local Occupation Matrices for each atom" << std::endl;
    this->write_occup_m(ucell, GlobalV::ofs_running, true, nspin, npol);

    // Write onsite.dm
    std::ofstream ofdftu;
    if (out_chg)
    {
        if (GlobalV::MY_RANK == 0)
        {
            ofdftu.open(global_out_dir + "onsite.dm");
        }
    }
    if (!ofdftu)
    {
        std::cout << " Plus_U::write_occup_m. Can't create file onsite.dm" << std::endl;
        exit(0);
    }
    this->write_occup_m(ucell, ofdftu, false, nspin, npol);
    ofdftu.close();

    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    GlobalV::ofs_running << " | # END DFT+U INFO    |" << std::endl;
    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>" << std::endl << std::endl;

    return;
}

// define the function calculate the eigenvalues of a matrix
std::vector<double> CalculateEigenvalues(std::vector<std::vector<double>>& A, int n);

void Plus_U::write_occup_m(const UnitCell& ucell,
                           std::ofstream& ofs,
                           bool diag,
                           int nspin,
                           int npol)
{
    ModuleBase::TITLE("Plus_U", "write_occup_m");

    if (GlobalV::MY_RANK != 0)
    {
        return;
    }

    for (int T = 0; T < ucell.ntype; T++)
    {
        if (!has_correlated_orbital(T))
        {
            continue;
        }
        const int NL = ucell.atoms[T].nwl + 1;
        const int LC = get_orbital_corr(T);

        for (int I = 0; I < ucell.atoms[T].na; I++)
        {
            const int iat = ucell.itia2iat(T, I);

            for (int l = 0; l < NL; l++)
            {
                if (l != get_orbital_corr(T))
                {
                    continue;
                }

                const int N = ucell.atoms[T].l_nchi[l];

                for (int n = 0; n < N; n++)
                {
                    // if(!Yukawa && n!=0) continue;
                    if (n != 0)
                    {
                        continue;
                    }

                    ofs << "\n Atom=" << iat+1;
                    ofs << " L=" << l;
                    ofs << " ORBITAL=" << n << std::endl;

                    if (nspin == 1 || nspin == 2)
                    {
                        double sum0[2];
                        for (int is = 0; is < 2; is++)
                        {
                            if (diag) // diagonalization for local occupation matrix and print the eigenvalues
                            {
                                std::vector<std::vector<double>> A(2 * l + 1, std::vector<double>(2 * l + 1));
                                for (int m0 = 0; m0 < 2 * l + 1; m0++)
                                {
                                    for (int m1 = 0; m1 < 2 * l + 1; m1++)
                                    {
                                        A[m0][m1] = locale[iat][l][n][is](m0, m1);
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
                            ofs << " spin=" << is+1 << std::endl;
                            for (int m0 = 0; m0 < 2 * l + 1; m0++)
                            {
                                for (int m1 = 0; m1 < 2 * l + 1; m1++)
                                {
                                    ofs << std::setw(12)
                                        << locale[iat][l][n][is](m0, m1);
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
                        if (diag) // diagonalization for local occupation matrix and print the eigenvalues
                        {         // output the eigenvalues for rho , mag_x, mag_y, mag_z
                            double sum0[4];
                            std::vector<std::vector<double>> A(2 * l + 1, std::vector<double>(2 * l + 1));
                            int index = 0;
                            for (int is = 0; is < 4; is++)
                            {
                                for (int m0 = 0; m0 < 2 * l + 1; m0++)
                                {
                                    for (int m1 = 0; m1 < 2 * l + 1; m1++)
                                    {
                                        A[m0][m1] = locale[iat][l][n][0].c[index];
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
                                                << locale[iat][l][n][0](m0_all, m1_all);
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

void Plus_U::read_occup_m(const UnitCell& ucell,
                          const std::string& fn,
                          const std::string& init_chg,
                          int nspin,
                          int npol)
{
    ModuleBase::TITLE("Plus_U", "read_occup_m");

    if (GlobalV::MY_RANK != 0)
    {
        return;
    }

    std::ifstream ifdftu(fn.c_str(), std::ios::in);

    if (!ifdftu)
    {
        if (omc > 0)
        {
            std::cout
                << "Plus_U::read_occup_m. Can not find the file initial_onsite.dm . Please check your initial_onsite.dm"
                << std::endl;
        }
        else
        {
            if (init_chg == "file")
            {
                std::cout << "Plus_U::read_occup_m. Can not find the file onsite.dm . Please do scf calculation first"
                          << std::endl;
            }
        }
        exit(0);
    }

    ifdftu.clear();
    ifdftu.seekg(0);

    char word[20];

    int T = 0;
    int iat = 0;
    int spin = 0;
    int L = 0;
    int zeta = 0;

    ifdftu.rdstate();

    while (ifdftu.good())
    {
        ifdftu >> word;
        if (ifdftu.eof())
        {
            break;
        }

        if (strcmp("Atom=", word) == 0)
        {
            ifdftu >> iat;
            iat -= 1;
            ifdftu >> word;

            if (strcmp("L=", word) != 0)
            {
                std::cout << "WRONG IN READING LOCAL OCCUPATION NUMBER MATRIX FROM Plus_U FILE" << std::endl;
                exit(0);
            }
            ifdftu >> L;
            ifdftu >> word;

            if (strcmp("ORBITAL=", word) != 0)
            {
                std::cout << "WRONG IN READING LOCAL OCCUPATION NUMBER MATRIX FROM Plus_U FILE" << std::endl;
                exit(0);
            }
            ifdftu >> zeta;
            ifdftu.ignore(150, '\n');

            T = ucell.iat2it[iat];
            const int NL = ucell.atoms[T].nwl + 1;
            const int LC = get_orbital_corr(T);

            for (int l = 0; l < NL; l++)
            {
                if (l != get_orbital_corr(T))
                {
                    continue;
                }

                if (nspin == 1 || nspin == 2)
                {
                    for (int is = 0; is < 2; is++)
                    {
                        ifdftu >> word;
                        if (strcmp("spin=", word) == 0)
                        {
                            ifdftu >> spin;
                            spin -= 1;
                            ifdftu.ignore(150, '\n');

                            double value = 0.0;
                            for (int m0 = 0; m0 < 2 * L + 1; m0++)
                            {
                                for (int m1 = 0; m1 < 2 * L + 1; m1++)
                                {
                                    ifdftu >> value;
                                    locale[iat][L][zeta][spin](m0, m1) = value;
                                }
                                ifdftu.ignore(150, '\n');
                            }
                        }
                        else
                        {
                            std::cout << "WRONG IN READING LOCAL OCCUPATION NUMBER MATRIX FROM Plus_U FILE"
                                      << std::endl;
                            exit(0);
                        }
                    }
                }
                else if (nspin == 4) // SOC
                {
                    double value = 0.0;
                    for (int m0 = 0; m0 < 2 * L + 1; m0++)
                    {
                        for (int ipol0 = 0; ipol0 < npol; ipol0++)
                        {
                            const int m0_all = m0 + (2 * L + 1) * ipol0;

                            for (int m1 = 0; m1 < 2 * L + 1; m1++)
                            {
                                for (int ipol1 = 0; ipol1 < npol; ipol1++)
                                {
                                    int m1_all = m1 + (2 * L + 1) * ipol1;
                                    ifdftu >> value;
                                    locale[iat][L][zeta][0](m0_all, m1_all) = value;
                                }
                            }
                            ifdftu.ignore(150, '\n');
                        }
                    }
                }
            }
        }
        else
        {
            std::cout << "WRONG IN READING LOCAL OCCUPATION NUMBER MATRIX FROM Plus_U FILE" << std::endl;
            exit(0);
        }

        ifdftu.rdstate();

        if (ifdftu.eof() != 0)
        {
            break;
        }
    }

    return;
}

void Plus_U::local_occup_bcast(const UnitCell& ucell,
                               int nspin,
                               int npol)
{
    ModuleBase::TITLE("Plus_U", "local_occup_bcast");

    for (int T = 0; T < ucell.ntype; T++)
    {
        if (!has_correlated_orbital(T))
        {
            continue;
        }

        for (int I = 0; I < ucell.atoms[T].na; I++)
        {
            const int iat = ucell.itia2iat(T, I);
            const int L = get_orbital_corr(T);

            for (int l = 0; l <= ucell.atoms[T].nwl; l++)
            {
                if (l != get_orbital_corr(T))
                {
                    continue;
                }

                for (int n = 0; n < ucell.atoms[T].l_nchi[l]; n++)
                {
                    // if(!Yukawa && n!=0) continue;
                    if (n != 0)
                    {
                        continue;
                    }

                    if (nspin == 1 || nspin == 2)
                    {
                        for (int spin = 0; spin < 2; spin++)
                        {
                            for (int m0 = 0; m0 < 2 * l + 1; m0++)
                            {
                                for (int m1 = 0; m1 < 2 * l + 1; m1++)
                                {
#ifdef __MPI
                                    MPI_Bcast(&locale[iat][l][n][spin](m0, m1), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
                                }
                            }
                        }
                    }
                    else if (nspin == 4) // SOC
                    {
                        for (int m0 = 0; m0 < 2 * L + 1; m0++)
                        {
                            for (int ipol0 = 0; ipol0 < npol; ipol0++)
                            {
                                const int m0_all = m0 + (2 * L + 1) * ipol0;

                                for (int m1 = 0; m1 < 2 * L + 1; m1++)
                                {
                                    for (int ipol1 = 0; ipol1 < npol; ipol1++)
                                    {
                                        int m1_all = m1 + (2 * L + 1) * ipol1;
#ifdef __MPI
                                        MPI_Bcast(&locale[iat][l][n][0](m0_all, m1_all),
                                                  1,
                                                  MPI_DOUBLE,
                                                  0,
                                                  MPI_COMM_WORLD);
#endif
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return;
}

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
