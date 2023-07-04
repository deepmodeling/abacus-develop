#include "dftu.h"
#include "module_base/timer.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

namespace ModuleDFTU
{

void DFTU::cal_VU_pot_mat_complex(const int spin, const bool newlocale, std::complex<double>* VU)
{
    ModuleBase::TITLE("DFTU", "cal_VU_pot_mat_complex");
    ModuleBase::GlobalFunc::ZEROS(VU, this->LM->ParaV->nloc);

    for (int it = 0; it < GlobalC::ucell.ntype; ++it)
    {
        if (INPUT.orbital_corr[it] == -1)
            continue;
        for (int ia = 0; ia < GlobalC::ucell.atoms[it].na; ia++)
        {
            const int iat = GlobalC::ucell.itia2iat(it, ia);
            for (int L = 0; L <= GlobalC::ucell.atoms[it].nwl; L++)
            {
                if (L != INPUT.orbital_corr[it])
                    continue;

                for (int n = 0; n < GlobalC::ucell.atoms[it].l_nchi[L]; n++)
                {
                    if (n != 0)
                        continue;

                    for (int m1 = 0; m1 < 2 * L + 1; m1++)
                    {
                        for (int ipol1 = 0; ipol1 < GlobalV::NPOL; ipol1++)
                        {
                            const int mu = this->LM->ParaV->global2local_row[this->iatlnmipol2iwt[iat][L][n][m1][ipol1]];
                            if (mu < 0)
                                continue;

                            for (int m2 = 0; m2 < 2 * L + 1; m2++)
                            {
                                for (int ipol2 = 0; ipol2 < GlobalV::NPOL; ipol2++)
                                {
                                    const int nu
                                        = this->LM->ParaV->global2local_col[this->iatlnmipol2iwt[iat][L][n][m2][ipol2]];
                                    if (nu < 0)
                                        continue;

                                    int m1_all = m1 + (2 * L + 1) * ipol1;
                                    int m2_all = m2 + (2 * L + 1) * ipol2;

                                    double val = get_onebody_eff_pot(it, iat, L, n, spin, m1_all, m2_all, newlocale);
                                    VU[nu * this->LM->ParaV->nrow + mu] = std::complex<double>(val, 0.0);
                                } // ipol2
                            } // m2
                        } // ipol1
                    } // m1
                } // n
            } // l
        } // ia
    } // it

    return;
}

void DFTU::cal_VU_pot_mat_real(const int spin, const bool newlocale, double* VU)
{
    ModuleBase::TITLE("DFTU", "cal_VU_pot_mat_real");
    ModuleBase::GlobalFunc::ZEROS(VU, this->LM->ParaV->nloc);

    for (int it = 0; it < GlobalC::ucell.ntype; ++it)
    {
        if (INPUT.orbital_corr[it] == -1)
            continue;
        for (int ia = 0; ia < GlobalC::ucell.atoms[it].na; ia++)
        {
            const int iat = GlobalC::ucell.itia2iat(it, ia);
            for (int L = 0; L <= GlobalC::ucell.atoms[it].nwl; L++)
            {
                if (L != INPUT.orbital_corr[it])
                    continue;

                for (int n = 0; n < GlobalC::ucell.atoms[it].l_nchi[L]; n++)
                {
                    if (n != 0)
                        continue;

                    for (int m1 = 0; m1 < 2 * L + 1; m1++)
                    {
                        for (int ipol1 = 0; ipol1 < GlobalV::NPOL; ipol1++)
                        {
                            const int mu = this->LM->ParaV->global2local_row[this->iatlnmipol2iwt[iat][L][n][m1][ipol1]];
                            if (mu < 0)
                                continue;
                            for (int m2 = 0; m2 < 2 * L + 1; m2++)
                            {
                                for (int ipol2 = 0; ipol2 < GlobalV::NPOL; ipol2++)
                                {
                                    const int nu
                                        = this->LM->ParaV->global2local_col[this->iatlnmipol2iwt[iat][L][n][m2][ipol2]];
                                    if (nu < 0)
                                        continue;

                                    int m1_all = m1 + (2 * L + 1) * ipol1;
                                    int m2_all = m2 + (2 * L + 1) * ipol2;

                                    VU[nu * this->LM->ParaV->nrow + mu]
                                        = this->get_onebody_eff_pot(it, iat, L, n, spin, m1_all, m2_all, newlocale);

                                } // ipol2
                            } // m2
                        } // ipol1
                    } // m1
                } // n
            } // l
        } // ia
    } // it

    return;
}

double DFTU::get_onebody_eff_pot(const int T,
                                 const int iat,
                                 const int L,
                                 const int N,
                                 const int spin,
                                 const int m0,
                                 const int m1,
                                 const bool newlocale)
{
    ModuleBase::TITLE("DFTU", "get_onebody_eff_pot");

    double VU = 0.0;

    switch (cal_type)
    {
    case 1: // rotationally invarient formalism and FLL double counting

        break;

    case 2: // rotationally invarient formalism and AMF double counting

        break;

    case 3: // simplified formalism and FLL double counting
        if (newlocale)
        {
            if (Yukawa)
            {
                if (m0 == m1)
                    VU = (this->U_Yukawa[T][L][N] - this->J_Yukawa[T][L][N])
                         * (0.5 - this->locale[iat][L][N][spin](m0, m1));
                else
                    VU = -(this->U_Yukawa[T][L][N] - this->J_Yukawa[T][L][N]) * this->locale[iat][L][N][spin](m0, m1);
            }
            else
            {
                if (m0 == m1)
                    VU = (this->U[T]) * (0.5 - this->locale[iat][L][N][spin](m0, m1));
                else
                    VU = -(this->U[T]) * this->locale[iat][L][N][spin](m0, m1);
            }
        }
        else
        {
            if (Yukawa)
            {
                if (m0 == m1)
                    VU = (this->U_Yukawa[T][L][N] - this->J_Yukawa[T][L][N])
                         * (0.5 - this->locale_save[iat][L][N][spin](m0, m1));
                else
                    VU = -(this->U_Yukawa[T][L][N] - this->J_Yukawa[T][L][N])
                         * this->locale_save[iat][L][N][spin](m0, m1);
            }
            else
            {
                if (m0 == m1)
                    VU = (this->U[T]) * (0.5 - this->locale_save[iat][L][N][spin](m0, m1));
                else
                    VU = -(this->U[T]) * this->locale_save[iat][L][N][spin](m0, m1);
            }
        }

        break;

    case 4: // simplified formalism and AMF double counting

        break;
    }

    return VU;
}
} // namespace ModuleDFTU