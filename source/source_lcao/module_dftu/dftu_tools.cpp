#include "dftu_lcao.h"

#ifdef __LCAO
void Plus_U::pot_onsite_complex(const int spin, const bool new_occ_mat, std::complex<double>* pot_onsite, const int npol)
{
    ModuleBase::TITLE("Plus_U", "pot_onsite_complex");
    ModuleBase::GlobalFunc::ZEROS(pot_onsite, this->paraV->nloc);

    for (int it = 0; it < this->ucell->ntype; ++it)
    {
        if (this->orbital_corr[it] == -1) 
        {
            continue;
        }
        for (int ia = 0; ia < this->ucell->atoms[it].na; ia++)
        {
            const int iat = this->ucell->itia2iat(it, ia);
            for (int L = 0; L <= this->ucell->atoms[it].nwl; L++)
            {
                if (L != this->orbital_corr[it]) 
                {
                    continue;
                }

                for (int n = 0; n < this->ucell->atoms[it].l_nchi[L]; n++)
                {
                    if (n != 0) 
                    {
                        continue;
                    }

                    for (int m1 = 0; m1 < 2 * L + 1; m1++)
                    {
                        for (int ipol1 = 0; ipol1 < npol; ipol1++)
                        {
                            const int mu = this->paraV->global2local_row(this->iatlnmipol2iwt[iat][L][n][m1][ipol1]);
                            if (mu < 0) 
                            {
                                continue;
                            }

                            for (int m2 = 0; m2 < 2 * L + 1; m2++)
                            {
                                for (int ipol2 = 0; ipol2 < npol; ipol2++)
                                {
                                    const int nu
                                        = this->paraV->global2local_col(this->iatlnmipol2iwt[iat][L][n][m2][ipol2]);
                                    if (nu < 0) 
                                    {
                                        continue;
                                    }
                                    int m1_all = m1 + (2 * L + 1) * ipol1;
                                    int m2_all = m2 + (2 * L + 1) * ipol2;
                                    double val = get_onebody_eff_pot(it, iat, L, n, spin, m1_all, m2_all, new_occ_mat);
                                    pot_onsite[nu * this->paraV->nrow + mu] = std::complex<double>(val, 0.0);
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

void Plus_U::pot_onsite_real(const int spin, const bool new_occ_mat, double* pot_onsite, const int npol)
{
    ModuleBase::TITLE("Plus_U", "pot_onsite_real");
    ModuleBase::GlobalFunc::ZEROS(pot_onsite, this->paraV->nloc);

    for (int it = 0; it < this->ucell->ntype; ++it)
    {
        if (this->orbital_corr[it] == -1) 
        {
            continue;
        }
        for (int ia = 0; ia < this->ucell->atoms[it].na; ia++)
        {
            const int iat = this->ucell->itia2iat(it, ia);
            for (int L = 0; L <= this->ucell->atoms[it].nwl; L++)
            {
                if (L != this->orbital_corr[it]) 
                {
                    continue;
                }

                for (int n = 0; n < this->ucell->atoms[it].l_nchi[L]; n++)
                {
                    if (n != 0) 
                    {
                        continue;
                    }
                    for (int m1 = 0; m1 < 2 * L + 1; m1++)
                    {
                        for (int ipol1 = 0; ipol1 < npol; ipol1++)
                        {
                            const int mu = this->paraV->global2local_row(this->iatlnmipol2iwt[iat][L][n][m1][ipol1]);
                            if (mu < 0) 
                            {
                                continue;
                            }
                            for (int m2 = 0; m2 < 2 * L + 1; m2++)
                            {
                                for (int ipol2 = 0; ipol2 < npol; ipol2++)
                                {
                                    const int nu
                                        = this->paraV->global2local_col(this->iatlnmipol2iwt[iat][L][n][m2][ipol2]);
                                    if (nu < 0) 
                                    {
                                        continue;
                                    }

                                    int m1_all = m1 + (2 * L + 1) * ipol1;
                                    int m2_all = m2 + (2 * L + 1) * ipol2;

                                    pot_onsite[nu * this->paraV->nrow + mu]
                                        = this->get_onebody_eff_pot(it, iat, L, n, spin, m1_all, m2_all, new_occ_mat);

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

double Plus_U::get_onebody_eff_pot(const int T,
                                 const int iat,
                                 const int L,
                                 const int N,
                                 const int spin,
                                 const int m0,
                                 const int m1,
                                 const bool new_occ_mat)
{
    ModuleBase::TITLE("Plus_U", "get_onebody_eff_pot");

    double pot_onsite = 0.0;

    switch (cal_type)
    {
    case 1: // rotationally invarient formalism and FLL double counting

        break;

    case 2: // rotationally invarient formalism and AMF double counting

        break;

    case 3: // simplified formalism and FLL double counting
        if (new_occ_mat)
        {
            if (use_yukawa_)
            {
                if (m0 == m1) 
                {
                    pot_onsite = (this->U_Yukawa[T][L][N] - this->J_Yukawa[T][L][N])
                         * (0.5 - this->occ_mat[iat][L][N][spin](m0, m1));
                } else {
                    pot_onsite = -(this->U_Yukawa[T][L][N] - this->J_Yukawa[T][L][N]) * this->occ_mat[iat][L][N][spin](m0, m1);
                }
            }
            else
            {
                if (m0 == m1) {
                    pot_onsite = (this->u_current[T]) * (0.5 - this->occ_mat[iat][L][N][spin](m0, m1));
                } else {
                    pot_onsite = -(this->u_current[T]) * this->occ_mat[iat][L][N][spin](m0, m1);
                }
            }
        }
        else
        {
            if (use_yukawa_)
            {
                if (m0 == m1) {
                    pot_onsite = (this->U_Yukawa[T][L][N] - this->J_Yukawa[T][L][N])
                         * (0.5 - this->occ_mat_save[iat][L][N][spin](m0, m1));
                } else {
                    pot_onsite = -(this->U_Yukawa[T][L][N] - this->J_Yukawa[T][L][N])
                         * this->occ_mat_save[iat][L][N][spin](m0, m1);
                }
            }
            else
            {
                if (m0 == m1) {
                    pot_onsite = (this->u_current[T]) * (0.5 - this->occ_mat_save[iat][L][N][spin](m0, m1));
                } else {
                    pot_onsite = -(this->u_current[T]) * this->occ_mat_save[iat][L][N][spin](m0, m1);
                }
            }
        }

        break;

    case 4: // simplified formalism and AMF double counting

        break;
    }

    return pot_onsite;
}
#endif
