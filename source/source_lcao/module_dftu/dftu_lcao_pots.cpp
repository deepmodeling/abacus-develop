#include "dftu_lcao.h"
#include "dftu_lcao_pots.h"

#ifdef __LCAO
void DFTU_LCAO::pot_onsite_complex(const Plus_U& dftu,
                                   const UnitCell& ucell,
                                   const Parallel_Orbitals* pv,
                                   const int spin,
                                   const bool new_occ_mat,
                                   std::complex<double>* pot_onsite,
                                   const int npol)
{
    ModuleBase::TITLE("DFTU_LCAO", "pot_onsite_complex");
    ModuleBase::GlobalFunc::ZEROS(pot_onsite, pv->nloc);

    const auto& iatlnmipol2iwt = dftu.get_iatlnmipol2iwt();

    for (int it = 0; it < ucell.ntype; ++it)
    {
        if (dftu.get_orbital_corr(it) == -1)
        {
            continue;
        }
        for (int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
            const int iat = ucell.itia2iat(it, ia);
            for (int L = 0; L <= ucell.atoms[it].nwl; L++)
            {
                if (L != dftu.get_orbital_corr(it))
                {
                    continue;
                }

                for (int n = 0; n < ucell.atoms[it].l_nchi[L]; n++)
                {
                    if (n != 0)
                    {
                        continue;
                    }

                    for (int m1 = 0; m1 < 2 * L + 1; m1++)
                    {
                        for (int ipol1 = 0; ipol1 < npol; ipol1++)
                        {
                            const int mu = pv->global2local_row(iatlnmipol2iwt[iat][L][n][m1][ipol1]);
                            if (mu < 0)
                            {
                                continue;
                            }

                            for (int m2 = 0; m2 < 2 * L + 1; m2++)
                            {
                                for (int ipol2 = 0; ipol2 < npol; ipol2++)
                                {
                                    const int nu
                                        = pv->global2local_col(iatlnmipol2iwt[iat][L][n][m2][ipol2]);
                                    if (nu < 0)
                                    {
                                        continue;
                                    }
                                    int m1_all = m1 + (2 * L + 1) * ipol1;
                                    int m2_all = m2 + (2 * L + 1) * ipol2;
                                    double val = get_onsite_pot(dftu, it, iat, L, n, spin,
                                                                m1_all, m2_all, new_occ_mat);
                                    pot_onsite[nu * pv->nrow + mu] = std::complex<double>(val, 0.0);
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

void DFTU_LCAO::pot_onsite_real(const Plus_U& dftu,
                                const UnitCell& ucell,
                                const Parallel_Orbitals* pv,
                                const int spin,
                                const bool new_occ_mat,
                                double* pot_onsite,
                                const int npol)
{
    ModuleBase::TITLE("DFTU_LCAO", "pot_onsite_real");
    ModuleBase::GlobalFunc::ZEROS(pot_onsite, pv->nloc);

    const auto& iatlnmipol2iwt = dftu.get_iatlnmipol2iwt();

    for (int it = 0; it < ucell.ntype; ++it)
    {
        if (dftu.get_orbital_corr(it) == -1)
        {
            continue;
        }
        for (int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
            const int iat = ucell.itia2iat(it, ia);
            for (int L = 0; L <= ucell.atoms[it].nwl; L++)
            {
                if (L != dftu.get_orbital_corr(it))
                {
                    continue;
                }

                for (int n = 0; n < ucell.atoms[it].l_nchi[L]; n++)
                {
                    if (n != 0)
                    {
                        continue;
                    }
                    for (int m1 = 0; m1 < 2 * L + 1; m1++)
                    {
                        for (int ipol1 = 0; ipol1 < npol; ipol1++)
                        {
                            const int mu = pv->global2local_row(iatlnmipol2iwt[iat][L][n][m1][ipol1]);
                            if (mu < 0)
                            {
                                continue;
                            }
                            for (int m2 = 0; m2 < 2 * L + 1; m2++)
                            {
                                for (int ipol2 = 0; ipol2 < npol; ipol2++)
                                {
                                    const int nu
                                        = pv->global2local_col(iatlnmipol2iwt[iat][L][n][m2][ipol2]);
                                    if (nu < 0)
                                    {
                                        continue;
                                    }

                                    int m1_all = m1 + (2 * L + 1) * ipol1;
                                    int m2_all = m2 + (2 * L + 1) * ipol2;

                                    pot_onsite[nu * pv->nrow + mu]
                                        = get_onsite_pot(dftu, it, iat, L, n, spin,
                                                          m1_all, m2_all, new_occ_mat);

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

double DFTU_LCAO::get_onsite_pot(const Plus_U& dftu,
                                 const int T,
                                 const int iat,
                                 const int L,
                                 const int N,
                                 const int spin,
                                 const int m0,
                                 const int m1,
                                 const bool new_occ_mat)
{
    ModuleBase::TITLE("DFTU_LCAO", "get_onsite_pot");

    double pot_onsite = 0.0;

    switch (dftu.get_cal_type())
    {
    case 1: // rotationally invarient formalism and FLL double counting
        break;

    case 2: // rotationally invarient formalism and AMF double counting
        break;

    case 3: // simplified formalism and FLL double counting
        if (new_occ_mat)
        {
            if (dftu.use_yukawa())
            {
                if (m0 == m1)
                {
                    pot_onsite = (dftu.get_U_Yukawa(T, L, N) - dftu.get_J_Yukawa(T, L, N))
                                 * (0.5 - dftu.get_occ_mat(iat, L, N, spin, m0, m1));
                }
                else
                {
                    pot_onsite = -(dftu.get_U_Yukawa(T, L, N) - dftu.get_J_Yukawa(T, L, N))
                                 * dftu.get_occ_mat(iat, L, N, spin, m0, m1);
                }
            }
            else
            {
                if (m0 == m1)
                {
                    pot_onsite = dftu.get_u_current(T)
                                 * (0.5 - dftu.get_occ_mat(iat, L, N, spin, m0, m1));
                }
                else
                {
                    pot_onsite = -dftu.get_u_current(T)
                                 * dftu.get_occ_mat(iat, L, N, spin, m0, m1);
                }
            }
        }
        else
        {
            if (dftu.use_yukawa())
            {
                if (m0 == m1)
                {
                    pot_onsite = (dftu.get_U_Yukawa(T, L, N) - dftu.get_J_Yukawa(T, L, N))
                                 * (0.5 - dftu.get_occ_mat_save(iat, L, N, spin, m0, m1));
                }
                else
                {
                    pot_onsite = -(dftu.get_U_Yukawa(T, L, N) - dftu.get_J_Yukawa(T, L, N))
                                 * dftu.get_occ_mat_save(iat, L, N, spin, m0, m1);
                }
            }
            else
            {
                if (m0 == m1)
                {
                    pot_onsite = dftu.get_u_current(T)
                                 * (0.5 - dftu.get_occ_mat_save(iat, L, N, spin, m0, m1));
                }
                else
                {
                    pot_onsite = -dftu.get_u_current(T)
                                 * dftu.get_occ_mat_save(iat, L, N, spin, m0, m1);
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
