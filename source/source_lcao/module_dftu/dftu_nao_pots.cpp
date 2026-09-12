#include "source_pw/module_pwdft/dftu_base.h"
#include "dftu_nao_pots.h"

#include "source_base/global_function.h"
#include "source_base/tool_title.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_cell/unitcell.h"

void DFTU_LCAO::pot_onsite_complex(const Plus_U_Base& dftu,
                                   const UnitCell& ucell,
                                   const Parallel_Orbitals* pv,
                                   const int spin,
                                   const bool new_occ_mat,
                                   std::complex<double>* pot_onsite,
                                   const int npol)
{
    ModuleBase::TITLE("DFTU_LCAO", "pot_onsite_complex");
    ModuleBase::GlobalFunc::ZEROS(pot_onsite, pv->nloc);

    const auto& iatlnmipol2iwt = dftu.occmat().iatlnmipol2iwt();

    for (int it = 0; it < ucell.ntype; ++it)
    {
        if (dftu.get_l_channel(it) == -1)
        {
            continue;
        }
        for (int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
            const int iat = ucell.itia2iat(it, ia);
            for (int L = 0; L <= ucell.atoms[it].nwl; L++)
            {
                if (L != dftu.get_l_channel(it))
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

void DFTU_LCAO::pot_onsite_real(const Plus_U_Base& dftu,
                                const UnitCell& ucell,
                                const Parallel_Orbitals* pv,
                                const int spin,
                                const bool new_occ_mat,
                                double* pot_onsite,
                                const int npol)
{
    ModuleBase::TITLE("DFTU_LCAO", "pot_onsite_real");
    ModuleBase::GlobalFunc::ZEROS(pot_onsite, pv->nloc);

    const auto& iatlnmipol2iwt = dftu.occmat().iatlnmipol2iwt();

    for (int it = 0; it < ucell.ntype; ++it)
    {
        if (dftu.get_l_channel(it) == -1)
        {
            continue;
        }
        for (int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
            const int iat = ucell.itia2iat(it, ia);
            for (int L = 0; L <= ucell.atoms[it].nwl; L++)
            {
                if (L != dftu.get_l_channel(it))
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

double DFTU_LCAO::get_onsite_pot(const Plus_U_Base& dftu,
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

    switch (dftu.get_form())
    {
    case Plus_U_Base::UForm::lich_fll: // Lichtenstein (rotationally invariant) + FLL DC
        break;

    case Plus_U_Base::UForm::lich_amf: // Lichtenstein (rotationally invariant) + AMF DC
        break;

    case Plus_U_Base::UForm::dud_fll: // Dudarev (simplified) + FLL DC
        if (new_occ_mat)
        {
            if (dftu.use_yukawa())
            {
                if (m0 == m1)
                {
                    pot_onsite = (dftu.yukawa().get_U(T, L, N) - dftu.yukawa().get_J(T, L, N))
                                 * (0.5 - dftu.occmat().get(iat, L, N, spin, m0, m1));
                }
                else
                {
                    pot_onsite = -(dftu.yukawa().get_U(T, L, N) - dftu.yukawa().get_J(T, L, N))
                                 * dftu.occmat().get(iat, L, N, spin, m0, m1);
                }
            }
            else
            {
                if (m0 == m1)
                {
                    pot_onsite = dftu.get_u_current(T)
                                 * (0.5 - dftu.occmat().get(iat, L, N, spin, m0, m1));
                }
                else
                {
                    pot_onsite = -dftu.get_u_current(T)
                                 * dftu.occmat().get(iat, L, N, spin, m0, m1);
                }
            }
        }
        else
        {
            if (dftu.use_yukawa())
            {
                if (m0 == m1)
                {
                    pot_onsite = (dftu.yukawa().get_U(T, L, N) - dftu.yukawa().get_J(T, L, N))
                                 * (0.5 - dftu.occmat().get_save(iat, L, N, spin, m0, m1));
                }
                else
                {
                    pot_onsite = -(dftu.yukawa().get_U(T, L, N) - dftu.yukawa().get_J(T, L, N))
                                 * dftu.occmat().get_save(iat, L, N, spin, m0, m1);
                }
            }
            else
            {
                if (m0 == m1)
                {
                    pot_onsite = dftu.get_u_current(T)
                                 * (0.5 - dftu.occmat().get_save(iat, L, N, spin, m0, m1));
                }
                else
                {
                    pot_onsite = -dftu.get_u_current(T)
                                 * dftu.occmat().get_save(iat, L, N, spin, m0, m1);
                }
            }
        }

        break;
    }

    return pot_onsite;
}
