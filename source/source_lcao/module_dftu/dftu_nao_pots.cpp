#include "source_pw/module_pwdft/dftu_base.h"
#include "dftu_nao_pots.h"

#include "source_base/global_function.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_cell/unitcell.h"

#include <algorithm>

/// On-site potential and Hubbard energy for one correlated shell:
///   pot_onsite(m,m') = U_eff * (0.5 * delta_{m,m'} - occ(m,m'))
///   EU = (U_eff / 2) * sum_{m,m'} occ(m,m') * (delta_{m,m'} - occ(m',m))
void DFTU_LCAO::cal_pot_onsite(const std::vector<double>& occ, const int m_size, const double u_value,
                               double* pot_onsite, double& eu)
{
    int spin_fold = occ.size() / m_size / m_size;
    if (spin_fold < 4) {
        for (int is = 0; is < spin_fold; ++is)
        {
            int start = is * m_size * m_size;
            for (int m1 = 0; m1 < m_size; m1++)
            {
                for (int m2 = 0; m2 < m_size; m2++)
                {
                    pot_onsite[start + m1 * m_size + m2] = u_value * (0.5 * (m1 == m2) - occ[start + m2 * m_size + m1]);
                    eu += u_value * 0.5 * occ[start + m2 * m_size + m1] * occ[start + m1 * m_size + m2];
                }
            }
        }
    } else
    {
        for (int m1 = 0; m1 < m_size; m1++)
        {
            for (int m2 = 0; m2 < m_size; m2++)
            {
                pot_onsite[m1 * m_size + m2] = u_value * (1.0 * (m1 == m2) - occ[m2 * m_size + m1]);
                eu += u_value * 0.25 * occ[m2 * m_size + m1] * occ[m1 * m_size + m2];
            }
        }
        for (int is = 1; is < spin_fold; ++is)
        {
            int start = is * m_size * m_size;
            for (int m1 = 0; m1 < m_size; m1++)
            {
                for (int m2 = 0; m2 < m_size; m2++)
                {
                    pot_onsite[start + m1 * m_size + m2] = u_value * (0 - occ[start + m2 * m_size + m1]);
                    eu += u_value * 0.25 * occ[start + m2 * m_size + m1] * occ[start + m1 * m_size + m2];
                }
            }
        }
    }
}

template <typename T>
void DFTU_LCAO::cal_pot_onsite(const Plus_U_Base& dftu,
                           const UnitCell& ucell,
                           const Parallel_Orbitals* pv,
                           const int spin,
                           const bool new_occ_mat,
                           T* pot_onsite)
{
    ModuleBase::TITLE("DFTU_LCAO", "pot_onsite");
    std::fill(pot_onsite, pot_onsite + pv->nloc, T(0));

    const int npol = dftu.occmat().npol();

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

                for (int m1 = 0; m1 < 2 * L + 1; m1++)
                {
                    for (int ipol1 = 0; ipol1 < npol; ipol1++)
                    {
                        const int mu = pv->global2local_row(dftu.occmat().corr_iwt(iat, L, m1, ipol1));
                        if (mu < 0)
                        {
                            continue;
                        }

                        for (int m2 = 0; m2 < 2 * L + 1; m2++)
                        {
                            for (int ipol2 = 0; ipol2 < npol; ipol2++)
                            {
                                const int nu
                                    = pv->global2local_col(dftu.occmat().corr_iwt(iat, L, m2, ipol2));
                                if (nu < 0)
                                {
                                    continue;
                                }
                                int m1_all = m1 + (2 * L + 1) * ipol1;
                                int m2_all = m2 + (2 * L + 1) * ipol2;
                                double val = get_onsite_pot(dftu, it, iat, L, spin,
                                                            m1_all, m2_all, new_occ_mat);
                                pot_onsite[nu * pv->nrow + mu] = static_cast<T>(val);
                            } // ipol2
                        } // m2
                    } // ipol1
                } // m1
            } // l
        } // ia
    } // it

    return;
}

// Explicit instantiation
template void DFTU_LCAO::cal_pot_onsite<double>(const Plus_U_Base& dftu,
                                                const UnitCell& ucell,
                                                const Parallel_Orbitals* pv,
                                                const int spin,
                                                const bool new_occ_mat,
                                                double* pot_onsite);
template void DFTU_LCAO::cal_pot_onsite<std::complex<double>>(const Plus_U_Base& dftu,
                                                              const UnitCell& ucell,
                                                              const Parallel_Orbitals* pv,
                                                              const int spin,
                                                              const bool new_occ_mat,
                                                              std::complex<double>* pot_onsite);

template <typename T>
void DFTU_LCAO::cal_pot_uterm(Plus_U_Base& dftu,
                              const UnitCell& ucell,
                              const Parallel_Orbitals* pv,
                              const int spin,
                              T* pot_uterm,
                              const T* sk)
{
    ModuleBase::TITLE("DFTU_LCAO", "cal_pot_uterm");
    ModuleBase::timer::start("DFTU_LCAO", "cal_pot_uterm");

    const int nlocal = pv->get_global_row_size();
    std::fill(pot_uterm, pot_uterm + pv->nloc, T(0));

    //=============================================================
    //   PART2: call pblas to calculate effective potential matrix
    //=============================================================
    const char transN = 'N';
    const char transT = 'T';
    const int one_int = 1;
    const T half = static_cast<T>(0.5);
    const T one = static_cast<T>(1.0);
    const T zero = static_cast<T>(0.0);

    std::vector<T> pot_onsite(pv->nloc);
    DFTU_LCAO::cal_pot_onsite(dftu, ucell, pv, spin, true, &pot_onsite[0]);

#ifdef __MPI
    ScalapackConnector::gemm(transN, transN,
            nlocal, nlocal, nlocal,
            half,
            ModuleBase::GlobalFunc::VECTOR_TO_PTR(pot_onsite), one_int, one_int, pv->desc,
            sk, one_int, one_int, pv->desc,
            zero,
            pot_uterm, one_int, one_int, pv->desc);
#endif

    for (int irc = 0; irc < pv->nloc; irc++)
    {
        pot_onsite[irc] = pot_uterm[irc];
    }

#ifdef __MPI
    ScalapackConnector::tranu(nlocal, nlocal,
            one,
            &pot_onsite[0], one_int, one_int, pv->desc,
            one,
            pot_uterm, one_int, one_int, pv->desc);
#endif

    ModuleBase::timer::end("DFTU_LCAO", "cal_pot_uterm");
    return;
}

// Explicit instantiation
template void DFTU_LCAO::cal_pot_uterm<double>(Plus_U_Base& dftu,
                                               const UnitCell& ucell,
                                               const Parallel_Orbitals* pv,
                                               const int spin,
                                               double* pot_uterm,
                                               const double* sk);
template void DFTU_LCAO::cal_pot_uterm<std::complex<double>>(Plus_U_Base& dftu,
                                                             const UnitCell& ucell,
                                                             const Parallel_Orbitals* pv,
                                                             const int spin,
                                                             std::complex<double>* pot_uterm,
                                                             const std::complex<double>* sk);

double DFTU_LCAO::get_onsite_pot(const Plus_U_Base& dftu,
                                 const int T,
                                 const int iat,
                                 const int L,
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
                    pot_onsite = (dftu.yukawa().get_U(T, L) - dftu.yukawa().get_J(T, L))
                                 * (0.5 - dftu.occmat().get(iat, L, spin, m0, m1));
                }
                else
                {
                    pot_onsite = -(dftu.yukawa().get_U(T, L) - dftu.yukawa().get_J(T, L))
                                 * dftu.occmat().get(iat, L, spin, m0, m1);
                }
            }
            else
            {
                if (m0 == m1)
                {
                    pot_onsite = dftu.get_u_current(T)
                                 * (0.5 - dftu.occmat().get(iat, L, spin, m0, m1));
                }
                else
                {
                    pot_onsite = -dftu.get_u_current(T)
                                 * dftu.occmat().get(iat, L, spin, m0, m1);
                }
            }
        }
        else
        {
            if (dftu.use_yukawa())
            {
                if (m0 == m1)
                {
                    pot_onsite = (dftu.yukawa().get_U(T, L) - dftu.yukawa().get_J(T, L))
                                 * (0.5 - dftu.occmat().get_save(iat, L, spin, m0, m1));
                }
                else
                {
                    pot_onsite = -(dftu.yukawa().get_U(T, L) - dftu.yukawa().get_J(T, L))
                                 * dftu.occmat().get_save(iat, L, spin, m0, m1);
                }
            }
            else
            {
                if (m0 == m1)
                {
                    pot_onsite = dftu.get_u_current(T)
                                 * (0.5 - dftu.occmat().get_save(iat, L, spin, m0, m1));
                }
                else
                {
                    pot_onsite = -dftu.get_u_current(T)
                                 * dftu.occmat().get_save(iat, L, spin, m0, m1);
                }
            }
        }

        break;
    }

    return pot_onsite;
}
