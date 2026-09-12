#include "source_pw/module_pwdft/dftu_base.h"
#include "dftu_nao_energy.h"
#include "dftu_nao_pots.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_cell/unitcell.h"

namespace DFTU_LCAO {

/**
 * @brief Accumulate the DFT+U energy term (0.5 * U * (n - n^2)) for one
 *        (T, iat, l, n=0) channel in the collinear case (nspin=1 or 2).
 *        Returns the per-atom contribution to energy_u.
 */
double calc_energy_u_collinear(const Plus_U_Base& dftu,
                               int T,
                               int iat,
                               int l,
                               int n)
{
    double energy_u_local = 0.0;
    const int m_tot = 2 * l + 1;
    for (int spin = 0; spin < 2; spin++)
    {
        double nm_trace = 0.0;
        double nm2_trace = 0.0;

        for (int m0 = 0; m0 < m_tot; m0++)
        {
            nm_trace += dftu.occmat().get(iat, l, n, spin, m0, m0);
            for (int m1 = 0; m1 < m_tot; m1++)
            {
                nm2_trace += dftu.occmat().get(iat, l, n, spin, m0, m1)
                             * dftu.occmat().get(iat, l, n, spin, m1, m0);
            }
        }
        if (dftu.use_yukawa())
        {
            energy_u_local += 0.5 * (dftu.yukawa().get_U(T, l, n) - dftu.yukawa().get_J(T, l, n))
                              * (nm_trace - nm2_trace);
        }
        else
        {
            energy_u_local += 0.5 * dftu.get_u_current(T) * (nm_trace - nm2_trace);
        }
    }
    return energy_u_local;
}

/**
 * @brief Accumulate the DFT+U energy term for one (T, iat, l, n=0) channel
 *        in the noncollinear case (nspin=4). Returns the per-atom
 *        contribution to energy_u.
 */
double calc_energy_u_noncollinear(const Plus_U_Base& dftu,
                                  int T,
                                  int iat,
                                  int l,
                                  int n,
                                  int npol)
{
    double energy_u_local = 0.0;
    const int m_tot = 2 * l + 1;
    double nm_trace = 0.0;
    double nm2_trace = 0.0;

    for (int m0 = 0; m0 < m_tot; m0++)
    {
        for (int ipol0 = 0; ipol0 < npol; ipol0++)
        {
            const int m0_all = m0 + m_tot * ipol0;
            nm_trace += dftu.occmat().get(iat, l, n, 0, m0_all, m0_all);

            for (int m1 = 0; m1 < m_tot; m1++)
            {
                for (int ipol1 = 0; ipol1 < npol; ipol1++)
                {
                    const int m1_all = m1 + m_tot * ipol1;
                    nm2_trace += dftu.occmat().get(iat, l, n, 0, m0_all, m1_all)
                                 * dftu.occmat().get(iat, l, n, 0, m1_all, m0_all);
                }
            }
        }
    }
    if (dftu.use_yukawa())
    {
        energy_u_local += 0.5 * (dftu.yukawa().get_U(T, l, n) - dftu.yukawa().get_J(T, l, n))
                          * (nm_trace - nm2_trace);
    }
    else
    {
        energy_u_local += 0.5 * dftu.get_u_current(T) * (nm_trace - nm2_trace);
    }
    return energy_u_local;
}

/**
 * @brief Accumulate the double-counting correction energy_dc for one
 *        (T, iat, l, n=0) channel by summing onsite_pot * occ over the
 *        (m1, ipol1, m2, ipol2) grid. Dispatches on nspin to choose the
 *        spin loop count. Returns the per-atom contribution to energy_dc.
 */
double calc_energy_dc_block(const Plus_U_Base& dftu,
                            int T,
                            int iat,
                            int l,
                            int n,
                            int nspin,
                            int npol)
{
    double energy_dc_local = 0.0;
    const int m_tot = 2 * l + 1;
    for (int m1 = 0; m1 < m_tot; m1++)
    {
        for (int ipol1 = 0; ipol1 < npol; ipol1++)
        {
            const int m1_all = m1 + ipol1 * m_tot;
            for (int m2 = 0; m2 < m_tot; m2++)
            {
                for (int ipol2 = 0; ipol2 < npol; ipol2++)
                {
                    const int m2_all = m2 + ipol2 * m_tot;

                    if (nspin == 1 || nspin == 2)
                    {
                        for (int is = 0; is < 2; is++)
                        {
                            const double pot_onsite = get_onsite_pot(dftu, T, iat, l, n, is, m1_all, m2_all, false);
                            energy_dc_local += pot_onsite * dftu.occmat().get(iat, l, n, is, m1_all, m2_all);
                        }
                    }
                    else if (nspin == 4)
                    {
                        const double pot_onsite = get_onsite_pot(dftu, T, iat, l, n, 0, m1_all, m2_all, false);
                        energy_dc_local += pot_onsite * dftu.occmat().get(iat, l, n, 0, m1_all, m2_all);
                    }
                }
            }
        }
    }
    return energy_dc_local;
}

/**
 * @brief DFT+U energy correction with the double-counting term subtracted.
 *
 * Computes energy_u from occ_mat and the onsite potential, then writes the
 * result back to dftu via set_energy. The spin channel count is sourced by
 * the caller from the input parameter to avoid a PARAM read here.
 */
void cal_energy_correction(Plus_U_Base& dftu,
                           const UnitCell& ucell,
                           int nspin,
                           int npol)
{
    ModuleBase::TITLE("DFTU_LCAO", "cal_energy_correction");
    ModuleBase::timer::start("DFTU_LCAO", "cal_energy_correction");
    if (!dftu.is_occmat_ready())
    {
        ModuleBase::timer::end("DFTU_LCAO", "cal_energy_correction");
        return;
    }

    double energy_u = 0.0;
    double energy_dc = 0.0;

    for (int T = 0; T < ucell.ntype; T++)
    {
        const int NL = ucell.atoms[T].nwl + 1;
        const int LC = dftu.get_l_channel(T);
        if (LC == -1)
        {
            continue;
        }

        for (int I = 0; I < ucell.atoms[T].na; I++)
        {
            const int iat = ucell.itia2iat(T, I);
            for (int l = 0; l < NL; l++)
            {
                if (l != dftu.get_l_channel(T))
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

                    // part 1: U-term contribution
                    if (nspin == 1 || nspin == 2)
                    {
                        energy_u += calc_energy_u_collinear(dftu, T, iat, l, n);
                    }
                    else if (nspin == 4)
                    {
                        energy_u += calc_energy_u_noncollinear(dftu, T, iat, l, n, npol);
                    }

                    // part 2: double-counting correction
                    energy_dc += calc_energy_dc_block(dftu, T, iat, l, n, nspin, npol);
                } // end n
            }     // end L
        }         // end I
    }             // end T

    // substract the double counting energy_dc included in band energy eband
    energy_u -= energy_dc;

    dftu.set_energy(energy_u);

    ModuleBase::timer::end("DFTU_LCAO", "cal_energy_correction");
    return;
}

} // namespace DFTU_LCAO
