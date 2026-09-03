#include "dftu_nao.h"
#include "dftu_nao_energy.h"
#include "dftu_nao_pots.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_cell/unitcell.h"
#include "source_io/module_parameter/parameter.h"

#ifdef __LCAO
void DFTU_LCAO::cal_energy_correction(Plus_U& dftu, const UnitCell& ucell)
{
    ModuleBase::TITLE("DFTU_LCAO", "cal_energy_correction");
    ModuleBase::timer::start("DFTU_LCAO", "cal_energy_correction");
    if (!dftu.is_occ_mat_initialized())
    {
        ModuleBase::timer::end("DFTU_LCAO", "cal_energy_correction");
        return;
    }

    // mohan update 20251106
    double energy_u = 0.0;

    double energy_dc = 0.0;

    // read from the global PARAM.inp.nspin instead of a Plus_U member;
    // the member indirection is being removed during the refactor
    const int nspin = PARAM.inp.nspin;
    const int npol = ucell.get_npol();

    for (int T = 0; T < ucell.ntype; T++)
    {
        const int NL = ucell.atoms[T].nwl + 1;
        const int LC = dftu.get_orbital_corr(T);
        for (int I = 0; I < ucell.atoms[T].na; I++)
        {
            if (LC == -1)
            {
                continue;
            }

            const int iat = ucell.itia2iat(T, I);
            const int L = dftu.get_orbital_corr(T);

            for (int l = 0; l < NL; l++)
            {
                if (l != dftu.get_orbital_corr(T))
                {
                    continue;
                }

                const int N = ucell.atoms[T].l_nchi[l];

                const int m_tot = 2 * l + 1;

                // part 1: calculate the DFT+U energy correction
                for (int n = 0; n < N; n++)
                {
                    if (n != 0)
                    {
                        continue;
                    }

                    if (nspin == 1 || nspin == 2)
                    {
                        for (int spin = 0; spin < 2; spin++)
                        {
                            double nm_trace = 0.0;
                            double nm2_trace = 0.0;

                            for (int m0 = 0; m0 < 2 * l + 1; m0++)
                            {
                                nm_trace += dftu.occmat().get(iat, l, n, spin, m0, m0);
                                for (int m1 = 0; m1 < 2 * l + 1; m1++)
                                {
                                    nm2_trace += dftu.occmat().get(iat, l, n, spin, m0, m1)
                                                 * dftu.occmat().get(iat, l, n, spin, m1, m0);
                                }
                            }
                            if (dftu.use_yukawa())
                            {
                                energy_u += 0.5 * (dftu.yukawa().get_U(T, l, n) - dftu.yukawa().get_J(T, l, n))
                                            * (nm_trace - nm2_trace);
                            }
                            else
                            {
                                energy_u += 0.5 * dftu.get_u_current(T) * (nm_trace - nm2_trace);
                            }
                        }
                    }
                    else if (nspin == 4)
                    {
                        double nm_trace = 0.0;
                        double nm2_trace = 0.0;

                        for (int m0 = 0; m0 < 2 * l + 1; m0++)
                        {
                            for (int ipol0 = 0; ipol0 < npol; ipol0++)
                            {
                                const int m0_all = m0 + (2 * l + 1) * ipol0;
                                nm_trace += dftu.occmat().get(iat, l, n, 0, m0_all, m0_all);

                                for (int m1 = 0; m1 < 2 * l + 1; m1++)
                                {
                                    for (int ipol1 = 0; ipol1 < npol; ipol1++)
                                    {
                                        int m1_all = m1 + (2 * l + 1) * ipol1;

                                        nm2_trace += dftu.occmat().get(iat, l, n, 0, m0_all, m1_all)
                                                     * dftu.occmat().get(iat, l, n, 0, m1_all, m0_all);
                                    }
                                }
                            }
                        }
                        if (dftu.use_yukawa())
                        {
                            energy_u += 0.5 * (dftu.yukawa().get_U(T, l, n) - dftu.yukawa().get_J(T, l, n))
                              * (nm_trace - nm2_trace);
                        }
                        else
                        {
                            energy_u += 0.5 * dftu.get_u_current(T) * (nm_trace - nm2_trace);
                        }
                    }

                    for (int m1 = 0; m1 < 2 * l + 1; m1++)
                    {
                        for (int ipol1 = 0; ipol1 < npol; ipol1++)
                        {
                            const int m1_all = m1 + ipol1 * (2 * l + 1);
                            for (int m2 = 0; m2 < 2 * l + 1; m2++)
                            {
                                for (int ipol2 = 0; ipol2 < npol; ipol2++)
                                {
                                    const int m2_all = m2 + ipol2 * (2 * l + 1);

                                    if (nspin == 1 || nspin == 2)
                                    {
                                        for (int is = 0; is < 2; is++)
                                        {
                                            double pot_onsite = 0.0;
                                            pot_onsite = get_onsite_pot(dftu, T, iat, l, n, is, m1_all, m2_all, false);
                                            energy_dc += pot_onsite * dftu.occmat().get(iat, l, n, is, m1_all, m2_all);
                                        }
                                    }
                                    else if (nspin == 4)
                                    {
                                        double pot_onsite = 0.0;
                                        pot_onsite = get_onsite_pot(dftu, T, iat, l, n, 0, m1_all, m2_all, false);
                                        energy_dc += pot_onsite * dftu.occmat().get(iat, l, n, 0, m1_all, m2_all);
                                    }
                                }
                            }
                        }
                    }
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
#endif
