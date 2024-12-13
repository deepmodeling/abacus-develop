#include "module_base/timer.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_domain.h"
#include "module_parameter/parameter.h"

namespace LCAO_domain
{

void single_derivative(ForceStressArrays& fsr,
                       const LCAO_Orbitals& orb,
                       const TwoCenterBundle& two_center_bundle,
                       const Parallel_Orbitals& pv,
                       const UnitCell& ucell,
                       const int nspin,
                       const bool cal_stress,
                       const int iw1_all,
                       const int iw2_all,
                       const int m1,
                       const int m2,
                       const char& dtype,
                       const int T1,
                       const int L1,
                       const int N1,
                       const int T2,
                       const int L2,
                       const int N2,
                       const ModuleBase::Vector3<double>& dtau,
                       const ModuleBase::Vector3<double>& tau1,
                       const ModuleBase::Vector3<double>& tau2,
                       const int npol,
                       const int jj,
                       const int jj0,
                       const int kk,
                       const int kk0,
                       int& nnr,
                       int& total_nnr,
                       double* olm // output value
)
{

    const bool gamma_only_local = PARAM.globalv.gamma_only_local;

    // convert m (0,1,...2l) to M (-l, -l+1, ..., l-1, l)
    const int M1 = (m1 % 2 == 0) ? -m1 / 2 : (m1 + 1) / 2;
    const int M2 = (m2 % 2 == 0) ? -m2 / 2 : (m2 + 1) / 2;
    switch (dtype)
    {
    case 'S':
        two_center_bundle.overlap_orb->calculate(T1, L1, N1, M1, T2, L2, N2, M2, dtau * ucell.lat0, nullptr, olm);
        break;
    case 'T':
        two_center_bundle.kinetic_orb->calculate(T1, L1, N1, M1, T2, L2, N2, M2, dtau * ucell.lat0, nullptr, olm);
        break;
    default: // not supposed to happen
        ModuleBase::WARNING_QUIT("LCAO_domain::build_ST_new", "dtype must be S or T");
    }

    // condition 7: gamma only or multiple k
    if (gamma_only_local)
    {
        LCAO_domain::set_force(pv,
                               iw1_all,
                               iw2_all,
                               olm[0],
                               olm[1],
                               olm[2],
                               dtype,
                               fsr.DSloc_x,
                               fsr.DSloc_y,
                               fsr.DSloc_z,
                               fsr.DHloc_fixed_x,
                               fsr.DHloc_fixed_y,
                               fsr.DHloc_fixed_z);

        if (cal_stress)
        {
            LCAO_domain::set_stress(pv,
                                    iw1_all,
                                    iw2_all,
                                    olm[0],
                                    olm[1],
                                    olm[2],
                                    dtype,
                                    dtau,
                                    fsr.DSloc_11,
                                    fsr.DSloc_12,
                                    fsr.DSloc_13,
                                    fsr.DSloc_22,
                                    fsr.DSloc_23,
                                    fsr.DSloc_33,
                                    fsr.DHloc_fixed_11,
                                    fsr.DHloc_fixed_12,
                                    fsr.DHloc_fixed_13,
                                    fsr.DHloc_fixed_22,
                                    fsr.DHloc_fixed_23,
                                    fsr.DHloc_fixed_33);
        } // end stress
    }     // end gamma_only
    else  // condition 7, multiple k-points algorithm
    {
        // condition 8, S or T
        if (dtype == 'S')
        {
            // condition 9, nspin
            if (nspin == 1 || nspin == 2)
            {
                fsr.DSloc_Rx[nnr] = olm[0];
                fsr.DSloc_Ry[nnr] = olm[1];
                fsr.DSloc_Rz[nnr] = olm[2];
            }
            else if (nspin == 4)
            {
                int is = (jj - jj0 * npol) + (kk - kk0 * npol) * 2;
                if (is == 0) // is==3 is not needed in force calculation
                {
                    fsr.DSloc_Rx[nnr] = olm[0];
                    fsr.DSloc_Ry[nnr] = olm[1];
                    fsr.DSloc_Rz[nnr] = olm[2];
                }
                else
                {
                    fsr.DSloc_Rx[nnr] = 0.0;
                    fsr.DSloc_Ry[nnr] = 0.0;
                    fsr.DSloc_Rz[nnr] = 0.0;
                }
            }
            else
            {
                ModuleBase::WARNING_QUIT("LCAO_domain::build_ST_new", "nspin must be 1, 2 or 4");
            } // end condition 9, nspin

            if (cal_stress)
            {
                fsr.DH_r[nnr * 3] = dtau.x;
                fsr.DH_r[nnr * 3 + 1] = dtau.y;
                fsr.DH_r[nnr * 3 + 2] = dtau.z;
            }
        }
        else if (dtype == 'T') // condition 8, S or T
        {
            // condtion 9, nspin
            if (nspin == 1 || nspin == 2)
            {
                fsr.DHloc_fixedR_x[nnr] = olm[0];
                fsr.DHloc_fixedR_y[nnr] = olm[1];
                fsr.DHloc_fixedR_z[nnr] = olm[2];
                if (cal_stress)
                {
                    fsr.stvnl11[nnr] = olm[0] * dtau.x;
                    fsr.stvnl12[nnr] = olm[0] * dtau.y;
                    fsr.stvnl13[nnr] = olm[0] * dtau.z;
                    fsr.stvnl22[nnr] = olm[1] * dtau.y;
                    fsr.stvnl23[nnr] = olm[1] * dtau.z;
                    fsr.stvnl33[nnr] = olm[2] * dtau.z;
                }
            }
            else if (nspin == 4) // condition 9
            {
                const int is = (jj - jj0 * npol) + (kk - kk0 * npol) * 2;
                // condition 10, details of nspin 4
                if (is == 0) // is==3 is not needed in force calculation
                {
                    fsr.DHloc_fixedR_x[nnr] = olm[0];
                    fsr.DHloc_fixedR_y[nnr] = olm[1];
                    fsr.DHloc_fixedR_z[nnr] = olm[2];
                    if (cal_stress)
                    {
                        fsr.stvnl11[nnr] = olm[0] * dtau.x;
                        fsr.stvnl12[nnr] = olm[0] * dtau.y;
                        fsr.stvnl13[nnr] = olm[0] * dtau.z;
                        fsr.stvnl22[nnr] = olm[1] * dtau.y;
                        fsr.stvnl23[nnr] = olm[1] * dtau.z;
                        fsr.stvnl33[nnr] = olm[2] * dtau.z;
                    }
                }
                else if (is == 1 || is == 2 || is == 3)
                {
                    fsr.DHloc_fixedR_x[nnr] = 0.0;
                    fsr.DHloc_fixedR_y[nnr] = 0.0;
                    fsr.DHloc_fixedR_z[nnr] = 0.0;
                    if (cal_stress)
                    {
                        fsr.stvnl11[nnr] = 0.0;
                        fsr.stvnl12[nnr] = 0.0;
                        fsr.stvnl13[nnr] = 0.0;
                        fsr.stvnl22[nnr] = 0.0;
                        fsr.stvnl23[nnr] = 0.0;
                        fsr.stvnl33[nnr] = 0.0;
                    }
                }
                else
                {
                    ModuleBase::WARNING_QUIT("LCAO_domain::build_ST_new", "is must be 0, 1, 2, 3");
                } // end condition 10, details of spin 4
            }
            else
            {
                ModuleBase::WARNING_QUIT("LCAO_domain::build_ST_new", "nspin must be 1, 2 or 4");
            } // end condition 9, nspin
        }     // end condition 8, S or T
        ++total_nnr;
        ++nnr;
    } // end condition 7, gamma or multiple k
}

void single_overlap(const LCAO_Orbitals& orb,
                    const TwoCenterBundle& two_center_bundle,
                    const Parallel_Orbitals& pv,
                    const UnitCell& ucell,
                    const int nspin,
                    const bool cal_stress,
                    const int iw1_all,
                    const int iw2_all,
                    const int m1,
                    const int m2,
                    const char& dtype,
                    const int T1,
                    const int L1,
                    const int N1,
                    const int T2,
                    const int L2,
                    const int N2,
                    const ModuleBase::Vector3<double>& dtau,
                    const ModuleBase::Vector3<double>& tau1,
                    const ModuleBase::Vector3<double>& tau2,
                    const int npol,
                    const int jj,
                    const int jj0,
                    const int kk,
                    const int kk0,
                    int& nnr,       // output value
                    int& total_nnr, // output value
                    double* olm,    // output value
                    double* HSloc   // output value
)
{
    const bool gamma_only_local = PARAM.globalv.gamma_only_local;

    // convert m (0,1,...2l) to M (-l, -l+1, ..., l-1, l)
    const int M1 = (m1 % 2 == 0) ? -m1 / 2 : (m1 + 1) / 2;
    const int M2 = (m2 % 2 == 0) ? -m2 / 2 : (m2 + 1) / 2;

    switch (dtype)
    {
    case 'S':
        two_center_bundle.overlap_orb->calculate(T1, L1, N1, M1, T2, L2, N2, M2, dtau * ucell.lat0, olm);
        break;
    case 'T':
        two_center_bundle.kinetic_orb->calculate(T1, L1, N1, M1, T2, L2, N2, M2, dtau * ucell.lat0, olm);
        break;
    default: // not supposed to happen
        ModuleBase::WARNING_QUIT("LCAO_domain::build_ST_new", "dtype must be S or T");
    }

    // When NSPIN == 4 , only diagonal term is calculated for T or S Operators
    // use olm1 to store the diagonal term with complex data type.
    std::complex<double> olm1[4];

    if (nspin == 4)
    {
        olm1[0] = std::complex<double>(olm[0], 0.0);
        olm1[1] = ModuleBase::ZERO;
        olm1[2] = ModuleBase::ZERO;
        olm1[3] = std::complex<double>(olm[0], 0.0);
    }

    // condition 7, gamma only or multiple k-points
    if (gamma_only_local)
    {
        // mohan add 2010-06-29
        // set the value in Hloc and Sloc
        // according to global2local_row and global2local_col
        // the last paramete: 1 for Sloc, 2 for Hloc
        // and 3 for Hloc_fixed.
        LCAO_domain::set_mat2d(iw1_all, iw2_all, olm[0], pv, HSloc);
    }
    else // condition 7, multiple k-points algorithm
    {
        // condition 8, S or T
        if (dtype == 'S')
        {
            // condition 9, nspin
            if (nspin == 1 || nspin == 2)
            {
                HSloc[nnr] = olm[0];
            }
            else
            {
                ModuleBase::WARNING_QUIT("LCAO_domain::build_ST_new", "nspin must be 1, 2 or 4");
            }
        }
        else if (dtype == 'T') // condition 8, S or T
        {
            // condition 9, nspin
            if (nspin == 1 || nspin == 2)
            {
                HSloc[nnr] = olm[0]; // <phi|kin|d phi>
            }
            else if (nspin == 4)
            { // only has diagonal term here.
            }
            else
            {
                ModuleBase::WARNING_QUIT("LCAO_domain::build_ST_new", "nspin must be 1, 2 or 4");
            }
        } // end condition 8, S or T
        ++total_nnr;
        ++nnr;
    } // end condition 7, gamma point or multiple k-points
}

void build_ST_new(ForceStressArrays& fsr,
                  const char& dtype,
                  const bool& calc_deri,
                  const bool& cal_stress,
                  const UnitCell& ucell,
                  const LCAO_Orbitals& orb,
                  const Parallel_Orbitals& pv,
                  const TwoCenterBundle& two_center_bundle,
                  const Grid_Driver* GridD,
                  double* HSloc,
                  bool cal_syns,
                  double dmax)
{
    ModuleBase::TITLE("LCAO_domain", "build_ST_new");
    ModuleBase::timer::tick("LCAO_domain", "build_ST_new");

    const int nspin = PARAM.inp.nspin;
    const int npol = PARAM.globalv.npol;
    const bool gamma_only_local = PARAM.globalv.gamma_only_local;

    int total_nnr = 0;
#ifdef _OPENMP
#pragma omp parallel reduction(+ : total_nnr)
    {
#endif
        // array to store data
        double olm[3] = {0.0, 0.0, 0.0};

        //\sum{T} e**{ikT} <\phi_{ia}|d\phi_{k\beta}(T)>
        ModuleBase::Vector3<double> tau1, tau2, dtau;
        ModuleBase::Vector3<double> dtau1, dtau2, tau0;

#ifdef _OPENMP
// use schedule(dynamic) for load balancing because adj_num is various
#pragma omp for schedule(dynamic)
#endif
        for (int iat1 = 0; iat1 < ucell.nat; iat1++) // loop 1, iat1
        {
            const int T1 = ucell.iat2it[iat1];
            const Atom* atom1 = &ucell.atoms[T1];
            const int I1 = ucell.iat2ia[iat1];

            tau1 = atom1->tau[I1];

            // GridD->Find_atom(tau1);
            AdjacentAtomInfo adjs;
            GridD->Find_atom(ucell, tau1, T1, I1, &adjs);
            // Record_adj.for_2d() may not called in some case
            int nnr = pv.nlocstart ? pv.nlocstart[iat1] : 0;

            if (cal_syns)
            {
                for (int k = 0; k < 3; k++)
                {
                    tau1[k] = tau1[k] - atom1->vel[I1][k] * PARAM.mdp.md_dt / ModuleBase::AU_to_FS / ucell.lat0;
                }
            }

            // loop 2, ad
            for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
            {
                const int T2 = adjs.ntype[ad];
                const int I2 = adjs.natom[ad];
                Atom* atom2 = &ucell.atoms[T2];
                tau2 = adjs.adjacent_tau[ad];
                dtau = tau2 - tau1;
                double distance = dtau.norm() * ucell.lat0;
                double rcut = orb.Phi[T1].getRcut() + orb.Phi[T2].getRcut();

                // condition 3, distance
                if (distance < rcut)
                {
                    int iw1_all = ucell.itiaiw2iwt(T1, I1, 0); // iw1_all = combined index (it, ia, iw)

                    // loop 4, jj
                    for (int jj = 0; jj < atom1->nw * npol; ++jj)
                    {
                        const int jj0 = jj / npol;
                        const int L1 = atom1->iw2l[jj0];
                        const int N1 = atom1->iw2n[jj0];
                        const int m1 = atom1->iw2m[jj0];

                        int iw2_all = ucell.itiaiw2iwt(T2, I2, 0); // zhengdy-soc

                        // loop 5, kk
                        for (int kk = 0; kk < atom2->nw * npol; ++kk)
                        {
                            const int kk0 = kk / npol;
                            const int L2 = atom2->iw2l[kk0];
                            const int N2 = atom2->iw2n[kk0];
                            const int m2 = atom2->iw2m[kk0];

                            // mohan add 2010-06-29
                            // this is in fact the same as in build_Nonlocal_mu,
                            // the difference is that here we use {L,N,m} for ccycle,
                            // build_Nonlocal_mu use atom.nw for cycle.
                            // so, here we use ParaO::in_this_processor,
                            // in build_Non... use global2local_row
                            // and global2local_col directly,
                            if (!pv.in_this_processor(iw1_all, iw2_all))
                            {
                                ++iw2_all;
                                continue;
                            }

                            olm[0] = 0.0;
                            olm[1] = 0.0;
                            olm[2] = 0.0;

                            // condition 6, not calculate the derivative
                            if (!calc_deri)
                            {
                                single_overlap(orb,
                                               two_center_bundle,
                                               pv,
                                               ucell,
                                               nspin,
                                               cal_stress,
                                               iw1_all,
                                               iw2_all,
                                               m1,
                                               m2,
                                               dtype,
                                               T1,
                                               L1,
                                               N1,
                                               T2,
                                               L2,
                                               N2,
                                               dtau,
                                               tau1,
                                               tau2,
                                               npol,
                                               jj,
                                               jj0,
                                               kk,
                                               kk0,
                                               nnr,
                                               total_nnr,
                                               olm,
                                               HSloc);
                            }
                            else // condition 6, calculate the derivative
                            {
                                single_derivative(fsr,
                                                  orb,
                                                  two_center_bundle,
                                                  pv,
                                                  ucell,
                                                  nspin,
                                                  cal_stress,
                                                  iw1_all,
                                                  iw2_all,
                                                  m1,
                                                  m2,
                                                  dtype,
                                                  T1,
                                                  L1,
                                                  N1,
                                                  T2,
                                                  L2,
                                                  N2,
                                                  dtau,
                                                  tau1,
                                                  tau2,
                                                  npol,
                                                  jj,
                                                  jj0,
                                                  kk,
                                                  kk0,
                                                  nnr,
                                                  total_nnr,
                                                  olm);
                            } // end condition 6, calc_deri
                            ++iw2_all;
                        } // end loop 5, kk
                        ++iw1_all;
                    } // end loop 4, jj
                }     // condition 3, distance
                else if (distance >= rcut && (!gamma_only_local))
                {
                    int start1 = ucell.itiaiw2iwt(T1, I1, 0);
                    int start2 = ucell.itiaiw2iwt(T2, I2, 0);

                    bool is_adj = false;
                    for (int ad0 = 0; ad0 < adjs.adj_num + 1; ++ad0)
                    {
                        const int T0 = adjs.ntype[ad0];
                        tau0 = adjs.adjacent_tau[ad0];
                        dtau1 = tau0 - tau1;
                        double distance1 = dtau1.norm() * ucell.lat0;
                        double rcut1 = orb.Phi[T1].getRcut() + ucell.infoNL.Beta[T0].get_rcut_max();
                        dtau2 = tau0 - tau2;
                        double distance2 = dtau2.norm() * ucell.lat0;
                        double rcut2 = orb.Phi[T2].getRcut() + ucell.infoNL.Beta[T0].get_rcut_max();
                        if (distance1 < rcut1 && distance2 < rcut2)
                        {
                            is_adj = true;
                            break;
                        }
                    } // ad0

                    if (is_adj)
                    {
                        for (int jj = 0; jj < atom1->nw * npol; ++jj)
                        {
                            const int mu = pv.global2local_row(start1 + jj);
                            if (mu < 0) {
                                continue;
}
                            for (int kk = 0; kk < atom2->nw * npol; ++kk)
                            {
                                const int nu = pv.global2local_col(start2 + kk);
                                if (nu < 0) {
                                    continue;
}
                                ++total_nnr;
                                ++nnr;
                            } // kk
                        }     // jj
                    }         // is_adj
                }             // distance, end condition 3
            }                 // end loop 2, ad
        }                     // end loop 1, iat1

#ifdef _OPENMP
    }
#endif

    if (!gamma_only_local)
    {
        if (total_nnr != pv.nnr)
        {
            std::cout << " nnr=" << total_nnr << " LNNR.nnr=" << pv.nnr << std::endl;
            GlobalV::ofs_running << " nnr=" << total_nnr << " LNNR.nnr=" << pv.nnr << std::endl;
            ModuleBase::WARNING_QUIT("LCAO_domain::build_ST_new", "nnr != LNNR.nnr");
        }
    }

    ModuleBase::timer::tick("LCAO_domain", "build_ST_new");
    return;
}

} // namespace LCAO_domain
