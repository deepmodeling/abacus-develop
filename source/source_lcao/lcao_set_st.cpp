#include "source_base/timer.h"
#include "source_lcao/lcao_domain.h"
#include "source_io/module_parameter/parameter.h"

namespace LCAO_domain
{

// write the dS/dH-fixed contributions of one S-matrix element at slot nnr.
// defined statically here so the compiler can inline it at the hot inner-loop
// call site of single_deriv.
static void set_deriv_s(ForceStressArrays& fsr,
                        const int nspin,
                        const int is,
                        const int nnr,
                        const double* olm,
                        const ModuleBase::Vector3<double>& dtau,
                        const bool cal_stress)
{
    // write DSloc_R* only when allocated (skipped in cal_dS where only DHloc_fixedR_* is used)
    const bool write_dsloc_r = !fsr.DSloc_Rx.empty();
    // DHloc_fixedR_* and DH_r must be pre-allocated by the caller for the derivative path
    const bool write_dhloc = !fsr.DHloc_fixedR_x.empty();
    const bool write_dhr = cal_stress && !fsr.DH_r.empty();
    // condition 9, nspin
    if (nspin == 1 || nspin == 2)
    {
        if (write_dsloc_r)
        {
            fsr.DSloc_Rx[nnr] = olm[0];
            fsr.DSloc_Ry[nnr] = olm[1];
            fsr.DSloc_Rz[nnr] = olm[2];
        }
        if (write_dhloc)
        {
            fsr.DHloc_fixedR_x[nnr] = olm[0];
            fsr.DHloc_fixedR_y[nnr] = olm[1];
            fsr.DHloc_fixedR_z[nnr] = olm[2];
        }
    }
    else if (nspin == 4)
    {
        const double v0 = (is == 0) ? olm[0] : 0.0; // is==3 is not needed in force calculation
        const double v1 = (is == 0) ? olm[1] : 0.0;
        const double v2 = (is == 0) ? olm[2] : 0.0;
        if (write_dsloc_r)
        {
            fsr.DSloc_Rx[nnr] = v0;
            fsr.DSloc_Ry[nnr] = v1;
            fsr.DSloc_Rz[nnr] = v2;
        }
        if (write_dhloc)
        {
            fsr.DHloc_fixedR_x[nnr] = v0;
            fsr.DHloc_fixedR_y[nnr] = v1;
            fsr.DHloc_fixedR_z[nnr] = v2;
        }
    }
    else
    {
        ModuleBase::WARNING_QUIT("LCAO_domain::set_deriv_s", "nspin must be 1, 2 or 4");
    } // end condition 9, nspin

    if (write_dhr)
    {
        fsr.DH_r[nnr * 3] = dtau.x;
        fsr.DH_r[nnr * 3 + 1] = dtau.y;
        fsr.DH_r[nnr * 3 + 2] = dtau.z;
    }
}

// write the dH-fixed and stress contributions of one T-matrix element at slot nnr.
// kept in this translation unit so the compiler can inline it at the hot
// inner-loop call site of single_deriv.
static void set_deriv_t(ForceStressArrays& fsr,
                        const int nspin,
                        const int is,
                        const int nnr,
                        const double* olm,
                        const ModuleBase::Vector3<double>& dtau,
                        const bool cal_stress)
{
    // DHloc_fixedR_* and stress arrays must be pre-allocated by the caller for the derivative path
    const bool write_dhloc = !fsr.DHloc_fixedR_x.empty();
    const bool write_stress = cal_stress && !fsr.stvnl11.empty();
    // condition 9, nspin
    if (nspin == 1 || nspin == 2)
    {
        if (write_dhloc)
        {
            fsr.DHloc_fixedR_x[nnr] = olm[0];
            fsr.DHloc_fixedR_y[nnr] = olm[1];
            fsr.DHloc_fixedR_z[nnr] = olm[2];
        }
        if (write_stress)
        {
            fsr.stvnl11[nnr] = olm[0] * dtau.x;
            fsr.stvnl12[nnr] = olm[0] * dtau.y;
            fsr.stvnl13[nnr] = olm[0] * dtau.z;
            fsr.stvnl22[nnr] = olm[1] * dtau.y;
            fsr.stvnl23[nnr] = olm[1] * dtau.z;
            fsr.stvnl33[nnr] = olm[2] * dtau.z;
        }
    }
    else if (nspin == 4)
    {
        // condition 10, details of nspin 4
        if (is == 0) // is==3 is not needed in force calculation
        {
            if (write_dhloc)
            {
                fsr.DHloc_fixedR_x[nnr] = olm[0];
                fsr.DHloc_fixedR_y[nnr] = olm[1];
                fsr.DHloc_fixedR_z[nnr] = olm[2];
            }
            if (write_stress)
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
            if (write_dhloc)
            {
                fsr.DHloc_fixedR_x[nnr] = 0.0;
                fsr.DHloc_fixedR_y[nnr] = 0.0;
                fsr.DHloc_fixedR_z[nnr] = 0.0;
            }
            if (write_stress)
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
            ModuleBase::WARNING_QUIT("LCAO_domain::set_deriv_t", "is must be 0, 1, 2, 3");
        } // end condition 10, details of spin 4
    }
    else
    {
        ModuleBase::WARNING_QUIT("LCAO_domain::set_deriv_t", "nspin must be 1, 2 or 4");
    } // end condition 9, nspin
}

void single_deriv(const ST_env& env,
                  const ST_elem& e,
                  ForceStressArrays& fsr,
                  int& nnr,
                  int& total_nnr,
                  double* olm // output value
)
{

    const bool gamma_only_local = env.gamma_only_local;
    const int nspin = env.nspin;
    const int npol = env.npol;
    const bool cal_stress = env.cal_stress;
    const int iw1_all = e.iw1_all;
    const int iw2_all = e.iw2_all;
    const char dtype = e.dtype;
    const int m1 = e.m1;
    const int m2 = e.m2;
    const int t1 = e.t1;
    const int l1 = e.l1;
    const int n1 = e.n1;
    const int t2 = e.t2;
    const int l2 = e.l2;
    const int n2 = e.n2;
    const ModuleBase::Vector3<double>& dtau = e.dtau;
    const int jj = e.jj;
    const int jj0 = e.jj0;
    const int kk = e.kk;
    const int kk0 = e.kk0;

    // convert m (0,1,...2l) to mm (-l, -l+1, ..., l-1, l)
    const int mm1 = (m1 % 2 == 0) ? -m1 / 2 : (m1 + 1) / 2;
    const int mm2 = (m2 % 2 == 0) ? -m2 / 2 : (m2 + 1) / 2;
    switch (dtype)
    {
    case 'S':
        env.two_center_bundle.overlap_orb->calculate(t1, l1, n1, mm1, t2, l2, n2, mm2, dtau * env.ucell.lat0, nullptr, olm);
        break;
    case 'T':
        env.two_center_bundle.kinetic_orb->calculate(t1, l1, n1, mm1, t2, l2, n2, mm2, dtau * env.ucell.lat0, nullptr, olm);
        break;
    default: // not supposed to happen
        ModuleBase::WARNING_QUIT("LCAO_domain::single_deriv", "dtype must be S or T");
    }

    // condition 7: gamma only or multiple k
    if (gamma_only_local)
    {
        LCAO_domain::set_force(env.pv,
                               iw1_all,
                               iw2_all,
                               olm[0],
                               olm[1],
                               olm[2],
                               dtype,
                               fsr.DSloc_x.data(),
                               fsr.DSloc_y.data(),
                               fsr.DSloc_z.data(),
                               fsr.DHloc_fixed_x.data(),
                               fsr.DHloc_fixed_y.data(),
                               fsr.DHloc_fixed_z.data());
    }     // end gamma_only
    else  // condition 7, multiple k-points algorithm
    {
        const int is = (jj - jj0 * npol) + (kk - kk0 * npol) * 2;
        // condition 8, S or T
        if (dtype == 'S')
        {
            set_deriv_s(fsr, nspin, is, nnr, olm, dtau, cal_stress);
        }
        else if (dtype == 'T')
        {
            set_deriv_t(fsr, nspin, is, nnr, olm, dtau, cal_stress);
        }     // end condition 8, S or T
        ++total_nnr;
        ++nnr;
    } // end condition 7, gamma or multiple k
}

void single_overlap(const ST_env& env,
                    const ST_elem& e,
                    int& nnr,       // output value
                    int& total_nnr, // output value
                    double* olm,    // output value
                    double* HSloc   // output value
)
{
    const bool gamma_only_local = env.gamma_only_local;
    const int nspin = env.nspin;
    const int iw1_all = e.iw1_all;
    const int iw2_all = e.iw2_all;
    const char dtype = e.dtype;
    const int m1 = e.m1;
    const int m2 = e.m2;
    const int t1 = e.t1;
    const int l1 = e.l1;
    const int n1 = e.n1;
    const int t2 = e.t2;
    const int l2 = e.l2;
    const int n2 = e.n2;
    const ModuleBase::Vector3<double>& dtau = e.dtau;

    // convert m (0,1,...2l) to mm (-l, -l+1, ..., l-1, l)
    const int mm1 = (m1 % 2 == 0) ? -m1 / 2 : (m1 + 1) / 2;
    const int mm2 = (m2 % 2 == 0) ? -m2 / 2 : (m2 + 1) / 2;

    switch (dtype)
    {
    case 'S':
        env.two_center_bundle.overlap_orb->calculate(t1, l1, n1, mm1, t2, l2, n2, mm2, dtau * env.ucell.lat0, olm);
        break;
    case 'T':
        env.two_center_bundle.kinetic_orb->calculate(t1, l1, n1, mm1, t2, l2, n2, mm2, dtau * env.ucell.lat0, olm);
        break;
    default: // not supposed to happen
        ModuleBase::WARNING_QUIT("LCAO_domain::single_overlap", "dtype must be S or T");
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
        LCAO_domain::set_mat2d(iw1_all, iw2_all, olm[0], env.pv, HSloc);
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
                ModuleBase::WARNING_QUIT("LCAO_domain::single_overlap", "nspin must be 1, 2 or 4");
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
                ModuleBase::WARNING_QUIT("LCAO_domain::single_overlap", "nspin must be 1, 2 or 4");
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
    ModuleBase::timer::start("LCAO_domain", "build_ST_new");

    const int nspin = PARAM.inp.nspin;
    const int npol = PARAM.globalv.npol;
    const bool gamma_only_local = PARAM.globalv.gamma_only_local;

    // derivative path must provide the target buffers
    if (calc_deri && !gamma_only_local)
    {
        if (fsr.DHloc_fixedR_x.empty() || fsr.DHloc_fixedR_y.empty() || fsr.DHloc_fixedR_z.empty())
        {
            ModuleBase::WARNING_QUIT("LCAO_domain::build_ST_new",
                "DHloc_fixedR_x/y/z must be allocated when calc_deri=true in multi-k mode");
        }
        if (cal_stress && fsr.DH_r.empty())
        {
            ModuleBase::WARNING_QUIT("LCAO_domain::build_ST_new",
                "DH_r must be allocated when calc_deri=true and cal_stress=true in multi-k mode");
        }
    }

    // read-only environment shared by every element of this build
    const ST_env env{orb, two_center_bundle, pv, ucell, nspin, npol, cal_stress, gamma_only_local};

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
            const int t1 = ucell.iat2it[iat1];
            const Atom* atom1 = &ucell.atoms[t1];
            const int i1 = ucell.iat2ia[iat1];

            tau1 = atom1->tau[i1];

            // GridD->Find_atom(tau1);
            AdjacentAtomInfo adjs;
            GridD->Find_atom(ucell, tau1, t1, i1, &adjs);
            // Record_adj.for_2d() may not called in some case
            int nnr = 0;
            if (!pv.nlocstart.empty())
            {
                nnr = pv.nlocstart[iat1];
            }

            if (cal_syns)
            {
                for (int k = 0; k < 3; k++)
                {
                    tau1[k] = tau1[k] - atom1->vel[i1][k] * PARAM.mdp.md_dt / ModuleBase::AU_to_FS / ucell.lat0;
                }
            }

            // loop 2, ad
            for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
            {
                const int t2 = adjs.ntype[ad];
                const int i2 = adjs.natom[ad];
                Atom* atom2 = &ucell.atoms[t2];
                tau2 = adjs.adjacent_tau[ad];
                dtau = tau2 - tau1;
                double distance = dtau.norm() * ucell.lat0;
                double rcut = orb.Phi[t1].getRcut() + orb.Phi[t2].getRcut();

                // condition 3, distance
                if (distance < rcut)
                {
                    int iw1_all = ucell.itiaiw2iwt(t1, i1, 0); // iw1_all = combined index (it, ia, iw)

                    // loop 4, jj
                    for (int jj = 0; jj < atom1->nw * npol; ++jj)
                    {
                        const int jj0 = jj / npol;
                        const int l1 = atom1->iw2l[jj0];
                        const int n1 = atom1->iw2n[jj0];
                        const int m1 = atom1->iw2m[jj0];

                        int iw2_all = ucell.itiaiw2iwt(t2, i2, 0); // zhengdy-soc

                        // loop 5, kk
                        for (int kk = 0; kk < atom2->nw * npol; ++kk)
                        {
                            const int kk0 = kk / npol;
                            const int l2 = atom2->iw2l[kk0];
                            const int n2 = atom2->iw2n[kk0];
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

                            const ST_elem elem{dtype,
                                               iw1_all,
                                               iw2_all,
                                               m1,
                                               m2,
                                               t1,
                                               l1,
                                               n1,
                                               t2,
                                               l2,
                                               n2,
                                               dtau,
                                               jj,
                                               jj0,
                                               kk,
                                               kk0};

                            // condition 6, not calculate the derivative
                            if (!calc_deri)
                            {
                                single_overlap(env, elem, nnr, total_nnr, olm, HSloc);
                            }
                            else // condition 6, calculate the derivative
                            {
                                single_deriv(env, elem, fsr, nnr, total_nnr, olm);
                            } // end condition 6, calc_deri
                            ++iw2_all;
                        } // end loop 5, kk
                        ++iw1_all;
                    } // end loop 4, jj
                }     // condition 3, distance
                else if (distance >= rcut && (!gamma_only_local))
                {
                    int start1 = ucell.itiaiw2iwt(t1, i1, 0);
                    int start2 = ucell.itiaiw2iwt(t2, i2, 0);

                    bool is_adj = false;
                    for (int ad0 = 0; ad0 < adjs.adj_num + 1; ++ad0)
                    {
                        const int t0 = adjs.ntype[ad0];
                        tau0 = adjs.adjacent_tau[ad0];
                        dtau1 = tau0 - tau1;
                        double distance1 = dtau1.norm() * ucell.lat0;
                        double rcut1 = orb.Phi[t1].getRcut() + ucell.infoNL->get_rcut_max(t0);
                        dtau2 = tau0 - tau2;
                        double distance2 = dtau2.norm() * ucell.lat0;
                        double rcut2 = orb.Phi[t2].getRcut() + ucell.infoNL->get_rcut_max(t0);
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
                            if (mu < 0)
                            {
                                continue;
                            }
                            for (int kk = 0; kk < atom2->nw * npol; ++kk)
                            {
                                const int nu = pv.global2local_col(start2 + kk);
                                if (nu < 0)
                                {
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

    ModuleBase::timer::end("LCAO_domain", "build_ST_new");
    return;
}

} // namespace LCAO_domain
