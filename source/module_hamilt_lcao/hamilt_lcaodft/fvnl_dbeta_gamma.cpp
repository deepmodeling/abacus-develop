#include "FORCE.h"
#include "module_base/timer.h"
#include "module_parameter/parameter.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hamilt_lcao/hamilt_lcaodft/stress_tools.h"

#include <unordered_map>

template <>
void Force_LCAO<double>::cal_fvnl_dbeta(const elecstate::DensityMatrix<double, double>* dm,
                                        const Parallel_Orbitals& pv,
                                        const UnitCell& ucell,
                                        const LCAO_Orbitals& orb,
                                        const TwoCenterIntegrator& intor_orb_beta,
                                        Grid_Driver& gd,
                                        const bool isforce,
                                        const bool isstress,
                                        ModuleBase::matrix& fvnl_dbeta,
                                        ModuleBase::matrix& svnl_dbeta)
{
    ModuleBase::TITLE("Force_LCAO", "cal_fvnl_dbeta");
    ModuleBase::timer::tick("Force_LCAO", "cal_fvnl_dbeta");

// use schedule(dynamic) for load balancing because adj_num is various
#pragma omp parallel
{ 
    ModuleBase::matrix local_svnl_dbeta(3, 3);
    #pragma omp for schedule(dynamic)
    for (int iat = 0; iat < ucell.nat; iat++)
    {
        const int it = ucell.iat2it[iat];
        const int ia = ucell.iat2ia[iat];
        const int T0 = it;
        const int I0 = ia;
        double r0[3]{0};
        double r1[3]{0};

        const ModuleBase::Vector3<double> tau0 = ucell.atoms[it].tau[ia];
        // find ajacent atom of atom ia
        // gd.Find_atom( ucell.atoms[it].tau[ia] );
        AdjacentAtomInfo adjs;
        gd.Find_atom(ucell, tau0, it, ia, &adjs);
        const double Rcut_Beta = ucell.infoNL.Beta[it].get_rcut_max();

        std::vector<std::unordered_map<int, std::vector<std::vector<double>>>> nlm_tot;
        nlm_tot.resize(adjs.adj_num + 1); // this saves <psi_i|beta> and <psi_i|\nabla|beta>

        // Step 1 : Calculate and save <psi_i|beta> and <psi_i|\nabla|beta>
        for (int ad1 = 0; ad1 < adjs.adj_num + 1; ad1++)
        {
            const int T1 = adjs.ntype[ad1];
            const Atom* atom1 = &ucell.atoms[T1];
            const int I1 = adjs.natom[ad1];
            const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
            const ModuleBase::Vector3<double> tau1 = adjs.adjacent_tau[ad1];
            const double Rcut_AO1 = orb.Phi[T1].getRcut();

            nlm_tot[ad1].clear();

            const double dist1 = (tau1 - tau0).norm() * ucell.lat0;
            if (dist1 > Rcut_Beta + Rcut_AO1)
            {
                continue;
            }

            for (int iw1 = 0; iw1 < ucell.atoms[T1].nw; ++iw1)
            {
                const int iw1_all = start1 + iw1;
                const int iw1_local = pv.global2local_row(iw1_all);
                const int iw2_local = pv.global2local_col(iw1_all);

                if (iw1_local < 0 && iw2_local < 0)
                {
                    continue;
                }

                std::vector<std::vector<double>> nlm;

                int L1 = atom1->iw2l[iw1];
                int N1 = atom1->iw2n[iw1];
                int m1 = atom1->iw2m[iw1];

                // convert m (0,1,...2l) to M (-l, -l+1, ..., l-1, l)
                int M1 = (m1 % 2 == 0) ? -m1 / 2 : (m1 + 1) / 2;

                ModuleBase::Vector3<double> dtau = ucell.atoms[T0].tau[I0] - tau1;

                intor_orb_beta.snap(T1, L1, N1, M1, T0, dtau * ucell.lat0, true, nlm);

                assert(nlm.size() == 4);
                nlm_tot[ad1].insert({iw1, nlm});
            }
        } // ad

        // Step 2 : sum to get beta<psi_i|beta><beta|\nabla|psi_j>
        for (int ad1 = 0; ad1 < adjs.adj_num + 1; ++ad1)
        {
            const int T1 = adjs.ntype[ad1];
            const Atom* atom1 = &ucell.atoms[T1];
            const int I1 = adjs.natom[ad1];
            const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
            const ModuleBase::Vector3<double> tau1 = adjs.adjacent_tau[ad1];
            const double Rcut_AO1 = orb.Phi[T1].getRcut();

            for (int ad2 = 0; ad2 < adjs.adj_num + 1; ad2++)
            {
                const int T2 = adjs.ntype[ad2];
                const Atom* atom2 = &ucell.atoms[T2];
                const int I2 = adjs.natom[ad2];
                const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
                const ModuleBase::Vector3<double> tau2 = adjs.adjacent_tau[ad2];
                const double Rcut_AO2 = orb.Phi[T2].getRcut();

                const double dist1 = (tau1 - tau0).norm() * ucell.lat0;
                const double dist2 = (tau2 - tau0).norm() * ucell.lat0;
                if (isstress)
                {
                    r1[0] = (tau1.x - tau0.x);
                    r1[1] = (tau1.y - tau0.y);
                    r1[2] = (tau1.z - tau0.z);
                    r0[0] = (tau2.x - tau0.x);
                    r0[1] = (tau2.y - tau0.y);
                    r0[2] = (tau2.z - tau0.z);
                }

                if (dist1 > Rcut_Beta + Rcut_AO1 || dist2 > Rcut_Beta + Rcut_AO2)
                {
                    continue;
                }

                for (int iw1 = 0; iw1 < ucell.atoms[T1].nw; ++iw1)
                {
                    const int iw1_all = start1 + iw1;
                    const int iw1_local = pv.global2local_row(iw1_all);
                    if (iw1_local < 0)
                    {
                        continue;
                    }
                    for (int iw2 = 0; iw2 < ucell.atoms[T2].nw; ++iw2)

                    {
                        const int iw2_all = start2 + iw2;
                        const int iw2_local = pv.global2local_col(iw2_all);

                        if (iw2_local < 0)
                        {
                            continue;
                        }

                        double nlm[3] = {0, 0, 0};
                        double nlm_t[3] = {0, 0, 0}; // transpose

                        const double* nlm1 = nlm_tot[ad1][iw1][0].data();
                        const double* nlm2[3];
                        for (int i = 0; i < 3; i++)
                        {
                            nlm2[i] = nlm_tot[ad2][iw2][i + 1].data();
                        }

                        const double* tmp_d = nullptr;
                        for (int no = 0; no < ucell.atoms[T0].ncpp.non_zero_count_soc[0]; no++)
                        {
                            const int p1 = ucell.atoms[T0].ncpp.index1_soc[0][no];
                            const int p2 = ucell.atoms[T0].ncpp.index2_soc[0][no];
                            ucell.atoms[T0].ncpp.get_d(0, p1, p2, tmp_d);
                            for (int ir = 0; ir < 3; ir++)
                            {
                                nlm[ir] += nlm2[ir][p2] * nlm1[p1] * (*tmp_d);
                            }
                        }

                        if (isstress)
                        {
                            nlm1 = nlm_tot[ad2][iw2][0].data();
                            for (int i = 0; i < 3; i++)
                            {
                                nlm2[i] = nlm_tot[ad1][iw1][i + 1].data();
                            }

                            const double* tmp_d = nullptr;
                            for (int no = 0; no < ucell.atoms[T0].ncpp.non_zero_count_soc[0]; no++)
                            {
                                const int p1 = ucell.atoms[T0].ncpp.index1_soc[0][no];
                                const int p2 = ucell.atoms[T0].ncpp.index2_soc[0][no];
                                ucell.atoms[T0].ncpp.get_d(0, p1, p2, tmp_d);
                                for (int ir = 0; ir < 3; ir++)
                                {
                                    nlm_t[ir] += nlm2[ir][p1] * nlm1[p2] * (*tmp_d);
                                }
                            }
                        }
                        // dbeta is minus, that's consistent.
                        // only one projector for each atom force.

                        double sum = 0.0;
                        for (int is = 0; is < PARAM.inp.nspin; ++is)
                        {
                            // sum += dm2d[is](iw2_local, iw1_local);
                            sum += dm->get_DMK(is + 1, 0, iw2_local, iw1_local);
                        }
                        sum *= 2.0;

                        if (isforce)
                        {
                            fvnl_dbeta(iat, 0) -= sum * nlm[0];
                            fvnl_dbeta(iat, 1) -= sum * nlm[1];
                            fvnl_dbeta(iat, 2) -= sum * nlm[2];
                        }

                        if (isstress)
                        {
                            for (int ipol = 0; ipol < 3; ipol++)
                            {
                                for (int jpol = ipol; jpol < 3; jpol++)
                                {
                                    local_svnl_dbeta(ipol, jpol)
                                    += sum / 2.0 * (nlm[ipol] * r0[jpol] + nlm_t[ipol] * r1[jpol]);
                                }
                            }
                        }
                    } // iw2
                }     // iw1
            }         // ad2
        }             // ad1
    }                 // iat

    // sum up local_svnl_dbeta to svnl_dbeta
    if (isstress)
    {
        #pragma omp critical
        {
            for (int ipol = 0; ipol < 3; ipol++)
            {
                for (int jpol = ipol; jpol < 3; jpol++)
                {
                    svnl_dbeta(ipol, jpol) += local_svnl_dbeta(ipol, jpol);
                }
            }
        }
    }
}
    if (isstress)
    {
        StressTools::stress_fill(ucell.lat0, ucell.omega, svnl_dbeta);
    }

    ModuleBase::timer::tick("Force_LCAO", "cal_fvnl_dbeta");
}
