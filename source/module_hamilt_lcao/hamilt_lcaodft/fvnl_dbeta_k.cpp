#include "FORCE.h"
#include "module_base/memory.h"
#include "module_parameter/parameter.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_base/tool_threading.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_elecstate/cal_dm.h"
#include "module_elecstate/elecstate_lcao.h"
#include "module_elecstate/module_dm/cal_dm_psi.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/write_HS.h"
#include "module_hamilt_lcao/hamilt_lcaodft/stress_tools.h"

#include <map>
#include <unordered_map>

#ifdef __DEEPKS
#include "module_hamilt_lcao/module_deepks/LCAO_deepks.h"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

typedef std::tuple<int, int, int, int> key_tuple;

// must consider three-center H matrix.
template <>
void Force_LCAO<std::complex<double>>::cal_fvnl_dbeta(const elecstate::DensityMatrix<std::complex<double>, double>* dm,
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

    const int nspin = PARAM.inp.nspin;
    const int nspin_DMR = (nspin == 2) ? 2 : 1;
    const int npol = PARAM.globalv.npol;

    // Data structure for storing <psi|beta>, for a detailed description
    // check out the same data structure in build_Nonlocal_mu_new
    std::vector<std::map<key_tuple, std::unordered_map<int, std::vector<std::vector<double>>>>> nlm_tot;

    nlm_tot.resize(ucell.nat);

#ifdef _OPENMP
// use schedule(dynamic) for load balancing because adj_num is various
#pragma omp parallel for schedule(dynamic)
#endif
    for (int iat = 0; iat < ucell.nat; iat++)
    {

        const int it = ucell.iat2it[iat];
        const int ia = ucell.iat2ia[iat];

        // Step 1 : generate <psi|beta>
        // type of atom; distance; atomic basis; projectors

        const double Rcut_Beta = ucell.infoNL.Beta[it].get_rcut_max();
        const ModuleBase::Vector3<double> tau = ucell.atoms[it].tau[ia];
        AdjacentAtomInfo adjs;
        gd.Find_atom(ucell, tau, it, ia, &adjs);

        nlm_tot[iat].clear();

        for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
        {
            const int T1 = adjs.ntype[ad];
            const int I1 = adjs.natom[ad];
            const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
            const double Rcut_AO1 = orb.Phi[T1].getRcut();

            const ModuleBase::Vector3<double>& tau1 = adjs.adjacent_tau[ad];
            const Atom* atom1 = &ucell.atoms[T1];
            const int nw1_tot = atom1->nw * npol;

            const ModuleBase::Vector3<double> dtau = tau1 - tau;
            const double dist1 = dtau.norm2() * ucell.lat0 * ucell.lat0;
            if (dist1 > (Rcut_Beta + Rcut_AO1) * (Rcut_Beta + Rcut_AO1))
            {
                continue;
            }

            std::unordered_map<int, std::vector<std::vector<double>>> nlm_cur;
            nlm_cur.clear();

            for (int iw1 = 0; iw1 < nw1_tot; ++iw1)
            {
                const int iw1_all = start1 + iw1;
                const int iw1_local = pv.global2local_row(iw1_all);
                const int iw2_local = pv.global2local_col(iw1_all);
                if (iw1_local < 0 && iw2_local < 0)
                {
                    continue;
                }
                const int iw1_0 = iw1 / npol;
                std::vector<std::vector<double>> nlm;
                int L1 = atom1->iw2l[iw1_0];
                int N1 = atom1->iw2n[iw1_0];
                int m1 = atom1->iw2m[iw1_0];

                // convert m (0,1,...2l) to M (-l, -l+1, ..., l-1, l)
                int M1 = (m1 % 2 == 0) ? -m1 / 2 : (m1 + 1) / 2;

                ModuleBase::Vector3<double> dtau = tau - tau1;
                intor_orb_beta.snap(T1, L1, N1, M1, it, dtau * ucell.lat0, true, nlm);

                nlm_cur.insert({iw1_all, nlm});
            } // end iw
            const int iat1 = ucell.itia2iat(T1, I1);
            const int rx1 = adjs.box[ad].x;
            const int ry1 = adjs.box[ad].y;
            const int rz1 = adjs.box[ad].z;
            key_tuple key_1(iat1, rx1, ry1, rz1);
            nlm_tot[iat][key_1] = nlm_cur;
        } // end ad
    }

    //=======================================================
    // Step2:
    // calculate sum_(L0,M0) beta<psi_i|beta><beta|psi_j>
    // and accumulate the value to Hloc_fixedR(i,j)
    //=======================================================
    int total_nnr = 0;
#ifdef _OPENMP
#pragma omp parallel reduction(+ : total_nnr)
    {
        ModuleBase::matrix local_svnl_dbeta(3, 3);
        const int num_threads = omp_get_num_threads();
#else
    ModuleBase::matrix& local_svnl_dbeta = svnl_dbeta;
#endif

        ModuleBase::Vector3<double> tau1;
        ModuleBase::Vector3<double> tau2;
        ModuleBase::Vector3<double> dtau;
        ModuleBase::Vector3<double> tau0;
        ModuleBase::Vector3<double> dtau1;
        ModuleBase::Vector3<double> dtau2;

        double rcut = 0.0;
        double distance = 0.0;

        double rcut1 = 0.0;
        double rcut2 = 0.0;
        double distance1 = 0.0;
        double distance2 = 0.0;

#ifdef _OPENMP
// use schedule(dynamic) for load balancing because adj_num is various
#pragma omp for schedule(dynamic)
#endif
        for (int iat1 = 0; iat1 < ucell.nat; iat1++)
        {
            const int T1 = ucell.iat2it[iat1];
            const Atom* atom1 = &ucell.atoms[T1];

            {
                const int I1 = ucell.iat2ia[iat1];
                tau1 = atom1->tau[I1];
                AdjacentAtomInfo adjs;
                gd.Find_atom(ucell, tau1, T1, I1, &adjs);
                const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
                int nnr = pv.nlocstart[iat1];

                /*
                    !!!!!!!!!!!!!!!!
                    This optimization is also improving the performance of single thread.
                    Making memory access more linearly in the core loop
                */
                bool iat_recorded = false;
                bool force_updated = false;
                // record iat of adjs
                std::vector<int> adj_iat;
                // record fvnl_dbeta diff of adjs
                std::vector<double> adj_fvnl_dbeta;
                if (isforce)
                {
                    adj_iat.resize(adjs.adj_num + 1);
                    adj_fvnl_dbeta.resize((adjs.adj_num + 1) * 3, 0.0);
                }

                for (int ad2 = 0; ad2 < adjs.adj_num + 1; ++ad2)
                {
                    const int T2 = adjs.ntype[ad2];
                    const Atom* atom2 = &ucell.atoms[T2];

                    const int I2 = adjs.natom[ad2];
                    const int iat2 = ucell.itia2iat(T2, I2);
                    const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
                    tau2 = adjs.adjacent_tau[ad2];

                    const int rx2 = adjs.box[ad2].x;
                    const int ry2 = adjs.box[ad2].y;
                    const int rz2 = adjs.box[ad2].z;

                    dtau = tau2 - tau1;
                    distance = dtau.norm2() * (ucell.lat0 * ucell.lat0);
                    // this rcut is in order to make nnr consistent
                    // with other matrix.
                    rcut = (orb.Phi[T1].getRcut() + orb.Phi[T2].getRcut()) * (orb.Phi[T1].getRcut() + orb.Phi[T2].getRcut());

                    // check if this a adjacent atoms.
                    bool is_adj = false;
                    if (distance < rcut)
                    {
                        is_adj = true;
                    }
                    else if (distance >= rcut)
                    {
                        for (int ad0 = 0; ad0 < adjs.adj_num + 1; ++ad0)
                        {
                            const int T0 = adjs.ntype[ad0];
                            if (ucell.infoNL.nproj[T0] == 0)
                            {
                                continue;
                            }
                            const int I0 = adjs.natom[ad0];
                            // const int iat0 = ucell.itia2iat(T0, I0);
                            // const int start0 = ucell.itiaiw2iwt(T0, I0, 0);

                            tau0 = adjs.adjacent_tau[ad0];
                            dtau1 = tau0 - tau1;
                            distance1 = dtau1.norm() * ucell.lat0;
                            rcut1 = orb.Phi[T1].getRcut() + ucell.infoNL.Beta[T0].get_rcut_max();

                            dtau2 = tau0 - tau2;
                            distance2 = dtau2.norm() * ucell.lat0;
                            rcut2 = orb.Phi[T2].getRcut() + ucell.infoNL.Beta[T0].get_rcut_max();

                            if (distance1 < rcut1 && distance2 < rcut2)
                            {
                                is_adj = true;
                                break;
                            }
                        }
                    }

                    if (is_adj)
                    {
                        // basematrix and its data pointer
                        if (pv.get_row_size(iat1) <= 0 || pv.get_col_size(iat2) <= 0)
                        {
                            continue;
                        }
                        std::vector<double*> tmp_matrix_ptr;
                        for (int is = 0; is < nspin_DMR; ++is)
                        {
                            auto* tmp_base_matrix = dm->get_DMR_pointer(is + 1)->find_matrix(iat1, iat2, rx2, ry2, rz2);
                            tmp_matrix_ptr.push_back(tmp_base_matrix->get_pointer());
                        }

                        for (int ad0 = 0; ad0 < adjs.adj_num + 1; ++ad0)
                        {
                            const int T0 = adjs.ntype[ad0];
                            const int I0 = adjs.natom[ad0];
                            const int iat = ucell.itia2iat(T0, I0);
                            if (!iat_recorded && isforce)
                            {
                                adj_iat[ad0] = iat;
                            }

                            // mohan add 2010-12-19
                            if (ucell.infoNL.nproj[T0] == 0)
                            {
                                continue;
                            }

                            // const int I0 = gd.getNatom(ad0);
                            // const int start0 = ucell.itiaiw2iwt(T0, I0, 0);
                            tau0 = adjs.adjacent_tau[ad0];

                            dtau1 = tau0 - tau1;
                            dtau2 = tau0 - tau2;
                            const double distance1 = dtau1.norm2() * ucell.lat0 * ucell.lat0;
                            const double distance2 = dtau2.norm2() * ucell.lat0 * ucell.lat0;

                            // seems a bug here!! mohan 2011-06-17
                            rcut1 = (orb.Phi[T1].getRcut() + ucell.infoNL.Beta[T0].get_rcut_max()) * (orb.Phi[T1].getRcut() + ucell.infoNL.Beta[T0].get_rcut_max());
                            rcut2 = (orb.Phi[T2].getRcut() + ucell.infoNL.Beta[T0].get_rcut_max()) * (orb.Phi[T2].getRcut() + ucell.infoNL.Beta[T0].get_rcut_max());

                            double r0[3];
                            double r1[3];
                            r1[0] = (tau1.x - tau0.x);
                            r1[1] = (tau1.y - tau0.y);
                            r1[2] = (tau1.z - tau0.z);
                            r0[0] = (tau2.x - tau0.x);
                            r0[1] = (tau2.y - tau0.y);
                            r0[2] = (tau2.z - tau0.z);

                            if (distance1 >= rcut1 || distance2 >= rcut2)
                            {
                                continue;
                            }

                            const int rx0 = adjs.box[ad0].x;
                            const int ry0 = adjs.box[ad0].y;
                            const int rz0 = adjs.box[ad0].z;
                            key_tuple key1(iat1, -rx0, -ry0, -rz0);
                            key_tuple key2(iat2, rx2 - rx0, ry2 - ry0, rz2 - rz0);

                            int nnr_inner = 0;
                            int dhsize = pv.get_row_size(iat1) * pv.get_col_size(iat2);
                            std::vector<std::vector<double>> dhvnl(3, std::vector<double>(dhsize, 0.0));
                            std::vector<std::vector<double>> dhvnl1(3, std::vector<double>(dhsize, 0.0));
                            for (int j = 0; j < atom1->nw * npol; j++)
                            {
                                const int j0 = j / npol; // added by zhengdy-soc
                                const int iw1_all = start1 + j;
                                const int mu = pv.global2local_row(iw1_all);
                                if (mu < 0)
                                {
                                    continue;
                                }

                                for (int k = 0; k < atom2->nw * npol; k++)
                                {
                                    const int k0 = k / npol;
                                    const int iw2_all = start2 + k;
                                    const int nu = pv.global2local_col(iw2_all);
                                    if (nu < 0)
                                    {
                                        continue;
                                    }

                                    if (nspin == 4)
                                    {
                                        const double* nlm_1 = nlm_tot[iat][key2][iw2_all][0].data();
                                        const double* nlm_2[3];
                                        for (int i = 0; i < 3; i++)
                                        {
                                            nlm_2[i] = nlm_tot[iat][key1][iw1_all][i + 1].data();
                                        }
                                        int is0 = (j - j0 * npol) + (k - k0 * npol) * 2;
                                        for (int no = 0; no < ucell.atoms[T0].ncpp.non_zero_count_soc[is0]; no++)
                                        {
                                            const int p1 = ucell.atoms[T0].ncpp.index1_soc[is0][no];
                                            const int p2 = ucell.atoms[T0].ncpp.index2_soc[is0][no];
                                            for (int ir = 0; ir < 3; ir++)
                                            {
                                                if (is0 == 0)
                                                {
                                                    dhvnl[ir][nnr_inner]
                                                        += nlm_2[ir][p1] * nlm_1[p2]
                                                           * (ucell.atoms[T0].ncpp.d_so(0, p2, p1).real()
                                                              + ucell.atoms[T0].ncpp.d_so(3, p2, p1).real())
                                                           * 0.5;
                                                }
                                                else if (is0 == 1)
                                                {
                                                    dhvnl[ir][nnr_inner]
                                                        += nlm_2[ir][p1] * nlm_1[p2]
                                                           * (ucell.atoms[T0].ncpp.d_so(1, p2, p1).real()
                                                              + ucell.atoms[T0].ncpp.d_so(2, p2, p1).real())
                                                           * 0.5;
                                                }
                                                else if (is0 == 2)
                                                {
                                                    dhvnl[ir][nnr_inner]
                                                        += nlm_2[ir][p1] * nlm_1[p2]
                                                           * (-ucell.atoms[T0].ncpp.d_so(1, p2, p1).imag()
                                                              + ucell.atoms[T0].ncpp.d_so(2, p2, p1).imag())
                                                           * 0.5;
                                                }
                                                else if (is0 == 3)
                                                {
                                                    dhvnl[ir][nnr_inner]
                                                        += nlm_2[ir][p1] * nlm_1[p2]
                                                           * (ucell.atoms[T0].ncpp.d_so(0, p2, p1).real()
                                                              - ucell.atoms[T0].ncpp.d_so(3, p2, p1).real())
                                                           * 0.5;
                                                }
                                            }
                                        }
                                        if (isstress)
                                        {
                                            const double* nlm_1 = nlm_tot[iat][key1][iw1_all][0].data();
                                            const double* nlm_2[3];
                                            for (int i = 0; i < 3; i++)
                                            {
                                                nlm_2[i] = nlm_tot[iat][key2][iw2_all][i + 1].data();
                                            }
                                            int is0 = (j - j0 * npol) + (k - k0 * npol) * 2;
                                            for (int no = 0; no < ucell.atoms[T0].ncpp.non_zero_count_soc[is0]; no++)
                                            {
                                                const int p1 = ucell.atoms[T0].ncpp.index1_soc[is0][no];
                                                const int p2 = ucell.atoms[T0].ncpp.index2_soc[is0][no];
                                                for (int ir = 0; ir < 3; ir++)
                                                {
                                                    if (is0 == 0)
                                                    {
                                                        dhvnl1[ir][nnr_inner]
                                                            += nlm_2[ir][p1] * nlm_1[p2]
                                                               * (ucell.atoms[T0].ncpp.d_so(0, p2, p1).real()
                                                                  + ucell.atoms[T0].ncpp.d_so(3, p2, p1).real())
                                                               * 0.5;
                                                    }
                                                    else if (is0 == 1)
                                                    {
                                                        dhvnl1[ir][nnr_inner]
                                                            += nlm_2[ir][p1] * nlm_1[p2]
                                                               * (ucell.atoms[T0].ncpp.d_so(1, p2, p1).real()
                                                                  + ucell.atoms[T0].ncpp.d_so(2, p2, p1).real())
                                                               * 0.5;
                                                    }
                                                    else if (is0 == 2)
                                                    {
                                                        dhvnl1[ir][nnr_inner]
                                                            += nlm_2[ir][p1] * nlm_1[p2]
                                                               * (-ucell.atoms[T0].ncpp.d_so(1, p2, p1).imag()
                                                                  + ucell.atoms[T0].ncpp.d_so(2, p2, p1).imag())
                                                               * 0.5;
                                                    }
                                                    else if (is0 == 3)
                                                    {
                                                        dhvnl1[ir][nnr_inner]
                                                            += nlm_2[ir][p1] * nlm_1[p2]
                                                               * (ucell.atoms[T0].ncpp.d_so(0, p2, p1).real()
                                                                  - ucell.atoms[T0].ncpp.d_so(3, p2, p1).real())
                                                               * 0.5;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    else if (nspin == 1 || nspin == 2)
                                    {
                                        // const Atom* atom0 = &ucell.atoms[T0];
                                        double nlm[3] = {0, 0, 0};
                                        const double* nlm_1 = nlm_tot[iat][key2][iw2_all][0].data();
                                        const double* nlm_2[3];
                                        for (int i = 0; i < 3; i++)
                                        {
                                            nlm_2[i] = nlm_tot[iat][key1][iw1_all][i + 1].data();
                                        }

                                        const double* tmp_d = nullptr;
                                        for (int no = 0; no < ucell.atoms[T0].ncpp.non_zero_count_soc[0]; no++)
                                        {
                                            const int p1 = ucell.atoms[T0].ncpp.index1_soc[0][no];
                                            const int p2 = ucell.atoms[T0].ncpp.index2_soc[0][no];
                                            ucell.atoms[T0].ncpp.get_d(0, p1, p2, tmp_d);
                                            for (int ir = 0; ir < 3; ir++)
                                            {
                                                nlm[ir] += nlm_2[ir][p2] * nlm_1[p1] * (*tmp_d);
                                            }
                                        }

                                        double nlm1[3] = {0, 0, 0};
                                        if (isstress)
                                        {
                                            const double* nlm_1 = nlm_tot[iat][key1][iw1_all][0].data();
                                            const double* nlm_2[3];
                                            for (int i = 0; i < 3; i++)
                                            {
                                                nlm_2[i] = nlm_tot[iat][key2][iw2_all][i + 1].data();
                                            }

                                            const double* tmp_d = nullptr;
                                            for (int no = 0; no < ucell.atoms[T0].ncpp.non_zero_count_soc[0]; no++)
                                            {
                                                const int p1 = ucell.atoms[T0].ncpp.index1_soc[0][no];
                                                const int p2 = ucell.atoms[T0].ncpp.index2_soc[0][no];
                                                ucell.atoms[T0].ncpp.get_d(0, p1, p2, tmp_d);
                                                for (int ir = 0; ir < 3; ir++)
                                                {
                                                    nlm1[ir] += nlm_2[ir][p1] * nlm_1[p2] * (*tmp_d);
                                                }
                                            }
                                        }

                                        /// only one projector for each atom force, but another projector for stress
                                        force_updated = true;
                                        // get DMR
                                        double dm2d1 = 0.0;
                                        for (int is = 0; is < nspin_DMR; ++is)
                                        {
                                            dm2d1 += tmp_matrix_ptr[is][nnr_inner];
                                        }
                                        double dm2d2 = 2.0 * dm2d1;
                                        //
                                        for (int jpol = 0; jpol < 3; jpol++)
                                        {
                                            if (isforce)
                                            {
                                                adj_fvnl_dbeta[ad0 * 3 + jpol] -= dm2d2 * nlm[jpol];
                                            }
                                            if (isstress)
                                            {
                                                for (int ipol = jpol; ipol < 3; ipol++)
                                                {
                                                    local_svnl_dbeta(jpol, ipol)
                                                        += dm2d1 * (nlm[jpol] * r1[ipol] + nlm1[jpol] * r0[ipol]);
                                                }
                                            }
                                        }
                                    }
                                    else
                                    {
                                        ModuleBase::WARNING_QUIT("Force_LCAO_k::cal_fvnl_dbeta_k",
                                                                 "nspin must be 1, 2 or 4");
                                    }
                                    //}
                                    nnr_inner++;
                                } // k
                            }     // j
                            if (nspin == 4)
                            {
                                force_updated = true;
                                if (isforce)
                                {
                                    for (int ir = 0; ir < dhsize; ir++)
                                    {
                                        adj_fvnl_dbeta[ad0 * 3 + 0] -= 2.0 * tmp_matrix_ptr[0][ir] * dhvnl[0][ir];
                                        adj_fvnl_dbeta[ad0 * 3 + 1] -= 2.0 * tmp_matrix_ptr[0][ir] * dhvnl[1][ir];
                                        adj_fvnl_dbeta[ad0 * 3 + 2] -= 2.0 * tmp_matrix_ptr[0][ir] * dhvnl[2][ir];
                                    }
                                }
                                if (isstress)
                                {
                                    for (int ir = 0; ir < dhsize; ir++)
                                    {
                                        for (int jpol = 0; jpol < 3; jpol++)
                                        {
                                            if (isstress)
                                            {
                                                for (int ipol = jpol; ipol < 3; ipol++)
                                                {
                                                    local_svnl_dbeta(jpol, ipol)
                                                        += tmp_matrix_ptr[0][ir]
                                                           * (dhvnl[jpol][ir] * r1[ipol] + dhvnl1[jpol][ir] * r0[ipol]);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        } // ad0

                        // outer circle : accumulate nnr
                        for (int j = 0; j < atom1->nw * npol; j++)
                        {
                            const int j0 = j / npol; // added by zhengdy-soc
                            const int iw1_all = start1 + j;
                            const int mu = pv.global2local_row(iw1_all);
                            if (mu < 0)
                            {
                                continue;
                            }

                            // fix a serious bug: atom2[T2] -> atom2
                            // mohan 2010-12-20
                            for (int k = 0; k < atom2->nw * npol; k++)
                            {
                                const int k0 = k / npol;
                                const int iw2_all = start2 + k;
                                const int nu = pv.global2local_col(iw2_all);
                                if (nu < 0)
                                {
                                    continue;
                                }
                                total_nnr++;
                                nnr++;
                            }
                        }
                        iat_recorded = true;
                    } // is_adj
                }     // ad2

                // sum the diff to fvnl_dbeta
                if (force_updated && isforce)
                {
#ifdef _OPENMP
                    if (num_threads > 1)
                    {
                        for (int ad0 = 0; ad0 < adjs.adj_num + 1; ++ad0)
                        {
#pragma omp atomic
                            fvnl_dbeta(adj_iat[ad0], 0) += adj_fvnl_dbeta[ad0 * 3 + 0];
#pragma omp atomic
                            fvnl_dbeta(adj_iat[ad0], 1) += adj_fvnl_dbeta[ad0 * 3 + 1];
#pragma omp atomic
                            fvnl_dbeta(adj_iat[ad0], 2) += adj_fvnl_dbeta[ad0 * 3 + 2];
                        }
                    }
                    else
#endif
                    {
                        for (int ad0 = 0; ad0 < adjs.adj_num + 1; ++ad0)
                        {
                            fvnl_dbeta(adj_iat[ad0], 0) += adj_fvnl_dbeta[ad0 * 3 + 0];
                            fvnl_dbeta(adj_iat[ad0], 1) += adj_fvnl_dbeta[ad0 * 3 + 1];
                            fvnl_dbeta(adj_iat[ad0], 2) += adj_fvnl_dbeta[ad0 * 3 + 2];
                        }
                    }
                }
            } // I1
        }     // T1

#ifdef _OPENMP
        if (isstress)
        {
#pragma omp critical(cal_fvnl_dbeta_k_new_reduce)
            {
                for (int l = 0; l < 3; l++)
                {
                    for (int m = 0; m < 3; m++)
                    {
                        svnl_dbeta(l, m) += local_svnl_dbeta(l, m);
                    }
                }
            }
        }
    }
#endif

    assert(total_nnr == pv.nnr);

    if (isstress)
    {
        StressTools::stress_fill(ucell.lat0, ucell.omega, svnl_dbeta);
    }

    ModuleBase::timer::tick("Force_LCAO", "cal_fvnl_dbeta");
    return;
}
