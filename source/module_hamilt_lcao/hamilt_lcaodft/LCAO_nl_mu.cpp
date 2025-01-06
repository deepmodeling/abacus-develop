#include "module_base/timer.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_domain.h"
#include "module_parameter/parameter.h"

namespace LCAO_domain
{

typedef std::tuple<int, int, int, int> key_tuple;

#include "record_adj.h" //mohan add 2012-07-06

void build_Nonlocal_mu_new(const Parallel_Orbitals& pv,
                           ForceStressArrays& fsr,
                           double* NLloc,
                           const bool& calc_deri,
                           const UnitCell& ucell,
                           const LCAO_Orbitals& orb,
                           const TwoCenterIntegrator& intor_orb_beta,
                           const Grid_Driver* GridD)
{
    ModuleBase::TITLE("LCAO_domain", "vnl_mu_new");
    ModuleBase::timer::tick("LCAO_domain", "vnl_mu_new");

    const int nspin = PARAM.inp.nspin;
    const int npol = PARAM.globalv.npol;
    const bool gamma_only_local = PARAM.globalv.gamma_only_local;

    // < phi1 | beta > < beta | phi2 >
    // phi1 is within the unitcell.
    // while beta is in the supercell.
    // while phi2 is in the supercell.

    // Step 1 : generate <psi|beta>

    // This is the data structure for storing <psi|beta>
    // It is a 4 layer data structure
    // The outmost layer is std::vector with size being number of atoms in unit cell
    // The second layer is a map, the key being a combination of 4 number (iat, dRx, dRy, dRz)
    // which identifies a unique adjacent atom of the first atom
    // The third layer is an unordered map, with key being the index of atomic basis |psi>
    // The inner layer is a vector, each element representing a projector |beta>
    // It then either stores the number <psi|beta> (nlm_tot)
    // or a vector of 4, storing additionally <d/dx_i psi|beta> (nlm_tot1) x_i=x,y,z
    std::vector<std::map<key_tuple, std::unordered_map<int, std::vector<double>>>> nlm_tot;
    std::vector<std::map<key_tuple, std::unordered_map<int, std::vector<std::vector<double>>>>> nlm_tot1;

    if (!calc_deri)
    {
        nlm_tot.resize(ucell.nat);
    }
    else
    {
        nlm_tot1.resize(ucell.nat);
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int iat = 0; iat < ucell.nat; iat++)
    {
        const int it = ucell.iat2it[iat];
        const int ia = ucell.iat2ia[iat];

        const double Rcut_Beta = ucell.infoNL.Beta[it].get_rcut_max();
        const ModuleBase::Vector3<double> tau = ucell.atoms[it].tau[ia];
        AdjacentAtomInfo adjs;
        GridD->Find_atom(ucell, tau, it, ia, &adjs);

        if (!calc_deri)
        {
            nlm_tot[iat].clear();
        }
        else
        {
            nlm_tot1[iat].clear();
        }

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
            const double dist1 = dtau.norm2() * pow(ucell.lat0, 2);

            if (dist1 > pow(Rcut_Beta + Rcut_AO1, 2))
            {
                continue;
            }
            std::unordered_map<int, std::vector<double>> nlm_cur;
            std::unordered_map<int, std::vector<std::vector<double>>> nlm_cur1;

            if (!calc_deri)
            {
                nlm_cur.clear();
            }
            else
            {
                nlm_cur1.clear();
            }
            for (int iw1 = 0; iw1 < nw1_tot; ++iw1)
            {
                const int iw1_all = start1 + iw1;
                const int iw1_local = pv.global2local_row(iw1_all);
                const int iw2_local = pv.global2local_col(iw1_all);
                if (iw1_local < 0 && iw2_local < 0) {
                    continue;
}
                const int iw1_0 = iw1 / npol;
                std::vector<std::vector<double>> nlm;
                // nlm is a vector of vectors, but size of outer vector is only 1 here
                // If we are calculating force, we need also to store the gradient
                // and size of outer vector is then 4
                // inner loop : all projectors (L0,M0)

                int L1 = atom1->iw2l[iw1_0];
                int N1 = atom1->iw2n[iw1_0];
                int m1 = atom1->iw2m[iw1_0];

                // convert m (0,1,...2l) to M (-l, -l+1, ..., l-1, l)
                int M1 = (m1 % 2 == 0) ? -m1 / 2 : (m1 + 1) / 2;

                ModuleBase::Vector3<double> dtau = tau - tau1;
                intor_orb_beta.snap(T1, L1, N1, M1, it, dtau * ucell.lat0, calc_deri, nlm);

                if (!calc_deri)
                {
                    nlm_cur.insert({iw1_all, nlm[0]});
                }
                else
                {
                    nlm_cur1.insert({iw1_all, nlm});
                }
            } // end iw

            const int iat1 = ucell.itia2iat(T1, I1);
            const int rx1 = adjs.box[ad].x;
            const int ry1 = adjs.box[ad].y;
            const int rz1 = adjs.box[ad].z;
            key_tuple key_1(iat1, rx1, ry1, rz1);

            if (!calc_deri)
            {
                nlm_tot[iat][key_1] = nlm_cur;
            }
            else
            {
                nlm_tot1[iat][key_1] = nlm_cur1;
            }
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
#endif
        ModuleBase::Vector3<double> tau1, tau2, dtau;
        ModuleBase::Vector3<double> dtau1, dtau2, tau0;
        double distance = 0.0;
        double rcut = 0.0;
        double rcut1 = 0.0;
        double rcut2 = 0.0;

        //	Record_adj RA;
        //	RA.for_2d();

        // psi1
#ifdef _OPENMP
// use schedule(dynamic) for load balancing because adj_num is various
#pragma omp for schedule(dynamic)
#endif
        for (int iat1 = 0; iat1 < ucell.nat; iat1++)
        {
            const int T1 = ucell.iat2it[iat1];
            const Atom* atom1 = &ucell.atoms[T1];
            const int I1 = ucell.iat2ia[iat1];
            {
                // GridD->Find_atom( atom1->tau[I1] );
                AdjacentAtomInfo adjs;
                GridD->Find_atom(ucell, atom1->tau[I1], T1, I1, &adjs);
                const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
                // Record_adj.for_2d() may not called in some case
                int nnr = pv.nlocstart ? pv.nlocstart[iat1] : 0;
                tau1 = atom1->tau[I1];

                // psi2
                for (int ad2 = 0; ad2 < adjs.adj_num + 1; ++ad2)
                {
                    const int T2 = adjs.ntype[ad2];
                    const Atom* atom2 = &ucell.atoms[T2];

                    const int I2 = adjs.natom[ad2];
                    const int iat2 = ucell.itia2iat(T2, I2);
                    const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
                    tau2 = adjs.adjacent_tau[ad2];

                    bool is_adj = false;

                    const int rx2 = adjs.box[ad2].x;
                    const int ry2 = adjs.box[ad2].y;
                    const int rz2 = adjs.box[ad2].z;

                    dtau = tau2 - tau1;
                    distance = dtau.norm2() * pow(ucell.lat0, 2);
                    // this rcut is in order to make nnr consistent
                    // with other matrix.
                    rcut = pow(orb.Phi[T1].getRcut() + orb.Phi[T2].getRcut(), 2);
                    if (distance < rcut) {
                        is_adj = true;
                    } else if (distance >= rcut)
                    {
                        for (int ad0 = 0; ad0 < adjs.adj_num + 1; ++ad0)
                        {
                            const int T0 = adjs.ntype[ad0];

                            tau0 = adjs.adjacent_tau[ad0];
                            dtau1 = tau0 - tau1;
                            dtau2 = tau0 - tau2;

                            const double distance1 = dtau1.norm2() * pow(ucell.lat0, 2);
                            const double distance2 = dtau2.norm2() * pow(ucell.lat0, 2);

                            rcut1 = pow(orb.Phi[T1].getRcut() + ucell.infoNL.Beta[T0].get_rcut_max(), 2);
                            rcut2 = pow(orb.Phi[T2].getRcut() + ucell.infoNL.Beta[T0].get_rcut_max(), 2);

                            if (distance1 < rcut1 && distance2 < rcut2)
                            {
                                is_adj = true;
                                break;
                            }
                        }
                    }

                    if (is_adj)
                    {
                        // < psi1 | all projectors | psi2 >
                        // ----------------------------- enter the nnr increaing zone -------------------------
                        for (int ad0 = 0; ad0 < adjs.adj_num + 1; ++ad0)
                        {
                            const int T0 = adjs.ntype[ad0];
                            const int I0 = adjs.natom[ad0];
                            const int iat = ucell.itia2iat(T0, I0);

                            // mohan add 2010-12-19
                            if (ucell.infoNL.nproj[T0] == 0)
                            {
                                continue;
                            }

                            tau0 = adjs.adjacent_tau[ad0];

                            dtau1 = tau0 - tau1;
                            dtau2 = tau0 - tau2;
                            const double distance1 = dtau1.norm2() * pow(ucell.lat0, 2);
                            const double distance2 = dtau2.norm2() * pow(ucell.lat0, 2);

                            // seems a bug here!! mohan 2011-06-17
                            rcut1 = pow(orb.Phi[T1].getRcut() + ucell.infoNL.Beta[T0].get_rcut_max(), 2);
                            rcut2 = pow(orb.Phi[T2].getRcut() + ucell.infoNL.Beta[T0].get_rcut_max(), 2);

                            if (distance1 >= rcut1 || distance2 >= rcut2)
                            {
                                continue;
                            }
                            // const Atom* atom0 = &ucell.atoms[T0];
                            const int rx0 = adjs.box[ad0].x;
                            const int ry0 = adjs.box[ad0].y;
                            const int rz0 = adjs.box[ad0].z;
                            key_tuple key1(iat1, -rx0, -ry0, -rz0);
                            key_tuple key2(iat2, rx2 - rx0, ry2 - ry0, rz2 - rz0);

                            std::unordered_map<int, std::vector<double>>* nlm_cur1_e; // left hand side, for energy
                            std::unordered_map<int, std::vector<std::vector<double>>>* nlm_cur1_f; // lhs, for force
                            std::unordered_map<int, std::vector<double>>* nlm_cur2_e;              // rhs, for energy
                            std::unordered_map<int, std::vector<std::vector<double>>>* nlm_cur2_f; // rhs, for force

                            if (!calc_deri)
                            {
                                nlm_cur1_e = &nlm_tot[iat][key1];
                                nlm_cur2_e = &nlm_tot[iat][key2];
                            }
                            else
                            {
                                nlm_cur1_f = &nlm_tot1[iat][key1];
                                nlm_cur2_f = &nlm_tot1[iat][key2];
                            }

                            int nnr_inner = 0;

                            for (int j = 0; j < atom1->nw * npol; j++)
                            {
                                const int j0 = j / npol; // added by zhengdy-soc
                                const int iw1_all = start1 + j;
                                const int mu = pv.global2local_row(iw1_all);
                                if (mu < 0) {
                                    continue;
}

                                // fix a serious bug: atom2[T2] -> atom2
                                // mohan 2010-12-20
                                for (int k = 0; k < atom2->nw * npol; k++)
                                {
                                    const int k0 = k / npol;
                                    const int iw2_all = start2 + k;
                                    const int nu = pv.global2local_col(iw2_all);
                                    if (nu < 0) {
                                        continue;
}

                                    if (!calc_deri)
                                    {
                                        std::vector<double> nlm_1 = (*nlm_cur1_e)[iw1_all];
                                        std::vector<double> nlm_2 = (*nlm_cur2_e)[iw2_all];
                                        if (nspin == 2 || nspin == 1)
                                        {
                                            double nlm_tmp = 0.0;
                                            const double* tmp_d = nullptr;
                                            for (int no = 0; no < ucell.atoms[T0].ncpp.non_zero_count_soc[0]; no++)
                                            {
                                                const int p1 = ucell.atoms[T0].ncpp.index1_soc[0][no];
                                                const int p2 = ucell.atoms[T0].ncpp.index2_soc[0][no];
                                                ucell.atoms[T0].ncpp.get_d(0, p1, p2, tmp_d);
                                                nlm_tmp += nlm_2[p2] * nlm_1[p1] * (*tmp_d);
                                            }

                                            if (gamma_only_local)
                                            {
                                                // mohan add 2010-12-20
                                                if (nlm_tmp != 0.0)
                                                {
                                                    LCAO_domain::set_mat2d(iw1_all,
                                                                   iw2_all,
                                                                   nlm_tmp,
                                                                   pv,
                                                                   NLloc); // N stands for nonlocal.
                                                }
                                            }
                                            else
                                            {
                                                if (nlm_tmp != 0.0)
                                                {
                                                    NLloc[nnr + nnr_inner] += nlm_tmp;
                                                }
                                            }
                                        } // end nspin
                                    }     // calc_deri
                                    else  // calculate the derivative
                                    {
                                        if (nspin == 4)
                                        {
                                            std::vector<double> nlm_1 = (*nlm_cur2_f)[iw2_all][0];
                                            std::vector<std::vector<double>> nlm_2;
                                            nlm_2.resize(3);
                                            for (int i = 0; i < 3; i++)
                                            {
                                                nlm_2[i] = (*nlm_cur1_f)[iw1_all][i + 1];
                                            }
                                            std::complex<double> nlm[4][3] = {ModuleBase::ZERO};
                                            int is0 = (j - j0 * npol) + (k - k0 * npol) * 2;
                                            for (int no = 0; no < ucell.atoms[T0].ncpp.non_zero_count_soc[is0]; no++)
                                            {
                                                const int p1 = ucell.atoms[T0].ncpp.index1_soc[is0][no];
                                                const int p2 = ucell.atoms[T0].ncpp.index2_soc[is0][no];
                                                if (is0 == 0)
                                                {
                                                    fsr.DHloc_fixedR_x[nnr + nnr_inner]
                                                        += nlm_2[0][p1] * nlm_1[p2]
                                                           * (ucell.atoms[T0].ncpp.d_so(0, p2, p1).real()
                                                              + ucell.atoms[T0].ncpp.d_so(3, p2, p1).real())
                                                           * 0.5;
                                                    fsr.DHloc_fixedR_y[nnr + nnr_inner]
                                                        += nlm_2[1][p1] * nlm_1[p2]
                                                           * (ucell.atoms[T0].ncpp.d_so(0, p2, p1).real()
                                                              + ucell.atoms[T0].ncpp.d_so(3, p2, p1).real())
                                                           * 0.5;
                                                    fsr.DHloc_fixedR_z[nnr + nnr_inner]
                                                        += nlm_2[2][p1] * nlm_1[p2]
                                                           * (ucell.atoms[T0].ncpp.d_so(0, p2, p1).real()
                                                              + ucell.atoms[T0].ncpp.d_so(3, p2, p1).real())
                                                           * 0.5;
                                                }
                                                else if (is0 == 1)
                                                {
                                                    fsr.DHloc_fixedR_x[nnr + nnr_inner]
                                                        += nlm_2[0][p1] * nlm_1[p2]
                                                           * (ucell.atoms[T0].ncpp.d_so(1, p2, p1).real()
                                                              + ucell.atoms[T0].ncpp.d_so(2, p2, p1).real())
                                                           * 0.5;
                                                    fsr.DHloc_fixedR_y[nnr + nnr_inner]
                                                        += nlm_2[1][p1] * nlm_1[p2]
                                                           * (ucell.atoms[T0].ncpp.d_so(1, p2, p1).real()
                                                              + ucell.atoms[T0].ncpp.d_so(2, p2, p1).real())
                                                           * 0.5;
                                                    fsr.DHloc_fixedR_z[nnr + nnr_inner]
                                                        += nlm_2[2][p1] * nlm_1[p2]
                                                           * (ucell.atoms[T0].ncpp.d_so(1, p2, p1).real()
                                                              + ucell.atoms[T0].ncpp.d_so(2, p2, p1).real())
                                                           * 0.5;
                                                }
                                                else if (is0 == 2)
                                                {
                                                    fsr.DHloc_fixedR_x[nnr + nnr_inner]
                                                        += nlm_2[0][p1] * nlm_1[p2]
                                                           * (-ucell.atoms[T0].ncpp.d_so(1, p2, p1).imag()
                                                              + ucell.atoms[T0].ncpp.d_so(2, p2, p1).imag())
                                                           * 0.5;
                                                    fsr.DHloc_fixedR_y[nnr + nnr_inner]
                                                        += nlm_2[1][p1] * nlm_1[p2]
                                                           * (-ucell.atoms[T0].ncpp.d_so(1, p2, p1).imag()
                                                              + ucell.atoms[T0].ncpp.d_so(2, p2, p1).imag())
                                                           * 0.5;
                                                    fsr.DHloc_fixedR_z[nnr + nnr_inner]
                                                        += nlm_2[2][p1] * nlm_1[p2]
                                                           * (-ucell.atoms[T0].ncpp.d_so(1, p2, p1).imag()
                                                              + ucell.atoms[T0].ncpp.d_so(2, p2, p1).imag())
                                                           * 0.5;
                                                }
                                                else if (is0 == 3)
                                                {
                                                    fsr.DHloc_fixedR_x[nnr + nnr_inner]
                                                        += nlm_2[0][p1] * nlm_1[p2]
                                                           * (ucell.atoms[T0].ncpp.d_so(0, p2, p1).real()
                                                              - ucell.atoms[T0].ncpp.d_so(3, p2, p1).real())
                                                           * 0.5;
                                                    fsr.DHloc_fixedR_y[nnr + nnr_inner]
                                                        += nlm_2[1][p1] * nlm_1[p2]
                                                           * (ucell.atoms[T0].ncpp.d_so(0, p2, p1).real()
                                                              - ucell.atoms[T0].ncpp.d_so(3, p2, p1).real())
                                                           * 0.5;
                                                    fsr.DHloc_fixedR_z[nnr + nnr_inner]
                                                        += nlm_2[2][p1] * nlm_1[p2]
                                                           * (ucell.atoms[T0].ncpp.d_so(0, p2, p1).real()
                                                              - ucell.atoms[T0].ncpp.d_so(3, p2, p1).real())
                                                           * 0.5;
                                                }
                                            }
                                        }
                                        else if (nspin == 1 || nspin == 2)
                                        {
                                            if (gamma_only_local)
                                            {
                                                double nlm[3] = {0, 0, 0};

                                                // sum all projectors for one atom.
                                                std::vector<double> nlm_1 = (*nlm_cur1_f)[iw1_all][0];
                                                std::vector<std::vector<double>> nlm_2;
                                                nlm_2.resize(3);
                                                for (int i = 0; i < 3; i++)
                                                {
                                                    nlm_2[i] = (*nlm_cur2_f)[iw2_all][i + 1];
                                                }

                                                assert(nlm_1.size() == nlm_2[0].size());

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

                                                LCAO_domain::set_force(pv,
                                                                       iw1_all,
                                                                       iw2_all,
                                                                       nlm[0],
                                                                       nlm[1],
                                                                       nlm[2],
                                                                       'N',
                                                                       fsr.DSloc_x,
                                                                       fsr.DSloc_y,
                                                                       fsr.DSloc_z,
                                                                       fsr.DHloc_fixed_x,
                                                                       fsr.DHloc_fixed_y,
                                                                       fsr.DHloc_fixed_z);
                                            }
                                            else
                                            {
                                                // mohan change the order on 2011-06-17
                                                // origin: < psi1 | beta > < beta | dpsi2/dtau >
                                                // now: < psi1/dtau | beta > < beta | psi2 >
                                                double nlm[3] = {0, 0, 0};

                                                // sum all projectors for one atom.
                                                std::vector<double> nlm_1 = (*nlm_cur2_f)[iw2_all][0];
                                                std::vector<std::vector<double>> nlm_2;
                                                nlm_2.resize(3);
                                                for (int i = 0; i < 3; i++)
                                                {
                                                    nlm_2[i] = (*nlm_cur1_f)[iw1_all][i + 1];
                                                }

                                                assert(nlm_1.size() == nlm_2[0].size());

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

                                                fsr.DHloc_fixedR_x[nnr + nnr_inner] += nlm[0];
                                                fsr.DHloc_fixedR_y[nnr + nnr_inner] += nlm[1];
                                                fsr.DHloc_fixedR_z[nnr + nnr_inner] += nlm[2];
                                            }
                                        }
                                        else
                                        {
                                            ModuleBase::WARNING_QUIT("LCAO_domain::build_Nonlocal_mu_new",
                                                                     "nspin must be 1, 2 or 4");
                                        }
                                    } //! calc_deri
                                    nnr_inner++;
                                } // k
                            }     // j
                        }         // ad0

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
                    } // end is_adj
                }     // ad2
            }         // I1
        }             // T1
#ifdef _OPENMP
    }
#endif
    if (!gamma_only_local)
    {
        if (total_nnr != pv.nnr)
        {
            GlobalV::ofs_running << " nr=" << total_nnr << std::endl;
            GlobalV::ofs_running << " pv->nnr=" << pv.nnr << std::endl;
            ModuleBase::WARNING_QUIT("LCAO_domain::build_Nonlocal_mu_new", "nnr!=LNNR.nnr");
        }
    }

    ModuleBase::timer::tick("LCAO_domain", "vnl_mu_new");
    return;
}

} // namespace LCAO_domain
